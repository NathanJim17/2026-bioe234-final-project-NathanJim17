"""CRISPR Validation Planner — recommends a sequencing/assay strategy
for validating a CRISPR experiment based on the edit type, target locus,
and delivery method.
"""


class ValidationPlanner:
    """
    Description:
        Takes upstream CRISPR design decisions (edit type, target locus, and
        optionally the delivery method) and returns a structured validation
        checklist recommending sequencing strategies, molecular assays,
        experimental controls, and an estimated timeline.

    Input:
        edit_type       (str):           The category of CRISPR edit.
                                         One of: knockout, knock-in, base_edit,
                                         prime_edit, CRISPRa, CRISPRi,
                                         large_deletion.
        locus           (str):           Target gene symbol or genomic locus
                                         (e.g. "BRCA1", "TP53").
        delivery_method (str, optional): How CRISPR reagents are delivered.
                                         Accepts freeform strings from the
                                         upstream Delivery Strategy Advisor;
                                         normalized internally. Default:
                                         "unspecified".

    Output:
        dict: A validation plan containing:
              - edit_type, locus, delivery_method (echoed back, normalized)
              - recommended_sequencing (str)
              - validation_checklist (list[dict])
              - controls (list[dict])
              - estimated_timeline_days (int)
              - notes (list[str])

    Tests:
        - Case:
            Input: edit_type="knockout", locus="BRCA1"
            Expected Output: dict with recommended_sequencing="Sanger sequencing"
                             and >=3 checklist items
            Description: Typical knockout produces genomic + protein validation.
        - Case:
            Input: edit_type="CRISPRa", locus="TP53"
            Expected Output: dict with no genomic-category checklist items
            Description: Activation edits do not alter DNA; only transcript/
                         protein assays should appear.
        - Case:
            Input: edit_type="knock-in", locus="EMX1",
                   delivery_method="lentiviral transduction"
            Expected Output: dict where delivery_method is normalized to
                             "lentiviral" and controls include an
                             integration-related entry
            Description: Freeform delivery string is normalized; delivery-
                         specific controls are appended.
        - Case:
            Input: edit_type="knockout", locus="BRCA1",
                   delivery_method="some_unknown_method"
            Expected Output: dict with delivery_method="unspecified" and a
                             note mentioning the unrecognized input
            Description: Unrecognized delivery method falls back gracefully.
        - Case:
            Input: edit_type="invalid_type", locus="BRCA1"
            Expected Exception: ValueError
            Description: Unsupported edit_type must raise ValueError.
        - Case:
            Input: edit_type="knockout", locus=""
            Expected Exception: ValueError
            Description: Empty locus string must raise ValueError.
    """

    _edit_rules: dict[str, dict]
    _delivery_alias: dict[str, str]
    _delivery_controls: dict[str, list[dict]]
    _valid_edit_types: set[str]

    def initiate(self) -> None:
        """One-time setup: populate rule dictionaries for all edit types
        and delivery methods.  All instance state set here is treated as
        immutable during run()."""

        self._valid_edit_types = {
            "knockout", "knock-in", "base_edit", "prime_edit",
            "CRISPRa", "CRISPRi", "large_deletion",
        }

        self._edit_rules = _build_edit_rules()
        self._delivery_alias = _build_delivery_alias_map()
        self._delivery_controls = _build_delivery_controls()

    def run(
        self,
        edit_type: str,
        locus: str,
        delivery_method: str = "unspecified",
    ) -> dict:
        """Return a structured validation plan for the given CRISPR edit."""

        if not locus or not locus.strip():
            raise ValueError("locus must be a non-empty string")

        edit_key = edit_type.strip().lower()

        key_to_canonical: dict[str, str] = {
            et.lower(): et for et in self._valid_edit_types
        }
        if edit_key not in key_to_canonical:
            allowed = ", ".join(sorted(self._valid_edit_types))
            raise ValueError(
                f"Unsupported edit_type '{edit_type}'. "
                f"Must be one of: {allowed}"
            )
        canonical_edit = key_to_canonical[edit_key]

        normalized_delivery, delivery_note = self._normalize_delivery(
            delivery_method
        )

        rule = self._edit_rules[canonical_edit]

        locus_clean = locus.strip()

        checklist: list[dict] = []
        for i, item in enumerate(rule["checklist"]):
            entry = {
                "step": i + 1,
                "category": item["category"],
                "assay": item["assay"],
                "purpose": item["purpose"].replace("{locus}", locus_clean),
                "priority": item["priority"],
                "expected_result": item["expected_result"].replace(
                    "{locus}", locus_clean
                ),
            }
            checklist.append(entry)

        controls: list[dict] = list(rule["base_controls"])
        controls.extend(self._delivery_controls.get(normalized_delivery, []))

        notes: list[str] = list(rule["notes"])
        if delivery_note:
            notes.append(delivery_note)
        notes = [n.replace("{locus}", locus_clean) for n in notes]

        return {
            "edit_type": canonical_edit,
            "locus": locus_clean,
            "delivery_method": normalized_delivery,
            "recommended_sequencing": rule["recommended_sequencing"],
            "validation_checklist": checklist,
            "controls": controls,
            "estimated_timeline_days": rule["estimated_timeline_days"],
            "notes": notes,
        }

    def _normalize_delivery(self, raw: str) -> tuple[str, str]:
        """Map a freeform delivery string to a canonical category.

        Returns (canonical_name, optional_warning_note).
        """
        key = raw.strip().lower()
        if key in self._delivery_alias:
            return self._delivery_alias[key], ""
        if key in ("", "unspecified"):
            return "unspecified", ""
        return (
            "unspecified",
            f"Delivery method '{raw}' was not recognized and was treated "
            f"as 'unspecified'. Supported values include: RNP, lipofection, "
            f"plasmid, lentiviral, AAV.",
        )


# Rule-building helpers (called once by initiate)

def _build_edit_rules() -> dict[str, dict]:
    """Return the master rule dictionary keyed by canonical edit type."""

    return {
        "knockout": {
            "recommended_sequencing": "Sanger sequencing",
            "estimated_timeline_days": 14,
            "checklist": [
                {
                    "category": "genomic",
                    "assay": "T7 Endonuclease I (T7E1) mismatch assay",
                    "purpose": "Rapid screening for heteroduplex DNA "
                               "indicating successful indel formation",
                    "priority": "required",
                    "expected_result": "Cleavage bands on agarose gel at "
                                       "expected fragment sizes",
                },
                {
                    "category": "genomic",
                    "assay": "Sanger sequencing of target amplicon",
                    "purpose": "Confirm the exact indel sequence at the "
                               "cut site",
                    "priority": "required",
                    "expected_result": "Mixed chromatogram peaks downstream "
                                       "of cut site (bulk) or clean indel "
                                       "sequence (clonal)",
                },
                {
                    "category": "genomic",
                    "assay": "NGS amplicon sequencing",
                    "purpose": "Quantify editing efficiency and characterize "
                               "the indel spectrum",
                    "priority": "recommended",
                    "expected_result": "High percentage of reads containing "
                                       "frameshifting indels at the target",
                },
                {
                    "category": "protein",
                    "assay": "Western blot",
                    "purpose": "Confirm loss of target protein expression",
                    "priority": "recommended",
                    "expected_result": "Absent or significantly reduced band "
                                       "at expected molecular weight for "
                                       "{locus}",
                },
                {
                    "category": "transcript",
                    "assay": "RT-qPCR",
                    "purpose": "Assess transcript-level reduction via "
                               "nonsense-mediated decay",
                    "priority": "optional",
                    "expected_result": "Significant reduction in {locus} "
                                       "mRNA relative to control",
                },
            ],
            "base_controls": [
                {
                    "control_type": "negative",
                    "description": "Non-targeting gRNA control (same "
                                   "delivery, scrambled/non-targeting guide)",
                },
                {
                    "control_type": "negative",
                    "description": "Untreated wild-type cells",
                },
            ],
            "notes": [
                "Consider clonal isolation for homozygous knockout "
                "confirmation before functional studies.",
                "If antibody for {locus} is unavailable, RT-qPCR or "
                "functional assays can substitute for Western blot.",
            ],
        },

        "knock-in": {
            "recommended_sequencing": "Sanger sequencing across junctions",
            "estimated_timeline_days": 21,
            "checklist": [
                {
                    "category": "genomic",
                    "assay": "5-prime junction PCR",
                    "purpose": "Confirm correct integration at the 5-prime "
                               "boundary of the insert",
                    "priority": "required",
                    "expected_result": "Single band of expected size spanning "
                                       "the genomic–insert junction",
                },
                {
                    "category": "genomic",
                    "assay": "3-prime junction PCR",
                    "purpose": "Confirm correct integration at the 3-prime "
                               "boundary of the insert",
                    "priority": "required",
                    "expected_result": "Single band of expected size spanning "
                                       "the insert–genomic junction",
                },
                {
                    "category": "genomic",
                    "assay": "Sanger sequencing of full insert",
                    "purpose": "Verify the inserted sequence is intact and "
                               "mutation-free",
                    "priority": "required",
                    "expected_result": "Clean chromatogram matching the "
                                       "expected donor sequence",
                },
                {
                    "category": "genomic",
                    "assay": "Long-read sequencing or Southern blot",
                    "purpose": "Rule out concatemeric insertions, partial "
                               "integrations, or backbone incorporation",
                    "priority": "recommended",
                    "expected_result": "Single-copy, full-length integration "
                                       "at the target locus",
                },
                {
                    "category": "protein",
                    "assay": "Western blot or reporter assay",
                    "purpose": "Confirm expression of the knock-in product",
                    "priority": "recommended",
                    "expected_result": "Band or signal at the expected size/"
                                       "wavelength for the inserted gene "
                                       "product",
                },
            ],
            "base_controls": [
                {
                    "control_type": "negative",
                    "description": "Non-targeting gRNA + donor template "
                                   "control (rules out random integration)",
                },
                {
                    "control_type": "negative",
                    "description": "Untreated wild-type cells",
                },
                {
                    "control_type": "positive",
                    "description": "gRNA-only control without donor template "
                                   "(confirms cutting at the locus)",
                },
            ],
            "notes": [
                "Screen multiple clones — HDR efficiency is typically low "
                "(1-30%) depending on cell type.",
                "For large inserts (>2 kb), long-read sequencing (PacBio or "
                "ONT) is strongly recommended over Southern blot.",
            ],
        },

        "base_edit": {
            "recommended_sequencing": "NGS amplicon sequencing",
            "estimated_timeline_days": 10,
            "checklist": [
                {
                    "category": "genomic",
                    "assay": "Sanger sequencing of target amplicon",
                    "purpose": "Confirm the precise base conversion at the "
                               "target position",
                    "priority": "required",
                    "expected_result": "Clean single-peak change at the "
                                       "target nucleotide (clonal) or mixed "
                                       "peak (bulk)",
                },
                {
                    "category": "genomic",
                    "assay": "NGS amplicon sequencing",
                    "purpose": "Quantify editing efficiency and detect "
                               "bystander edits within the editing window",
                    "priority": "required",
                    "expected_result": "High on-target conversion rate with "
                                       "minimal bystander edits at adjacent "
                                       "positions",
                },
                {
                    "category": "genomic",
                    "assay": "Off-target site sequencing",
                    "purpose": "Check computationally predicted off-target "
                               "sites for unintended base conversions",
                    "priority": "recommended",
                    "expected_result": "No significant editing above "
                                       "background at top predicted "
                                       "off-target loci",
                },
                {
                    "category": "protein",
                    "assay": "Western blot",
                    "purpose": "Confirm altered protein if the edit is a "
                               "missense change",
                    "priority": "optional",
                    "expected_result": "Shifted band or altered abundance "
                                       "consistent with the amino acid change",
                },
            ],
            "base_controls": [
                {
                    "control_type": "negative",
                    "description": "Catalytically dead base editor control "
                                   "(same gRNA, dead deaminase)",
                },
                {
                    "control_type": "negative",
                    "description": "Untreated wild-type cells",
                },
            ],
            "notes": [
                "ABE (adenine) and CBE (cytosine) editors have different "
                "bystander profiles — check all C's or A's within the "
                "~4-8 nt editing window.",
                "Indel byproducts are rare with base editors but should "
                "be monitored via NGS.",
            ],
        },

        "prime_edit": {
            "recommended_sequencing": "NGS amplicon sequencing",
            "estimated_timeline_days": 10,
            "checklist": [
                {
                    "category": "genomic",
                    "assay": "Sanger sequencing of target amplicon",
                    "purpose": "Confirm the precise intended edit at the "
                               "target site",
                    "priority": "required",
                    "expected_result": "Clean sequence showing the designed "
                                       "edit (clonal) or mixed signal (bulk)",
                },
                {
                    "category": "genomic",
                    "assay": "NGS amplicon sequencing",
                    "purpose": "Quantify prime editing efficiency and "
                               "detect indel byproducts",
                    "priority": "required",
                    "expected_result": "High fraction of reads with the "
                                       "precise intended edit and low "
                                       "indel rate",
                },
                {
                    "category": "genomic",
                    "assay": "Scaffold incorporation check",
                    "purpose": "Detect unintended incorporation of the "
                               "pegRNA scaffold sequence at the target site",
                    "priority": "recommended",
                    "expected_result": "No reads containing pegRNA scaffold "
                                       "sequence at the edit locus",
                },
                {
                    "category": "genomic",
                    "assay": "Off-target site sequencing",
                    "purpose": "Assess off-target prime editing or nicking "
                               "at predicted sites",
                    "priority": "recommended",
                    "expected_result": "No significant editing above "
                                       "background at off-target loci",
                },
            ],
            "base_controls": [
                {
                    "control_type": "negative",
                    "description": "Non-targeting pegRNA control (same PE "
                                   "protein, scrambled pegRNA)",
                },
                {
                    "control_type": "negative",
                    "description": "Untreated wild-type cells",
                },
            ],
            "notes": [
                "Prime editing generates fewer indel byproducts than Cas9 "
                "nuclease cutting, but they should still be quantified.",
                "PE3 strategies (nicking guide) may increase efficiency but "
                "also increase indel rates — monitor both.",
            ],
        },

        "CRISPRa": {
            "recommended_sequencing": "RT-qPCR",
            "estimated_timeline_days": 7,
            "checklist": [
                {
                    "category": "transcript",
                    "assay": "RT-qPCR",
                    "purpose": "Measure fold-change in {locus} transcript "
                               "levels upon activation",
                    "priority": "required",
                    "expected_result": "Significant upregulation (typically "
                                       ">2-fold) of {locus} mRNA relative "
                                       "to non-targeting control",
                },
                {
                    "category": "protein",
                    "assay": "Western blot",
                    "purpose": "Confirm increased protein expression of "
                               "{locus}",
                    "priority": "recommended",
                    "expected_result": "Increased band intensity at expected "
                                       "molecular weight for {locus}",
                },
                {
                    "category": "transcript",
                    "assay": "RNA-seq",
                    "purpose": "Assess genome-wide transcriptional "
                               "specificity of the activation",
                    "priority": "optional",
                    "expected_result": "{locus} is the top differentially "
                                       "expressed gene with minimal "
                                       "off-target transcriptional changes",
                },
                {
                    "category": "functional",
                    "assay": "Phenotypic or functional assay",
                    "purpose": "Confirm downstream biological effect of "
                               "{locus} overexpression",
                    "priority": "optional",
                    "expected_result": "Measurable phenotypic change "
                                       "consistent with {locus} gain of "
                                       "function",
                },
            ],
            "base_controls": [
                {
                    "control_type": "negative",
                    "description": "Non-targeting gRNA with the same dCas9 "
                                   "activator system",
                },
                {
                    "control_type": "negative",
                    "description": "Untreated wild-type cells",
                },
            ],
            "notes": [
                "CRISPRa does not alter the DNA sequence — genomic "
                "sequencing is not necessary for validation.",
                "Multiple gRNAs tiling the promoter region of {locus} can "
                "be tested to maximize activation.",
            ],
        },

        "CRISPRi": {
            "recommended_sequencing": "RT-qPCR",
            "estimated_timeline_days": 7,
            "checklist": [
                {
                    "category": "transcript",
                    "assay": "RT-qPCR",
                    "purpose": "Measure reduction in {locus} transcript "
                               "levels upon repression",
                    "priority": "required",
                    "expected_result": "Significant downregulation (typically "
                                       ">50% reduction) of {locus} mRNA "
                                       "relative to non-targeting control",
                },
                {
                    "category": "protein",
                    "assay": "Western blot",
                    "purpose": "Confirm decreased protein expression of "
                               "{locus}",
                    "priority": "recommended",
                    "expected_result": "Reduced band intensity at expected "
                                       "molecular weight for {locus}",
                },
                {
                    "category": "transcript",
                    "assay": "RNA-seq",
                    "purpose": "Assess genome-wide transcriptional "
                               "specificity of the repression",
                    "priority": "optional",
                    "expected_result": "{locus} is among the top "
                                       "downregulated genes with minimal "
                                       "off-target transcriptional changes",
                },
                {
                    "category": "functional",
                    "assay": "Phenotypic or functional assay",
                    "purpose": "Confirm downstream biological effect of "
                               "{locus} repression",
                    "priority": "optional",
                    "expected_result": "Measurable phenotypic change "
                                       "consistent with {locus} loss of "
                                       "function",
                },
            ],
            "base_controls": [
                {
                    "control_type": "negative",
                    "description": "Non-targeting gRNA with the same "
                                   "dCas9-KRAB repressor system",
                },
                {
                    "control_type": "negative",
                    "description": "Untreated wild-type cells",
                },
            ],
            "notes": [
                "CRISPRi does not alter the DNA sequence — genomic "
                "sequencing is not necessary for validation.",
                "Knockdown is reversible upon gRNA removal; confirm "
                "sustained repression over the experimental time course.",
            ],
        },

        "large_deletion": {
            "recommended_sequencing": "Sanger sequencing of deletion junction",
            "estimated_timeline_days": 14,
            "checklist": [
                {
                    "category": "genomic",
                    "assay": "Deletion-spanning PCR",
                    "purpose": "Confirm the expected amplicon size shift "
                               "indicating successful deletion",
                    "priority": "required",
                    "expected_result": "Shorter PCR product consistent with "
                                       "the expected deletion size, alongside "
                                       "or instead of the wild-type band",
                },
                {
                    "category": "genomic",
                    "assay": "Sanger sequencing of breakpoint junction",
                    "purpose": "Confirm the exact deletion boundaries and "
                               "any microhomology-mediated repair junctions",
                    "priority": "required",
                    "expected_result": "Clean chromatogram showing the "
                                       "expected junction sequence",
                },
                {
                    "category": "genomic",
                    "assay": "qPCR or ddPCR within deleted region",
                    "purpose": "Quantitatively confirm copy-number loss "
                               "in the deleted segment",
                    "priority": "recommended",
                    "expected_result": "~50% reduction (heterozygous) or "
                                       "~100% loss (homozygous) of signal "
                                       "from an internal primer pair",
                },
                {
                    "category": "genomic",
                    "assay": "FISH or cytogenetic analysis",
                    "purpose": "Visualize large-scale deletions at the "
                               "chromosomal level",
                    "priority": "optional",
                    "expected_result": "Loss of fluorescent signal at the "
                                       "{locus} region on target chromosome",
                },
                {
                    "category": "protein",
                    "assay": "Western blot",
                    "purpose": "Confirm loss of protein if the deletion "
                               "disrupts a coding gene",
                    "priority": "optional",
                    "expected_result": "Absent or reduced band at expected "
                                       "molecular weight for {locus}",
                },
            ],
            "base_controls": [
                {
                    "control_type": "negative",
                    "description": "Single-gRNA-only controls (each guide "
                                   "individually, to rule out single-cut "
                                   "artifacts)",
                },
                {
                    "control_type": "negative",
                    "description": "Untreated wild-type cells",
                },
            ],
            "notes": [
                "Large deletions require two gRNAs flanking the region — "
                "validate each guide's cutting individually before "
                "combining.",
                "Inversions are a common byproduct of dual-gRNA strategies; "
                "design inversion-detecting primers as well.",
            ],
        },
    }


def _build_delivery_alias_map() -> dict[str, str]:
    """Return a lowercase-key -> canonical-name mapping for delivery methods.

    Covers common phrasings that the upstream Delivery Strategy Advisor
    or a human user might provide.
    """
    aliases: dict[str, str] = {}

    for term in ("rnp", "rnp electroporation", "ribonucleoprotein",
                 "ribonucleoprotein electroporation"):
        aliases[term] = "RNP"

    for term in ("lipofection", "lipofectamine", "lipid transfection",
                 "lipid nanoparticle", "lnp"):
        aliases[term] = "lipofection"

    for term in ("plasmid", "plasmid transfection", "chemical transfection",
                 "capo4", "calcium phosphate", "pei"):
        aliases[term] = "plasmid"

    for term in ("lentiviral", "lentivirus", "lenti",
                 "lentiviral transduction"):
        aliases[term] = "lentiviral"

    for term in ("aav", "aav transduction", "adeno-associated virus",
                 "adeno-associated viral"):
        aliases[term] = "AAV"

    for term in ("unspecified", ""):
        aliases[term] = "unspecified"

    return aliases


def _build_delivery_controls() -> dict[str, list[dict]]:
    """Return delivery-method-specific controls to append to any plan."""
    return {
        "RNP": [
            {
                "control_type": "delivery",
                "description": "Mock electroporation control (electroporation "
                               "without RNP cargo) to assess electroporation-"
                               "induced toxicity",
            },
        ],
        "lipofection": [
            {
                "control_type": "delivery",
                "description": "Transfection-reagent-only control "
                               "(lipofectamine without nucleic acid cargo) "
                               "to assess reagent cytotoxicity",
            },
        ],
        "plasmid": [
            {
                "control_type": "delivery",
                "description": "Empty-vector control (backbone plasmid "
                               "without gRNA/Cas9 insert)",
            },
            {
                "control_type": "delivery",
                "description": "Check for random plasmid backbone "
                               "integration by PCR with backbone-specific "
                               "primers",
            },
        ],
        "lentiviral": [
            {
                "control_type": "delivery",
                "description": "Untransduced cells to control for "
                               "lentiviral toxicity and insertional effects",
            },
            {
                "control_type": "delivery",
                "description": "Integration copy-number qPCR to confirm "
                               "low-copy proviral integration",
            },
        ],
        "AAV": [
            {
                "control_type": "delivery",
                "description": "ITR-junction PCR to check for unintended "
                               "AAV vector integration at the target or "
                               "off-target sites",
            },
            {
                "control_type": "delivery",
                "description": "Empty-capsid AAV control to assess "
                               "viral-entry-related effects",
            },
        ],
        "unspecified": [],
    }


# ---------------------------------------------------------------------------
# Module-level alias
#
#   from modules.seq_basics.tools.validation_planner import validation_planner
#
# Creates ONE shared instance and exposes run() as a plain callable.
# ---------------------------------------------------------------------------
_instance = ValidationPlanner()
_instance.initiate()
validation_planner = _instance.run


if __name__ == "__main__":
    import json as _json

    plan = validation_planner("knockout", "BRCA1")
    print(_json.dumps(plan, indent=2))

    plan2 = validation_planner("CRISPRa", "TP53", "lipofection")
    print(_json.dumps(plan2, indent=2))
