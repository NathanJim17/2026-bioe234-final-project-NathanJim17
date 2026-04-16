# Tool #5: CRISPR Validation Planner

**Objective:** Transform upstream CRISPR design decisions (edit type, target locus, and delivery method) into a structured, biologically rigorous validation checklist — recommending the exact sequencing strategies, molecular assays, experimental controls, and timeline a researcher needs to confirm their edit worked.

---

## The Core Functions

**`ValidationPlanner.initiate()`**: One-time setup that loads three immutable rule dictionaries into memory:
- **Edit rules** — Maps each of the 7 supported edit types to its validation strategy (checklist of assays, base controls, primary sequencing recommendation, estimated timeline, and contextual notes).
- **Delivery alias map** — A ~25-entry lookup table that normalizes freeform delivery method strings (e.g., `"RNP electroporation"`, `"lipofectamine"`, `"lenti"`) into 6 canonical categories. This allows seamless integration with the upstream Delivery Strategy Advisor tool.
- **Delivery controls** — Maps each canonical delivery category to its method-specific experimental controls (e.g., lentiviral requires integration copy-number qPCR; RNP requires a mock electroporation control).

**`ValidationPlanner.run(edit_type, locus, delivery_method)`**: The core logic. It validates inputs, normalizes the delivery string, looks up edit-type-specific rules, substitutes the locus name into all templated text fields, appends delivery-specific controls, and assembles the final validation plan dictionary.

---

## The Input

Gemini will parse the user's request and send the function a JSON object:

- **`edit_type`** (String, required): The category of CRISPR edit. One of: `knockout`, `knock-in`, `base_edit`, `prime_edit`, `CRISPRa`, `CRISPRi`, `large_deletion`.
- **`locus`** (String, required): The target gene symbol or genomic locus (e.g., `"BRCA1"`, `"TP53"`, `"EMX1"`).
- **`delivery_method`** (String, optional): How the CRISPR reagents are delivered. Accepts freeform strings from the upstream Delivery Strategy Advisor (e.g., `"RNP electroporation"`, `"lipofection"`, `"lentiviral transduction"`). Normalized internally. Defaults to `"unspecified"`.

---

## The Output

The function returns a structured JSON validation plan. This provides the data that allows Gemini to walk the user through a complete post-editing validation strategy.

```json
{
  "edit_type": "knockout",
  "locus": "BRCA1",
  "delivery_method": "RNP",
  "recommended_sequencing": "Sanger sequencing",
  "validation_checklist": [
    {
      "step": 1,
      "category": "genomic",
      "assay": "T7 Endonuclease I (T7E1) mismatch assay",
      "purpose": "Rapid screening for heteroduplex DNA indicating successful indel formation",
      "priority": "required",
      "expected_result": "Cleavage bands on agarose gel at expected fragment sizes"
    }
  ],
  "controls": [
    {
      "control_type": "negative",
      "description": "Non-targeting gRNA control"
    },
    {
      "control_type": "delivery",
      "description": "Mock electroporation control"
    }
  ],
  "estimated_timeline_days": 14,
  "notes": [
    "Consider clonal isolation for homozygous knockout confirmation."
  ]
}
```

---

## The Final Output (The "Deliverable")

The final output of this MCP tool is a **Validation Checklist and Assay Strategy Report**.

When the end-user asks Gemini, *"Design a CRISPR knockout of BRCA1 in HEK293 cells using SpCas9,"* the upstream tools handle target selection, guide design, off-target analysis, and delivery recommendation. This tool provides the data that allows Gemini to say:

> *"To validate your knockout, here is your recommended validation plan: Start with a T7E1 mismatch assay for rapid screening (Day 1-3), then confirm the exact indel by Sanger sequencing (Day 3-7). Run a Western blot against BRCA1 to verify protein loss (Day 7-14). Your controls should include a non-targeting gRNA and untreated wild-type cells. Since you're using RNP electroporation, also include a mock electroporation control. Estimated total timeline: ~14 days."*

---

## Edit-Type Validation Logic

Each edit type maps to a distinct validation strategy reflecting the biological mechanism:

| Edit Type | Primary Assay | Key Validation Focus | Timeline |
|---|---|---|---|
| **Knockout** | Sanger sequencing | T7E1 screen → Sanger → Western blot for protein loss | 14 days |
| **Knock-in** | Sanger across junctions | 5'/3' junction PCR → full-insert sequencing → expression check | 21 days |
| **Base edit** | NGS amplicon-seq | Precise conversion confirmation + bystander edit quantification | 10 days |
| **Prime edit** | NGS amplicon-seq | Precise edit + indel byproducts + scaffold incorporation check | 10 days |
| **CRISPRa** | RT-qPCR | Transcript upregulation → Western blot → RNA-seq specificity | 7 days |
| **CRISPRi** | RT-qPCR | Transcript downregulation → Western blot → RNA-seq specificity | 7 days |
| **Large deletion** | Sanger of junction | Deletion-spanning PCR → breakpoint sequencing → copy-number qPCR | 14 days |

---

## Delivery Method Normalization

The tool accepts freeform delivery strings from the upstream Delivery Strategy Advisor and normalizes them internally:

| Upstream Input Examples | Canonical Category | Delivery-Specific Controls |
|---|---|---|
| `"RNP electroporation"`, `"ribonucleoprotein"`, `"rnp"` | **RNP** | Mock electroporation (no cargo) |
| `"lipofection"`, `"lipofectamine"`, `"lipid transfection"` | **lipofection** | Reagent-only control (no nucleic acid) |
| `"plasmid transfection"`, `"plasmid"`, `"PEI"` | **plasmid** | Empty-vector control + backbone integration check |
| `"lentiviral transduction"`, `"lentivirus"`, `"lenti"` | **lentiviral** | Untransduced cells + integration copy-number qPCR |
| `"AAV transduction"`, `"adeno-associated virus"` | **AAV** | ITR-junction PCR + empty-capsid control |
| Unrecognized strings | **unspecified** | Generic negative control + warning note |

---

## Literature Basis for Validation Rules

The assay recommendations encoded in this tool are derived from the following foundational CRISPR literature:

### General CRISPR Validation Protocols
- **Ran, F.A. et al.** (2013). "Genome engineering using the CRISPR-Cas9 system." *Nature Protocols*, 8(11), 2281–2308. — Established the standard T7E1 → Sanger → Western blot validation workflow for Cas9-based knockouts and knock-ins.
- **Cong, L. et al.** (2013). "Multiplex Genome Engineering Using CRISPR/Cas Systems." *Science*, 339(6121), 819–823. — Original demonstration of CRISPR gene editing in mammalian cells; T7E1 and Sanger sequencing used as primary validation.
- **Mali, P. et al.** (2013). "RNA-Guided Human Genome Engineering via Cas9." *Science*, 339(6121), 823–826. — Parallel foundational paper establishing Sanger-based validation of Cas9 edits.

### NGS-Based Quantification
- **Clement, K. et al.** (2019). "CRISPResso2 provides accurate and rapid genome editing sequence analysis." *Nature Biotechnology*, 37(3), 224–226. — The standard tool for NGS amplicon analysis of CRISPR edits; basis for recommending NGS for efficiency quantification.

### Base Editing
- **Komor, A.C. et al.** (2016). "Programmable editing of a target base in genomic DNA without double-stranded DNA cleavage." *Nature*, 533(7603), 420–424. — Introduced cytosine base editors (CBE); established the need to check bystander edits within the ~4-8 nt editing window.
- **Gaudelli, N.M. et al.** (2017). "Programmable base editing of A·T to G·C in genomic DNA without DNA cleavage." *Nature*, 551(7681), 464–471. — Introduced adenine base editors (ABE); demonstrated position-dependent editing efficiency requiring NGS quantification.

### Prime Editing
- **Anzalone, A.V. et al.** (2019). "Search-and-replace genome editing without double-strand breaks or donor DNA." *Nature*, 576(7785), 149–157. — Introduced prime editing; documented scaffold incorporation as a failure mode requiring specific validation.

### CRISPRa / CRISPRi
- **Gilbert, L.A. et al.** (2014). "Genome-Scale CRISPR-Mediated Control of Gene Repression and Activation." *Cell*, 159(3), 647–661. — Established RT-qPCR and RNA-seq as the primary validation readouts for dCas9-based transcriptional modulation.
- **Konermann, S. et al.** (2015). "Genome-scale transcriptional activation by an engineered CRISPR-Cas9 complex." *Nature*, 517(7536), 583–588. — Demonstrated synergistic activation mediator (SAM) CRISPRa system; validated by RT-qPCR fold-change.

### Large Deletions
- **Canver, M.C. et al.** (2014). "Characterization of genomic deletion efficiency mediated by clustered regularly interspaced short palindromic repeats (CRISPR)/Cas9 nuclease system in mammalian cells." *Journal of Biological Chemistry*, 289(31), 21312–21324. — Characterized deletion-spanning PCR and inversion byproducts from dual-gRNA strategies.
- **Kraft, K. et al.** (2015). "Deletions, Inversions, Duplications: Engineering of Structural Variants using CRISPR/Cas in Mice." *Cell Reports*, 10(5), 833–839. — Documented inversions and duplications as common byproducts of large-deletion strategies requiring dedicated detection primers.

---

## Files

| Deliverable | Path |
|---|---|
| Python class | `modules/seq_basics/tools/validation_planner.py` |
| MCP JSON wrapper | `modules/seq_basics/tools/validation_planner.json` |
| Pytests (41 tests) | `tests/test_tools.py` |
| Test prompts | `modules/seq_basics/tools/prompts.json` |
