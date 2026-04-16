"""
Unit tests for the seq_basics example tools.
 
:: Each tool is now a class following the Python Function Object Pattern
(initiate / run).  Tests cover both the canonical class interface AND the
module-level alias  (for example: `reverse_complement = _instance.run`)  so that
direct imports continue to work for students who prefer that style.
"""

import pytest

from modules.seq_basics.tools.translate import translate
from modules.seq_basics.tools.reverse_complement import reverse_complement
from modules.seq_basics.tools.validation_planner import (
    ValidationPlanner,
    validation_planner,
)


def test_reverse_complement_basic():
    assert reverse_complement("ATGC") == "GCAT"


def test_reverse_complement_ambiguity_codes():
    # Should not error for supported IUPAC subset
    assert reverse_complement("ATRYSWKMN")


def test_translate_basic():
    assert translate("ATGGCT") == "MA"


def test_translate_frame_validation():
    with pytest.raises(ValueError):
        translate("ATGGCT", frame=0)
    with pytest.raises(ValueError):
        translate("ATGGCT", frame=4)


def test_translate_with_coordinates_and_frame():
    # sequence: A ATG GCT AAA
    # start=1 → ATGGCTAAA
    # frame=1 → ATG GCT AAA → M A K
    assert translate("AATGGCTAAA", start=1, end=None, frame=1) == "MAK"


# Validation Planner

class TestValidationPlannerClassInterface:
    """Verify the class can be instantiated and used independently."""

    def setup_method(self):
        self.vp = ValidationPlanner()
        self.vp.initiate()

    def test_initiate_populates_rules(self):
        assert len(self.vp._edit_rules) == 7
        assert len(self.vp._delivery_alias) > 0
        assert "knockout" in self.vp._edit_rules

    def test_run_returns_dict(self):
        result = self.vp.run("knockout", "BRCA1")
        assert isinstance(result, dict)


class TestValidationPlannerOutputStructure:
    """Every plan must contain the required top-level keys and nested shapes."""

    REQUIRED_KEYS = {
        "edit_type", "locus", "delivery_method",
        "recommended_sequencing", "validation_checklist",
        "controls", "estimated_timeline_days", "notes",
    }

    CHECKLIST_KEYS = {
        "step", "category", "assay", "purpose", "priority", "expected_result",
    }

    def test_all_top_level_keys_present(self):
        plan = validation_planner("knockout", "BRCA1")
        assert self.REQUIRED_KEYS.issubset(plan.keys())

    def test_checklist_item_structure(self):
        plan = validation_planner("knockout", "BRCA1")
        for item in plan["validation_checklist"]:
            assert self.CHECKLIST_KEYS.issubset(item.keys())

    def test_controls_item_structure(self):
        plan = validation_planner("knockout", "BRCA1")
        for ctrl in plan["controls"]:
            assert "control_type" in ctrl
            assert "description" in ctrl

    def test_steps_are_sequential(self):
        plan = validation_planner("knock-in", "EMX1")
        steps = [item["step"] for item in plan["validation_checklist"]]
        assert steps == list(range(1, len(steps) + 1))


class TestValidationPlannerEditTypes:
    """Each supported edit type produces a biologically appropriate plan."""

    def test_knockout_has_genomic_assays(self):
        plan = validation_planner("knockout", "BRCA1")
        categories = {c["category"] for c in plan["validation_checklist"]}
        assert "genomic" in categories

    def test_knockout_recommends_sanger(self):
        plan = validation_planner("knockout", "BRCA1")
        assert plan["recommended_sequencing"] == "Sanger sequencing"

    def test_knockin_has_junction_pcr(self):
        plan = validation_planner("knock-in", "EMX1")
        assays = [c["assay"] for c in plan["validation_checklist"]]
        junction_assays = [a for a in assays if "junction" in a.lower()]
        assert len(junction_assays) >= 2  # 5' and 3'

    def test_base_edit_checks_bystander(self):
        plan = validation_planner("base_edit", "TP53")
        purposes = " ".join(c["purpose"] for c in plan["validation_checklist"])
        assert "bystander" in purposes.lower()

    def test_prime_edit_checks_scaffold(self):
        plan = validation_planner("prime_edit", "HBB")
        assays = " ".join(c["assay"] for c in plan["validation_checklist"])
        assert "scaffold" in assays.lower()

    def test_crispra_has_no_genomic_assays(self):
        plan = validation_planner("CRISPRa", "TP53")
        categories = {c["category"] for c in plan["validation_checklist"]}
        assert "genomic" not in categories

    def test_crispri_has_no_genomic_assays(self):
        plan = validation_planner("CRISPRi", "MYC")
        categories = {c["category"] for c in plan["validation_checklist"]}
        assert "genomic" not in categories

    def test_crispra_recommends_rtqpcr(self):
        plan = validation_planner("CRISPRa", "TP53")
        assert plan["recommended_sequencing"] == "RT-qPCR"

    def test_large_deletion_has_spanning_pcr(self):
        plan = validation_planner("large_deletion", "BCL11A")
        assays = [c["assay"].lower() for c in plan["validation_checklist"]]
        assert any("deletion-spanning" in a or "spanning" in a for a in assays)

    def test_all_edit_types_produce_plans(self):
        for et in ("knockout", "knock-in", "base_edit", "prime_edit",
                    "CRISPRa", "CRISPRi", "large_deletion"):
            plan = validation_planner(et, "TEST_GENE")
            assert len(plan["validation_checklist"]) >= 3


class TestValidationPlannerDeliveryNormalization:
    """Freeform delivery strings are normalized to canonical categories."""

    @pytest.mark.parametrize("raw,expected", [
        ("RNP electroporation", "RNP"),
        ("rnp", "RNP"),
        ("ribonucleoprotein", "RNP"),
        ("lipofection", "lipofection"),
        ("lipofectamine", "lipofection"),
        ("lentiviral transduction", "lentiviral"),
        ("lenti", "lentiviral"),
        ("AAV transduction", "AAV"),
        ("adeno-associated virus", "AAV"),
        ("plasmid transfection", "plasmid"),
        ("plasmid", "plasmid"),
    ])
    def test_alias_normalization(self, raw, expected):
        plan = validation_planner("knockout", "BRCA1", raw)
        assert plan["delivery_method"] == expected

    def test_unknown_delivery_falls_back_to_unspecified(self):
        plan = validation_planner("knockout", "BRCA1", "carrier_pigeon")
        assert plan["delivery_method"] == "unspecified"

    def test_unknown_delivery_adds_warning_note(self):
        plan = validation_planner("knockout", "BRCA1", "carrier_pigeon")
        assert any("not recognized" in n for n in plan["notes"])

    def test_lentiviral_appends_integration_controls(self):
        plan = validation_planner("knockout", "BRCA1", "lentiviral")
        delivery_ctrls = [
            c for c in plan["controls"] if c["control_type"] == "delivery"
        ]
        assert len(delivery_ctrls) >= 1
        descriptions = " ".join(c["description"] for c in delivery_ctrls)
        assert "integration" in descriptions.lower()

    def test_rnp_appends_mock_electroporation_control(self):
        plan = validation_planner("knockout", "BRCA1", "RNP")
        delivery_ctrls = [
            c for c in plan["controls"] if c["control_type"] == "delivery"
        ]
        assert len(delivery_ctrls) == 1
        assert "electroporation" in delivery_ctrls[0]["description"].lower()

    def test_unspecified_appends_no_delivery_controls(self):
        plan = validation_planner("knockout", "BRCA1")
        delivery_ctrls = [
            c for c in plan["controls"] if c["control_type"] == "delivery"
        ]
        assert len(delivery_ctrls) == 0


class TestValidationPlannerLocusSubstitution:
    """The {locus} placeholder is replaced throughout the output."""

    def test_locus_in_notes(self):
        plan = validation_planner("knockout", "BRCA1")
        for note in plan["notes"]:
            assert "{locus}" not in note
        assert any("BRCA1" in n for n in plan["notes"])

    def test_locus_in_checklist(self):
        plan = validation_planner("knockout", "XYZGENE")
        combined = " ".join(
            c["purpose"] + c["expected_result"]
            for c in plan["validation_checklist"]
        )
        assert "{locus}" not in combined
        assert "XYZGENE" in combined

    def test_locus_whitespace_stripped(self):
        plan = validation_planner("knockout", "  BRCA1  ")
        assert plan["locus"] == "BRCA1"


class TestValidationPlannerEdgeCases:
    """Invalid inputs and edge cases."""

    def test_invalid_edit_type_raises(self):
        with pytest.raises(ValueError, match="Unsupported edit_type"):
            validation_planner("gene_drive", "BRCA1")

    def test_empty_locus_raises(self):
        with pytest.raises(ValueError, match="non-empty"):
            validation_planner("knockout", "")

    def test_whitespace_only_locus_raises(self):
        with pytest.raises(ValueError, match="non-empty"):
            validation_planner("knockout", "   ")

    def test_edit_type_case_insensitive(self):
        plan = validation_planner("KNOCKOUT", "BRCA1")
        assert plan["edit_type"] == "knockout"

    def test_edit_type_crispra_case_insensitive(self):
        plan = validation_planner("crispra", "TP53")
        assert plan["edit_type"] == "CRISPRa"

    def test_estimated_timeline_is_positive_int(self):
        plan = validation_planner("knockout", "BRCA1")
        assert isinstance(plan["estimated_timeline_days"], int)
        assert plan["estimated_timeline_days"] > 0
