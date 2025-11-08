"""
Unit Tests for Drug Panel Configuration

Tests panel loading, line filtering, and appropriateness calculations.
"""

import pytest
from ..services.panel_config import (
    get_panel_for_disease,
    get_drugs_for_line,
    get_line_metadata,
    calculate_line_appropriateness,
    OVARIAN_PANEL,
    BREAST_HER2_PANEL
)


def test_ovarian_panel_loaded():
    """Test that ovarian panel is loaded with drugs"""
    panel = get_panel_for_disease("ovarian_cancer")
    
    assert len(panel) > 0
    assert any(entry.drug_name == "carboplatin+paclitaxel" for entry in panel)
    assert any(entry.drug_name == "olaparib" for entry in panel)


def test_breast_her2_panel_loaded():
    """Test that breast HER2+ panel is loaded with drugs"""
    panel = get_panel_for_disease("breast_her2_positive")
    
    assert len(panel) > 0
    assert any(entry.drug_name == "trastuzumab+pertuzumab+paclitaxel" for entry in panel)
    assert any(entry.drug_name == "trastuzumab deruxtecan" for entry in panel)


def test_get_drugs_for_line_ovarian_l1():
    """Test getting first-line ovarian drugs"""
    drugs = get_drugs_for_line("ovarian_cancer", line=1)
    
    drug_names = [d.drug_name for d in drugs]
    assert "carboplatin+paclitaxel" in drug_names
    assert "olaparib" in drug_names  # First-line maintenance


def test_get_drugs_for_line_ovarian_l2():
    """Test getting second-line ovarian drugs"""
    drugs = get_drugs_for_line("ovarian_cancer", line=2)
    
    drug_names = [d.drug_name for d in drugs]
    assert "olaparib" in drug_names
    assert "niraparib" in drug_names
    assert "pegylated liposomal doxorubicin" in drug_names


def test_get_drugs_for_line_breast_l2():
    """Test getting second-line breast HER2+ drugs"""
    drugs = get_drugs_for_line("breast_her2_positive", line=2)
    
    drug_names = [d.drug_name for d in drugs]
    assert "trastuzumab deruxtecan" in drug_names
    assert "lapatinib+capecitabine" in drug_names


def test_get_line_metadata_platinum_l1():
    """Test getting line metadata for platinum first-line"""
    meta = get_line_metadata("carboplatin+paclitaxel", "ovarian_cancer", line=1)
    
    assert meta is not None
    assert meta.nccn_category == "1"
    assert meta.is_standard is True
    assert meta.sequencing_preference == 1


def test_get_line_metadata_olaparib_l2():
    """Test getting line metadata for olaparib second-line"""
    meta = get_line_metadata("olaparib", "ovarian_cancer", line=2)
    
    assert meta is not None
    assert meta.nccn_category == "1"
    assert meta.is_standard is True


def test_get_line_metadata_tdxd_l2():
    """Test getting line metadata for T-DXd second-line"""
    meta = get_line_metadata("trastuzumab deruxtecan", "breast_her2_positive", line=2)
    
    assert meta is not None
    assert meta.nccn_category == "1"
    assert meta.sequencing_preference == 1


def test_calculate_line_appropriateness_perfect():
    """Test appropriateness calculation for perfect match (NCCN Cat 1)"""
    score, rationale = calculate_line_appropriateness(
        "carboplatin+paclitaxel",
        "ovarian_cancer",
        current_line=1
    )
    
    assert score == 1.0
    assert "standard" in rationale.lower() or "first-line" in rationale.lower()


def test_calculate_line_appropriateness_nccn_2a():
    """Test appropriateness for NCCN Category 2A"""
    score, rationale = calculate_line_appropriateness(
        "gemcitabine",
        "ovarian_cancer",
        current_line=3
    )
    
    assert score == 0.9
    assert "2A" in rationale


def test_calculate_line_appropriateness_wrong_line():
    """Test appropriateness when drug not appropriate for line"""
    score, rationale = calculate_line_appropriateness(
        "topotecan",
        "ovarian_cancer",
        current_line=1  # Topotecan is for line 3+
    )
    
    assert score == 0.6  # Default for no guidance


def test_biomarker_requirements_olaparib():
    """Test that olaparib has biomarker requirements"""
    panel = get_panel_for_disease("ovarian_cancer")
    olaparib = next(e for e in panel if e.drug_name == "olaparib")
    
    assert olaparib.biomarker_requirements is not None
    assert "BRCA1_mutation" in olaparib.biomarker_requirements or "BRCA2_mutation" in olaparib.biomarker_requirements


def test_biomarker_requirements_tdxd():
    """Test that T-DXd requires HER2 amplification"""
    panel = get_panel_for_disease("breast_her2_positive")
    tdxd = next(e for e in panel if e.drug_name == "trastuzumab deruxtecan")
    
    assert tdxd.biomarker_requirements is not None
    assert "HER2_amplification" in tdxd.biomarker_requirements or "ERBB2_amplification" in tdxd.biomarker_requirements


def test_sequencing_preferences_ovarian_l1():
    """Test that first-line ovarian drugs have sequencing preferences"""
    drugs = get_drugs_for_line("ovarian_cancer", line=1)
    
    # Should have at least 2 drugs with sequencing preferences
    with_prefs = [d for d in drugs if any(lm.sequencing_preference for lm in d.line_metadata)]
    assert len(with_prefs) >= 2


def test_case_insensitive_drug_lookup():
    """Test that drug name lookup is case-insensitive"""
    meta1 = get_line_metadata("OLAPARIB", "ovarian_cancer", line=2)
    meta2 = get_line_metadata("olaparib", "ovarian_cancer", line=2)
    
    assert meta1 is not None
    assert meta2 is not None
    assert meta1.nccn_category == meta2.nccn_category


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


