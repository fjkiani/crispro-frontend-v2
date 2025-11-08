"""
Unit Tests for Cross-Resistance Logic

Tests the cross-resistance mapping and calculation functions.
"""

import pytest
from ..services.cross_resistance_map import (
    get_cross_resistance_risk,
    calculate_aggregate_cross_resistance,
    get_all_cross_resistant_classes
)


def test_platinum_parp_cross_resistance():
    """Test known platinum-PARP cross-resistance"""
    risk, rationale = get_cross_resistance_risk("carboplatin", "olaparib")
    
    assert risk == 0.4
    assert "DNA repair" in rationale
    assert rationale is not None


def test_parp_platinum_cross_resistance():
    """Test reverse direction (PARP→platinum)"""
    risk, rationale = get_cross_resistance_risk("olaparib", "carboplatin")
    
    assert risk == 0.4
    assert "DNA repair" in rationale


def test_no_cross_resistance():
    """Test drugs with no known cross-resistance"""
    risk, rationale = get_cross_resistance_risk("carboplatin", "bevacizumab")
    
    assert risk == 0.0
    assert rationale is None


def test_her2_cross_resistance():
    """Test HER2 agent cross-resistance"""
    risk, rationale = get_cross_resistance_risk("trastuzumab", "trastuzumab deruxtecan")
    
    assert risk == 0.3
    assert "HER2" in rationale


def test_unknown_drug():
    """Test behavior with unknown drug"""
    risk, rationale = get_cross_resistance_risk("unknown_drug", "olaparib")
    
    assert risk == 0.0
    assert rationale is None


def test_aggregate_cross_resistance():
    """Test aggregate calculation across multiple prior therapies"""
    prior = ["carboplatin", "paclitaxel"]
    candidate = "olaparib"
    
    max_risk, rationales = calculate_aggregate_cross_resistance(prior, candidate)
    
    assert max_risk == 0.4  # Max from carboplatin
    assert len(rationales) == 1  # Only carboplatin has cross-resistance
    assert "DNA repair" in rationales[0]


def test_multiple_cross_resistances():
    """Test multiple cross-resistance relationships"""
    prior = ["trastuzumab", "pertuzumab"]
    candidate = "trastuzumab deruxtecan"
    
    max_risk, rationales = calculate_aggregate_cross_resistance(prior, candidate)
    
    assert max_risk == 0.3
    assert len(rationales) >= 1


def test_get_all_cross_resistant_classes():
    """Test getting all cross-resistant classes"""
    classes = get_all_cross_resistant_classes("PARP_inhibitor")
    
    assert "platinum_agent" in classes


def test_case_insensitive_drug_names():
    """Test that drug names are case-insensitive"""
    risk1, _ = get_cross_resistance_risk("CARBOPLATIN", "olaparib")
    risk2, _ = get_cross_resistance_risk("carboplatin", "OLAPARIB")
    risk3, _ = get_cross_resistance_risk("Carboplatin", "Olaparib")
    
    assert risk1 == risk2 == risk3 == 0.4


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

