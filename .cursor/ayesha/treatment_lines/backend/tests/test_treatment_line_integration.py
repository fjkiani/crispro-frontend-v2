"""
Unit Tests for Treatment Line Integration

Tests the integration service that computes treatment line features for SAE.
"""

import pytest
import sys
import os

# Add treatment_lines to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

from backend.services.treatment_line_integration import (
    compute_treatment_line_features,
    modulate_confidence_with_treatment_line
)
from backend.schemas.treatment_history import TreatmentHistory


def test_compute_features_no_history():
    """Test feature computation with no treatment history"""
    features = compute_treatment_line_features(
        drug_name="olaparib",
        disease="ovarian_cancer",
        treatment_history=None
    )
    
    assert features["line_appropriateness"] == 0.6  # Default
    assert features["cross_resistance_risk"] == 0.0
    assert features["current_line"] == 1
    assert features["prior_therapies"] == []


def test_compute_features_ovarian_l2_case():
    """Test Ovarian L2 case: post-platinum → olaparib"""
    history = TreatmentHistory(
        current_line=2,
        prior_therapies=["carboplatin", "paclitaxel"]
    )
    
    features = compute_treatment_line_features(
        drug_name="olaparib",
        disease="ovarian_cancer",
        treatment_history=history
    )
    
    # Should have perfect line appropriateness (NCCN Cat 1 for L2)
    assert features["line_appropriateness"] == 1.0
    assert features["nccn_category"] == "1"
    
    # Should have 0.4 cross-resistance risk with platinum
    assert features["cross_resistance_risk"] == 0.4
    assert "DNA repair" in features["cross_resistance_rationale"]
    
    # Sequencing fitness = 1.0 × (1 - 0.4) = 0.6
    assert abs(features["sequencing_fitness"] - 0.6) < 0.01


def test_compute_features_breast_l3_tdxd():
    """Test breast HER2+ L3 post-T-DXd → tucatinib"""
    history = TreatmentHistory(
        current_line=3,
        prior_therapies=["trastuzumab+pertuzumab+paclitaxel", "trastuzumab deruxtecan"]
    )
    
    features = compute_treatment_line_features(
        drug_name="tucatinib+trastuzumab+capecitabine",
        disease="breast_her2_positive",
        treatment_history=history
    )
    
    # Should have perfect line appropriateness for L3
    assert features["line_appropriateness"] == 1.0
    assert features["nccn_category"] == "1"
    
    # Should have 0.2 cross-resistance with T-DXd
    assert features["cross_resistance_risk"] == 0.2
    
    # Sequencing fitness = 1.0 × (1 - 0.2) = 0.8
    assert abs(features["sequencing_fitness"] - 0.8) < 0.01


def test_compute_features_first_line():
    """Test first-line therapy (no prior therapies)"""
    history = TreatmentHistory(
        current_line=1,
        prior_therapies=[]
    )
    
    features = compute_treatment_line_features(
        drug_name="carboplatin+paclitaxel",
        disease="ovarian_cancer",
        treatment_history=history
    )
    
    # Perfect appropriateness for first-line platinum
    assert features["line_appropriateness"] == 1.0
    assert features["nccn_category"] == "1"
    
    # No cross-resistance (no prior therapies)
    assert features["cross_resistance_risk"] == 0.0
    
    # Sequencing fitness = 1.0 × (1 - 0.0) = 1.0
    assert features["sequencing_fitness"] == 1.0


def test_modulate_confidence_no_penalty():
    """Test confidence modulation with no cross-resistance"""
    features = {
        "cross_resistance_risk": 0.0,
        "line_appropriateness": 1.0
    }
    
    modulated, rationale = modulate_confidence_with_treatment_line(
        base_confidence=0.8,
        treatment_line_features=features
    )
    
    assert modulated == 0.8  # No change
    assert "No treatment line adjustments" in rationale


def test_modulate_confidence_with_penalty():
    """Test confidence modulation with cross-resistance penalty"""
    features = {
        "cross_resistance_risk": 0.4,  # 40% risk
        "line_appropriateness": 1.0
    }
    
    modulated, rationale = modulate_confidence_with_treatment_line(
        base_confidence=0.8,
        treatment_line_features=features
    )
    
    # Penalty = 0.4 × 0.2 = 0.08 (8%)
    # New confidence = 0.8 - 0.08 = 0.72
    assert abs(modulated - 0.72) < 0.01
    assert "Reduced by 8.0%" in rationale


def test_modulate_confidence_poor_line_fit():
    """Test confidence modulation with poor line appropriateness"""
    features = {
        "cross_resistance_risk": 0.0,
        "line_appropriateness": 0.6  # Poor fit
    }
    
    modulated, rationale = modulate_confidence_with_treatment_line(
        base_confidence=0.8,
        treatment_line_features=features
    )
    
    # No penalty from cross-resistance, but rationale mentions poor fit
    assert modulated == 0.8
    assert "Limited line appropriateness" in rationale


def test_modulate_confidence_combined_penalty():
    """Test confidence modulation with both cross-resistance and poor fit"""
    features = {
        "cross_resistance_risk": 0.5,  # 50% risk
        "line_appropriateness": 0.7   # Moderate fit
    }
    
    modulated, rationale = modulate_confidence_with_treatment_line(
        base_confidence=0.9,
        treatment_line_features=features
    )
    
    # Penalty = 0.5 × 0.2 = 0.1 (10%)
    # New confidence = 0.9 - 0.1 = 0.8
    assert abs(modulated - 0.8) < 0.01
    assert "Reduced by 10.0%" in rationale
    assert "Limited line appropriateness" in rationale


def test_confidence_floor():
    """Test that confidence never goes below 0.0"""
    features = {
        "cross_resistance_risk": 1.0,  # Maximum risk
        "line_appropriateness": 0.5
    }
    
    modulated, _ = modulate_confidence_with_treatment_line(
        base_confidence=0.1,
        treatment_line_features=features
    )
    
    # Should not go negative
    assert modulated >= 0.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])









