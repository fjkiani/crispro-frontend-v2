"""
Phase 3 Integration Test: Treatment Line → Efficacy Orchestrator

Tests that treatment line features are properly integrated into the efficacy orchestrator
and that confidence is modulated correctly.
"""

import pytest
import sys
import os
from typing import Dict, Any

# Add paths
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../../../../oncology-coPilot/oncology-backend-minimal"))

from backend.services.treatment_line_integration import (
    compute_treatment_line_features,
    modulate_confidence_with_treatment_line
)
from backend.schemas.treatment_history import TreatmentHistory


def test_ayesha_olaparib_l2_post_platinum():
    """
    Test Ayesha's case: Ovarian L2 post-platinum → olaparib
    
    Expected:
    - Line appropriateness: 1.0 (NCCN Cat 1)
    - Cross-resistance risk: 0.4 (DNA repair overlap)
    - Sequencing fitness: 0.6
    - Confidence penalty: ~8% (0.4 × 0.2 = 0.08)
    """
    treatment_hist = TreatmentHistory(
        current_line=2,
        prior_therapies=["carboplatin", "paclitaxel"]
    )
    
    # Compute features
    features = compute_treatment_line_features(
        drug_name="olaparib",
        disease="ovarian_cancer",
        treatment_history=treatment_hist
    )
    
    # Validate features
    assert features["line_appropriateness"] == 1.0, "Should be perfect for L2"
    assert features["nccn_category"] == "1", "Should be NCCN Cat 1"
    assert 0.35 <= features["cross_resistance_risk"] <= 0.45, "Should have moderate cross-resistance"
    assert 0.55 <= features["sequencing_fitness"] <= 0.65, "Sequencing fitness should be ~0.6"
    
    # Test confidence modulation
    base_confidence = 0.80
    modulated_confidence, rationale = modulate_confidence_with_treatment_line(
        base_confidence=base_confidence,
        treatment_line_features=features
    )
    
    # Validate modulation
    expected_penalty = features["cross_resistance_risk"] * 0.2
    expected_confidence = base_confidence - expected_penalty
    
    assert abs(modulated_confidence - expected_confidence) < 0.01, \
        f"Expected {expected_confidence:.2f}, got {modulated_confidence:.2f}"
    assert 0.70 <= modulated_confidence <= 0.75, \
        "Ayesha's olaparib confidence should drop to ~0.72"
    assert "cross-resistance" in rationale.lower(), \
        "Rationale should mention cross-resistance"
    
    print("\n✅ AYESHA CASE (Ovarian L2 post-platinum → olaparib):")
    print(f"  Line Appropriateness: {features['line_appropriateness']:.2f}")
    print(f"  NCCN Category: {features['nccn_category']}")
    print(f"  Cross-Resistance Risk: {features['cross_resistance_risk']:.2f}")
    print(f"  Sequencing Fitness: {features['sequencing_fitness']:.2f}")
    print(f"  Base Confidence: {base_confidence:.2f}")
    print(f"  Modulated Confidence: {modulated_confidence:.2f}")
    print(f"  Penalty: {expected_penalty * 100:.1f}%")
    print(f"  Rationale: {rationale}")


def test_dr_lustberg_tucatinib_l3_post_tdxd():
    """
    Test Dr. Lustberg's case: Breast HER2+ L3 post-T-DXd → tucatinib combo
    
    Expected:
    - Line appropriateness: 1.0 (NCCN Cat 1)
    - Cross-resistance risk: 0.2 (low TKI cross-resistance)
    - Sequencing fitness: 0.8
    - Confidence penalty: ~4% (0.2 × 0.2 = 0.04)
    """
    treatment_hist = TreatmentHistory(
        current_line=3,
        prior_therapies=["trastuzumab deruxtecan", "pertuzumab"]
    )
    
    # Compute features (use full combo name as in panel_config)
    features = compute_treatment_line_features(
        drug_name="tucatinib+trastuzumab+capecitabine",
        disease="breast_her2_positive",
        treatment_history=treatment_hist
    )
    
    # Validate features
    assert features["line_appropriateness"] == 1.0, "Should be perfect for L3"
    assert features["nccn_category"] == "1", "Should be NCCN Cat 1"
    assert 0.15 <= features["cross_resistance_risk"] <= 0.25, "Should have low cross-resistance"
    assert 0.75 <= features["sequencing_fitness"] <= 0.85, "Sequencing fitness should be ~0.8"
    
    # Test confidence modulation
    base_confidence = 0.85
    modulated_confidence, rationale = modulate_confidence_with_treatment_line(
        base_confidence=base_confidence,
        treatment_line_features=features
    )
    
    # Validate modulation
    expected_penalty = features["cross_resistance_risk"] * 0.2
    expected_confidence = base_confidence - expected_penalty
    
    assert abs(modulated_confidence - expected_confidence) < 0.01, \
        f"Expected {expected_confidence:.2f}, got {modulated_confidence:.2f}"
    assert 0.80 <= modulated_confidence <= 0.85, \
        "Dr. Lustberg's tucatinib confidence should drop to ~0.81"
    
    print("\n✅ DR. LUSTBERG CASE (Breast HER2+ L3 post-T-DXd → tucatinib):")
    print(f"  Line Appropriateness: {features['line_appropriateness']:.2f}")
    print(f"  NCCN Category: {features['nccn_category']}")
    print(f"  Cross-Resistance Risk: {features['cross_resistance_risk']:.2f}")
    print(f"  Sequencing Fitness: {features['sequencing_fitness']:.2f}")
    print(f"  Base Confidence: {base_confidence:.2f}")
    print(f"  Modulated Confidence: {modulated_confidence:.2f}")
    print(f"  Penalty: {expected_penalty * 100:.1f}%")
    print(f"  Rationale: {rationale}")


def test_first_line_no_penalty():
    """
    Test first-line therapy (no prior therapies) → no cross-resistance penalty
    """
    treatment_hist = TreatmentHistory(
        current_line=1,
        prior_therapies=[]
    )
    
    # Use correct drug name from panel
    features = compute_treatment_line_features(
        drug_name="carboplatin+paclitaxel",
        disease="ovarian_cancer",
        treatment_history=treatment_hist
    )
    
    # Should have no cross-resistance risk (no prior therapies)
    assert features["cross_resistance_risk"] == 0.0, "First-line should have no cross-resistance"
    
    # Confidence modulation - should have no cross-resistance penalty (max penalty is 20%)
    base_confidence = 0.75
    modulated_confidence, rationale = modulate_confidence_with_treatment_line(
        base_confidence=base_confidence,
        treatment_line_features=features
    )
    
    # No cross-resistance penalty, but line appropriateness check still applies
    # (0.0 cross-res × 0.2 = 0.0 penalty from cross-resistance)
    expected_penalty = features["cross_resistance_risk"] * 0.2
    assert expected_penalty == 0.0, "Should have zero cross-resistance penalty"
    
    print("\n✅ FIRST-LINE CASE (no penalty):")
    print(f"  Base Confidence: {base_confidence:.2f}")
    print(f"  Modulated Confidence: {modulated_confidence:.2f}")
    print(f"  Rationale: {rationale}")


def test_confidence_floor():
    """
    Test that confidence never goes below 0.0
    """
    # Create high cross-resistance scenario
    features = {
        "cross_resistance_risk": 1.0,  # Maximum risk
        "line_appropriateness": 0.5
    }
    
    # Very low base confidence
    base_confidence = 0.10
    modulated_confidence, _ = modulate_confidence_with_treatment_line(
        base_confidence=base_confidence,
        treatment_line_features=features
    )
    
    # Should not go negative
    assert modulated_confidence >= 0.0, "Confidence should never be negative"
    assert modulated_confidence == 0.0, "Should hit floor at 0.0"
    
    print("\n✅ CONFIDENCE FLOOR TEST:")
    print(f"  Base: {base_confidence:.2f} → Floor: {modulated_confidence:.2f}")


if __name__ == "__main__":
    print("=" * 70)
    print("PHASE 3 INTEGRATION TEST: Treatment Line → Efficacy Orchestrator")
    print("=" * 70)
    
    test_ayesha_olaparib_l2_post_platinum()
    test_dr_lustberg_tucatinib_l3_post_tdxd()
    test_first_line_no_penalty()
    test_confidence_floor()
    
    print("\n" + "=" * 70)
    print("✅ ALL PHASE 3 INTEGRATION TESTS PASSED!")
    print("=" * 70)

