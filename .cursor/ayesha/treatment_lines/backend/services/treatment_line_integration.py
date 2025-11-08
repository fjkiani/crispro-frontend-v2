"""
Treatment Line Integration Service

Computes treatment line context for SAE features and confidence modulation.
Integrates panel_config.py and cross_resistance_map.py.
"""

import sys
import os
from typing import Dict, Any, Optional, List

# Add treatment_lines to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

from backend.services.panel_config import calculate_line_appropriateness, get_line_metadata
from backend.services.cross_resistance_map import calculate_aggregate_cross_resistance
from backend.schemas.treatment_history import TreatmentHistory


def compute_treatment_line_features(
    drug_name: str,
    disease: str,
    treatment_history: Optional[TreatmentHistory] = None
) -> Dict[str, Any]:
    """
    Compute treatment line features for a specific drug.
    
    Args:
        drug_name: Drug/regimen name
        disease: Disease name (e.g., "ovarian_cancer", "breast_her2_positive")
        treatment_history: Patient's treatment history (optional)
    
    Returns:
        Dictionary with treatment line features:
        - line_appropriateness: float (0.0-1.0)
        - line_rationale: str
        - nccn_category: str
        - cross_resistance_risk: float (0.0-1.0)
        - cross_resistance_rationale: str
        - sequencing_fitness: float (0.0-1.0)
        - current_line: int
        - prior_therapies: List[str]
    """
    # Default values (no treatment history)
    if not treatment_history:
        return {
            "line_appropriateness": 0.6,
            "line_rationale": "No treatment history provided - using default appropriateness",
            "nccn_category": "unknown",
            "cross_resistance_risk": 0.0,
            "cross_resistance_rationale": "No prior therapies to assess cross-resistance",
            "sequencing_fitness": 0.6,
            "current_line": 1,  # Assume first-line
            "prior_therapies": []
        }
    
    current_line = treatment_history.current_line
    prior_therapies = treatment_history.prior_therapies
    
    # Compute line appropriateness
    line_appropriateness, line_rationale = calculate_line_appropriateness(
        drug_name=drug_name,
        disease=disease,
        current_line=current_line
    )
    
    # Get NCCN category
    line_meta = get_line_metadata(drug_name, disease, current_line)
    nccn_category = line_meta.nccn_category if line_meta else "unknown"
    
    # Compute cross-resistance risk
    cross_resistance_risk = 0.0
    cross_resistance_rationale = "No prior therapies"
    
    if prior_therapies:
        cross_resistance_risk, rationales = calculate_aggregate_cross_resistance(
            prior_therapies=prior_therapies,
            candidate_drug=drug_name
        )
        if rationales:
            cross_resistance_rationale = "; ".join(rationales)
        else:
            cross_resistance_rationale = "No known cross-resistance with prior therapies"
    
    # Compute sequencing fitness (combines appropriateness and cross-resistance)
    # Formula: fitness = line_appropriateness × (1 - cross_resistance_risk)
    sequencing_fitness = line_appropriateness * (1.0 - cross_resistance_risk)
    
    return {
        "line_appropriateness": line_appropriateness,
        "line_rationale": line_rationale,
        "nccn_category": nccn_category,
        "cross_resistance_risk": cross_resistance_risk,
        "cross_resistance_rationale": cross_resistance_rationale,
        "sequencing_fitness": sequencing_fitness,
        "current_line": current_line,
        "prior_therapies": prior_therapies
    }


def modulate_confidence_with_treatment_line(
    base_confidence: float,
    treatment_line_features: Dict[str, Any]
) -> tuple[float, str]:
    """
    Modulate confidence score based on treatment line context.
    
    Applies linear penalty for cross-resistance:
    confidence -= cross_resistance_risk × 0.2 (max -20%)
    
    Args:
        base_confidence: Base confidence score (0.0-1.0)
        treatment_line_features: Treatment line features from compute_treatment_line_features
    
    Returns:
        Tuple of (modulated_confidence, rationale)
    """
    cross_resistance_risk = treatment_line_features.get("cross_resistance_risk", 0.0)
    line_appropriateness = treatment_line_features.get("line_appropriateness", 1.0)
    
    # Apply cross-resistance penalty (max -20%)
    penalty = cross_resistance_risk * 0.2
    modulated_confidence = max(0.0, min(1.0, base_confidence - penalty))
    
    # Build rationale
    rationale_parts = []
    
    if penalty > 0.05:  # Significant penalty (>5%)
        rationale_parts.append(
            f"Reduced by {penalty * 100:.1f}% due to cross-resistance risk"
        )
    
    if line_appropriateness < 0.8:  # Not ideal for this line
        rationale_parts.append(
            f"Limited line appropriateness ({line_appropriateness:.2f}) affects confidence"
        )
    
    if not rationale_parts:
        rationale = "No treatment line adjustments applied"
    else:
        rationale = "; ".join(rationale_parts)
    
    return modulated_confidence, rationale

