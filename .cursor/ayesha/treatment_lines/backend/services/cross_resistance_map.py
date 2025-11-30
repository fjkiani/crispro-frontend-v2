"""
Cross-Resistance Mapping for Treatment Line Integration

Defines known cross-resistance relationships between drug classes based on MoA overlap.
"""

from typing import Dict, List, Optional
from .drug_class_map import get_drug_class


# Cross-resistance map (P0: Ovarian + Breast HER2+ only)
CROSS_RESISTANCE_MAP: Dict[str, Dict] = {
    # ===== OVARIAN CANCER =====
    "PARP_inhibitor": {
        "cross_resistant_with": ["platinum_agent"],
        "risk_level": 0.4,
        "rationale": "DNA repair pathway overlap - both target DNA damage response"
    },
    "platinum_agent": {
        "cross_resistant_with": ["PARP_inhibitor"],
        "risk_level": 0.4,
        "rationale": "DNA repair pathway overlap - platinum resistance often predicts PARP resistance"
    },
    
    # ===== BREAST HER2+ =====
    "T-DXd": {
        "cross_resistant_with": ["trastuzumab", "pertuzumab"],
        "risk_level": 0.3,
        "rationale": "HER2-targeted mechanism overlap - prior HER2 blockade may affect T-DXd efficacy"
    },
    "tucatinib_combo": {
        "cross_resistant_with": ["trastuzumab", "T-DXd", "lapatinib_combo"],
        "risk_level": 0.2,
        "rationale": "HER2 TKI resistance after prior HER2 blockade - cross-resistance possible but lower"
    },
    "neratinib_combo": {
        "cross_resistant_with": ["trastuzumab", "T-DXd", "lapatinib_combo"],
        "risk_level": 0.2,
        "rationale": "HER2 TKI resistance after prior HER2 blockade"
    },
    "lapatinib_combo": {
        "cross_resistant_with": ["trastuzumab", "neratinib_combo", "tucatinib_combo"],
        "risk_level": 0.25,
        "rationale": "HER2 TKI cross-resistance between lapatinib, neratinib, tucatinib"
    },
}


def get_cross_resistance_risk(
    prior_drug: str,
    candidate_drug: str
) -> tuple[float, Optional[str]]:
    """
    Calculate cross-resistance risk between a prior drug and a candidate drug.
    
    Args:
        prior_drug: Name of previously used drug
        candidate_drug: Name of candidate drug being evaluated
        
    Returns:
        Tuple of (risk_level, rationale)
        - risk_level: 0.0-1.0 indicating cross-resistance risk
        - rationale: Human-readable explanation or None
        
    Examples:
        >>> get_cross_resistance_risk("carboplatin", "olaparib")
        (0.4, "DNA repair pathway overlap - both target DNA damage response")
        
        >>> get_cross_resistance_risk("carboplatin", "bevacizumab")
        (0.0, None)
    """
    # Get drug classes
    prior_class = get_drug_class(prior_drug)
    candidate_class = get_drug_class(candidate_drug)
    
    if not prior_class or not candidate_class:
        # Drug not in our map - no known cross-resistance (P0)
        return 0.0, None
    
    # Check if candidate class has known cross-resistance with prior class
    cross_res_info = CROSS_RESISTANCE_MAP.get(candidate_class, {})
    cross_resistant_with = cross_res_info.get("cross_resistant_with", [])
    
    if prior_class in cross_resistant_with:
        risk_level = cross_res_info.get("risk_level", 0.4)
        rationale = cross_res_info.get("rationale", "Known cross-resistance relationship")
        return risk_level, rationale
    
    # No known cross-resistance
    return 0.0, None


def get_all_cross_resistant_classes(drug_class: str) -> List[str]:
    """
    Get all drug classes that have known cross-resistance with the given class.
    
    Args:
        drug_class: Drug class name
        
    Returns:
        List of cross-resistant drug classes
        
    Examples:
        >>> get_all_cross_resistant_classes("PARP_inhibitor")
        ["platinum_agent"]
    """
    cross_res_info = CROSS_RESISTANCE_MAP.get(drug_class, {})
    return cross_res_info.get("cross_resistant_with", [])


def calculate_aggregate_cross_resistance(
    prior_therapies: List[str],
    candidate_drug: str
) -> tuple[float, List[str]]:
    """
    Calculate aggregate cross-resistance risk across multiple prior therapies.
    
    Args:
        prior_therapies: List of prior drug names
        candidate_drug: Candidate drug being evaluated
        
    Returns:
        Tuple of (max_risk, rationales)
        - max_risk: Maximum cross-resistance risk encountered
        - rationales: List of rationale strings for each cross-resistance relationship
        
    Examples:
        >>> calculate_aggregate_cross_resistance(
        ...     ["carboplatin", "paclitaxel"],
        ...     "olaparib"
        ... )
        (0.4, ["DNA repair pathway overlap - both target DNA damage response"])
    """
    max_risk = 0.0
    rationales = []
    
    for prior_drug in prior_therapies:
        risk, rationale = get_cross_resistance_risk(prior_drug, candidate_drug)
        if risk > 0.0 and rationale:
            max_risk = max(max_risk, risk)
            if rationale not in rationales:  # Deduplicate
                rationales.append(rationale)
    
    return max_risk, rationales









