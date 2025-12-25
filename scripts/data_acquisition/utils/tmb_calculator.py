#!/usr/bin/env python3
"""
TMB Calculation Utilities
==========================
Modular utilities for calculating Tumor Mutational Burden (TMB).

Author: Agent
Date: January 28, 2025
"""

from typing import List, Dict, Optional

# Standard TCGA WES exome size (Mb)
TCGA_EXOME_SIZE_MB = 38.0

# Variant classifications that count as nonsynonymous
NONSYNONYMOUS_CLASSIFICATIONS = {
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Splice_Site",
    "Translation_Start_Site",
    "Nonstop_Mutation",
    "De_novo_Start_OutOfFrame",
    "De_novo_Start_InFrame"
}


def count_nonsynonymous_mutations(mutations: List[Dict]) -> int:
    """
    Count nonsynonymous mutations from mutation list.
    
    Args:
        mutations: List of mutation dicts with 'variant_classification' field
    
    Returns:
        Count of nonsynonymous mutations
    """
    if not mutations:
        return 0
    
    count = 0
    for mut in mutations:
        classification = mut.get("variant_classification", "").strip()
        if classification in NONSYNONYMOUS_CLASSIFICATIONS:
            count += 1
    
    return count


def calculate_tmb_from_mutations(mutations: List[Dict], exome_size_mb: float = TCGA_EXOME_SIZE_MB) -> float:
    """
    Calculate TMB from mutation list.
    
    Args:
        mutations: List of mutation dicts
        exome_size_mb: Exome size in megabases (default: 38 Mb for TCGA WES)
    
    Returns:
        TMB in mutations per megabase
    """
    nonsyn_count = count_nonsynonymous_mutations(mutations)
    if exome_size_mb <= 0:
        return 0.0
    return nonsyn_count / exome_size_mb


def extract_tmb_from_clinical(clinical_data: Dict[str, str]) -> Optional[float]:
    """
    Extract TMB value from clinical data dict.
    Looks for common TMB field names.
    
    Args:
        clinical_data: Dict of clinical attributes
    
    Returns:
        TMB value as float, or None if not found
    """
    tmb_field_names = [
        "TMB_SCORE",
        "TMB",
        "TMB_NONSYNONYMOUS",
        "MUTATION_COUNT",
        "TMB_(NONSYNONYMOUS)",
        "TMB_NONSYNONYMOUS_MUTATIONS"
    ]
    
    for field_name in tmb_field_names:
        if field_name in clinical_data:
            try:
                value = clinical_data[field_name]
                if value and str(value).strip().lower() not in ["", "na", "n/a", "null", "none"]:
                    return float(value)
            except (ValueError, TypeError):
                continue
    
    return None


def normalize_tmb_value(tmb: Optional[float]) -> Optional[float]:
    """
    Validate and normalize TMB value.
    
    Args:
        tmb: TMB value to validate
    
    Returns:
        Normalized TMB value, or None if invalid
    """
    if tmb is None:
        return None
    
    # Validate range (0-500 mut/Mb is reasonable)
    if tmb < 0 or tmb > 500:
        return None
    
    return round(tmb, 2)


def get_tmb_classification(tmb: Optional[float], threshold: float = 10.0) -> Optional[bool]:
    """
    Classify TMB as high (True) or low (False).
    
    Args:
        tmb: TMB value
        threshold: Threshold for TMB-H (default: 10 mut/Mb)
    
    Returns:
        True if TMB-H, False if TMB-L, None if TMB unavailable
    """
    if tmb is None:
        return None
    return tmb >= threshold

