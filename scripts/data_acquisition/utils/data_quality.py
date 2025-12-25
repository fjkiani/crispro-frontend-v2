#!/usr/bin/env python3
"""
Data Quality Validation Utilities
==================================
Modular utilities for validating and flagging data quality.

Author: Agent
Date: January 28, 2025
"""

from typing import List, Dict, Optional, Set

# MSI status normalization
MSI_HIGH_VALUES = {"msi-h", "msi high", "high", "1", "yes", "true"}
MSI_LOW_VALUES = {"msi-l", "msi low", "low"}
MSS_VALUES = {"mss", "microsatellite stable", "stable", "0", "no", "false"}


def normalize_msi_status(msi_value: Optional[str]) -> Optional[str]:
    """
    Normalize MSI status to standard values.
    
    Args:
        msi_value: Raw MSI status value from clinical data
    
    Returns:
        Normalized value: "MSI-H", "MSI-L", "MSS", or "Unknown"
    """
    if not msi_value:
        return "Unknown"
    
    msi_lower = str(msi_value).strip().lower()
    
    if msi_lower in MSI_HIGH_VALUES:
        return "MSI-H"
    elif msi_lower in MSI_LOW_VALUES:
        return "MSI-L"
    elif msi_lower in MSS_VALUES:
        return "MSS"
    else:
        return "Unknown"


def extract_msi_from_clinical(clinical_data: Dict[str, str]) -> Optional[str]:
    """
    Extract MSI status from clinical data dict.
    Looks for common MSI field names.
    
    Args:
        clinical_data: Dict of clinical attributes
    
    Returns:
        Normalized MSI status or None
    """
    msi_field_names = [
        "MSI_STATUS",
        "MSI",
        "MICROSATELLITE_INSTABILITY",
        "MSI_SCORE",
        "MSI_STATUS_SCORE"
    ]
    
    for field_name in msi_field_names:
        if field_name in clinical_data:
            value = clinical_data[field_name]
            if value and str(value).strip().lower() not in ["", "na", "n/a", "null", "none"]:
                return normalize_msi_status(value)
    
    return None


def generate_data_quality_flags(
    has_tmb: bool,
    has_msi: bool,
    has_mutations: bool,
    msi_status: Optional[str],
    tmb_value: Optional[float],
    has_mbd4: bool,
    has_tp53: bool,
    has_ddr: bool
) -> List[str]:
    """
    Generate data quality flags for a patient.
    
    Args:
        has_tmb: Whether TMB value is present
        has_msi: Whether MSI status is present
        has_mutations: Whether mutations were extracted
        msi_status: Normalized MSI status
        tmb_value: TMB value
        has_mbd4: Has MBD4 mutation
        has_tp53: Has TP53 mutation
        has_ddr: Has DDR pathway mutation
    
    Returns:
        List of data quality flag strings
    """
    flags = []
    
    if has_tmb:
        flags.append("has_tmb")
    if has_msi:
        flags.append("has_msi")
    if has_mutations:
        flags.append("has_mutations")
    
    if not has_tmb and not has_msi:
        flags.append("biomarkers_incomplete")
    
    if msi_status == "MSI-H":
        flags.append("msi_high")
    
    if tmb_value is not None and tmb_value >= 10.0:
        flags.append("tmb_high")
    
    if has_mbd4:
        flags.append("mbd4_mutant")
    if has_tp53:
        flags.append("tp53_mutant")
    if has_ddr:
        flags.append("ddr_mutant")
    
    return flags


def deduplicate_patients(patients: List[Dict], key: str = "patient_id") -> List[Dict]:
    """
    Deduplicate patients by patient ID, keeping first occurrence.
    
    Args:
        patients: List of patient dicts
        key: Key to use for deduplication (default: "patient_id")
    
    Returns:
        Deduplicated list of patients
    """
    seen = set()
    deduplicated = []
    
    for patient in patients:
        patient_id = patient.get(key)
        if patient_id and patient_id not in seen:
            seen.add(patient_id)
            deduplicated.append(patient)
    
    return deduplicated


def validate_patient_record(patient: Dict) -> Dict:
    """
    Validate and clean a patient record.
    
    Args:
        patient: Patient record dict
    
    Returns:
        Validated patient record
    """
    # Validate TMB
    tmb = patient.get("tmb_score")
    if tmb is not None:
        if tmb < 0 or tmb > 500:
            patient["tmb_score"] = None
            patient["tmb_calculated"] = None
    
    # Normalize MSI
    msi = patient.get("msi_status")
    if msi:
        patient["msi_status"] = normalize_msi_status(msi)
    
    return patient

