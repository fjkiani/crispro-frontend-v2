"""
Drug Class Mapping for Treatment Line Integration

Maps individual drug/regimen names to drug classes for cross-resistance and sequencing logic.
"""

from typing import Optional

# Drug class mapping (P0: Ovarian + Breast HER2+ only)
DRUG_CLASS_MAP = {
    # ===== OVARIAN CANCER =====
    "carboplatin": "platinum_agent",
    "cisplatin": "platinum_agent",
    "carboplatin+paclitaxel": "platinum_agent",  # combo still considered platinum class
    "cisplatin+paclitaxel": "platinum_agent",
    
    "olaparib": "PARP_inhibitor",
    "niraparib": "PARP_inhibitor",
    "rucaparib": "PARP_inhibitor",
    "talazoparib": "PARP_inhibitor",
    
    "bevacizumab": "bevacizumab_combo",
    "bevacizumab+carboplatin": "bevacizumab_combo",
    "bevacizumab+paclitaxel": "bevacizumab_combo",
    
    "topotecan": "topotecan",
    "gemcitabine": "gemcitabine",
    "pegylated liposomal doxorubicin": "anthracycline",
    "doxorubicin": "anthracycline",
    
    # ===== BREAST HER2+ =====
    "trastuzumab+pertuzumab+paclitaxel": "TP_taxane_combo",
    "trastuzumab+pertuzumab+docetaxel": "TP_taxane_combo",
    "trastuzumab+paclitaxel": "TP_taxane_combo",
    "trastuzumab+docetaxel": "TP_taxane_combo",
    "pertuzumab+trastuzumab+taxane": "TP_taxane_combo",
    
    "trastuzumab deruxtecan": "T-DXd",
    "tdxd": "T-DXd",
    "t-dxd": "T-DXd",
    "enhertu": "T-DXd",
    
    "tucatinib+trastuzumab+capecitabine": "tucatinib_combo",
    "tucatinib": "tucatinib_combo",
    
    "neratinib+capecitabine": "neratinib_combo",
    "neratinib": "neratinib_combo",
    
    "lapatinib+capecitabine": "lapatinib_combo",
    "lapatinib": "lapatinib_combo",
    
    # Generic HER2 agents
    "trastuzumab": "trastuzumab",
    "pertuzumab": "pertuzumab",
}


def get_drug_class(drug_name: str) -> Optional[str]:
    """
    Get the drug class for a given drug/regimen name.
    
    Args:
        drug_name: Drug or regimen name (case-insensitive)
        
    Returns:
        Drug class string or None if not found
        
    Examples:
        >>> get_drug_class("carboplatin")
        "platinum_agent"
        >>> get_drug_class("Olaparib")
        "PARP_inhibitor"
        >>> get_drug_class("unknown_drug")
        None
    """
    normalized = drug_name.lower().strip()
    return DRUG_CLASS_MAP.get(normalized)


def get_all_drugs_in_class(drug_class: str) -> list[str]:
    """
    Get all drugs in a given drug class.
    
    Args:
        drug_class: Drug class name (e.g., "platinum_agent")
        
    Returns:
        List of drug names in that class
        
    Examples:
        >>> get_all_drugs_in_class("platinum_agent")
        ["carboplatin", "cisplatin", ...]
    """
    return [
        drug for drug, cls in DRUG_CLASS_MAP.items()
        if cls == drug_class
    ]









