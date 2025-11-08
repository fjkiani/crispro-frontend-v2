"""
Disease-Specific Drug Panel Configuration with Treatment Line Metadata

Defines drug panels for ovarian and breast HER2+ cancer with:
- Treatment line appropriateness (1st/2nd/3rd/4th line)
- NCCN category per line
- Standard regimens
- Sequencing preferences
"""

from typing import Dict, List, Optional
from dataclasses import dataclass


@dataclass
class DrugLineMetadata:
    """Metadata for a drug at a specific treatment line"""
    line: int
    nccn_category: str  # "1", "2A", "2B", "3"
    is_standard: bool
    sequencing_preference: Optional[int] = None  # 1 = first choice, 2 = second choice, etc.
    rationale: Optional[str] = None


@dataclass
class DrugPanelEntry:
    """Complete drug panel entry with all treatment line metadata"""
    drug_name: str
    drug_class: str
    disease: str
    line_metadata: List[DrugLineMetadata]
    contraindications: Optional[List[str]] = None
    biomarker_requirements: Optional[List[str]] = None


# ===== OVARIAN CANCER PANEL =====

OVARIAN_PANEL: List[DrugPanelEntry] = [
    # First-line platinum-based
    DrugPanelEntry(
        drug_name="carboplatin+paclitaxel",
        drug_class="platinum_agent",
        disease="ovarian_cancer",
        line_metadata=[
            DrugLineMetadata(
                line=1,
                nccn_category="1",
                is_standard=True,
                sequencing_preference=1,
                rationale="Standard first-line therapy for advanced ovarian cancer"
            )
        ],
        biomarker_requirements=None
    ),
    
    DrugPanelEntry(
        drug_name="carboplatin+paclitaxel+bevacizumab",
        drug_class="bevacizumab_combo",
        disease="ovarian_cancer",
        line_metadata=[
            DrugLineMetadata(
                line=1,
                nccn_category="1",
                is_standard=True,
                sequencing_preference=2,
                rationale="First-line with anti-angiogenic therapy for high-risk disease"
            )
        ],
        biomarker_requirements=None
    ),
    
    # PARP inhibitors - maintenance and treatment
    DrugPanelEntry(
        drug_name="olaparib",
        drug_class="PARP_inhibitor",
        disease="ovarian_cancer",
        line_metadata=[
            DrugLineMetadata(
                line=1,
                nccn_category="1",
                is_standard=True,
                sequencing_preference=1,
                rationale="First-line maintenance for BRCA-mutated ovarian cancer"
            ),
            DrugLineMetadata(
                line=2,
                nccn_category="1",
                is_standard=True,
                sequencing_preference=1,
                rationale="Platinum-sensitive recurrence with BRCA mutation"
            )
        ],
        biomarker_requirements=["BRCA1_mutation", "BRCA2_mutation", "HRD_positive"]
    ),
    
    DrugPanelEntry(
        drug_name="niraparib",
        drug_class="PARP_inhibitor",
        disease="ovarian_cancer",
        line_metadata=[
            DrugLineMetadata(
                line=1,
                nccn_category="1",
                is_standard=True,
                sequencing_preference=2,
                rationale="First-line maintenance regardless of biomarker status"
            ),
            DrugLineMetadata(
                line=2,
                nccn_category="1",
                is_standard=True,
                sequencing_preference=2,
                rationale="Platinum-sensitive recurrence"
            )
        ],
        biomarker_requirements=None  # Can be used without biomarker
    ),
    
    DrugPanelEntry(
        drug_name="rucaparib",
        drug_class="PARP_inhibitor",
        disease="ovarian_cancer",
        line_metadata=[
            DrugLineMetadata(
                line=1,
                nccn_category="1",
                is_standard=True,
                sequencing_preference=3,
                rationale="First-line maintenance for BRCA-mutated or HRD-positive"
            ),
            DrugLineMetadata(
                line=2,
                nccn_category="1",
                is_standard=True,
                sequencing_preference=3,
                rationale="Platinum-sensitive recurrence with BRCA/HRD"
            )
        ],
        biomarker_requirements=["BRCA1_mutation", "BRCA2_mutation", "HRD_positive"]
    ),
    
    # Later-line options
    DrugPanelEntry(
        drug_name="topotecan",
        drug_class="topotecan",
        disease="ovarian_cancer",
        line_metadata=[
            DrugLineMetadata(
                line=3,
                nccn_category="1",
                is_standard=True,
                sequencing_preference=1,
                rationale="Standard third-line option for platinum-resistant disease"
            ),
            DrugLineMetadata(
                line=4,
                nccn_category="2A",
                is_standard=True,
                sequencing_preference=1,
                rationale="Subsequent therapy option"
            )
        ],
        biomarker_requirements=None
    ),
    
    DrugPanelEntry(
        drug_name="gemcitabine",
        drug_class="gemcitabine",
        disease="ovarian_cancer",
        line_metadata=[
            DrugLineMetadata(
                line=3,
                nccn_category="2A",
                is_standard=True,
                sequencing_preference=2,
                rationale="Platinum-resistant disease option"
            )
        ],
        biomarker_requirements=None
    ),
    
    DrugPanelEntry(
        drug_name="pegylated liposomal doxorubicin",
        drug_class="anthracycline",
        disease="ovarian_cancer",
        line_metadata=[
            DrugLineMetadata(
                line=2,
                nccn_category="1",
                is_standard=True,
                sequencing_preference=2,
                rationale="Platinum-sensitive recurrence alternative"
            ),
            DrugLineMetadata(
                line=3,
                nccn_category="2A",
                is_standard=True,
                sequencing_preference=3,
                rationale="Later-line option"
            )
        ],
        biomarker_requirements=None
    ),
]


# ===== BREAST HER2+ PANEL =====

BREAST_HER2_PANEL: List[DrugPanelEntry] = [
    # First-line standard
    DrugPanelEntry(
        drug_name="trastuzumab+pertuzumab+paclitaxel",
        drug_class="TP_taxane_combo",
        disease="breast_her2_positive",
        line_metadata=[
            DrugLineMetadata(
                line=1,
                nccn_category="1",
                is_standard=True,
                sequencing_preference=1,
                rationale="Standard first-line dual HER2 blockade for metastatic HER2+ breast cancer"
            )
        ],
        biomarker_requirements=["HER2_amplification", "ERBB2_amplification"]
    ),
    
    DrugPanelEntry(
        drug_name="trastuzumab+pertuzumab+docetaxel",
        drug_class="TP_taxane_combo",
        disease="breast_her2_positive",
        line_metadata=[
            DrugLineMetadata(
                line=1,
                nccn_category="1",
                is_standard=True,
                sequencing_preference=2,
                rationale="Alternative first-line dual HER2 blockade with docetaxel"
            )
        ],
        biomarker_requirements=["HER2_amplification", "ERBB2_amplification"]
    ),
    
    # T-DXd - Second-line standard
    DrugPanelEntry(
        drug_name="trastuzumab deruxtecan",
        drug_class="T-DXd",
        disease="breast_her2_positive",
        line_metadata=[
            DrugLineMetadata(
                line=2,
                nccn_category="1",
                is_standard=True,
                sequencing_preference=1,
                rationale="Preferred second-line after trastuzumab-based therapy (DESTINY-Breast03)"
            ),
            DrugLineMetadata(
                line=3,
                nccn_category="1",
                is_standard=True,
                sequencing_preference=1,
                rationale="Can be used in later lines if not given earlier"
            )
        ],
        biomarker_requirements=["HER2_amplification", "ERBB2_amplification"],
        contraindications=["severe_ILD_history"]
    ),
    
    # Tucatinib combo - Later-line TKI
    DrugPanelEntry(
        drug_name="tucatinib+trastuzumab+capecitabine",
        drug_class="tucatinib_combo",
        disease="breast_her2_positive",
        line_metadata=[
            DrugLineMetadata(
                line=3,
                nccn_category="1",
                is_standard=True,
                sequencing_preference=1,
                rationale="Preferred third-line, especially with brain metastases (HER2CLIMB)"
            ),
            DrugLineMetadata(
                line=4,
                nccn_category="1",
                is_standard=True,
                sequencing_preference=1,
                rationale="Subsequent therapy option"
            )
        ],
        biomarker_requirements=["HER2_amplification", "ERBB2_amplification"]
    ),
    
    # Neratinib combo
    DrugPanelEntry(
        drug_name="neratinib+capecitabine",
        drug_class="neratinib_combo",
        disease="breast_her2_positive",
        line_metadata=[
            DrugLineMetadata(
                line=3,
                nccn_category="2A",
                is_standard=True,
                sequencing_preference=2,
                rationale="Alternative TKI option after ≥2 prior HER2-directed regimens (NALA)"
            ),
            DrugLineMetadata(
                line=4,
                nccn_category="2A",
                is_standard=True,
                sequencing_preference=2,
                rationale="Later-line TKI option"
            )
        ],
        biomarker_requirements=["HER2_amplification", "ERBB2_amplification"]
    ),
    
    # Lapatinib combo
    DrugPanelEntry(
        drug_name="lapatinib+capecitabine",
        drug_class="lapatinib_combo",
        disease="breast_her2_positive",
        line_metadata=[
            DrugLineMetadata(
                line=2,
                nccn_category="1",
                is_standard=True,
                sequencing_preference=3,
                rationale="Alternative second-line TKI option"
            ),
            DrugLineMetadata(
                line=3,
                nccn_category="2A",
                is_standard=True,
                sequencing_preference=3,
                rationale="Later-line if other TKIs not available"
            )
        ],
        biomarker_requirements=["HER2_amplification", "ERBB2_amplification"]
    ),
    
    # Trastuzumab monotherapy (maintenance/later line)
    DrugPanelEntry(
        drug_name="trastuzumab+chemotherapy",
        drug_class="trastuzumab",
        disease="breast_her2_positive",
        line_metadata=[
            DrugLineMetadata(
                line=2,
                nccn_category="1",
                is_standard=True,
                sequencing_preference=4,
                rationale="Trastuzumab with alternative chemotherapy backbone"
            ),
            DrugLineMetadata(
                line=3,
                nccn_category="2A",
                is_standard=True,
                sequencing_preference=4,
                rationale="Continue HER2 blockade with different chemotherapy"
            )
        ],
        biomarker_requirements=["HER2_amplification", "ERBB2_amplification"]
    ),
]


# ===== PANEL ACCESS FUNCTIONS =====

def get_panel_for_disease(disease: str) -> List[DrugPanelEntry]:
    """
    Get the complete drug panel for a disease.
    
    Args:
        disease: Disease name ("ovarian_cancer", "breast_her2_positive")
        
    Returns:
        List of DrugPanelEntry objects
    """
    if disease == "ovarian_cancer":
        return OVARIAN_PANEL
    elif disease == "breast_her2_positive":
        return BREAST_HER2_PANEL
    else:
        return []


def get_drugs_for_line(disease: str, line: int) -> List[DrugPanelEntry]:
    """
    Get drugs appropriate for a specific treatment line.
    
    Args:
        disease: Disease name
        line: Treatment line number
        
    Returns:
        Filtered list of drugs appropriate for that line
    """
    panel = get_panel_for_disease(disease)
    
    result = []
    for entry in panel:
        # Check if drug has metadata for this line
        for line_meta in entry.line_metadata:
            if line_meta.line == line:
                result.append(entry)
                break
    
    return result


def get_line_metadata(drug_name: str, disease: str, line: int) -> Optional[DrugLineMetadata]:
    """
    Get line-specific metadata for a drug.
    
    Args:
        drug_name: Drug/regimen name
        disease: Disease name
        line: Treatment line
        
    Returns:
        DrugLineMetadata or None if not found
    """
    panel = get_panel_for_disease(disease)
    
    for entry in panel:
        if entry.drug_name.lower() == drug_name.lower():
            for line_meta in entry.line_metadata:
                if line_meta.line == line:
                    return line_meta
    
    return None


def calculate_line_appropriateness(
    drug_name: str,
    disease: str,
    current_line: int
) -> tuple[float, str]:
    """
    Calculate how appropriate a drug is for a given treatment line.
    
    Args:
        drug_name: Drug/regimen name
        disease: Disease name
        current_line: Current treatment line
        
    Returns:
        Tuple of (appropriateness_score, rationale)
        - appropriateness_score: 0.0-1.0
        - rationale: Explanation string
    """
    line_meta = get_line_metadata(drug_name, disease, current_line)
    
    if not line_meta:
        # Drug not in panel for this line
        return 0.6, f"No specific line guidance for {drug_name} at line {current_line}"
    
    # Perfect appropriateness if standard for this line
    if line_meta.is_standard and line_meta.nccn_category == "1":
        return 1.0, line_meta.rationale or f"Standard NCCN Category 1 for line {current_line}"
    
    # Good appropriateness for NCCN 2A
    if line_meta.is_standard and line_meta.nccn_category == "2A":
        return 0.9, line_meta.rationale or f"NCCN Category 2A for line {current_line}"
    
    # Moderate for NCCN 2B
    if line_meta.nccn_category == "2B":
        return 0.75, line_meta.rationale or f"NCCN Category 2B for line {current_line}"
    
    # Lower for NCCN 3
    if line_meta.nccn_category == "3":
        return 0.6, line_meta.rationale or f"NCCN Category 3 (less evidence) for line {current_line}"
    
    return 0.8, "Appropriate option for this line"

