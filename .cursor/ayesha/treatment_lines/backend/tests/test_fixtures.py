"""
Test Fixtures for Treatment Line Integration

Provides test cases for the 6 required smoke tests.
"""

from typing import Dict, Any

# ===== OVARIAN CANCER TEST CASES =====

OVARIAN_L1_CASE: Dict[str, Any] = {
    "name": "Ovarian 1st Line",
    "disease": "ovarian_cancer",
    "treatment_history": {
        "current_line": 1,
        "prior_therapies": [],
        "outcomes": None
    },
    "variants": [
        {"gene": "BRCA2", "hgvs_p": "c.5946delT", "consequence": "frameshift"}
    ],
    "expected": {
        "top_drug": "carboplatin",
        "top_drug_class": "platinum_agent",
        "line_appropriateness": 1.0,
        "cross_resistance_risk": 0.0,
        "sequencing_fitness": 1.0,
        "nccn_category": "1"
    }
}

OVARIAN_L2_POST_PLATINUM_CASE: Dict[str, Any] = {
    "name": "Ovarian 2nd Line Post-Platinum (Ayesha's Case)",
    "disease": "ovarian_cancer",
    "treatment_history": {
        "current_line": 2,
        "prior_therapies": ["carboplatin", "paclitaxel"],
        "outcomes": [
            {"therapy": "carboplatin+paclitaxel", "pfs_months": 8, "response": "partial"}
        ]
    },
    "variants": [
        {"gene": "BRCA2", "hgvs_p": "c.5946delT", "consequence": "frameshift"}
    ],
    "expected": {
        "top_drug": "olaparib",
        "top_drug_class": "PARP_inhibitor",
        "line_appropriateness": 1.0,
        "cross_resistance_risk": 0.4,  # DNA repair overlap with platinum
        "sequencing_fitness": 0.8,  # Reduced due to cross-resistance
        "nccn_category": "1"
    }
}

# ===== BREAST HER2+ TEST CASES =====

BREAST_HER2_L1_CASE: Dict[str, Any] = {
    "name": "Breast HER2+ 1st Line",
    "disease": "breast_her2_positive",
    "treatment_history": {
        "current_line": 1,
        "prior_therapies": [],
        "outcomes": None
    },
    "variants": [
        {"gene": "ERBB2", "hgvs_p": "amplification", "consequence": "amplification"}
    ],
    "expected": {
        "top_drug": "trastuzumab+pertuzumab+paclitaxel",
        "top_drug_class": "TP_taxane_combo",
        "line_appropriateness": 1.0,
        "cross_resistance_risk": 0.0,
        "sequencing_fitness": 1.0,
        "nccn_category": "1"
    }
}

BREAST_HER2_L2_POST_TRASTUZUMAB_CASE: Dict[str, Any] = {
    "name": "Breast HER2+ 2nd Line Post-Trastuzumab",
    "disease": "breast_her2_positive",
    "treatment_history": {
        "current_line": 2,
        "prior_therapies": ["trastuzumab+pertuzumab+paclitaxel"],
        "outcomes": [
            {"therapy": "trastuzumab+pertuzumab+paclitaxel", "pfs_months": 12, "response": "complete"}
        ]
    },
    "variants": [
        {"gene": "ERBB2", "hgvs_p": "amplification", "consequence": "amplification"}
    ],
    "expected": {
        "top_drug": "trastuzumab deruxtecan",
        "top_drug_class": "T-DXd",
        "line_appropriateness": 1.0,
        "cross_resistance_risk": 0.3,  # HER2 mechanism overlap
        "sequencing_fitness": 0.85,
        "nccn_category": "1"
    }
}

BREAST_HER2_L3_POST_TDXD_CASE: Dict[str, Any] = {
    "name": "Breast HER2+ 3rd Line Post-T-DXd (Dr. Lustberg's Case)",
    "disease": "breast_her2_positive",
    "treatment_history": {
        "current_line": 3,
        "prior_therapies": [
            "trastuzumab+pertuzumab+paclitaxel",
            "trastuzumab deruxtecan"
        ],
        "outcomes": [
            {"therapy": "trastuzumab+pertuzumab+paclitaxel", "pfs_months": 12, "response": "complete"},
            {"therapy": "trastuzumab deruxtecan", "pfs_months": 3, "response": "progression"}
        ]
    },
    "variants": [
        {"gene": "ERBB2", "hgvs_p": "amplification", "consequence": "amplification"}
    ],
    "expected": {
        "top_drug": "tucatinib+trastuzumab+capecitabine",
        "top_drug_class": "tucatinib_combo",
        "line_appropriateness": 1.0,
        "cross_resistance_risk": 0.2,  # Lower TKI cross-resistance
        "sequencing_fitness": 0.9,  # Good sequencing fitness
        "nccn_category": "1"
    }
}

# ===== EDGE CASE =====

EDGE_CASE_L4_EXHAUSTED: Dict[str, Any] = {
    "name": "Edge Case: 4th Line (Options Exhausted)",
    "disease": "ovarian_cancer",
    "treatment_history": {
        "current_line": 4,
        "prior_therapies": [
            "carboplatin+paclitaxel",
            "olaparib",
            "topotecan"
        ],
        "outcomes": [
            {"therapy": "carboplatin+paclitaxel", "pfs_months": 8, "response": "partial"},
            {"therapy": "olaparib", "pfs_months": 5, "response": "stable"},
            {"therapy": "topotecan", "pfs_months": 3, "response": "progression"}
        ]
    },
    "variants": [
        {"gene": "BRCA2", "hgvs_p": "c.5946delT", "consequence": "frameshift"}
    ],
    "expected": {
        "top_recommendation": "clinical_trial",
        "line_appropriateness": 0.6,  # Limited options remaining
        "cross_resistance_risk": 0.5,  # High due to multiple prior therapies
        "sequencing_fitness": 0.5,  # Lower due to exhausted options
        "show_trial_flag": True
    }
}

# All test cases for easy iteration
ALL_TEST_CASES = [
    OVARIAN_L1_CASE,
    OVARIAN_L2_POST_PLATINUM_CASE,
    BREAST_HER2_L1_CASE,
    BREAST_HER2_L2_POST_TRASTUZUMAB_CASE,
    BREAST_HER2_L3_POST_TDXD_CASE,
    EDGE_CASE_L4_EXHAUSTED
]


