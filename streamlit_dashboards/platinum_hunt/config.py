"""
Configuration for Platinum Response Dashboard
"""

from pathlib import Path
from typing import Dict

# Data file paths (relative to project root)
PROJECT_ROOT = Path(__file__).parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "validation"

# Data files
JR2_DATA_FILE = DATA_DIR / "tcga_ov_platinum_response_labels.json"
ZO_DATA_FILE = DATA_DIR / "tcga_ov_full_validation_dataset.json"

# Color scheme
COLORS = {
    "primary": "#1f77b4",
    "secondary": "#ff7f0e",
    "success": "#2ca02c",
    "danger": "#d62728",
    "warning": "#ff7f0e",
    "info": "#17a2b8",
    "sensitive": "#2ca02c",  # Green
    "resistant": "#d62728",  # Red
    "refractory": "#ff7f0e",  # Orange
    "background": "#f8f9fa",
    "text": "#212529",
}

# Response type labels
RESPONSE_LABELS = {
    "sensitive": "Sensitive",
    "resistant": "Resistant",
    "refractory": "Refractory",
}

# Response colors for charts
RESPONSE_COLORS = {
    "sensitive": COLORS["sensitive"],
    "resistant": COLORS["resistant"],
    "refractory": COLORS["refractory"],
}

# Validation thresholds
VALIDATION_THRESHOLDS = {
    "min_n": 40,
    "target_n": 100,
    "current_n": 161,  # Will be updated from data
}

# Page configuration
PAGE_CONFIG = {
    "page_title": "Platinum Response Data Hunt",
    "page_icon": "⚔️",
    "layout": "wide",
    "initial_sidebar_state": "expanded",
}


