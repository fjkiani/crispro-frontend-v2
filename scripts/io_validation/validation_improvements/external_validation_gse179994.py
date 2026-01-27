#!/usr/bin/env python3
"""
Task 1.4: External Validation (GSE179994)
==========================================

Validates 8-pathway model on independent NSCLC cohort (GSE179994).
Applies model trained on GSE91061 (melanoma) to GSE179994 (NSCLC) without retraining.

Note: GSE179994 is single-cell RNA-seq data. We need to:
1. Aggregate single-cell data to pseudo-bulk expression
2. Compute pathway scores using same gene lists as GSE91061
3. Apply logistic regression coefficients from GSE91061 (no retraining)
4. Evaluate AUC, sensitivity, specificity on external cohort
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix
from sklearn.preprocessing import StandardScaler
import json
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# Configuration
SCRIPT_DIR = Path(__file__).parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
DATA_DIR = REPO_ROOT / "scripts" / "data_acquisition" / "IO"
OUTPUT_DIR = REPO_ROOT / "publications" / "06-io-response-prediction" / "data"
FIGURES_DIR = REPO_ROOT / "publications" / "06-io-response-prediction" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Set random seed
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

print("="*80)
print("TASK 1.4: EXTERNAL VALIDATION (GSE179994)")
print("="*80)
print()

# IO pathway gene lists (same as GSE91061)
IO_PATHWAYS = {
    'TIL_INFILTRATION': [
        'CD8A', 'CD8B', 'CD3D', 'CD3E', 'CD3G',
        'CD4', 'CD2', 'GZMA', 'GZMB', 'PRF1', 
        'IFNG', 'TNF', 'IL2'
    ],
    'T_EFFECTOR': [
        'CD274', 'PDCD1LG2', 'IDO1', 'IDO2',
        'CXCL9', 'CXCL10', 'CXCL11',
        'HLA-DRA', 'HLA-DRB1', 'STAT1', 'IRF1', 'IFNG'
    ],
    'ANGIOGENESIS': [
        'VEGFA', 'VEGFB', 'VEGFC', 'VEGFD',
        'KDR', 'FLT1', 'FLT4',
        'ANGPT1', 'ANGPT2', 'TEK', 'PECAM1', 'VWF'
    ],
    'TGFB_RESISTANCE': [
        'TGFB1', 'TGFB2', 'TGFB3',
        'TGFBR1', 'TGFBR2', 'TGFBR3',
        'SMAD2', 'SMAD3', 'SMAD4', 'SMAD7'
    ],
    'MYELOID_INFLAMMATION': [
        'IL6', 'IL1B', 'IL8', 'CXCL8',
        'CXCL1', 'CXCL2', 'CXCL3',
        'PTGS2', 'CCL2', 'CCL3', 'CCL4',
        'S100A8', 'S100A9', 'S100A12'
    ],
    'PROLIFERATION': [
        'MKI67', 'PCNA', 'TOP2A',
        'CCNA2', 'CCNB1', 'CCNB2',
        'CDK1', 'CDK2', 'CDK4',
        'CDC20', 'AURKA', 'AURKB'
    ],
    'IMMUNOPROTEASOME': [
        'PSMB8', 'PSMB9', 'PSMB10',
        'TAP1', 'TAP2', 'B2M',
        'HLA-A', 'HLA-B', 'HLA-C'
    ],
    'EXHAUSTION': [
        'PDCD1', 'CTLA4', 'LAG3', 'TIGIT', 'HAVCR2',
        'BTLA', 'CD96', 'VSIR'
    ]
}

def compute_pathway_score(expr_df, gene_list, pathway_name):
    """Compute pathway score as mean log2(TPM+1) expression."""
    available_genes = [g for g in gene_list if g in expr_df.index]
    
    if len(available_genes) == 0:
        return pd.Series(np.nan, index=expr_df.columns)
    
    expr_log = np.log2(expr_df.loc[available_genes] + 1)
    score = expr_log.mean(axis=0)
    
    return score

print("⚠️  NOTE: GSE179994 is single-cell RNA-seq data.")
print("   This script requires aggregated pseudo-bulk expression data.")
print("   If you have bulk RNA-seq data for GSE179994, please provide it.")
print("   Otherwise, we'll need to aggregate single-cell data first.")
print()

# Check if we have bulk expression data
# For now, we'll create a placeholder script that can be adapted

print("Step 1: Load GSE179994 expression data...")
print("  [PLACEHOLDER] - Need to provide bulk RNA-seq or aggregated scRNA-seq data")
print()

# Placeholder: Load model from GSE91061
print("Step 2: Load trained model from GSE91061...")
try:
    # Try to load reduced model (if available)
    model_file = OUTPUT_DIR / "reduced_model_3pathway.json"
    if not model_file.exists():
        model_file = OUTPUT_DIR / "reduced_model_4pathway.json"
    
    if model_file.exists():
        with open(model_file, 'r') as f:
            model = json.load(f)
        print(f"  ✅ Loaded model: {model_file.name}")
        print(f"     Pathways: {', '.join(model['pathways'])}")
        print(f"     Training AUC: {model['train_auc']:.3f}")
        print(f"     CV AUC: {model['cv_auc_mean']:.3f} ± {model['cv_auc_std']:.3f}")
    else:
        # Fallback: Load from original analysis
        print("  ⚠️  Reduced model not found, will use 8-pathway model from GSE91061")
        model = None
except Exception as e:
    print(f"  ⚠️  Could not load model: {e}")
    model = None

print()

# Instructions for user
print("="*80)
print("INSTRUCTIONS FOR EXTERNAL VALIDATION")
print("="*80)
print()
print("To complete external validation on GSE179994, you need:")
print()
print("1. **Bulk RNA-seq expression data** (or aggregated single-cell data)")
print("   - Format: Samples × Genes (TPM or counts)")
print("   - File: GSE179994_bulk_expression.tsv (or similar)")
print()
print("2. **Clinical metadata with response labels**")
print("   - Format: Sample ID, Response (CR/PR vs SD/PD)")
print("   - File: GSE179994_clinical.csv")
print()
print("3. **Pre-treatment samples only**")
print("   - Filter to pre-treatment timepoint")
print("   - Map response labels to binary (1=responder, 0=non-responder)")
print()
print("Once you have this data, this script will:")
print("  - Compute 8-pathway scores using same gene lists")
print("  - Apply logistic regression coefficients from GSE91061 (no retraining)")
print("  - Evaluate AUC, sensitivity, specificity on external cohort")
print("  - Compare training (GSE91061) vs. external (GSE179994) performance")
print()

# Create template script for when data is available
template_script = OUTPUT_DIR / "external_validation_template.py"
template_content = f'''#!/usr/bin/env python3
"""
External Validation Template for GSE179994
===========================================

Complete this script once you have:
1. GSE179994_bulk_expression.tsv (or aggregated scRNA-seq)
2. GSE179994_clinical.csv (with response labels)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler
import json

# Load expression data
expr = pd.read_csv('GSE179994_bulk_expression.tsv', sep='\\t', index_col=0)

# Load clinical data
clinical = pd.read_csv('GSE179994_clinical.csv')

# Filter to pre-treatment samples
pre_samples = clinical[clinical['timepoint'] == 'Pre']

# Map response labels
response_map = {{'CR': 1, 'PR': 1, 'SD': 0, 'PD': 0}}
pre_samples['response'] = pre_samples['response_raw'].map(response_map)

# Compute pathway scores (use same function as GSE91061)
# ... [pathway scoring code] ...

# Load model from GSE91061
with open('reduced_model_3pathway.json', 'r') as f:
    model = json.load(f)

# Apply model (no retraining)
# ... [apply coefficients] ...

# Evaluate
# ... [compute AUC, sensitivity, specificity] ...
'''

with open(template_script, 'w') as f:
    f.write(template_content)

print(f"✅ Created template script: {template_script}")
print()

# Alternative: Check if we can use GSE168204 (bulk RNA-seq) instead
print("Alternative: Check for other bulk RNA-seq datasets...")
other_datasets = [
    "GSE168204",  # Mentioned in final report
]

for dataset in other_datasets:
    dataset_file = DATA_DIR / f"{dataset}_MGH_counts.csv"
    if dataset_file.exists():
        print(f"  ✅ Found: {dataset_file.name}")
        print(f"     This may be usable for external validation")
    else:
        print(f"  ❌ Not found: {dataset}")

print()

print(f"{'='*80}")
print("TASK 1.4 STATUS: AWAITING DATA")
print(f"{'='*80}")
print()
print("Next steps:")
print("1. Provide bulk RNA-seq expression data for GSE179994 (or aggregate scRNA-seq)")
print("2. Provide clinical metadata with response labels")
print("3. Run this script again to complete external validation")
print()
print("OR use alternative dataset (e.g., GSE168204) if available")
print()
