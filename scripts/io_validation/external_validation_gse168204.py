#!/usr/bin/env python3
"""
EXTERNAL VALIDATION OF IO 2-PATHWAY MODEL ON GSE168204
======================================================

Validates the 2-pathway model (TIL_INFILTRATION + EXHAUSTION) trained on GSE91061
on an independent melanoma cohort (GSE168204, MGH cohort).

This is REAL external validation - no retraining, just apply model.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
import json
from datetime import datetime

# Paths
SCRIPT_DIR = Path(__file__).parent
IO_DIR = Path("/Users/fahadkiani/Desktop/development/crispr-assistant-main/scripts/data_acquisition/IO")
OUTPUT_DIR = Path("/Users/fahadkiani/Desktop/development/crispr-assistant-main/publications/06-io-response-prediction/data")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("EXTERNAL VALIDATION: GSE168204 (MGH Melanoma)")
print("Model trained on: GSE91061 (Riaz et al. Cell 2017)")
print("=" * 80)
print()

# Pathway gene lists (same as training)
IO_PATHWAYS = {
    'TIL_INFILTRATION': [
        'CD8A', 'CD8B', 'CD3D', 'CD3E', 'CD3G',
        'CD4', 'CD2', 'GZMA', 'GZMB', 'PRF1', 
        'IFNG', 'TNF', 'IL2'
    ],
    'EXHAUSTION': [
        'PDCD1', 'CTLA4', 'LAG3', 'TIGIT', 'HAVCR2',
        'BTLA', 'CD96', 'VSIR'
    ]
}

# Model coefficients (from GSE91061 training)
# Using 2-pathway model: TIL_INFILTRATION + EXHAUSTION
MODEL_COEFFICIENTS = {
    'intercept': -1.388,
    'TIL_INFILTRATION': 0.807,
    'EXHAUSTION': 0.633
}

def compute_pathway_score(expr_df, gene_list):
    """Compute pathway score as mean log2(TPM+1) expression."""
    available_genes = [g for g in gene_list if g in expr_df.index]
    if len(available_genes) == 0:
        return pd.Series(np.nan, index=expr_df.columns)
    expr_log = np.log2(expr_df.loc[available_genes] + 1)
    return expr_log.mean(axis=0)

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

# Step 1: Load GSE168204 expression data
print("📋 Step 1: Loading GSE168204 expression data...")
expr_file = IO_DIR / "GSE168204_MGH_counts.csv"
expr_df = pd.read_csv(expr_file, index_col=0)

# Clean up - remove empty columns and non-numeric columns
expr_df = expr_df.loc[:, ~expr_df.columns.str.contains('Unnamed')]

# Remove 'gene' column if present (duplicate index)
if 'gene' in expr_df.columns:
    expr_df = expr_df.drop(columns=['gene'])

# Convert all columns to numeric, coercing errors to NaN
for col in expr_df.columns:
    expr_df[col] = pd.to_numeric(expr_df[col], errors='coerce')

# Drop rows/cols with all NaN
expr_df = expr_df.dropna(axis=1, how='all')
expr_df = expr_df.dropna(axis=0, how='all')

print(f"✅ Loaded {expr_df.shape[1]} samples, {expr_df.shape[0]} genes")

# Step 2: Parse sample metadata from column names
print("\n🔍 Step 2: Parsing sample metadata...")

# Sample names contain patient info, dates, etc.
# e.g., "MGH200_FFPE_031116", "MGHIPIPD1001_041814"
# We need to identify responders vs non-responders

# Based on GSE168204 study design:
# MGH cohort samples are labeled with treatment info
# Let me check the sample names
sample_names = expr_df.columns.tolist()
print(f"  Sample names: {sample_names[:5]}...")

# Step 3: Download/create clinical labels
# For now, we'll try to infer from the published paper or GEO metadata
print("\n📊 Step 3: Creating clinical labels from GEO metadata...")

# From the GSE168204 paper (Jiang et al.), we know:
# - MGH cohort had responders (R) and non-responders (NR)
# - Pre-treatment and on-treatment samples
# 
# Since we don't have explicit labels, let's try to:
# 1. Compute pathway scores
# 2. Apply the model
# 3. Check if predictions make biological sense

# Compute pathway scores
print("\n🧬 Step 4: Computing pathway scores...")
pathway_scores = pd.DataFrame(index=expr_df.columns)

for pathway_name, gene_list in IO_PATHWAYS.items():
    score = compute_pathway_score(expr_df, gene_list)
    pathway_scores[pathway_name] = score
    
    # Count available genes
    available = [g for g in gene_list if g in expr_df.index]
    print(f"  {pathway_name}: {len(available)}/{len(gene_list)} genes available")

# Drop samples with missing scores
pathway_scores = pathway_scores.dropna()
print(f"  Samples with complete scores: {len(pathway_scores)}")

# Standardize (using training set statistics - should have these from GSE91061)
# For now, standardize within this cohort (suboptimal but workable)
print("\n📐 Step 5: Standardizing pathway scores...")
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
pathway_scores_std = pd.DataFrame(
    scaler.fit_transform(pathway_scores),
    index=pathway_scores.index,
    columns=pathway_scores.columns
)

# Apply model
print("\n🎯 Step 6: Applying 2-pathway model (no retraining)...")
linear_pred = (
    MODEL_COEFFICIENTS['intercept'] +
    MODEL_COEFFICIENTS['TIL_INFILTRATION'] * pathway_scores_std['TIL_INFILTRATION'] +
    MODEL_COEFFICIENTS['EXHAUSTION'] * pathway_scores_std['EXHAUSTION']
)

predicted_proba = sigmoid(linear_pred)
pathway_scores['predicted_proba'] = predicted_proba.values
pathway_scores['predicted_response'] = (predicted_proba > 0.5).astype(int).values

print(f"  Predictions:")
print(f"    Predicted responders: {pathway_scores['predicted_response'].sum()}")
print(f"    Predicted non-responders: {(1 - pathway_scores['predicted_response']).sum()}")

# Step 7: Try to get actual response labels from GEO
print("\n📥 Step 7: Fetching clinical metadata from GEO...")

# Download series matrix to extract sample characteristics
import subprocess

geo_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/GSE168204/matrix/GSE168204_series_matrix.txt.gz"
series_matrix = IO_DIR / "GSE168204_series_matrix.txt"

if not series_matrix.exists():
    print("  Downloading series matrix from GEO...")
    try:
        subprocess.run([
            "curl", "-o", str(IO_DIR / "GSE168204_series_matrix.txt.gz"), geo_url
        ], check=True, capture_output=True)
        subprocess.run([
            "gunzip", str(IO_DIR / "GSE168204_series_matrix.txt.gz")
        ], check=True, capture_output=True)
        print("  ✅ Downloaded")
    except Exception as e:
        print(f"  ⚠️ Download failed: {e}")

# Try to parse series matrix for response labels
if series_matrix.exists():
    print("  Parsing series matrix for clinical labels...")
    
    with open(series_matrix, 'r', errors='ignore') as f:
        lines = f.readlines()
    
    # Look for sample characteristics
    characteristics = {}
    sample_ids = []
    
    for line in lines:
        if line.startswith('!Sample_geo_accession'):
            sample_ids = line.strip().split('\t')[1:]
        elif line.startswith('!Sample_characteristics_ch1'):
            parts = line.strip().split('\t')[1:]
            # Parse characteristic
            key_val = parts[0].split(': ') if ': ' in parts[0] else [parts[0], '']
            if len(key_val) == 2:
                key = key_val[0]
                if key not in characteristics:
                    characteristics[key] = []
                for p in parts:
                    if ': ' in p:
                        characteristics[key].append(p.split(': ')[1])
                    else:
                        characteristics[key].append(p)
    
    print(f"    Found characteristics: {list(characteristics.keys())}")
    
    # Check if response data is available
    if 'response' in str(characteristics).lower():
        print("    ✅ Response labels found!")
    else:
        print("    ⚠️ Response labels not explicitly in series matrix")
        print("    Will need to cross-reference with paper supplement")

# Step 8: Summary without ground truth (for now)
print("\n" + "=" * 80)
print("EXTERNAL VALIDATION SUMMARY (Pending Ground Truth Labels)")
print("=" * 80)

print(f"""
Dataset: GSE168204 (MGH Melanoma Cohort)
Samples: {len(pathway_scores)}
Model: 2-pathway (TIL_INFILTRATION + EXHAUSTION)

Pathway Score Statistics (External Cohort):
  TIL_INFILTRATION: {pathway_scores['TIL_INFILTRATION'].mean():.3f} ± {pathway_scores['TIL_INFILTRATION'].std():.3f}
  EXHAUSTION: {pathway_scores['EXHAUSTION'].mean():.3f} ± {pathway_scores['EXHAUSTION'].std():.3f}

Model Predictions:
  Predicted Responders: {pathway_scores['predicted_response'].sum()} ({pathway_scores['predicted_response'].mean()*100:.1f}%)
  Predicted Non-responders: {(1-pathway_scores['predicted_response']).sum()} ({(1-pathway_scores['predicted_response'].mean())*100:.1f}%)

NEXT STEPS:
1. Obtain actual response labels from GSE168204 supplementary data
2. Compute external AUC, sensitivity, specificity
3. Update manuscript with external validation results
""")

# Save results
results = {
    "timestamp": datetime.now().isoformat(),
    "dataset": "GSE168204",
    "n_samples": len(pathway_scores),
    "model": "2-pathway (TIL_INFILTRATION + EXHAUSTION)",
    "coefficients": MODEL_COEFFICIENTS,
    "pathway_stats": {
        "TIL_INFILTRATION": {
            "mean": float(pathway_scores['TIL_INFILTRATION'].mean()),
            "std": float(pathway_scores['TIL_INFILTRATION'].std())
        },
        "EXHAUSTION": {
            "mean": float(pathway_scores['EXHAUSTION'].mean()),
            "std": float(pathway_scores['EXHAUSTION'].std())
        }
    },
    "predicted_responders": int(pathway_scores['predicted_response'].sum()),
    "status": "PENDING_GROUND_TRUTH"
}

results_file = OUTPUT_DIR / "external_validation_gse168204.json"
with open(results_file, 'w') as f:
    json.dump(results, f, indent=2)

print(f"✅ Results saved: {results_file}")

# Save predictions for later matching with ground truth
predictions_file = OUTPUT_DIR / "external_validation_predictions_gse168204.csv"
pathway_scores.to_csv(predictions_file)
print(f"✅ Predictions saved: {predictions_file}")

print("=" * 80)
