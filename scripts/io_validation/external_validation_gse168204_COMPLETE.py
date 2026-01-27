#!/usr/bin/env python3
"""
COMPLETE EXTERNAL VALIDATION OF IO 2-PATHWAY MODEL ON GSE168204
================================================================

Uses ground truth response labels from Nature Communications supplementary data.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix
from sklearn.preprocessing import StandardScaler
import json
from datetime import datetime

# Paths
IO_DIR = Path("/Users/fahadkiani/Desktop/development/crispr-assistant-main/scripts/data_acquisition/IO")
OUTPUT_DIR = Path("/Users/fahadkiani/Desktop/development/crispr-assistant-main/publications/06-io-response-prediction/data")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("EXTERNAL VALIDATION: GSE168204 (MGH Melanoma) - WITH GROUND TRUTH")
print("=" * 80)
print()

# Pathway gene lists
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

# Step 1: Load expression data
print("📋 Step 1: Loading GSE168204 expression data...")
expr_file = IO_DIR / "GSE168204_MGH_counts.csv"
expr_df = pd.read_csv(expr_file, index_col=0)
expr_df = expr_df.loc[:, ~expr_df.columns.str.contains('Unnamed')]
if 'gene' in expr_df.columns:
    expr_df = expr_df.drop(columns=['gene'])
for col in expr_df.columns:
    expr_df[col] = pd.to_numeric(expr_df[col], errors='coerce')
expr_df = expr_df.dropna(axis=1, how='all').dropna(axis=0, how='all')
print(f"✅ Loaded {expr_df.shape[1]} samples, {expr_df.shape[0]} genes")

# Step 2: Load clinical data with response labels
print("\n📊 Step 2: Loading clinical data with response labels...")
clinical_file = IO_DIR / "GSE168204_MGH_clinical.csv"
clinical_df = pd.read_csv(clinical_file)

# Clean up - remove header row and filter PRE-treatment only for fair comparison
clinical_df = clinical_df[clinical_df['Patient_ID'] != 'Patient ID']
clinical_df = clinical_df[clinical_df['Sex'] == 'PRE']  # PRE-treatment samples only

# Map response
clinical_df['response_binary'] = (clinical_df['Resp_NoResp'] == 'Response').astype(int)

print(f"  PRE-treatment samples: {len(clinical_df)}")
print(f"  Responders: {clinical_df['response_binary'].sum()}")
print(f"  Non-responders: {(1-clinical_df['response_binary']).sum()}")

# Step 3: Match samples to expression data
print("\n🔗 Step 3: Matching samples...")

# Expression sample names contain patient IDs (e.g., "MGH200_FFPE_031116")
# Clinical has patient IDs (e.g., "200")

# Create mapping
sample_to_patient = {}
for sample in expr_df.columns:
    # Extract patient number from sample name
    parts = sample.replace('MGH', '').split('_')
    patient_num = parts[0].replace('IPIPD', '').replace('PDL', '').replace('BI', '')
    # Remove any letters
    patient_num = ''.join(c for c in patient_num if c.isdigit())
    sample_to_patient[sample] = patient_num

# Match
matched_samples = []
for sample, patient_id in sample_to_patient.items():
    clinical_match = clinical_df[clinical_df['Patient_ID'] == patient_id]
    if len(clinical_match) > 0:
        response = clinical_match.iloc[0]['response_binary']
        matched_samples.append({
            'sample': sample,
            'patient_id': patient_id,
            'response': response
        })

matched_df = pd.DataFrame(matched_samples)
print(f"  Matched samples: {len(matched_df)}")
print(f"  Responders: {matched_df['response'].sum()}")
print(f"  Non-responders: {(1-matched_df['response']).sum()}")

# Step 4: Compute pathway scores
print("\n🧬 Step 4: Computing pathway scores...")
pathway_scores = pd.DataFrame(index=matched_df['sample'])

for pathway_name, gene_list in IO_PATHWAYS.items():
    score = compute_pathway_score(expr_df[matched_df['sample'].tolist()], gene_list)
    pathway_scores[pathway_name] = score.values
    available = [g for g in gene_list if g in expr_df.index]
    print(f"  {pathway_name}: {len(available)}/{len(gene_list)} genes")

pathway_scores = pathway_scores.dropna()
matched_df = matched_df[matched_df['sample'].isin(pathway_scores.index)].reset_index(drop=True)
pathway_scores = pathway_scores.loc[matched_df['sample']]

print(f"  Samples with complete scores: {len(pathway_scores)}")

# Step 5: Standardize and apply model
print("\n🎯 Step 5: Applying 2-pathway model (NO RETRAINING)...")

scaler = StandardScaler()
pathway_scores_std = pd.DataFrame(
    scaler.fit_transform(pathway_scores),
    index=pathway_scores.index,
    columns=pathway_scores.columns
)

linear_pred = (
    MODEL_COEFFICIENTS['intercept'] +
    MODEL_COEFFICIENTS['TIL_INFILTRATION'] * pathway_scores_std['TIL_INFILTRATION'].values +
    MODEL_COEFFICIENTS['EXHAUSTION'] * pathway_scores_std['EXHAUSTION'].values
)

predicted_proba = sigmoid(linear_pred)
matched_df['predicted_proba'] = predicted_proba
matched_df['predicted_response'] = (predicted_proba > 0.5).astype(int)

# Step 6: Compute metrics
print("\n📈 Step 6: Computing external validation metrics...")

y_true = matched_df['response'].values
y_pred_proba = matched_df['predicted_proba'].values
y_pred = matched_df['predicted_response'].values

# AUC
if len(np.unique(y_true)) > 1:
    auc = roc_auc_score(y_true, y_pred_proba)
    print(f"\n  ✅ EXTERNAL AUC: {auc:.3f}")
    
    # Confusion matrix
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0
    
    print(f"  Sensitivity: {sensitivity:.3f}")
    print(f"  Specificity: {specificity:.3f}")
    print(f"  PPV: {ppv:.3f}")
    print(f"  NPV: {npv:.3f}")
    
    print(f"\n  Confusion Matrix:")
    print(f"    TP={tp}, FP={fp}")
    print(f"    FN={fn}, TN={tn}")
else:
    auc = None
    print("  ⚠️ Cannot compute AUC (only one class present)")

# Step 7: Save results
print("\n💾 Step 7: Saving results...")

results = {
    "timestamp": datetime.now().isoformat(),
    "dataset": "GSE168204 (MGH Melanoma)",
    "validation_type": "EXTERNAL VALIDATION",
    "training_dataset": "GSE91061 (Riaz et al. Cell 2017)",
    "n_samples": len(matched_df),
    "n_responders": int(matched_df['response'].sum()),
    "n_non_responders": int((1-matched_df['response']).sum()),
    "model": "2-pathway (TIL_INFILTRATION + EXHAUSTION)",
    "external_auc": float(auc) if auc else None,
    "sensitivity": float(sensitivity) if auc else None,
    "specificity": float(specificity) if auc else None,
    "ppv": float(ppv) if auc else None,
    "npv": float(npv) if auc else None,
    "status": "COMPLETE"
}

results_file = OUTPUT_DIR / "external_validation_gse168204_COMPLETE.json"
with open(results_file, 'w') as f:
    json.dump(results, f, indent=2)

predictions_file = OUTPUT_DIR / "external_validation_predictions_gse168204_COMPLETE.csv"
matched_df.to_csv(predictions_file, index=False)

print(f"✅ Results: {results_file}")
print(f"✅ Predictions: {predictions_file}")

# Final summary
print("\n" + "=" * 80)
print("EXTERNAL VALIDATION COMPLETE")
print("=" * 80)
print(f"""
TRAINING DATA: GSE91061 (Riaz et al. Cell 2017, n=105)
  - Training AUC: 0.621
  - Nested CV AUC: 0.601 ± 0.071
  - Held-out test AUC: 0.806

EXTERNAL VALIDATION DATA: GSE168204 (MGH Melanoma, n={len(matched_df)})
  - Responders: {int(matched_df['response'].sum())}
  - Non-responders: {int((1-matched_df['response']).sum())}
  - EXTERNAL AUC: {auc:.3f if auc else 'N/A'}
  - Sensitivity: {sensitivity:.3f if auc else 'N/A'}
  - Specificity: {specificity:.3f if auc else 'N/A'}

CONCLUSION:
{'✅ EXTERNAL VALIDATION SUCCESSFUL' if auc and auc > 0.5 else '⚠️ EXTERNAL VALIDATION INCONCLUSIVE'}
{'Model maintains predictive power on independent cohort!' if auc and auc > 0.6 else 'Model requires refinement for external validation.'}
""")
print("=" * 80)
