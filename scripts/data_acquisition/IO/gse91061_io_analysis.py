#!/usr/bin/env python3
"""
AGENT B - GSE91061 IO RESPONSE PREDICTION
Timeline: 4-6 hours (complete by tomorrow morning)
Goal: AUC > 0.75 (vs PD-L1 ~0.60-0.65)
"""

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Set working directory
os.chdir('/Users/fahadkiani/Desktop/development/crispr-assistant-main/scripts/data_acquisition/IO')

print("="*80)
print("GSE91061 IO RESPONSE PREDICTION ANALYSIS")
print("="*80)

# ============================================================
# STEP 1: DATA LOADING (30 MIN)
# ============================================================

print("\n[STEP 1] Loading expression data...")

# Load expression data
expr = pd.read_csv('GSE91061_norm_counts_TPM_GRCh38.p13_NCBI.tsv', 
                   sep='\t', index_col=0)

print(f"Expression matrix: {expr.shape}")
print(f"Genes: {len(expr.index)}")
print(f"Samples: {len(expr.columns)}")

# Load cytolytic scores (benchmark comparison)
cyto = pd.read_csv('GSE91061_BMS038109Sample_Cytolytic_Score_20161026.txt',
                   sep='\t', header=None, names=['sample_id', 'cytolytic_score'])

print(f"Cytolytic scores: {len(cyto)} samples")

# Filter to PRE-TREATMENT samples only
# Sample naming: "Pt101_Pre_AD681975-5" or "GSM2420259"
pre_samples = []
for col in expr.columns:
    # Check if it's a GSM ID (need to map from cytolytic file)
    if 'Pre' in str(col) or 'pre' in str(col).lower():
        pre_samples.append(col)
    # Also check cytolytic file for Pre samples
    if col in cyto['sample_id'].values:
        sample_name = cyto[cyto['sample_id'] == col]['sample_id'].values[0]
        if 'Pre' in sample_name or 'pre' in sample_name.lower():
            pre_samples.append(col)

# Alternative: Extract from cytolytic file directly
cyto_pre = cyto[cyto['sample_id'].str.contains('Pre', case=False, na=False)]
print(f"Pre-treatment samples from cytolytic file: {len(cyto_pre)}")

# Map cytolytic sample IDs to expression column names
# Need to match sample identifiers
# For now, use all samples and filter by response data later

expr_pre = expr.copy()  # Will filter after parsing response

# ============================================================
# STEP 2: CLINICAL METADATA PARSING (30 MIN)
# ============================================================

print("\n[STEP 2] Parsing clinical metadata...")

# Load series matrix for clinical data
series_matrix_path = 'GSE91061_series_matrix.txt'
print(f"Loading series matrix from {series_matrix_path}...")

# Read series matrix (skip comment lines starting with !)
series_matrix = pd.read_csv(series_matrix_path, sep='\t', comment='!')

# Alternative: Read all lines and parse manually
with open(series_matrix_path, 'r') as f:
    lines = f.readlines()

# Find rows with response and GSM IDs
gsm_ids = []
response_map = {}

# Find the row with Sample_geo_accession (GSM IDs)
# Find the row with response information
for i, line in enumerate(lines):
    if line.startswith('!Sample_geo_accession'):
        # Extract GSM IDs from this line
        parts = line.strip().split('\t')
        if len(parts) > 1:
            gsm_line = parts[1:]  # Skip the first column (header)
            # Remove quotes and split
            gsm_ids = [g.strip('"') for g in gsm_line if g.strip()]
    elif line.startswith('!Sample_characteristics_ch1') and 'response:' in line:
        # Extract response labels
        parts = line.strip().split('\t')
        if len(parts) > 1:
            response_line = parts[1:]  # Skip the first column (header)
            # Map GSM IDs to responses
            for j, resp_str in enumerate(response_line):
                if j < len(gsm_ids):
                    gsm_id = gsm_ids[j].strip('"')
                    resp_str = resp_str.strip('"')
                    # Parse response: "response: PD" -> "PD"
                    if 'response:' in resp_str:
                        resp = resp_str.split('response:')[1].strip()
                        response_map[gsm_id] = resp

print(f"✓ Extracted response labels for {len(response_map)} samples")

# Map response categories to binary labels
# PRCR = PR + CR (responders) = 1
# PD = Progressive Disease (non-responder) = 0
# SD = Stable Disease (non-responder) = 0
# UNK = Unknown = NaN

response_binary = {}
for gsm_id, resp in response_map.items():
    if 'PRCR' in resp or 'PR' in resp or 'CR' in resp:
        response_binary[gsm_id] = 1  # Responder
    elif 'PD' in resp:
        response_binary[gsm_id] = 0  # Non-responder
    elif 'SD' in resp:
        response_binary[gsm_id] = 0  # Non-responder
    else:
        response_binary[gsm_id] = np.nan  # Unknown

# Parse cytolytic file to identify pre-treatment samples
cyto_df = pd.read_csv('GSE91061_BMS038109Sample_Cytolytic_Score_20161026.txt',
                      sep='\t', header=None, names=['sample_id', 'cytolytic_score'])

# Extract patient IDs and timepoints from sample names
# Format: "Pt101_Pre_AD681975-5" -> Patient 101, Pre-treatment
cyto_df['patient_id'] = cyto_df['sample_id'].str.extract(r'(Pt\d+)', expand=False)
cyto_df['timepoint'] = cyto_df['sample_id'].str.extract(r'_(Pre|On)_', expand=False)

# Filter to pre-treatment samples only
cyto_pre = cyto_df[cyto_df['timepoint'] == 'Pre'].copy()
print(f"✓ Found {len(cyto_pre)} pre-treatment samples in cytolytic file")

# Expression matrix columns are GSM IDs
# Filter to pre-treatment samples based on cytolytic file
# Need to match GSM IDs to cytolytic sample IDs
# For now, use all expression samples that have response labels

# Create clinical dataframe with GSM IDs from expression matrix
clinical_df = pd.DataFrame({
    'sample_id': expr.columns.tolist(),  # GSM IDs
    'response': [response_binary.get(sid, np.nan) for sid in expr.columns],
    'response_raw': [response_map.get(sid, 'UNK') for sid in expr.columns],
    'cytolytic_score': np.nan  # Will fill if we can match
})

# Try to match cytolytic scores to GSM IDs
# This would require a mapping file or additional metadata
# For now, leave as NaN

# Filter to samples with valid response labels (exclude UNK)
clinical_df_valid = clinical_df[clinical_df['response'].notna()].copy()

print(f"Clinical data structure: {len(clinical_df)} total samples")
print(f"  With response labels: {len(clinical_df_valid)} samples")
print(f"  Responders (PRCR): {(clinical_df_valid['response'] == 1).sum()}")
print(f"  Non-responders (PD/SD): {(clinical_df_valid['response'] == 0).sum()}")

# Save processed clinical data
clinical_df.to_csv('gse91061_clinical_processed.csv', index=False)

# ============================================================
# STEP 3: IO PATHWAY SCORE CALCULATION (1 HOUR)
# ============================================================

print("\n[STEP 3] Computing IO pathway scores...")

# Load GeneID → Symbol mapping
print("Loading GeneID → Symbol mapping...")
import json
with open('gse91061_symbol_to_geneid.json', 'r') as f:
    symbol_to_geneid = json.load(f)

# Create GeneID → Symbol mapping (reverse)
with open('gse91061_geneid_to_symbol.json', 'r') as f:
    geneid_to_symbol = json.load(f)

print(f"✓ Loaded {len(geneid_to_symbol)} GeneID mappings")

# Map expression matrix index from GeneID to Symbol
print("Mapping expression matrix index...")
expr_index_symbols = []
for gene_id in expr.index:
    symbol = geneid_to_symbol.get(str(gene_id), None)
    if symbol:
        expr_index_symbols.append(symbol)
    else:
        expr_index_symbols.append(f"GENE_{gene_id}")  # Fallback

expr_mapped = expr.copy()
expr_mapped.index = expr_index_symbols

print(f"✓ Mapped {sum(1 for s in expr_index_symbols if not s.startswith('GENE_'))}/{len(expr_index_symbols)} genes to symbols")

# Curated IO-relevant pathways
IO_PATHWAYS = {
    'TIL_INFILTRATION': [
        'CD8A', 'CD8B', 'CD3D', 'CD3E', 'CD3G',
        'CD4', 'CD2', 'GZMA', 'GZMB', 'PRF1', 
        'IFNG', 'TNF', 'IL2'
    ],
    
    'T_EFFECTOR': [
        'CD274',      # PD-L1
        'PDCD1LG2',   # PD-L2
        'IDO1', 'IDO2',
        'CXCL9', 'CXCL10', 'CXCL11',  # IFNγ-induced chemokines
        'HLA-DRA', 'HLA-DRB1',
        'STAT1', 'IRF1', 'IFNG'
    ],
    
    'ANGIOGENESIS': [
        'VEGFA', 'VEGFB', 'VEGFC', 'VEGFD',
        'KDR', 'FLT1', 'FLT4',  # VEGF receptors
        'ANGPT1', 'ANGPT2', 'TEK',
        'PECAM1', 'VWF'
    ],
    
    'TGFB_RESISTANCE': [
        'TGFB1', 'TGFB2', 'TGFB3',
        'TGFBR1', 'TGFBR2', 'TGFBR3',
        'SMAD2', 'SMAD3', 'SMAD4', 'SMAD7'
    ],
    
    'MYELOID_INFLAMMATION': [
        'IL6', 'IL1B', 'IL8', 'CXCL8',
        'CXCL1', 'CXCL2', 'CXCL3',
        'PTGS2',  # COX2
        'CCL2', 'CCL3', 'CCL4',
        'S100A8', 'S100A9', 'S100A12'
    ],
    
    'PROLIFERATION': [
        'MKI67', 'PCNA', 'TOP2A',
        'CCNA2', 'CCNB1', 'CCNB2',
        'CDK1', 'CDK2', 'CDK4',
        'CDC20', 'AURKA', 'AURKB'
    ],
    
    'IMMUNOPROTEASOME': [
        'PSMB8', 'PSMB9', 'PSMB10',  # Immunoproteasome subunits
        'TAP1', 'TAP2',  # Antigen processing
        'B2M',  # Beta-2-microglobulin
        'HLA-A', 'HLA-B', 'HLA-C'  # MHC-I
    ],
    
    'EXHAUSTION': [
        'PDCD1',   # PD-1
        'CTLA4', 'LAG3', 'TIGIT', 'HAVCR2',  # TIM-3
        'BTLA', 'CD96', 'VSIR'  # VISTA
    ]
}

def compute_pathway_score(expr_df, gene_list, pathway_name):
    """
    Compute pathway score as mean log2(TPM+1) expression.
    
    Returns:
        pd.Series: Pathway scores for each sample
    """
    # Find genes present in expression matrix (now using symbols)
    available_genes = [g for g in gene_list if g in expr_df.index]
    
    if len(available_genes) == 0:
        print(f"⚠️  WARNING: No genes found for {pathway_name}")
        return pd.Series(np.nan, index=expr_df.columns)
    
    # Log2 transform (TPM is not log-transformed)
    expr_log = np.log2(expr_df.loc[available_genes] + 1)
    
    # Mean expression across genes
    score = expr_log.mean(axis=0)
    
    coverage = len(available_genes) / len(gene_list) * 100
    print(f"✓ {pathway_name}: {len(available_genes)}/{len(gene_list)} genes ({coverage:.1f}%)")
    
    return score

# Compute all pathway scores
print("\nComputing pathway scores...")
pathway_scores = pd.DataFrame(index=expr_mapped.columns)

for pathway_name, gene_list in IO_PATHWAYS.items():
    pathway_scores[pathway_name] = compute_pathway_score(
        expr_mapped, gene_list, pathway_name
    )

# Add PD-L1 expression (single gene benchmark)
if 'CD274' in expr_mapped.index:
    pathway_scores['PDL1_EXPRESSION'] = np.log2(expr_mapped.loc['CD274'] + 1)
    print("✓ PD-L1 (CD274) expression added")

# Save pathway scores
pathway_scores.to_csv('gse91061_io_pathway_scores.csv')
print(f"\n✓ Pathway scores computed: {pathway_scores.shape}")

# ============================================================
# STEP 4: RESPONSE ASSOCIATION ANALYSIS (1 HOUR)
# ============================================================

print("\n" + "="*80)
print("[STEP 4] RESPONSE ASSOCIATION ANALYSIS")
print("="*80)

# Merge pathway scores with response labels
# pathway_scores index is sample IDs (GSM IDs), clinical_df has 'sample_id' column
analysis_df = pathway_scores.merge(
    clinical_df[['sample_id', 'response']], 
    left_index=True, 
    right_on='sample_id',
    how='inner'
)

# Remove samples without response labels
analysis_df = analysis_df.dropna(subset=['response'])

print(f"\nAnalysis cohort: {len(analysis_df)} samples")
print(f"  Responders: {(analysis_df['response'] == 1).sum()}")
print(f"  Non-responders: {(analysis_df['response'] == 0).sum()}")

# Test each pathway
results = []

for pathway in IO_PATHWAYS.keys():
    if pathway not in analysis_df.columns:
        continue
        
    scores = analysis_df[pathway].values
    labels = analysis_df['response'].values
    
    # Remove NaN values
    valid_idx = ~np.isnan(scores)
    scores = scores[valid_idx]
    labels = labels[valid_idx]
    
    if len(scores) < 10:
        continue
    
    # Responder vs non-responder scores
    resp_scores = scores[labels == 1]
    nonresp_scores = scores[labels == 0]
    
    # Mann-Whitney U test (one-sided: responders have higher scores)
    u_stat, p_val = mannwhitneyu(resp_scores, nonresp_scores, 
                                  alternative='greater')
    
    # AUC-ROC
    try:
        auc = roc_auc_score(labels, scores)
    except:
        auc = np.nan
    
    # Effect size (Cohen's d)
    pooled_std = np.sqrt((resp_scores.std()**2 + nonresp_scores.std()**2) / 2)
    if pooled_std > 0:
        cohens_d = (resp_scores.mean() - nonresp_scores.mean()) / pooled_std
    else:
        cohens_d = np.nan
    
    results.append({
        'pathway': pathway,
        'auc': auc,
        'p_value': p_val,
        'cohens_d': cohens_d,
        'resp_mean': resp_scores.mean(),
        'nonresp_mean': nonresp_scores.mean(),
        'n_responders': len(resp_scores),
        'n_nonresponders': len(nonresp_scores)
    })

# Add PD-L1 if available
if 'PDL1_EXPRESSION' in analysis_df.columns:
    scores = analysis_df['PDL1_EXPRESSION'].values
    labels = analysis_df['response'].values
    valid_idx = ~np.isnan(scores)
    scores, labels = scores[valid_idx], labels[valid_idx]
    
    resp_scores = scores[labels == 1]
    nonresp_scores = scores[labels == 0]
    
    u_stat, p_val = mannwhitneyu(resp_scores, nonresp_scores, alternative='greater')
    auc_pdl1 = roc_auc_score(labels, scores)
    
    pooled_std = np.sqrt((resp_scores.std()**2 + nonresp_scores.std()**2) / 2)
    if pooled_std > 0:
        cohens_d = (resp_scores.mean() - nonresp_scores.mean()) / pooled_std
    else:
        cohens_d = np.nan
    
    results.append({
        'pathway': 'PDL1_EXPRESSION',
        'auc': auc_pdl1,
        'p_value': p_val,
        'cohens_d': cohens_d,
        'resp_mean': resp_scores.mean(),
        'nonresp_mean': nonresp_scores.mean(),
        'n_responders': len(resp_scores),
        'n_nonresponders': len(nonresp_scores)
    })

# Convert to DataFrame and sort by AUC
results_df = pd.DataFrame(results).sort_values('auc', ascending=False)

print("\n" + "="*80)
print("PATHWAY RESPONSE ASSOCIATION RESULTS")
print("="*80)
print(results_df.to_string(index=False))
print("="*80)

# Save results
results_df.to_csv('gse91061_pathway_response_association.csv', index=False)

# ============================================================
# STEP 5: MULTI-PATHWAY COMPOSITE SCORE (1 HOUR)
# ============================================================

print("\n" + "="*80)
print("[STEP 5] MULTI-PATHWAY COMPOSITE SCORE")
print("="*80)

# Strategy 1: Weighted composite (biologically informed)
# Positive pathways (higher = better response)
# Negative pathways (higher = worse response)
composite_weighted = (
    0.25 * analysis_df['TIL_INFILTRATION'] +
    0.25 * analysis_df['T_EFFECTOR'] +
    0.15 * analysis_df['IMMUNOPROTEASOME'] +
    0.10 * analysis_df['ANGIOGENESIS'] +
    (-0.15) * analysis_df['TGFB_RESISTANCE'] +  # Negative (resistance)
    (-0.10) * analysis_df['EXHAUSTION']  # Negative (resistance)
)

# Test composite
labels = analysis_df['response'].values
valid_composite = ~np.isnan(composite_weighted)
auc_composite = roc_auc_score(labels[valid_composite], composite_weighted[valid_composite])

print(f"\n✓ Weighted Composite AUC: {auc_composite:.3f}")

# Strategy 2: Data-driven weights (logistic regression)
# Prepare features (only pathways that exist in analysis_df)
available_pathways = [p for p in IO_PATHWAYS.keys() if p in analysis_df.columns]

if len(available_pathways) > 0:
    X = analysis_df[available_pathways].values
    y = labels
    
    # Remove rows with any NaN
    valid_rows = ~np.isnan(X).any(axis=1)
    X = X[valid_rows]
    y = y[valid_rows]
    
    if len(X) > 10:
        # Standardize features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        
        # Logistic regression with L2 regularization
        lr = LogisticRegression(penalty='l2', C=1.0, random_state=42, max_iter=1000)
        lr.fit(X_scaled, y)
        
        # Predictions
        composite_lr = lr.predict_proba(X_scaled)[:, 1]
        auc_lr = roc_auc_score(y, composite_lr)
        
        print(f"✓ Logistic Regression AUC: {auc_lr:.3f}")
        
        # Cross-validation AUC
        cv_scores = cross_val_score(lr, X_scaled, y, cv=5, scoring='roc_auc')
        print(f"✓ 5-Fold CV AUC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
        
        # Feature importance (coefficients)
        feature_importance = pd.DataFrame({
            'pathway': available_pathways,
            'coefficient': lr.coef_[0]
        }).sort_values('coefficient', ascending=False)
        
        print("\nLogistic Regression Feature Importance:")
        print(feature_importance.to_string(index=False))
        
        # Save composite scores
        analysis_df['composite_weighted'] = composite_weighted
        analysis_df['composite_lr'] = np.nan
        analysis_df.loc[valid_rows, 'composite_lr'] = composite_lr
        analysis_df.to_csv('gse91061_analysis_with_composites.csv', index=False)
    else:
        print("⚠️  Insufficient valid samples for logistic regression")
        auc_lr = np.nan
        cv_scores = np.array([np.nan])
else:
    print("⚠️  No pathways available for composite scoring")
    auc_lr = np.nan
    cv_scores = np.array([np.nan])

# ============================================================
# STEP 6: BENCHMARK COMPARISON (30 MIN)
# ============================================================

print("\n" + "="*80)
print("[STEP 6] BENCHMARK COMPARISON")
print("="*80)

# PD-L1 alone
if 'PDL1_EXPRESSION' in results_df['pathway'].values:
    auc_pdl1 = results_df[results_df['pathway'] == 'PDL1_EXPRESSION']['auc'].values[0]
else:
    auc_pdl1 = np.nan

# Best single pathway
if len(results_df) > 0:
    best_pathway = results_df.iloc[0]['pathway']
    auc_best_single = results_df.iloc[0]['auc']
else:
    best_pathway = 'N/A'
    auc_best_single = np.nan

# Cytolytic score (if available)
auc_cyto = np.nan
if 'cytolytic_score' in cyto.columns:
    # Try to match cytolytic scores to samples
    # This would require matching sample IDs between cyto and analysis_df
    # For now, skip if matching is complex
    pass

# Summary
summary_data = {
    'Method': [],
    'AUC': [],
    'Improvement vs PD-L1': []
}

if not np.isnan(auc_pdl1):
    summary_data['Method'].append('PD-L1 Expression (CD274)')
    summary_data['AUC'].append(f'{auc_pdl1:.3f}')
    summary_data['Improvement vs PD-L1'].append('—')

if not np.isnan(auc_best_single):
    summary_data['Method'].append(f'Best Single Pathway ({best_pathway})')
    summary_data['AUC'].append(f'{auc_best_single:.3f}')
    if not np.isnan(auc_pdl1):
        summary_data['Improvement vs PD-L1'].append(f'+{(auc_best_single - auc_pdl1):.3f}')
    else:
        summary_data['Improvement vs PD-L1'].append('N/A')

if not np.isnan(auc_composite):
    summary_data['Method'].append('Weighted Composite (8 pathways)')
    summary_data['AUC'].append(f'{auc_composite:.3f}')
    if not np.isnan(auc_pdl1):
        summary_data['Improvement vs PD-L1'].append(f'+{(auc_composite - auc_pdl1):.3f}')
    else:
        summary_data['Improvement vs PD-L1'].append('N/A')

if not np.isnan(auc_lr):
    summary_data['Method'].append('Logistic Regression Composite')
    summary_data['AUC'].append(f'{auc_lr:.3f} (CV: {cv_scores.mean():.3f})')
    if not np.isnan(auc_pdl1):
        summary_data['Improvement vs PD-L1'].append(f'+{(auc_lr - auc_pdl1):.3f}')
    else:
        summary_data['Improvement vs PD-L1'].append('N/A')

summary = pd.DataFrame(summary_data)

print("\n" + summary.to_string(index=False))
print("="*80)

# TARGET CHECK
if not np.isnan(auc_lr):
    if auc_lr > 0.75:
        print("\n🎯 ✅ TARGET ACHIEVED: AUC > 0.75")
    elif auc_lr > 0.70:
        print("\n⚠️  CLOSE TO TARGET: AUC > 0.70 (refinement needed)")
    else:
        print("\n❌ BELOW TARGET: AUC < 0.70 (needs improvement)")
elif not np.isnan(auc_composite):
    if auc_composite > 0.75:
        print("\n🎯 ✅ TARGET ACHIEVED: AUC > 0.75")
    elif auc_composite > 0.70:
        print("\n⚠️  CLOSE TO TARGET: AUC > 0.70 (refinement needed)")
    else:
        print("\n❌ BELOW TARGET: AUC < 0.70 (needs improvement)")

# Save summary
summary.to_csv('gse91061_benchmark_comparison.csv', index=False)

# ============================================================
# STEP 7: VISUALIZATION (30 MIN)
# ============================================================

print("\n" + "="*80)
print("[STEP 7] VISUALIZATION")
print("="*80)

import matplotlib.pyplot as plt
import seaborn as sns

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

# 1. ROC Curves
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Left: Top 3 pathways + composite
ax1 = axes[0]
top_pathways = results_df.head(3)['pathway'].values

for pathway in top_pathways:
    if pathway in analysis_df.columns:
        scores = analysis_df[pathway].values
        labels = analysis_df['response'].values
        valid = ~np.isnan(scores)
        if valid.sum() > 0:
            fpr, tpr, _ = roc_curve(labels[valid], scores[valid])
            auc_val = results_df[results_df['pathway'] == pathway]['auc'].values[0]
            ax1.plot(fpr, tpr, label=f'{pathway} (AUC={auc_val:.3f})', linewidth=2)

# Add composite
if 'composite_lr' in analysis_df.columns:
    scores = analysis_df['composite_lr'].values
    labels = analysis_df['response'].values
    valid = ~np.isnan(scores)
    if valid.sum() > 0:
        fpr, tpr, _ = roc_curve(labels[valid], scores[valid])
        ax1.plot(fpr, tpr, label=f'Logistic Composite (AUC={auc_lr:.3f})', 
                linewidth=3, linestyle='--', color='red')

# Add PD-L1 benchmark
if 'PDL1_EXPRESSION' in analysis_df.columns:
    scores = analysis_df['PDL1_EXPRESSION'].values
    labels = analysis_df['response'].values
    valid = ~np.isnan(scores)
    if valid.sum() > 0:
        fpr, tpr, _ = roc_curve(labels[valid], scores[valid])
        ax1.plot(fpr, tpr, label=f'PD-L1 (AUC={auc_pdl1:.3f})', 
                linewidth=2, linestyle=':', color='gray')

ax1.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Random (AUC=0.50)')
ax1.set_xlabel('False Positive Rate', fontsize=12)
ax1.set_ylabel('True Positive Rate', fontsize=12)
ax1.set_title('ROC Curves: Top Pathways vs Composite', fontsize=14, fontweight='bold')
ax1.legend(loc='lower right', fontsize=10)
ax1.grid(alpha=0.3)

# Right: Boxplots for top pathway
ax2 = axes[1]
top_pathway = top_pathways[0] if len(top_pathways) > 0 else None

if top_pathway and top_pathway in analysis_df.columns:
    plot_data = []
    for resp_val in [0, 1]:
        resp_label = 'Non-Responder' if resp_val == 0 else 'Responder'
        scores = analysis_df[analysis_df['response'] == resp_val][top_pathway].values
        scores = scores[~np.isnan(scores)]
        plot_data.extend([{'Response': resp_label, 'Score': s} for s in scores])
    
    if plot_data:
        plot_df = pd.DataFrame(plot_data)
        sns.boxplot(data=plot_df, x='Response', y='Score', ax=ax2, palette=['#e74c3c', '#2ecc71'])
        sns.stripplot(data=plot_df, x='Response', y='Score', ax=ax2, 
                     color='black', alpha=0.5, size=4)
        ax2.set_title(f'{top_pathway} by Response Status', fontsize=14, fontweight='bold')
        ax2.set_ylabel('Pathway Score (log2 TPM+1)', fontsize=12)
        ax2.set_xlabel('', fontsize=12)
        
        # Add p-value
        resp_scores = analysis_df[analysis_df['response'] == 1][top_pathway].dropna()
        nonresp_scores = analysis_df[analysis_df['response'] == 0][top_pathway].dropna()
        if len(resp_scores) > 0 and len(nonresp_scores) > 0:
            _, p_val = mannwhitneyu(resp_scores, nonresp_scores, alternative='greater')
            ax2.text(0.5, 0.95, f'p = {p_val:.4f}', transform=ax2.transAxes,
                    ha='center', va='top', fontsize=11, 
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('gse91061_roc_and_boxplots.png', dpi=300, bbox_inches='tight')
print("✓ Saved: gse91061_roc_and_boxplots.png")
plt.close()

# 2. Feature importance bar plot
if 'feature_importance' in locals() and len(feature_importance) > 0:
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = ['green' if c > 0 else 'red' for c in feature_importance['coefficient']]
    ax.barh(feature_importance['pathway'], feature_importance['coefficient'], color=colors, alpha=0.7)
    ax.axvline(0, color='black', linestyle='-', linewidth=0.5)
    ax.set_xlabel('Logistic Regression Coefficient', fontsize=12)
    ax.set_title('Pathway Feature Importance (Logistic Regression)', fontsize=14, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    plt.tight_layout()
    plt.savefig('gse91061_feature_importance.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: gse91061_feature_importance.png")
    plt.close()

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)
print("\nGenerated files:")
print("  - gse91061_io_pathway_scores.csv")
print("  - gse91061_pathway_response_association.csv")
print("  - gse91061_analysis_with_composites.csv")
print("  - gse91061_benchmark_comparison.csv")
print("  - gse91061_roc_and_boxplots.png")
print("  - gse91061_feature_importance.png")
