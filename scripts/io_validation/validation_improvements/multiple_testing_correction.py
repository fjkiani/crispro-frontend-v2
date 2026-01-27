#!/usr/bin/env python3
"""
Task 1.3A: Multiple Testing Correction
=======================================

Applies Benjamini-Hochberg FDR correction and Bonferroni correction
to 8 pathway tests to control for multiple testing.

Expected false positives: 8 tests × 0.05 = 0.4 false discoveries
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Try to import statsmodels, fallback to manual implementation
try:
    from statsmodels.stats.multitest import multipletests
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False
    print("⚠️  statsmodels not available, using manual FDR/Bonferroni correction")

# Configuration
SCRIPT_DIR = Path(__file__).parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
DATA_DIR = REPO_ROOT / "scripts" / "data_acquisition" / "IO"
OUTPUT_DIR = REPO_ROOT / "publications" / "06-io-response-prediction" / "data"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("="*80)
print("TASK 1.3A: MULTIPLE TESTING CORRECTION")
print("="*80)
print()

# Load pathway statistics
print("Loading pathway statistics...")
pathway_stats = pd.read_csv(DATA_DIR / "gse91061_pathway_response_association.csv")

# Filter to pathways (exclude PDL1_EXPRESSION for now, we'll add it back)
pathway_stats_filtered = pathway_stats[
    pathway_stats['pathway'] != 'PDL1_EXPRESSION'
].copy()

print(f"Testing {len(pathway_stats_filtered)} pathways")
print()

# Extract p-values
p_values = pathway_stats_filtered['p_value'].values
pathways = pathway_stats_filtered['pathway'].values

print("Original p-values:")
for pathway, p_val in zip(pathways, p_values):
    sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
    print(f"  {pathway:25s}: {p_val:.4f} {sig}")
print()

# Benjamini-Hochberg FDR correction
if HAS_STATSMODELS:
    rejected_fdr, p_adjusted_fdr, _, _ = multipletests(
        p_values, 
        alpha=0.05, 
        method='fdr_bh'
    )
    
    # Bonferroni correction
    rejected_bonf, p_adjusted_bonf, _, _ = multipletests(
        p_values, 
        alpha=0.05, 
        method='bonferroni'
    )
else:
    # Manual FDR correction (Benjamini-Hochberg)
    n = len(p_values)
    sorted_indices = np.argsort(p_values)
    sorted_p = p_values[sorted_indices]
    
    p_adjusted_fdr = np.zeros(n)
    for i in range(n):
        rank = i + 1
        p_adjusted_fdr[sorted_indices[i]] = sorted_p[i] * n / rank
    
    # Cap at 1.0
    p_adjusted_fdr = np.minimum(p_adjusted_fdr, 1.0)
    rejected_fdr = p_adjusted_fdr < 0.05
    
    # Manual Bonferroni correction
    p_adjusted_bonf = p_values * n
    p_adjusted_bonf = np.minimum(p_adjusted_bonf, 1.0)
    rejected_bonf = p_adjusted_bonf < 0.05

# Bonferroni threshold
bonferroni_threshold = 0.05 / len(p_values)

# Add corrected p-values to dataframe
pathway_stats_filtered['p_value_original'] = p_values
pathway_stats_filtered['fdr_q_value'] = p_adjusted_fdr
pathway_stats_filtered['bonferroni_p_value'] = p_adjusted_bonf
pathway_stats_filtered['fdr_significant'] = rejected_fdr
pathway_stats_filtered['bonferroni_significant'] = rejected_bonf

# Sort by original p-value
pathway_stats_filtered = pathway_stats_filtered.sort_values('p_value_original')

print(f"{'='*80}")
print("MULTIPLE TESTING CORRECTION RESULTS")
print(f"{'='*80}")
print()
print(f"Bonferroni threshold: p < {bonferroni_threshold:.4f} (0.05 / {len(p_values)})")
print(f"FDR threshold: q < 0.05")
print()

print("Pathway | Original p | FDR q-value | Bonferroni p | FDR sig | Bonf sig")
print("-" * 80)

for _, row in pathway_stats_filtered.iterrows():
    pathway = row['pathway']
    p_orig = row['p_value_original']
    fdr_q = row['fdr_q_value']
    bonf_p = row['bonferroni_p_value']
    fdr_sig = "✓" if row['fdr_significant'] else "✗"
    bonf_sig = "✓" if row['bonferroni_significant'] else "✗"
    
    print(f"{pathway:25s} | {p_orig:8.4f} | {fdr_q:10.4f} | {bonf_p:12.4f} | {fdr_sig:7s} | {bonf_sig:8s}")

print()

# Summary
n_fdr_sig = pathway_stats_filtered['fdr_significant'].sum()
n_bonf_sig = pathway_stats_filtered['bonferroni_significant'].sum()

print(f"Summary:")
print(f"  Pathways significant after FDR correction: {n_fdr_sig}/{len(pathway_stats_filtered)}")
print(f"  Pathways significant after Bonferroni correction: {n_bonf_sig}/{len(pathway_stats_filtered)}")
print()

# Expected false positives
expected_fp = len(p_values) * 0.05
print(f"Expected false positives (uncorrected): {expected_fp:.1f}")
print(f"Actual significant (FDR): {n_fdr_sig}")
print(f"Actual significant (Bonferroni): {n_bonf_sig}")
print()

# Save results
pathway_stats_filtered.to_csv(OUTPUT_DIR / "pathway_multiple_testing_corrected.csv", index=False)
print(f"✅ Saved results: {OUTPUT_DIR / 'pathway_multiple_testing_corrected.csv'}")

# Create publication-ready table
table_publication = pathway_stats_filtered[[
    'pathway', 'auc', 'p_value_original', 'fdr_q_value', 'bonferroni_p_value',
    'cohens_d', 'resp_mean', 'nonresp_mean'
]].copy()

table_publication = table_publication.rename(columns={
    'pathway': 'Pathway',
    'auc': 'AUC',
    'p_value_original': 'p-value',
    'fdr_q_value': 'FDR q-value',
    'bonferroni_p_value': 'Bonferroni p-value',
    'cohens_d': "Cohen's d",
    'resp_mean': 'Responder Mean',
    'nonresp_mean': 'Non-Responder Mean'
})

# Format p-values
table_publication['p-value'] = table_publication['p-value'].apply(
    lambda x: f"{x:.4f}" if x >= 0.0001 else f"{x:.2e}"
)
table_publication['FDR q-value'] = table_publication['FDR q-value'].apply(
    lambda x: f"{x:.4f}" if x >= 0.0001 else f"{x:.2e}"
)
table_publication['Bonferroni p-value'] = table_publication['Bonferroni p-value'].apply(
    lambda x: f"{x:.4f}" if x >= 0.0001 else f"{x:.2e}"
)

table_publication.to_csv(OUTPUT_DIR / "table1_pathway_performance_corrected.csv", index=False)
print(f"✅ Saved publication table: {OUTPUT_DIR / 'table1_pathway_performance_corrected.csv'}")

print(f"\n{'='*80}")
print("TASK 1.3A COMPLETE")
print(f"{'='*80}")
