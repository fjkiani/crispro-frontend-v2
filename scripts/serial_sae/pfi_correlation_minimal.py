#!/usr/bin/env python3
"""
PFI Correlation Analysis - Minimal Version (Statistics Only)
Avoids matplotlib blocking issues
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu
from datetime import datetime
import sys

RESULTS_DIR = Path("data/serial_sae/gse165897/results")

# PFI data (COMPLETE - all 11 GSE165897 patients)
# Source: Table 1 from Science Advances supplementary materials
# PFI in days, converted to months (days / 30.44)
# Resistance: R = PFI < 180 days (< 6 months), S = PFI ≥ 180 days (≥ 6 months)
PFI_DATA = {
    'patient_id': ['EOC1005', 'EOC136', 'EOC153', 'EOC227', 'EOC3', 
                   'EOC349', 'EOC372', 'EOC443', 'EOC540', 'EOC733', 'EOC87'],
    'pfi_days': [65, 520, 393, 230, 14, 36, 460, 177, 126, 83, 30],
    'pfi_months': [2.14, 17.08, 12.91, 7.56, 0.46, 1.18, 15.11, 5.81, 4.14, 2.73, 0.99],
    'resistant': ['R', 'S', 'S', 'R', 'R', 'R', 'S', 'R', 'R', 'R', 'R'],
    'stage': ['IVA', 'IVA', 'IVA', 'IVA', 'IVA', 'IVB', 'IIIC', 'IVA', 'IIIC', 'IVA', 'IIIC'],
    'age': [73, 64, 78, 74, 67, 67, 68, 54, 62, 72, 62]
}

def cohens_d(group1, group2):
    """Calculate Cohen's d effect size."""
    n1, n2 = len(group1), len(group2)
    if n1 == 0 or n2 == 0:
        return 0.0
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    if pooled_std == 0:
        return 0.0
    return (np.mean(group1) - np.mean(group2)) / pooled_std

def main():
    print("="*80)
    print("PFI CORRELATION ANALYSIS (MINIMAL - STATISTICS ONLY)")
    print("="*80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    sys.stdout.flush()
    
    # Load data
    kinetics_df = pd.read_csv(RESULTS_DIR / "pathway_kinetics.csv")
    print(f"✅ Loaded {len(kinetics_df)} patients\n")
    sys.stdout.flush()
    
    # Merge PFI
    pfi_df = pd.DataFrame(PFI_DATA)
    merged_df = kinetics_df.merge(pfi_df, on='patient_id', how='left')
    evaluable_df = merged_df.dropna(subset=['pfi_months']).copy()
    
    print(f"✅ Evaluable patients: n={len(evaluable_df)}")
    print(f"   Resistant: {len(evaluable_df[evaluable_df['resistant']=='R'])}")
    print(f"   Sensitive: {len(evaluable_df[evaluable_df['resistant']=='S'])}\n")
    
    if len(evaluable_df) < 5:
        print("⚠️  WARNING: Small sample size (n<5) limits statistical power")
        print("   Results should be interpreted as exploratory/pilot findings\n")
    sys.stdout.flush()
    
    # Correlations
    print("="*80)
    print("SPEARMAN CORRELATIONS")
    print("="*80)
    sys.stdout.flush()
    
    pathway_deltas = ['ddr_delta', 'mapk_delta', 'pi3k_delta', 'vegf_delta']
    correlation_results = []
    
    for pathway in pathway_deltas:
        valid_data = evaluable_df[[pathway, 'pfi_months']].dropna()
        if len(valid_data) >= 3:
            r, p = spearmanr(valid_data['pfi_months'], valid_data[pathway])
            correlation_results.append({
                'pathway': pathway.replace('_delta', ''),
                'spearman_r': r,
                'spearman_p': p,
                'n': len(valid_data)
            })
            print(f"{pathway.upper()}: r={r:.3f}, p={p:.4f} (n={len(valid_data)})")
            if p < 0.05:
                print(f"  ✅ Significant (p < 0.05)")
            elif p < 0.1:
                print(f"  ⚠️  Trend (p < 0.1)")
            sys.stdout.flush()
    
    correlation_df = pd.DataFrame(correlation_results)
    correlation_df.to_csv(RESULTS_DIR / "pfi_correlations.csv", index=False)
    print(f"\n💾 Saved: {RESULTS_DIR / 'pfi_correlations.csv'}\n")
    sys.stdout.flush()
    
    # Group comparison
    print("="*80)
    print("MANN-WHITNEY U TEST")
    print("="*80)
    sys.stdout.flush()
    
    resistant_df = evaluable_df[evaluable_df['resistant'] == 'R']
    sensitive_df = evaluable_df[evaluable_df['resistant'] == 'S']
    
    comparison_results = []
    
    for pathway in pathway_deltas:
        resistant_vals = resistant_df[pathway].dropna().values
        sensitive_vals = sensitive_df[pathway].dropna().values
        
        if len(resistant_vals) > 0 and len(sensitive_vals) > 0:
            u_stat, u_p = mannwhitneyu(resistant_vals, sensitive_vals, alternative='two-sided')
            d = cohens_d(resistant_vals, sensitive_vals)
            
            comparison_results.append({
                'pathway': pathway.replace('_delta', ''),
                'resistant_mean': np.mean(resistant_vals),
                'resistant_n': len(resistant_vals),
                'sensitive_mean': np.mean(sensitive_vals),
                'sensitive_n': len(sensitive_vals),
                'mannwhitney_p': u_p,
                'cohens_d': d
            })
            
            print(f"{pathway.upper()}:")
            print(f"  Resistant: {np.mean(resistant_vals):.4f} (n={len(resistant_vals)})")
            print(f"  Sensitive: {np.mean(sensitive_vals):.4f} (n={len(sensitive_vals)})")
            print(f"  Mann-Whitney U: p={u_p:.4f}, d={d:.3f}")
            sys.stdout.flush()
    
    comparison_df = pd.DataFrame(comparison_results)
    comparison_df.to_csv(RESULTS_DIR / "pfi_group_comparison.csv", index=False)
    print(f"\n💾 Saved: {RESULTS_DIR / 'pfi_group_comparison.csv'}\n")
    sys.stdout.flush()
    
    # Key findings
    print("="*80)
    print("KEY FINDINGS")
    print("="*80)
    sys.stdout.flush()
    
    vegf_row = correlation_df[correlation_df['pathway'] == 'vegf']
    if not vegf_row.empty:
        r = vegf_row.iloc[0]['spearman_r']
        p = vegf_row.iloc[0]['spearman_p']
        print(f"\n🎯 PRIMARY HYPOTHESIS (VEGF Activation → Short PFI):")
        print(f"   r = {r:.3f}, p = {p:.4f}")
        if r < -0.4 and p < 0.1:
            print(f"   ✅ STRONG SIGNAL: Hypothesis supported!")
        elif r < 0:
            print(f"   ⚠️  Direction correct but needs larger cohort")
        sys.stdout.flush()
    
    ddr_row = correlation_df[correlation_df['pathway'] == 'ddr']
    if not ddr_row.empty:
        r = ddr_row.iloc[0]['spearman_r']
        p = ddr_row.iloc[0]['spearman_p']
        print(f"\n🎯 SECONDARY HYPOTHESIS (DDR Suppression → Short PFI):")
        print(f"   r = {r:.3f}, p = {p:.4f}")
        if r > 0.4 and p < 0.1:
            print(f"   ✅ STRONG SIGNAL: Hypothesis supported!")
        elif r > 0:
            print(f"   ⚠️  Direction correct but needs larger cohort")
        sys.stdout.flush()
    
    print("\n✅ Analysis complete!\n")
    sys.stdout.flush()

if __name__ == "__main__":
    main()
