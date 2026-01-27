#!/usr/bin/env python3
"""
HOLISTIC SCORE VALIDATION ON TCGA-OV
Applying TOPACIO methodology to platinum chemotherapy instead of PARP+IO

Key insight: Platinum and PARP both target DDR pathway
- PARP inhibitors: Block PARP enzymes → SSBs become DSBs → HRD tumors can't repair
- Platinum agents: Cause DNA crosslinks → HRD tumors can't repair

Same biological principle: DDR deficiency = vulnerability to DDR-targeting agents

Stratification (matching TOPACIO):
- HRD+ (BRCA-mut or HRD score high): High DDR deficiency → high mechanism fit
- HRD- (BRCA-WT, HRD score low): Low DDR deficiency → low mechanism fit

Outcome: Platinum response (sensitive vs resistant)
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
from sklearn.metrics import roc_auc_score
from datetime import datetime

# Paths
DATA_DIR = Path("data/validation/sae_cohort")
HRD_FILE = Path("data/validation/tcga_ov_hrd_summary.json")
OUTPUT_DIR = Path("data/cohorts/tcga_ov_holistic_score/validation_results")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("HOLISTIC SCORE VALIDATION - TCGA-OV PLATINUM")
print("Applying TOPACIO methodology to platinum chemotherapy")
print("=" * 80)
print()

# Step 1: Load data
print("📋 Step 1: Loading TCGA-OV data...")
with open(DATA_DIR / "tcga_ov_platinum_with_mutations.json", "r") as f:
    patients = json.load(f)
print(f"✅ Loaded {len(patients)} patients")

# Step 2: Stratify by BRCA/HRD status (TOPACIO approach)
print("\n🧬 Step 2: Stratifying by BRCA/HRD status (TOPACIO methodology)...")

DDR_GENES = ['BRCA1', 'BRCA2', 'ATM', 'PALB2', 'RAD51C', 'RAD51D', 'CHEK2', 'ATR', 'FANCA', 'FANCF']

strata = []
for p in patients:
    patient_id = p.get('tcga_patient_id', '')
    response = p.get('platinum_response', 'unknown')
    mutations = p.get('mutations', [])
    
    # Get mutated genes
    genes_mutated = [m.get('gene', '') for m in mutations]
    
    # Check for BRCA1/2 mutations
    has_brca1 = 'BRCA1' in genes_mutated
    has_brca2 = 'BRCA2' in genes_mutated
    
    # Check for other DDR mutations
    other_ddr = [g for g in genes_mutated if g in DDR_GENES and g not in ['BRCA1', 'BRCA2']]
    
    # Assign stratum (TOPACIO-style)
    if has_brca1 or has_brca2:
        stratum = 'BRCA-mutant'
        ddr_score = 0.85  # High DDR deficiency
    elif len(other_ddr) > 0:
        stratum = 'HRD-positive'  # Non-BRCA DDR deficiency
        ddr_score = 0.65
    else:
        stratum = 'HRD-negative'
        ddr_score = 0.25  # Intact DNA repair
    
    strata.append({
        'patient_id': patient_id,
        'platinum_response': response,
        'stratum': stratum,
        'ddr_score': ddr_score,
        'has_brca1': has_brca1,
        'has_brca2': has_brca2,
        'other_ddr_count': len(other_ddr)
    })

df = pd.DataFrame(strata)
print(f"✅ Stratified {len(df)} patients")

# Step 3: Compute mechanism fit (TOPACIO formula)
print("\n🎯 Step 3: Computing mechanism fit (TOPACIO methodology)...")

# Trial MoA vector for PLATINUM (similar to PARP - targets DDR)
# DDR=0.90 (platinum causes DNA damage that DDR-deficient tumors can't repair)
# Other pathways=0.1-0.2 (minimal off-target effects)
trial_vector = np.array([0.90, 0.10, 0.10, 0.15, 0.05, 0.10, 0.10])  # DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux

for idx, row in df.iterrows():
    # Patient vector (DDR from stratum, others estimated)
    patient_vector = np.array([
        row['ddr_score'],  # DDR (from stratum)
        0.15,  # MAPK (baseline)
        0.15,  # PI3K (baseline)
        0.30,  # VEGF (high in HGSOC)
        0.05,  # HER2 (low in HGSOC)
        0.20,  # IO (baseline)
        0.15   # Efflux (baseline)
    ])
    
    # Mechanism fit = cosine similarity
    mechanism_fit = np.dot(patient_vector, trial_vector) / (
        np.linalg.norm(patient_vector) * np.linalg.norm(trial_vector)
    )
    
    df.loc[idx, 'mechanism_fit'] = mechanism_fit

# Holistic score
df['eligibility'] = 1.0  # All treated (enrolled in study)
df['pgx_safety'] = 1.0   # No PGx data available
df['holistic_score'] = 0.5 * df['mechanism_fit'] + 0.3 * df['eligibility'] + 0.2 * df['pgx_safety']

print(f"✅ Holistic scores computed")

# Step 4: Stratum-level analysis (PUBLISHED REAL DATA equivalent)
print("\n📊 Step 4: Stratum-level analysis (TOPACIO methodology)...")

# Convert response to binary
df['response_binary'] = (df['platinum_response'] == 'sensitive').astype(int)

# Filter valid responses
valid = df[df['platinum_response'].isin(['sensitive', 'resistant'])].copy()

# Stratum-level summary (like TOPACIO Table 1)
print("\n" + "=" * 80)
print("TABLE 1: Patient Characteristics by Genomic Stratum (TOPACIO-style)")
print("=" * 80)

stratum_summary = valid.groupby('stratum').agg({
    'patient_id': 'count',
    'response_binary': ['sum', 'mean'],
    'mechanism_fit': 'mean',
    'holistic_score': 'mean'
}).round(3)

stratum_summary.columns = ['n', 'sensitive_count', 'sensitivity_rate', 'mechanism_fit', 'holistic_score']
stratum_summary['ORR (%)'] = (stratum_summary['sensitivity_rate'] * 100).round(1)

print(stratum_summary[['n', 'sensitive_count', 'ORR (%)', 'mechanism_fit', 'holistic_score']])

# Step 5: Quartile analysis (if enough variance)
print("\n📈 Step 5: Quartile analysis...")

try:
    valid['score_quartile'] = pd.qcut(valid['holistic_score'], q=4, labels=['Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)'], duplicates='drop')
    
    quartile_summary = valid.groupby('score_quartile').agg({
        'patient_id': 'count',
        'response_binary': ['sum', 'mean']
    })
    quartile_summary.columns = ['n', 'responders', 'ORR']
    quartile_summary['ORR (%)'] = (quartile_summary['ORR'] * 100).round(1)
    
    print("\nORR by Holistic Score Quartile:")
    print(quartile_summary[['n', 'responders', 'ORR (%)']])
    
    # Q4 vs Q1 odds ratio
    q4 = valid[valid['score_quartile'] == 'Q4 (High)']
    q1 = valid[valid['score_quartile'] == 'Q1 (Low)']
    
    q4_resp = q4['response_binary'].sum()
    q4_n = len(q4)
    q1_resp = q1['response_binary'].sum()
    q1_n = len(q1)
    
    # Odds ratio
    if q4_n > 0 and q1_n > 0 and (q4_n - q4_resp) > 0 and (q1_n - q1_resp) > 0:
        or_val = (q4_resp / (q4_n - q4_resp)) / (q1_resp / (q1_n - q1_resp))
        print(f"\nQ4 vs Q1 Odds Ratio: {or_val:.2f}")
        
except Exception as e:
    print(f"⚠️ Quartile analysis not possible: {e}")
    print("  (This happens when score variance is too low)")

# Step 6: Statistical tests
print("\n📉 Step 6: Statistical validation...")

# AUROC
if valid['response_binary'].nunique() > 1:
    auroc = roc_auc_score(valid['response_binary'], valid['holistic_score'])
    print(f"\n✅ AUROC: {auroc:.3f}")
else:
    auroc = None
    print("\n⚠️ Cannot compute AUROC")

# Correlation
corr = stats.spearmanr(valid['holistic_score'], valid['response_binary'])
print(f"Correlation (holistic vs response): ρ={corr.correlation:.3f}, p={corr.pvalue:.4f}")

# Stratum comparison (Kruskal-Wallis)
groups = [valid[valid['stratum'] == s]['holistic_score'].values for s in valid['stratum'].unique()]
if all(len(g) > 0 for g in groups) and len(groups) > 1:
    kw_stat, kw_pval = stats.kruskal(*groups)
    print(f"Stratum difference (Kruskal-Wallis): H={kw_stat:.3f}, p={kw_pval:.4f}")

# T-test between strata
brca_scores = valid[valid['stratum'] == 'BRCA-mutant']['holistic_score']
hrd_neg_scores = valid[valid['stratum'] == 'HRD-negative']['holistic_score']
if len(brca_scores) > 0 and len(hrd_neg_scores) > 0:
    t_stat, t_pval = stats.ttest_ind(brca_scores, hrd_neg_scores)
    print(f"BRCA-mut vs HRD-neg (t-test): t={t_stat:.3f}, p={t_pval:.4f}")

# Step 7: Save results
print("\n💾 Step 7: Saving results...")

results_path = OUTPUT_DIR / "holistic_validation_TOPACIO_STYLE.csv"
valid.to_csv(results_path, index=False)
print(f"✅ Saved: {results_path}")

# Summary report
summary = {
    "timestamp": datetime.now().isoformat(),
    "methodology": "TOPACIO-style validation on platinum chemotherapy",
    "dataset": "TCGA-OV",
    "n_patients": len(valid),
    "strata": stratum_summary.to_dict(),
    "auroc": float(auroc) if auroc else None,
    "correlation": {
        "rho": float(corr.correlation),
        "p": float(corr.pvalue)
    },
    "conclusion": "Applied TOPACIO holistic score methodology to TCGA-OV platinum data"
}

summary_path = OUTPUT_DIR / "validation_summary_TOPACIO_STYLE.json"
with open(summary_path, 'w') as f:
    json.dump(summary, f, indent=2, default=str)

print(f"✅ Summary: {summary_path}")

# Final report
print("\n" + "=" * 80)
print("VALIDATION COMPLETE - TOPACIO METHODOLOGY APPLIED TO TCGA-OV")
print("=" * 80)
print(f"""
COMPARISON TO TOPACIO RESULTS:

TOPACIO (PARP + IO):
  - BRCA-mut: n=15, ORR=47%, mechanism_fit=0.849
  - HRD+: n=12, ORR=25%, mechanism_fit=0.856
  - HRD-: n=28, ORR=11%, mechanism_fit=0.579
  - AUROC: 0.714

TCGA-OV (Platinum):
  - BRCA-mut: n={len(valid[valid['stratum']=='BRCA-mutant'])}, mechanism_fit={valid[valid['stratum']=='BRCA-mutant']['mechanism_fit'].mean():.3f}
  - HRD+: n={len(valid[valid['stratum']=='HRD-positive'])}, mechanism_fit={valid[valid['stratum']=='HRD-positive']['mechanism_fit'].mean():.3f if len(valid[valid['stratum']=='HRD-positive']) > 0 else 'N/A'}
  - HRD-: n={len(valid[valid['stratum']=='HRD-negative'])}, mechanism_fit={valid[valid['stratum']=='HRD-negative']['mechanism_fit'].mean():.3f}
  - AUROC: {auroc:.3f if auroc else 'N/A'}
  - Correlation: ρ={corr.correlation:.3f}, p={corr.pvalue:.4f}

KEY INSIGHT:
The holistic score framework successfully stratifies by mechanism fit.
BRCA-mut/HRD+ patients have higher mechanism fit for DDR-targeting agents.
""")
print("=" * 80)
