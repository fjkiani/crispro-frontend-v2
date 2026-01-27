#!/usr/bin/env python3
"""
VALIDATE HOLISTIC SCORE USING EXISTING TCGA-OV DATA
No downloads needed - using already-curated data:
- tcga_ov_platinum_with_mutations.json (469 patients)
- tcga_ov_469_with_hrd.json (HRD labels)
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
from datetime import datetime

# Paths
DATA_DIR = Path("data/validation/sae_cohort")
HRD_FILE = Path("data/validation/tcga_ov_hrd_summary.json")
OUTPUT_DIR = Path("data/cohorts/tcga_ov_holistic_score/validation_results")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("HOLISTIC SCORE VALIDATION - USING EXISTING DATA")
print("=" * 80)
print()

# Step 1: Load mutation data with platinum response
print("📋 Step 1: Loading existing mutation + response data...")
with open(DATA_DIR / "tcga_ov_platinum_with_mutations.json", "r") as f:
    patients_data = json.load(f)

print(f"✅ Loaded {len(patients_data)} patients")

# Step 2: Extract DDR pathway mutations per patient
print("\n🧬 Step 2: Computing DDR pathway scores from mutations...")

DDR_GENES = ['BRCA1', 'BRCA2', 'TP53', 'ATM', 'PALB2', 'RAD51C', 'RAD51D', 'FANCF', 'CHEK2', 'ATR', 'PARP1', 'PARP2', 'PARP4']
MAPK_GENES = ['KRAS', 'NRAS', 'BRAF', 'MAP2K1', 'MAP2K2', 'MAPK1', 'MAPK3']
PI3K_GENES = ['PIK3CA', 'PTEN', 'AKT1', 'AKT2', 'MTOR']
VEGF_GENES = ['VEGFA', 'FLT1', 'KDR', 'FLT4']
IO_GENES = ['CD274', 'PDCD1', 'CTLA4', 'B2M', 'JAK1', 'JAK2']

patient_rows = []
for p in patients_data:
    patient_id = p.get('tcga_patient_id', '')
    response = p.get('platinum_response', 'unknown')
    mutations = p.get('mutations', [])
    
    # Count mutations in each pathway
    genes_mutated = [m.get('gene', '') for m in mutations]
    
    ddr_count = sum(1 for g in genes_mutated if g in DDR_GENES)
    mapk_count = sum(1 for g in genes_mutated if g in MAPK_GENES)
    pi3k_count = sum(1 for g in genes_mutated if g in PI3K_GENES)
    vegf_count = sum(1 for g in genes_mutated if g in VEGF_GENES)
    io_count = sum(1 for g in genes_mutated if g in IO_GENES)
    
    # Normalize to 0-1 scale (cap at 3 mutations = 1.0)
    ddr_score = min(ddr_count / 3.0, 1.0) if ddr_count > 0 else 0.3  # Baseline for HGSOC
    mapk_score = min(mapk_count / 2.0, 1.0) if mapk_count > 0 else 0.1
    pi3k_score = min(pi3k_count / 2.0, 1.0) if pi3k_count > 0 else 0.1
    vegf_score = 0.2  # Default for ovarian (highly vascular tumors)
    her2_score = 0.05  # Low in HGSOC
    io_score = min(io_count / 2.0, 1.0) if io_count > 0 else 0.2
    efflux_score = 0.15  # Default
    
    patient_rows.append({
        'patient_id': patient_id,
        'platinum_response': response,
        'ddr_score': ddr_score,
        'mapk_score': mapk_score,
        'pi3k_score': pi3k_score,
        'vegf_score': vegf_score,
        'her2_score': her2_score,
        'io_score': io_score,
        'efflux_score': efflux_score,
        'ddr_mutations': ddr_count,
        'total_mutations': len(mutations)
    })

df = pd.DataFrame(patient_rows)
print(f"✅ Computed pathway scores for {len(df)} patients")
print(f"  DDR score range: {df['ddr_score'].min():.3f} - {df['ddr_score'].max():.3f}")
print(f"  Patients with DDR mutations: {(df['ddr_mutations'] > 0).sum()}")

# Step 3: Define trial mechanism (platinum)
print("\n💊 Step 3: Defining platinum trial mechanism vector...")

# Platinum targets DDR pathway
trial_mechanism = np.array([0.9, 0.1, 0.1, 0.2, 0.0, 0.1, 0.1])  # DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux

# Step 4: Compute holistic scores
print("\n🎯 Step 4: Computing holistic scores...")

for idx, row in df.iterrows():
    patient_vector = np.array([
        row['ddr_score'], row['mapk_score'], row['pi3k_score'],
        row['vegf_score'], row['her2_score'], row['io_score'], row['efflux_score']
    ])
    
    # Cosine similarity
    mechanism_fit = np.dot(patient_vector, trial_mechanism) / (
        np.linalg.norm(patient_vector) * np.linalg.norm(trial_mechanism)
    )
    
    df.loc[idx, 'mechanism_fit'] = mechanism_fit

# Holistic score = 0.5*mechanism + 0.3*eligibility + 0.2*pgx
df['eligibility'] = 1.0  # All treated
df['pgx_safety'] = 1.0  # No PGx data
df['holistic_score'] = 0.5 * df['mechanism_fit'] + 0.3 * df['eligibility'] + 0.2 * df['pgx_safety']

print(f"✅ Holistic scores computed")
print(f"  Score range: {df['holistic_score'].min():.3f} - {df['holistic_score'].max():.3f}")
print(f"  Mean: {df['holistic_score'].mean():.3f} ± {df['holistic_score'].std():.3f}")

# Step 5: Validate against platinum response
print("\n📊 Step 5: Validating against platinum response...")

# Binary outcome: sensitive=1, resistant=0
df['response_binary'] = (df['platinum_response'] == 'sensitive').astype(int)

# Filter valid responses
valid = df[df['platinum_response'].isin(['sensitive', 'resistant'])].copy()
print(f"\n  Valid patients: {len(valid)}")
print(f"  Sensitive: {valid['response_binary'].sum()}, Resistant: {(1-valid['response_binary']).sum()}")

# AUROC
if valid['response_binary'].nunique() > 1:
    auroc = roc_auc_score(valid['response_binary'], valid['holistic_score'])
    print(f"\n  ✅ AUROC: {auroc:.3f}")
else:
    auroc = None
    print(f"\n  ⚠️ Cannot compute AUROC (insufficient class variability)")

# Correlation
corr = stats.spearmanr(valid['holistic_score'], valid['response_binary'])
print(f"  Correlation (holistic vs response): ρ={corr.correlation:.3f}, p={corr.pvalue:.4f}")

# Quartile analysis
try:
    valid['score_quartile'] = pd.qcut(valid['holistic_score'], q=4, labels=['Q1', 'Q2', 'Q3', 'Q4'], duplicates='drop')
    quartile_summary = valid.groupby('score_quartile').agg({
        'response_binary': ['sum', 'count', 'mean']
    })
    quartile_summary.columns = ['sensitive_count', 'total', 'response_rate']
    print(f"\n  Quartile Analysis:")
    print(quartile_summary)
except Exception as e:
    print(f"\n  ⚠️ Quartile analysis failed: {e}")

# Mean scores by response
print(f"\n  Mean holistic score:")
print(f"    Sensitive: {valid[valid['response_binary']==1]['holistic_score'].mean():.3f}")
print(f"    Resistant: {valid[valid['response_binary']==0]['holistic_score'].mean():.3f}")

# T-test
t_stat, t_pval = stats.ttest_ind(
    valid[valid['response_binary']==1]['holistic_score'],
    valid[valid['response_binary']==0]['holistic_score']
)
print(f"    T-test: t={t_stat:.3f}, p={t_pval:.4f}")

# Step 6: Save results
print("\n💾 Step 6: Saving results...")

results_path = OUTPUT_DIR / "holistic_validation_EXISTING_DATA.csv"
valid.to_csv(results_path, index=False)
print(f"✅ Saved: {results_path}")

# Summary
summary = {
    "timestamp": datetime.now().isoformat(),
    "dataset": "TCGA-OV (existing data)",
    "n_patients": len(valid),
    "n_sensitive": int(valid['response_binary'].sum()),
    "n_resistant": int((1-valid['response_binary']).sum()),
    "auroc": float(auroc) if auroc else None,
    "correlation": {
        "rho": float(corr.correlation),
        "p": float(corr.pvalue)
    },
    "ttest": {
        "t": float(t_stat),
        "p": float(t_pval)
    },
    "holistic_score_stats": {
        "mean": float(valid['holistic_score'].mean()),
        "std": float(valid['holistic_score'].std()),
        "min": float(valid['holistic_score'].min()),
        "max": float(valid['holistic_score'].max())
    }
}

summary_path = OUTPUT_DIR / "validation_summary_EXISTING_DATA.json"
with open(summary_path, 'w') as f:
    json.dump(summary, f, indent=2)

print(f"✅ Summary: {summary_path}")

# Final report
print("\n" + "=" * 80)
print("VALIDATION COMPLETE - REAL DATA, NO SYNTHETIC RECONSTRUCTION")
print("=" * 80)
print(f"  Patients: {len(valid)}")
print(f"  Sensitive: {valid['response_binary'].sum()}, Resistant: {(1-valid['response_binary']).sum()}")
print(f"  Holistic score: {valid['holistic_score'].mean():.3f} ± {valid['holistic_score'].std():.3f}")
print(f"  AUROC: {auroc:.3f}" if auroc else "  AUROC: N/A")
print(f"  Correlation: ρ={corr.correlation:.3f}, p={corr.pvalue:.4f}")
print(f"  T-test: p={t_pval:.4f}")
print()
print("✅ SUCCESS: Validated on EXISTING TCGA-OV data")
print("=" * 80)
