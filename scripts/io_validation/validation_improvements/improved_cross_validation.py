#!/usr/bin/env python3
"""
Task 1.3C: Improved Cross-Validation
====================================

Compares 5-fold CV, LOOCV, and bootstrap validation for small sample size (n=51).
5-fold CV on n=51 yields test folds with ~5 responders each, which is insufficient.

Methods:
1. Leave-One-Out Cross-Validation (LOOCV) - 51 iterations, 50 train / 1 test
2. Bootstrap Validation - 1000 iterations, resample with replacement
3. 5-fold CV (baseline) - 5 folds, ~40 train / 11 test
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, StratifiedKFold, LeaveOneOut
from sklearn.metrics import roc_auc_score
import json
from tqdm import tqdm

# Configuration
SCRIPT_DIR = Path(__file__).parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
DATA_DIR = REPO_ROOT / "scripts" / "data_acquisition" / "IO"
OUTPUT_DIR = REPO_ROOT / "publications" / "06-io-response-prediction" / "data"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Set random seed for reproducibility
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Bootstrap parameters
N_BOOTSTRAP = 1000

print("="*80)
print("TASK 1.3C: IMPROVED CROSS-VALIDATION")
print("="*80)
print()

# Load data
print("Loading GSE91061 data...")
df = pd.read_csv(DATA_DIR / "gse91061_analysis_with_composites.csv")

# Pathway columns
pathway_cols = ['TIL_INFILTRATION', 'T_EFFECTOR', 'ANGIOGENESIS', 'TGFB_RESISTANCE',
                'MYELOID_INFLAMMATION', 'PROLIFERATION', 'IMMUNOPROTEASOME', 'EXHAUSTION']

# Prepare features and labels
X = df[pathway_cols].values
y = df['response'].values

# Remove rows with NaN
valid_rows = ~np.isnan(X).any(axis=1)
X = X[valid_rows]
y = y[valid_rows]

print(f"Valid samples: {len(X)}")
print(f"  Responders: {y.sum()}")
print(f"  Non-responders: {len(y) - y.sum()}")
print()

# Standardize features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

results = {}

# Method 1: 5-fold CV (baseline)
print("Method 1: 5-fold Cross-Validation...")
cv_5fold = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_SEED)
lr_5fold = LogisticRegression(penalty='l2', C=1.0, random_state=RANDOM_SEED, max_iter=1000)

cv_5fold_scores = cross_val_score(lr_5fold, X_scaled, y, cv=cv_5fold, scoring='roc_auc')

results['5fold_cv'] = {
    'mean': float(cv_5fold_scores.mean()),
    'std': float(cv_5fold_scores.std()),
    'min': float(cv_5fold_scores.min()),
    'max': float(cv_5fold_scores.max()),
    'scores': cv_5fold_scores.tolist(),
    'n_folds': 5,
    'n_test_per_fold': len(X) // 5
}

print(f"  5-fold CV AUC: {cv_5fold_scores.mean():.3f} ± {cv_5fold_scores.std():.3f}")
print(f"  Range: [{cv_5fold_scores.min():.3f}, {cv_5fold_scores.max():.3f}]")
print(f"  Test samples per fold: ~{len(X) // 5}")
print()

# Method 2: Leave-One-Out CV (LOOCV)
print("Method 2: Leave-One-Out Cross-Validation...")
print("  (This may take a few minutes...)")
cv_loocv = LeaveOneOut()
lr_loocv = LogisticRegression(penalty='l2', C=1.0, random_state=RANDOM_SEED, max_iter=1000)

cv_loocv_scores = cross_val_score(lr_loocv, X_scaled, y, cv=cv_loocv, scoring='roc_auc')

results['loocv'] = {
    'mean': float(cv_loocv_scores.mean()),
    'std': float(cv_loocv_scores.std()),
    'min': float(cv_loocv_scores.min()),
    'max': float(cv_loocv_scores.max()),
    'scores': cv_loocv_scores.tolist(),
    'n_folds': len(X),
    'n_test_per_fold': 1
}

print(f"  LOOCV AUC: {cv_loocv_scores.mean():.3f} ± {cv_loocv_scores.std():.3f}")
print(f"  Range: [{cv_loocv_scores.min():.3f}, {cv_loocv_scores.max():.3f}]")
print(f"  Test samples per fold: 1")
print()

# Method 3: Bootstrap Validation
print(f"Method 3: Bootstrap Validation ({N_BOOTSTRAP} iterations)...")
print("  (This may take a few minutes...)")

bootstrap_scores = []

for i in tqdm(range(N_BOOTSTRAP), desc="Bootstrap"):
    # Resample with replacement
    n_samples = len(X)
    bootstrap_indices = np.random.choice(n_samples, size=n_samples, replace=True)
    
    # Get out-of-bag indices
    all_indices = set(range(n_samples))
    bootstrap_set = set(bootstrap_indices)
    oob_indices = np.array(list(all_indices - bootstrap_set))
    
    # Skip if no out-of-bag samples
    if len(oob_indices) == 0:
        continue
    
    # Bootstrap training set
    X_boot = X_scaled[bootstrap_indices]
    y_boot = y[bootstrap_indices]
    
    # Out-of-bag test set
    X_oob = X_scaled[oob_indices]
    y_oob = y[oob_indices]
    
    # Train model
    lr_boot = LogisticRegression(
        penalty='l2', 
        C=1.0, 
        random_state=RANDOM_SEED + i, 
        max_iter=1000
    )
    lr_boot.fit(X_boot, y_boot)
    
    # Evaluate on out-of-bag samples
    y_pred_oob = lr_boot.predict_proba(X_oob)[:, 1]
    oob_auc = roc_auc_score(y_oob, y_pred_oob)
    
    bootstrap_scores.append(oob_auc)

bootstrap_scores = np.array(bootstrap_scores)

results['bootstrap'] = {
    'mean': float(bootstrap_scores.mean()),
    'std': float(bootstrap_scores.std()),
    'min': float(bootstrap_scores.min()),
    'max': float(bootstrap_scores.max()),
    'scores': bootstrap_scores.tolist(),
    'n_iterations': int(len(bootstrap_scores)),
    'n_test_per_iteration': 'variable (OOB)'
}

print(f"  Bootstrap AUC: {bootstrap_scores.mean():.3f} ± {bootstrap_scores.std():.3f}")
print(f"  Range: [{bootstrap_scores.min():.3f}, {bootstrap_scores.max():.3f}]")
print(f"  Iterations: {len(bootstrap_scores)}")
print()

# Comparison
print(f"{'='*80}")
print("CROSS-VALIDATION COMPARISON")
print(f"{'='*80}")
print()

comparison_df = pd.DataFrame({
    'Method': ['5-fold CV', 'LOOCV', 'Bootstrap'],
    'Mean AUC': [
        results['5fold_cv']['mean'],
        results['loocv']['mean'],
        results['bootstrap']['mean']
    ],
    'Std AUC': [
        results['5fold_cv']['std'],
        results['loocv']['std'],
        results['bootstrap']['std']
    ],
    'Min AUC': [
        results['5fold_cv']['min'],
        results['loocv']['min'],
        results['bootstrap']['min']
    ],
    'Max AUC': [
        results['5fold_cv']['max'],
        results['loocv']['max'],
        results['bootstrap']['max']
    ],
    'Coefficient of Variation': [
        results['5fold_cv']['std'] / results['5fold_cv']['mean'] * 100,
        results['loocv']['std'] / results['loocv']['mean'] * 100,
        results['bootstrap']['std'] / results['bootstrap']['mean'] * 100
    ]
})

print(comparison_df.to_string(index=False))
print()

# Interpretation
print("Interpretation:")
print()

# Find most stable method (lowest CV)
cv_values = [
    results['5fold_cv']['std'] / results['5fold_cv']['mean'],
    results['loocv']['std'] / results['loocv']['mean'],
    results['bootstrap']['std'] / results['bootstrap']['mean']
]
most_stable_idx = np.argmin(cv_values)
methods = ['5-fold CV', 'LOOCV', 'Bootstrap']
most_stable = methods[most_stable_idx]

print(f"  Most stable method: {most_stable} (lowest coefficient of variation)")
print()

if results['loocv']['std'] < results['5fold_cv']['std']:
    print("  ✅ LOOCV has lower variance than 5-fold CV (recommended for n=51)")
else:
    print("  ⚠️  5-fold CV has lower variance, but test folds are small (~5 responders)")

if results['bootstrap']['std'] < results['5fold_cv']['std']:
    print("  ✅ Bootstrap has lower variance than 5-fold CV (recommended for n=51)")
print()

# Save results
output_file = OUTPUT_DIR / "improved_cv_comparison.json"
with open(output_file, 'w') as f:
    json.dump(results, f, indent=2)
print(f"✅ Saved results: {output_file}")

comparison_df.to_csv(OUTPUT_DIR / "improved_cv_comparison.csv", index=False)
print(f"✅ Saved comparison: {OUTPUT_DIR / 'improved_cv_comparison.csv'}")

print(f"\n{'='*80}")
print("TASK 1.3C COMPLETE")
print(f"{'='*80}")
