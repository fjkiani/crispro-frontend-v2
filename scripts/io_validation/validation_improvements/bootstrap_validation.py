#!/usr/bin/env python3
"""
Task 1.1B: Bootstrap Validation (Optimism Correction)
======================================================

Performs bootstrap validation to correct for overfitting.
Computes optimism-corrected AUC to get realistic performance estimate.

Method:
1. Resample with replacement (n=51)
2. Train model on resampled data
3. Evaluate on original data (optimistic)
4. Evaluate on out-of-bag samples (realistic)
5. Optimism = optimistic - realistic
6. Corrected AUC = training AUC - mean optimism
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
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
print("TASK 1.1B: BOOTSTRAP VALIDATION (OPTIMISM CORRECTION)")
print("="*80)
print()

# Load data
print("Loading GSE91061 data...")
df = pd.read_csv(DATA_DIR / "gse91061_analysis_with_composites.csv")

# Pathway columns (use 8-pathway model for now, can be reduced later)
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

# Train full model on all data (for comparison)
lr_full = LogisticRegression(
    penalty='l2', 
    C=1.0, 
    random_state=RANDOM_SEED, 
    max_iter=1000
)
lr_full.fit(X_scaled, y)
y_pred_proba_full = lr_full.predict_proba(X_scaled)[:, 1]
train_auc_full = roc_auc_score(y, y_pred_proba_full)

print(f"Full model (all data):")
print(f"  Training AUC: {train_auc_full:.3f}")
print()

# Bootstrap validation
print(f"Running bootstrap validation ({N_BOOTSTRAP} iterations)...")
print("This may take a few minutes...")
print()

optimistic_aucs = []
realistic_aucs = []
optimisms = []

for i in tqdm(range(N_BOOTSTRAP), desc="Bootstrap"):
    # Resample with replacement
    n_samples = len(X)
    bootstrap_indices = np.random.choice(n_samples, size=n_samples, replace=True)
    
    # Get out-of-bag indices (samples not in bootstrap)
    all_indices = set(range(n_samples))
    bootstrap_set = set(bootstrap_indices)
    oob_indices = np.array(list(all_indices - bootstrap_set))
    
    # Skip if no out-of-bag samples
    if len(oob_indices) == 0:
        continue
    
    # Bootstrap training set
    X_boot = X_scaled[bootstrap_indices]
    y_boot = y[bootstrap_indices]
    
    # Train model on bootstrap sample
    lr_boot = LogisticRegression(
        penalty='l2', 
        C=1.0, 
        random_state=RANDOM_SEED + i,  # Vary seed slightly
        max_iter=1000
    )
    lr_boot.fit(X_boot, y_boot)
    
    # Optimistic: Evaluate on original data (all samples)
    y_pred_optimistic = lr_boot.predict_proba(X_scaled)[:, 1]
    optimistic_auc = roc_auc_score(y, y_pred_optimistic)
    
    # Realistic: Evaluate on out-of-bag samples
    if len(oob_indices) > 0:
        X_oob = X_scaled[oob_indices]
        y_oob = y[oob_indices]
        y_pred_realistic = lr_boot.predict_proba(X_oob)[:, 1]
        realistic_auc = roc_auc_score(y_oob, y_pred_realistic)
        
        # Optimism
        optimism = optimistic_auc - realistic_auc
        
        optimistic_aucs.append(optimistic_auc)
        realistic_aucs.append(realistic_auc)
        optimisms.append(optimism)

# Compute statistics
optimistic_auc_mean = np.mean(optimistic_aucs)
optimistic_auc_std = np.std(optimistic_aucs)
realistic_auc_mean = np.mean(realistic_aucs)
realistic_auc_std = np.std(realistic_aucs)
optimism_mean = np.mean(optimisms)
optimism_std = np.std(optimisms)

# Corrected AUC
corrected_auc = train_auc_full - optimism_mean

# Bootstrap CI for corrected AUC (percentile method)
optimism_sorted = np.sort(optimisms)
ci_lower_idx = int(0.025 * len(optimism_sorted))
ci_upper_idx = int(0.975 * len(optimism_sorted))
optimism_lower = optimism_sorted[ci_lower_idx]
optimism_upper = optimism_sorted[ci_upper_idx]

corrected_auc_lower = train_auc_full - optimism_upper  # Upper optimism = lower corrected AUC
corrected_auc_upper = train_auc_full - optimism_lower  # Lower optimism = upper corrected AUC

# Results
print(f"\n{'='*80}")
print("BOOTSTRAP VALIDATION RESULTS")
print(f"{'='*80}")
print()
print(f"Training AUC (full model): {train_auc_full:.3f}")
print()
print(f"Bootstrap Statistics (n={len(optimisms)} iterations):")
print(f"  Optimistic AUC: {optimistic_auc_mean:.3f} ± {optimistic_auc_std:.3f}")
print(f"  Realistic AUC:  {realistic_auc_mean:.3f} ± {realistic_auc_std:.3f}")
print(f"  Optimism:       {optimism_mean:.3f} ± {optimism_std:.3f}")
print(f"    95% CI:      [{optimism_lower:.3f}, {optimism_upper:.3f}]")
print()
print(f"Optimism-Corrected AUC:")
print(f"  Corrected AUC:  {corrected_auc:.3f}")
print(f"  95% CI:        [{corrected_auc_lower:.3f}, {corrected_auc_upper:.3f}]")
print()
print(f"Overfitting Assessment:")
print(f"  Training AUC:  {train_auc_full:.3f}")
print(f"  Corrected AUC: {corrected_auc:.3f}")
print(f"  Overfitting:   {train_auc_full - corrected_auc:.3f} ({(train_auc_full - corrected_auc)/train_auc_full*100:.1f}%)")
print()

# Save results
results = {
    'training_auc': float(train_auc_full),
    'bootstrap_n_iterations': int(len(optimisms)),
    'optimistic_auc_mean': float(optimistic_auc_mean),
    'optimistic_auc_std': float(optimistic_auc_std),
    'realistic_auc_mean': float(realistic_auc_mean),
    'realistic_auc_std': float(realistic_auc_std),
    'optimism_mean': float(optimism_mean),
    'optimism_std': float(optimism_std),
    'optimism_ci_lower': float(optimism_lower),
    'optimism_ci_upper': float(optimism_upper),
    'corrected_auc': float(corrected_auc),
    'corrected_auc_ci_lower': float(corrected_auc_lower),
    'corrected_auc_ci_upper': float(corrected_auc_upper),
    'overfitting_amount': float(train_auc_full - corrected_auc),
    'overfitting_pct': float((train_auc_full - corrected_auc) / train_auc_full * 100)
}

output_file = OUTPUT_DIR / "bootstrap_validation_results.json"
with open(output_file, 'w') as f:
    json.dump(results, f, indent=2)
print(f"✅ Saved results: {output_file}")

# Save detailed bootstrap iterations
bootstrap_df = pd.DataFrame({
    'iteration': range(len(optimisms)),
    'optimistic_auc': optimistic_aucs,
    'realistic_auc': realistic_aucs,
    'optimism': optimisms
})
bootstrap_df.to_csv(OUTPUT_DIR / "bootstrap_iterations.csv", index=False)
print(f"✅ Saved iterations: {OUTPUT_DIR / 'bootstrap_iterations.csv'}")

print(f"\n{'='*80}")
print("TASK 1.1B COMPLETE")
print(f"{'='*80}")
