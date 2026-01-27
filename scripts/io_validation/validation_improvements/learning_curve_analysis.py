#!/usr/bin/env python3
"""
Task 1.1C: Learning Curve Analysis
===================================

Assesses sample size adequacy by training models on increasing sample sizes.
Shows how performance changes as more data is added.

Method:
- Train models on n=20, 30, 40, 51 samples
- Evaluate on held-out test set
- Plot AUC vs. sample size
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
import json

# Configuration
SCRIPT_DIR = Path(__file__).parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
DATA_DIR = REPO_ROOT / "scripts" / "data_acquisition" / "IO"
OUTPUT_DIR = REPO_ROOT / "publications" / "06-io-response-prediction" / "data"
FIGURES_DIR = REPO_ROOT / "publications" / "06-io-response-prediction" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Set random seed for reproducibility
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Sample sizes to test
SAMPLE_SIZES = [20, 30, 40, 51]
N_REPEATS = 50  # Repeat each sample size multiple times for stability

print("="*80)
print("TASK 1.1C: LEARNING CURVE ANALYSIS")
print("="*80)
print()

# Load data
print("Loading GSE91061 data...")
df = pd.read_csv(DATA_DIR / "gse91061_analysis_with_composites.csv")

# Pathway columns (use 8-pathway model)
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

# Split into train/test (80/20) - test set is held out
X_train_full, X_test, y_train_full, y_test = train_test_split(
    X, y, test_size=0.2, random_state=RANDOM_SEED, stratify=y
)

print(f"Train set: {len(X_train_full)} samples")
print(f"Test set: {len(X_test)} samples")
print()

# Standardize using full training set
scaler = StandardScaler()
X_train_full_scaled = scaler.fit_transform(X_train_full)
X_test_scaled = scaler.transform(X_test)

# Learning curve
print(f"Computing learning curve (n={N_REPEATS} repeats per sample size)...")
print("This may take a few minutes...")
print()

results = []

for n_train in SAMPLE_SIZES:
    if n_train > len(X_train_full):
        continue
    
    print(f"Sample size: {n_train}...")
    
    train_aucs = []
    test_aucs = []
    
    for repeat in range(N_REPEATS):
        # Sample n_train samples from training set
        indices = np.random.choice(len(X_train_full), size=n_train, replace=False)
        X_train_subset = X_train_full_scaled[indices]
        y_train_subset = y_train_full[indices]
        
        # Train model
        lr = LogisticRegression(
            penalty='l2', 
            C=1.0, 
            random_state=RANDOM_SEED + repeat, 
            max_iter=1000
        )
        lr.fit(X_train_subset, y_train_subset)
        
        # Evaluate on training subset
        y_pred_train = lr.predict_proba(X_train_subset)[:, 1]
        train_auc = roc_auc_score(y_train_subset, y_pred_train)
        
        # Evaluate on held-out test set
        y_pred_test = lr.predict_proba(X_test_scaled)[:, 1]
        test_auc = roc_auc_score(y_test, y_pred_test)
        
        train_aucs.append(train_auc)
        test_aucs.append(test_auc)
    
    # Statistics
    train_auc_mean = np.mean(train_aucs)
    train_auc_std = np.std(train_aucs)
    test_auc_mean = np.mean(test_aucs)
    test_auc_std = np.std(test_aucs)
    
    results.append({
        'n_train': n_train,
        'train_auc_mean': train_auc_mean,
        'train_auc_std': train_auc_std,
        'test_auc_mean': test_auc_mean,
        'test_auc_std': test_auc_std,
        'gap': train_auc_mean - test_auc_mean,
        'gap_pct': (train_auc_mean - test_auc_mean) / train_auc_mean * 100
    })
    
    print(f"  Train AUC: {train_auc_mean:.3f} ± {train_auc_std:.3f}")
    print(f"  Test AUC:  {test_auc_mean:.3f} ± {test_auc_std:.3f}")
    print(f"  Gap:       {train_auc_mean - test_auc_mean:.3f} ({(train_auc_mean - test_auc_mean)/train_auc_mean*100:.1f}%)")
    print()

# Results summary
results_df = pd.DataFrame(results)
print(f"{'='*80}")
print("LEARNING CURVE RESULTS")
print(f"{'='*80}")
print(results_df.to_string(index=False))
print()

# Save results
results_df.to_csv(OUTPUT_DIR / "learning_curve_results.csv", index=False)
print(f"✅ Saved results: {OUTPUT_DIR / 'learning_curve_results.csv'}")

# Plot learning curve
fig, ax = plt.subplots(figsize=(10, 6))

n_train_sizes = results_df['n_train'].values
train_means = results_df['train_auc_mean'].values
train_stds = results_df['train_auc_std'].values
test_means = results_df['test_auc_mean'].values
test_stds = results_df['test_auc_std'].values

ax.plot(n_train_sizes, train_means, 'o-', label='Training AUC', linewidth=2, markersize=8)
ax.fill_between(n_train_sizes, train_means - train_stds, train_means + train_stds, alpha=0.2)

ax.plot(n_train_sizes, test_means, 's-', label='Test AUC', linewidth=2, markersize=8)
ax.fill_between(n_train_sizes, test_means - test_stds, test_means + test_stds, alpha=0.2)

ax.axhline(y=0.75, color='gray', linestyle='--', linewidth=1, alpha=0.5, label='Target (AUC=0.75)')
ax.axhline(y=0.70, color='gray', linestyle=':', linewidth=1, alpha=0.5, label='Acceptable (AUC=0.70)')

ax.set_xlabel('Training Sample Size', fontsize=12, fontweight='bold')
ax.set_ylabel('AUC', fontsize=12, fontweight='bold')
ax.set_title('Learning Curve: 8-Pathway Model Performance vs Sample Size', 
             fontsize=14, fontweight='bold')
ax.legend(loc='lower right', fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_ylim([0.5, 1.0])

plt.tight_layout()
plt.savefig(FIGURES_DIR / "learning_curve.png", dpi=300, bbox_inches='tight')
plt.savefig(FIGURES_DIR / "learning_curve.pdf", bbox_inches='tight')
print(f"✅ Saved figure: {FIGURES_DIR / 'learning_curve.png'}")

# Interpretation
print(f"{'='*80}")
print("INTERPRETATION")
print(f"{'='*80}")
print()

final_result = results_df.iloc[-1]
print(f"Final model (n={final_result['n_train']}):")
print(f"  Test AUC: {final_result['test_auc_mean']:.3f} ± {final_result['test_auc_std']:.3f}")
print(f"  Train-Test Gap: {final_result['gap']:.3f} ({final_result['gap_pct']:.1f}%)")
print()

if final_result['test_auc_mean'] >= 0.75:
    print("✅ Test AUC ≥ 0.75: Performance is acceptable")
elif final_result['test_auc_mean'] >= 0.70:
    print("⚠️  Test AUC 0.70-0.75: Performance is marginal, consider reducing model complexity")
else:
    print("❌ Test AUC < 0.70: Performance is insufficient, reduce model complexity")

if final_result['gap_pct'] < 10:
    print("✅ Train-test gap < 10%: Model is not overfitting significantly")
elif final_result['gap_pct'] < 20:
    print("⚠️  Train-test gap 10-20%: Moderate overfitting, consider regularization")
else:
    print("❌ Train-test gap > 20%: Significant overfitting, reduce model complexity")

# Check if performance plateaus
if len(results_df) >= 2:
    last_two = results_df.tail(2)
    test_auc_diff = last_two.iloc[1]['test_auc_mean'] - last_two.iloc[0]['test_auc_mean']
    if abs(test_auc_diff) < 0.01:
        print("✅ Performance plateaus: Current sample size (n=51) appears adequate")
    else:
        print(f"⚠️  Performance still improving: Adding more samples may help (ΔAUC = {test_auc_diff:.3f})")

print(f"\n{'='*80}")
print("TASK 1.1C COMPLETE")
print(f"{'='*80}")
