#!/usr/bin/env python3
"""
Task 1.3B: Calibration Analysis
================================

Assesses prediction calibration by comparing predicted vs. observed response rates.
High AUC doesn't guarantee good calibration.

Method:
1. Bin predicted probabilities into 5 quantile groups
2. Compute observed response rate per bin
3. Plot predicted vs. observed (calibration plot)
4. Compute Hosmer-Lemeshow test for calibration goodness-of-fit
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from scipy import stats
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

# Number of bins for calibration
N_BINS = 5

print("="*80)
print("TASK 1.3B: CALIBRATION ANALYSIS")
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

# Train logistic regression model
lr = LogisticRegression(
    penalty='l2', 
    C=1.0, 
    random_state=RANDOM_SEED, 
    max_iter=1000
)
lr.fit(X_scaled, y)

# Get predicted probabilities
y_pred_proba = lr.predict_proba(X_scaled)[:, 1]

# AUC
auc = roc_auc_score(y, y_pred_proba)
print(f"Model AUC: {auc:.3f}")
print()

# Calibration analysis
print(f"Computing calibration (n={N_BINS} bins)...")

# Bin predicted probabilities into quantiles
bin_edges = np.percentile(y_pred_proba, np.linspace(0, 100, N_BINS + 1))
bin_edges[0] = 0.0  # Ensure first bin starts at 0
bin_edges[-1] = 1.0  # Ensure last bin ends at 1

# Assign samples to bins
bin_indices = np.digitize(y_pred_proba, bin_edges) - 1
bin_indices = np.clip(bin_indices, 0, N_BINS - 1)  # Clip to valid range

# Compute calibration statistics per bin
calibration_data = []

for i in range(N_BINS):
    bin_mask = bin_indices == i
    n_in_bin = bin_mask.sum()
    
    if n_in_bin == 0:
        continue
    
    # Predicted probability (mean in bin)
    pred_mean = y_pred_proba[bin_mask].mean()
    pred_std = y_pred_proba[bin_mask].std()
    
    # Observed response rate
    obs_mean = y[bin_mask].mean()
    obs_std = y[bin_mask].std()
    
    # Standard error for observed rate
    obs_se = np.sqrt(obs_mean * (1 - obs_mean) / n_in_bin) if n_in_bin > 0 else 0
    
    calibration_data.append({
        'bin': i + 1,
        'n_samples': n_in_bin,
        'predicted_mean': pred_mean,
        'predicted_std': pred_std,
        'observed_mean': obs_mean,
        'observed_std': obs_std,
        'observed_se': obs_se,
        'calibration_error': abs(pred_mean - obs_mean)
    })

calibration_df = pd.DataFrame(calibration_data)

print()
print("Calibration Statistics:")
print(calibration_df.to_string(index=False))
print()

# Hosmer-Lemeshow test
# H0: Model is well-calibrated (predicted ≈ observed)
# Test statistic: χ² = Σ[(observed - expected)² / (expected * (1 - expected) / n)]

hl_statistic = 0
for _, row in calibration_df.iterrows():
    n = row['n_samples']
    observed = row['observed_mean']
    expected = row['predicted_mean']
    
    if n > 0 and expected > 0 and expected < 1:
        numerator = (observed - expected) ** 2
        denominator = expected * (1 - expected) / n
        if denominator > 0:
            hl_statistic += numerator / denominator

# Degrees of freedom: n_bins - 2 (for logistic regression)
df = len(calibration_df) - 2
hl_p_value = 1 - stats.chi2.cdf(hl_statistic, df)

print(f"Hosmer-Lemeshow Test:")
print(f"  χ² statistic: {hl_statistic:.3f}")
print(f"  Degrees of freedom: {df}")
print(f"  p-value: {hl_p_value:.4f}")
print()

if hl_p_value > 0.05:
    print("✅ p > 0.05: Model is well-calibrated (fail to reject H0)")
else:
    print("❌ p ≤ 0.05: Model is poorly calibrated (reject H0)")

# Mean calibration error
mean_calibration_error = calibration_df['calibration_error'].mean()
print(f"\nMean Calibration Error: {mean_calibration_error:.3f}")
print()

# Save results
results = {
    'auc': float(auc),
    'n_bins': int(N_BINS),
    'mean_calibration_error': float(mean_calibration_error),
    'hosmer_lemeshow_statistic': float(hl_statistic),
    'hosmer_lemeshow_df': int(df),
    'hosmer_lemeshow_p_value': float(hl_p_value),
    'well_calibrated': bool(hl_p_value > 0.05),
    'calibration_bins': calibration_df.to_dict('records')
}

output_file = OUTPUT_DIR / "calibration_analysis_results.json"
with open(output_file, 'w') as f:
    json.dump(results, f, indent=2)
print(f"✅ Saved results: {output_file}")

calibration_df.to_csv(OUTPUT_DIR / "calibration_bins.csv", index=False)
print(f"✅ Saved bins: {OUTPUT_DIR / 'calibration_bins.csv'}")

# Plot calibration curve
fig, ax = plt.subplots(figsize=(8, 8))

# Plot perfect calibration line
ax.plot([0, 1], [0, 1], 'k--', linewidth=2, label='Perfect Calibration', alpha=0.5)

# Plot calibration curve
predicted = calibration_df['predicted_mean'].values
observed = calibration_df['observed_mean'].values
observed_se = calibration_df['observed_se'].values

ax.errorbar(predicted, observed, yerr=observed_se, 
            fmt='o-', linewidth=2, markersize=8, capsize=5,
            label='Model Calibration', color='blue')

# Add sample size annotations
for _, row in calibration_df.iterrows():
    ax.annotate(f"n={row['n_samples']}", 
                xy=(row['predicted_mean'], row['observed_mean']),
                xytext=(5, 5), textcoords='offset points', fontsize=9, alpha=0.7)

ax.set_xlabel('Predicted Probability', fontsize=12, fontweight='bold')
ax.set_ylabel('Observed Response Rate', fontsize=12, fontweight='bold')
ax.set_title('Calibration Plot: Predicted vs. Observed Response Rates', 
             fontsize=14, fontweight='bold')
ax.legend(loc='upper left', fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])

# Add text with statistics
textstr = f'AUC: {auc:.3f}\nMean Error: {mean_calibration_error:.3f}\nH-L p: {hl_p_value:.4f}'
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig(FIGURES_DIR / "calibration_plot.png", dpi=300, bbox_inches='tight')
plt.savefig(FIGURES_DIR / "calibration_plot.pdf", bbox_inches='tight')
print(f"✅ Saved figure: {FIGURES_DIR / 'calibration_plot.png'}")

print(f"\n{'='*80}")
print("TASK 1.3B COMPLETE")
print(f"{'='*80}")
