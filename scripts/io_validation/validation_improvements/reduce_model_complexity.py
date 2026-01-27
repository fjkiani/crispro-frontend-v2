#!/usr/bin/env python3
"""
Task 1.1A: Reduce Model Complexity
===================================

Addresses overfitting by reducing 8-pathway model to 3-4 pathways.
EPV ratio: 23 responders / 8 features = 2.9 (below minimum 10)
Target: 23 responders / 3-4 features = 7.7-5.8 EPV (better, though still below 10)

Selects top pathways by AUC and p-value significance.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve
import json

# Configuration
# Script is at: scripts/io_validation/validation_improvements/
SCRIPT_DIR = Path(__file__).parent  # validation_improvements/
REPO_ROOT = SCRIPT_DIR.parent.parent.parent  # repo root
DATA_DIR = REPO_ROOT / "scripts" / "data_acquisition" / "IO"
OUTPUT_DIR = REPO_ROOT / "publications" / "06-io-response-prediction" / "data"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Set random seed for reproducibility
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

print("="*80)
print("TASK 1.1A: REDUCE MODEL COMPLEXITY")
print("="*80)
print()

# Load data
print("Loading GSE91061 data...")
df = pd.read_csv(DATA_DIR / "gse91061_analysis_with_composites.csv")
pathway_stats = pd.read_csv(DATA_DIR / "gse91061_pathway_response_association.csv")

print(f"Loaded {len(df)} samples")
print(f"  Responders: {(df['response'] == 1).sum()}")
print(f"  Non-responders: {(df['response'] == 0).sum()}")
print()

# Pathway columns (8 pathways)
pathway_cols = ['TIL_INFILTRATION', 'T_EFFECTOR', 'ANGIOGENESIS', 'TGFB_RESISTANCE',
                'MYELOID_INFLAMMATION', 'PROLIFERATION', 'IMMUNOPROTEASOME', 'EXHAUSTION']

# Response labels
y = df['response'].values
n_responders = y.sum()
n_nonresponders = len(y) - n_responders

print(f"Current Model: 8 pathways")
print(f"  EPV ratio: {n_responders} responders / 8 features = {n_responders/8:.1f} EPV")
print(f"  Status: ⚠️  Below minimum 10 EPV (overfitting risk)")
print()

# Select top pathways by AUC and significance
print("Selecting top pathways...")
pathway_stats_sorted = pathway_stats.sort_values('auc', ascending=False)

# Filter to pathways with p < 0.05 (or p < 0.10 for borderline)
significant_pathways = pathway_stats_sorted[
    (pathway_stats_sorted['p_value'] < 0.10) & 
    (pathway_stats_sorted['pathway'] != 'PDL1_EXPRESSION')
]

print(f"Significant pathways (p < 0.10): {len(significant_pathways)}")
print(significant_pathways[['pathway', 'auc', 'p_value']].to_string(index=False))
print()

# Option 1: Top 3 pathways (EXHAUSTION, TIL_INFILTRATION, T_EFFECTOR)
top3_pathways = significant_pathways.head(3)['pathway'].tolist()
print(f"Option 1: Top 3 pathways")
print(f"  Pathways: {', '.join(top3_pathways)}")
print(f"  EPV ratio: {n_responders} responders / 3 features = {n_responders/3:.1f} EPV")
print()

# Option 2: Top 4 pathways (add ANGIOGENESIS if significant)
top4_pathways = significant_pathways.head(4)['pathway'].tolist()
print(f"Option 2: Top 4 pathways")
print(f"  Pathways: {', '.join(top4_pathways)}")
print(f"  EPV ratio: {n_responders} responders / 4 features = {n_responders/4:.1f} EPV")
print()

# Train models and compare performance
results = []

for n_pathways, pathway_list in [(3, top3_pathways), (4, top4_pathways), (8, pathway_cols)]:
    print(f"\n{'='*80}")
    print(f"Training {n_pathways}-pathway model...")
    print(f"{'='*80}")
    
    # Prepare features
    X = df[pathway_list].values
    
    # Remove rows with NaN
    valid_rows = ~np.isnan(X).any(axis=1)
    X = X[valid_rows]
    y_valid = y[valid_rows]
    
    print(f"Valid samples: {len(X)} (removed {len(y) - len(X)} with missing data)")
    
    if len(X) < 10:
        print(f"⚠️  Insufficient samples for {n_pathways}-pathway model")
        continue
    
    # Standardize features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Train logistic regression
    lr = LogisticRegression(
        penalty='l2', 
        C=1.0, 
        random_state=RANDOM_SEED, 
        max_iter=1000
    )
    lr.fit(X_scaled, y_valid)
    
    # Training AUC
    y_pred_proba = lr.predict_proba(X_scaled)[:, 1]
    train_auc = roc_auc_score(y_valid, y_pred_proba)
    
    # 5-fold CV
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_SEED)
    cv_scores = cross_val_score(lr, X_scaled, y_valid, cv=cv, scoring='roc_auc')
    cv_mean = cv_scores.mean()
    cv_std = cv_scores.std()
    
    # EPV ratio
    epv_ratio = y_valid.sum() / n_pathways
    
    # Store results
    results.append({
        'n_pathways': n_pathways,
        'pathways': ', '.join(pathway_list),
        'n_samples': len(X),
        'n_responders': y_valid.sum(),
        'epv_ratio': epv_ratio,
        'train_auc': train_auc,
        'cv_auc_mean': cv_mean,
        'cv_auc_std': cv_std,
        'cv_auc_lower': cv_mean - 1.96 * cv_std,
        'cv_auc_upper': cv_mean + 1.96 * cv_std,
        'auc_drop': train_auc - cv_mean,
        'auc_drop_pct': (train_auc - cv_mean) / train_auc * 100
    })
    
    print(f"  Training AUC: {train_auc:.3f}")
    print(f"  CV AUC: {cv_mean:.3f} ± {cv_std:.3f}")
    print(f"  AUC drop: {train_auc - cv_mean:.3f} ({(train_auc - cv_mean)/train_auc*100:.1f}%)")
    print(f"  EPV ratio: {epv_ratio:.1f}")
    
    # Save coefficients
    coef_df = pd.DataFrame({
        'pathway': pathway_list,
        'coefficient': lr.coef_[0],
        'abs_coefficient': np.abs(lr.coef_[0])
    }).sort_values('abs_coefficient', ascending=False)
    
    print(f"\n  Coefficients:")
    print(coef_df.to_string(index=False))
    
    # Save model
    model_output = {
        'n_pathways': n_pathways,
        'pathways': pathway_list,
        'coefficients': lr.coef_[0].tolist(),
        'intercept': float(lr.intercept_[0]),
        'scaler_mean': scaler.mean_.tolist(),
        'scaler_scale': scaler.scale_.tolist(),
        'train_auc': float(train_auc),
        'cv_auc_mean': float(cv_mean),
        'cv_auc_std': float(cv_std),
        'epv_ratio': float(epv_ratio),
        'n_samples': int(len(X)),
        'n_responders': int(y_valid.sum())
    }
    
    output_file = OUTPUT_DIR / f"reduced_model_{n_pathways}pathway.json"
    with open(output_file, 'w') as f:
        json.dump(model_output, f, indent=2)
    print(f"\n  ✅ Saved model: {output_file}")

# Compare models
print(f"\n{'='*80}")
print("MODEL COMPARISON")
print(f"{'='*80}")

results_df = pd.DataFrame(results)
print(results_df.to_string(index=False))

# Save comparison
results_df.to_csv(OUTPUT_DIR / "reduced_model_comparison.csv", index=False)
print(f"\n✅ Saved comparison: {OUTPUT_DIR / 'reduced_model_comparison.csv'}")

# Recommendation
print(f"\n{'='*80}")
print("RECOMMENDATION")
print(f"{'='*80}")

best_model = results_df.loc[results_df['cv_auc_mean'].idxmax()]
print(f"Best model: {best_model['n_pathways']}-pathway model")
print(f"  CV AUC: {best_model['cv_auc_mean']:.3f} ± {best_model['cv_auc_std']:.3f}")
print(f"  EPV ratio: {best_model['epv_ratio']:.1f}")
print(f"  AUC drop: {best_model['auc_drop']:.3f} ({best_model['auc_drop_pct']:.1f}%)")
print()

if best_model['epv_ratio'] < 10:
    print("⚠️  EPV ratio still below 10, but improved from 2.9")
    print("   Consider this a proof-of-concept model requiring external validation")
else:
    print("✅ EPV ratio ≥ 10, model complexity is acceptable")

print(f"\n{'='*80}")
print("TASK 1.1A COMPLETE")
print(f"{'='*80}")
