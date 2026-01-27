#!/usr/bin/env python3
"""
Validate BRCA Model with Nested Cross-Validation
=================================================

FIXES CRITICAL ISSUE: Feature selection leakage
- Feature selection performed INSIDE each CV fold
- Prevents data leakage
- Provides unbiased performance estimate
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Any, Tuple
from datetime import datetime
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from scipy import stats
# FDR correction (simple implementation)
def fdr_correction(p_values, alpha=0.05):
    """Simple FDR correction (Benjamini-Hochberg)."""
    if len(p_values) == 0:
        return p_values
    
    p_array = np.array(p_values)
    n = len(p_array)
    
    # Sort p-values
    sorted_indices = np.argsort(p_array)
    sorted_p = p_array[sorted_indices]
    
    # Compute adjusted p-values
    adjusted_p = np.zeros(n)
    for i in range(n-1, -1, -1):
        if i == n-1:
            adjusted_p[sorted_indices[i]] = sorted_p[i]
        else:
            adjusted_p[sorted_indices[i]] = min(
                adjusted_p[sorted_indices[i+1]],
                sorted_p[i] * n / (i + 1)
            )
    
    return adjusted_p

# Paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
BRCA_CHECKPOINT = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "checkpoints" / "BRCA_TCGA_TRUE_SAE_cohort.json"
RECURRENCE_LABELS = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "brca_recurrence_labels.json"
OUTPUT_DIR = PROJECT_ROOT / ".cursor" / "MOAT" / "SAE_INTELLIGENCE" / "BRCA_VALIDATION"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUTPUT_DIR / "brca_nested_cv_validation.json"

# Statistical thresholds
P_VALUE_THRESHOLD = 0.05
COHENS_D_THRESHOLD = 0.5
MAX_FEATURES = 10  # Reduced from 17 to prevent overfitting


def build_feature_matrix(patients: List[Dict]) -> Tuple[np.ndarray, List[str]]:
    """Build feature matrix from BRCA SAE features."""
    n_patients = len(patients)
    n_features = 32768
    
    feature_matrix = np.zeros((n_patients, n_features))
    patient_ids = []
    
    for i, patient in enumerate(patients):
        patient_id = patient.get("patient_id")
        patient_ids.append(patient_id)
        
        variants = patient.get("variants", [])
        for variant in variants:
            top_features = variant.get("top_features", [])
            for tf in top_features:
                idx = tf.get("index")
                val = tf.get("value", 0.0)
                if idx is not None and 0 <= idx < n_features:
                    feature_matrix[i, idx] += abs(float(val))
        
        if len(variants) > 0:
            feature_matrix[i, :] = feature_matrix[i, :] / len(variants)
    
    return feature_matrix, patient_ids


def select_features_nested(
    X_train: np.ndarray,
    y_train: np.ndarray,
    p_threshold: float = 0.05,
    cohen_d_threshold: float = 0.5,
    max_features: int = 10
) -> List[int]:
    """Select features on training fold only (prevents leakage)."""
    n_features = X_train.shape[1]
    
    # Compute correlations
    pearson_r = np.zeros(n_features)
    pearson_p = np.ones(n_features)
    
    for i in range(n_features):
        feature_values = X_train[:, i]
        if np.std(feature_values) == 0:
            continue
        r, p = stats.pearsonr(feature_values, y_train)
        pearson_r[i] = r
        pearson_p[i] = p
    
    # Compute Cohen's d
    cohen_d = np.zeros(n_features)
    recurrence_mask = y_train.astype(bool)
    
    for i in range(n_features):
        feature_values = X_train[:, i]
        recurrence_values = feature_values[recurrence_mask]
        no_recurrence_values = feature_values[~recurrence_mask]
        
        if len(recurrence_values) == 0 or len(no_recurrence_values) == 0:
            continue
        
        mean_diff = np.mean(recurrence_values) - np.mean(no_recurrence_values)
        pooled_std = np.sqrt((np.var(recurrence_values) + np.var(no_recurrence_values)) / 2)
        
        if pooled_std == 0:
            continue
        
        cohen_d[i] = abs(mean_diff / pooled_std)
    
    # Apply FDR correction
    valid_mask = pearson_p < 1.0
    if np.sum(valid_mask) > 0:
        valid_p = pearson_p[valid_mask]
        p_corrected = fdr_correction(valid_p, alpha=p_threshold)
        # Map back to full array
        p_corrected_full = np.ones(n_features)
        p_corrected_full[valid_mask] = p_corrected
    else:
        p_corrected_full = pearson_p
    
    # Select features
    significant = []
    for i in range(n_features):
        if p_corrected_full[i] < p_threshold and cohen_d[i] >= cohen_d_threshold:
            significant.append((i, cohen_d[i], pearson_r[i], p_corrected_full[i]))
    
    # Sort by Cohen's d, take top N
    significant.sort(key=lambda x: x[1], reverse=True)
    selected_indices = [f[0] for f in significant[:max_features]]
    
    return selected_indices


def nested_cv_validation() -> Dict[str, Any]:
    """Run nested cross-validation to prevent feature selection leakage."""
    print("=" * 80)
    print("🔬 NESTED CV VALIDATION - Fixing Feature Selection Leakage")
    print("=" * 80)
    print()
    
    # Load data
    print("📥 Loading data...")
    with open(BRCA_CHECKPOINT, 'r') as f:
        brca_sae = json.load(f)
    
    patients_dict = brca_sae.get("data", {})
    
    with open(RECURRENCE_LABELS, 'r') as f:
        labels_data = json.load(f)
    
    outcomes = labels_data.get("outcomes", {})
    
    # Prepare data
    patients_list = []
    outcome_labels = []
    
    for patient_id, patient_data in patients_dict.items():
        outcome = outcomes.get(patient_id, {})
        recurrence = outcome.get("recurrence")
        
        if recurrence is None:
            continue
        
        patient_data["patient_id"] = patient_id
        patients_list.append(patient_data)
        outcome_labels.append(1 if recurrence else 0)
    
    print(f"   ✅ Loaded {len(patients_list)} patients")
    print(f"   Recurrence: {sum(outcome_labels)} True, {len(outcome_labels) - sum(outcome_labels)} False")
    
    # Build full feature matrix
    print()
    print("🔍 Building feature matrix...")
    feature_matrix, _ = build_feature_matrix(patients_list)
    X = feature_matrix
    y = np.array(outcome_labels)
    
    print(f"   ✅ Feature matrix: {X.shape}")
    
    # Nested CV: Outer loop for model evaluation, inner loop for feature selection
    print()
    print("📊 Running nested cross-validation...")
    print("   Outer CV: 5-fold (model evaluation)")
    print("   Inner CV: Feature selection on training fold only")
    print()
    
    outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_scores = []
    fold_results = []
    
    for fold_idx, (train_idx, test_idx) in enumerate(outer_cv.split(X, y)):
        print(f"   Fold {fold_idx + 1}/5...")
        
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        
        # Feature selection on training fold only
        selected_features = select_features_nested(
            X_train, y_train,
            p_threshold=P_VALUE_THRESHOLD,
            cohen_d_threshold=COHENS_D_THRESHOLD,
            max_features=MAX_FEATURES
        )
        
        if len(selected_features) == 0:
            print(f"      ⚠️  No features selected, skipping fold")
            continue
        
        # Train model on selected features
        X_train_selected = X_train[:, selected_features]
        X_test_selected = X_test[:, selected_features]
        
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train_selected)
        X_test_scaled = scaler.transform(X_test_selected)
        
        model = LogisticRegression(
            class_weight='balanced',
            max_iter=1000,
            random_state=42,
            solver='lbfgs'
        )
        
        model.fit(X_train_scaled, y_train)
        y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
        
        # Compute AUROC
        if len(set(y_test)) >= 2:
            auroc = roc_auc_score(y_test, y_pred_proba)
            cv_scores.append(auroc)
            fold_results.append({
                "fold": fold_idx + 1,
                "n_features": len(selected_features),
                "features": selected_features,
                "auroc": float(auroc),
                "n_train": len(y_train),
                "n_test": len(y_test),
                "n_pos_train": int(sum(y_train)),
                "n_pos_test": int(sum(y_test))
            })
            print(f"      ✅ AUROC: {auroc:.3f} ({len(selected_features)} features)")
        else:
            print(f"      ⚠️  Insufficient outcome diversity in test fold")
    
    # Summary
    print()
    print("=" * 80)
    print("📊 NESTED CV RESULTS")
    print("=" * 80)
    
    if len(cv_scores) > 0:
        mean_auroc = np.mean(cv_scores)
        std_auroc = np.std(cv_scores)
        
        print(f"CV AUROC: {mean_auroc:.3f} ± {std_auroc:.3f}")
        print(f"CV scores: {[f'{s:.3f}' for s in cv_scores]}")
        print(f"Folds completed: {len(cv_scores)}/5")
        
        # Compare to original (with leakage)
        original_cv = 0.844
        print()
        print("Comparison to original (with leakage):")
        print(f"  Original CV AUROC: {original_cv:.3f}")
        print(f"  Nested CV AUROC: {mean_auroc:.3f}")
        print(f"  Difference: {mean_auroc - original_cv:+.3f}")
        
        if mean_auroc < original_cv - 0.05:
            print("  ⚠️  Significant drop - original was inflated by leakage")
        elif mean_auroc > original_cv - 0.02:
            print("  ✅ Similar performance - leakage had minimal impact")
        else:
            print("  ⚠️  Some impact from leakage")
    else:
        print("❌ No valid folds completed")
        mean_auroc = 0.0
        std_auroc = 0.0
    
    # Save results
    results = {
        "validation_date": datetime.now().isoformat(),
        "method": "nested_cross_validation",
        "outer_cv_folds": 5,
        "feature_selection": "inside_cv_fold",
        "fdr_correction": True,
        "max_features": MAX_FEATURES,
        "performance": {
            "cv_auroc_mean": float(mean_auroc),
            "cv_auroc_std": float(std_auroc),
            "cv_scores": [float(s) for s in cv_scores],
            "n_folds": len(cv_scores)
        },
        "fold_results": fold_results,
        "comparison": {
            "original_cv_auroc": 0.844,
            "nested_cv_auroc": float(mean_auroc),
            "difference": float(mean_auroc - 0.844)
        }
    }
    
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(results, f, indent=2)
    
    print()
    print(f"💾 Results saved to: {OUTPUT_FILE}")
    print("=" * 80)
    
    return results


if __name__ == "__main__":
    try:
        result = nested_cv_validation()
        print("\n✅ NESTED CV VALIDATION COMPLETE")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ VALIDATION FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
