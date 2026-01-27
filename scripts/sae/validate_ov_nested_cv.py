#!/usr/bin/env python3
"""
Validate Ovarian SAE Model with Nested Cross-Validation
========================================================

CRITICAL: Prevents data leakage by performing feature selection
INSIDE each CV fold (not before CV).

Based on BRCA validation fixes - applies same safeguards to ovarian.
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

# Paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
TIER3_CHECKPOINT = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "checkpoints" / "Tier3_validation_cohort.json"
OUTPUT_DIR = PROJECT_ROOT / ".cursor" / "MOAT" / "SAE_INTELLIGENCE" / "OVARIAN_VALIDATION"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUTPUT_DIR / "ov_nested_cv_validation.json"

# Statistical thresholds
P_VALUE_THRESHOLD = 0.05
COHENS_D_THRESHOLD = 0.5
MAX_FEATURES = 29  # Match original: 9 diamonds + 20 top features


def fdr_correction(p_values: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    """FDR correction using Benjamini-Hochberg method."""
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


def build_feature_matrix(patients: List[Dict]) -> Tuple[np.ndarray, List[str]]:
    """
    Build feature matrix from ovarian SAE features.
    
    Aggregation: Mean across all variants per patient (matches manuscript).
    """
    n_patients = len(patients)
    n_features = 32768  # Full SAE feature space
    
    feature_matrix = np.zeros((n_patients, n_features))
    patient_ids = []
    
    for i, patient in enumerate(patients):
        patient_id = patient.get("patient_id") or patient.get("id")
        patient_ids.append(patient_id)
        
        variants = patient.get("variants", [])
        if not variants:
            continue
        
        # Aggregate features across variants (mean aggregation)
        variant_features = []
        for variant in variants:
            top_features = variant.get("top_features", [])
            for tf in top_features:
                idx = tf.get("index")
                val = tf.get("value", 0.0)
                if idx is not None and 0 <= idx < n_features:
                    variant_features.append((idx, abs(float(val))))
        
        # Mean aggregation: sum all feature values, then divide by variant count
        for idx, val in variant_features:
            feature_matrix[i, idx] += val
        
        if len(variants) > 0:
            feature_matrix[i, :] = feature_matrix[i, :] / len(variants)
    
    return feature_matrix, patient_ids


def select_features_nested(
    X_train: np.ndarray,
    y_train: np.ndarray,
    p_threshold: float = 0.05,
    cohen_d_threshold: float = 0.5,
    max_features: int = 29
) -> List[int]:
    """
    Select features on training fold only (prevents leakage).
    
    Method:
    1. Compute Mann-Whitney U test (resistant vs sensitive)
    2. Compute Cohen's d (effect size)
    3. Apply FDR correction
    4. Select features with p < threshold AND Cohen's d >= threshold
    5. Rank by Cohen's d, take top N
    """
    n_features = X_train.shape[1]
    
    # Compute statistics on training fold ONLY
    p_values = np.ones(n_features)
    cohen_d_values = np.zeros(n_features)
    
    resistant_mask = y_train.astype(bool)
    sensitive_mask = ~resistant_mask
    
    for i in range(n_features):
        feature_values = X_train[:, i]
        resistant_vals = feature_values[resistant_mask]
        sensitive_vals = feature_values[sensitive_mask]
        
        if len(resistant_vals) == 0 or len(sensitive_vals) == 0:
            continue
        
        # Mann-Whitney U test (one-sided: resistant > sensitive)
        try:
            u_stat, p_val = stats.mannwhitneyu(
                resistant_vals, sensitive_vals, alternative='greater'
            )
            p_values[i] = p_val
        except:
            p_values[i] = 1.0
        
        # Cohen's d
        mean_resistant = np.mean(resistant_vals)
        mean_sensitive = np.mean(sensitive_vals)
        std_resistant = np.std(resistant_vals)
        std_sensitive = np.std(sensitive_vals)
        
        pooled_std = np.sqrt((std_resistant**2 + std_sensitive**2) / 2)
        if pooled_std > 0:
            cohen_d_values[i] = abs((mean_resistant - mean_sensitive) / pooled_std)
    
    # Apply FDR correction
    valid_mask = p_values < 1.0
    if np.sum(valid_mask) > 0:
        valid_p = p_values[valid_mask]
        p_corrected = fdr_correction(valid_p, alpha=p_threshold)
        # Map back to full array
        p_corrected_full = np.ones(n_features)
        p_corrected_full[valid_mask] = p_corrected
    else:
        p_corrected_full = p_values
    
    # Select features: FDR-corrected p < threshold AND Cohen's d >= threshold
    significant = []
    for i in range(n_features):
        if p_corrected_full[i] < p_threshold and cohen_d_values[i] >= cohen_d_threshold:
            significant.append((i, cohen_d_values[i], p_corrected_full[i]))
    
    # Sort by Cohen's d (descending), take top N
    significant.sort(key=lambda x: x[1], reverse=True)
    selected_indices = [f[0] for f in significant[:max_features]]
    
    return selected_indices


def nested_cv_validation() -> Dict[str, Any]:
    """
    Run nested cross-validation to prevent feature selection leakage.
    
    CRITICAL SAFEGUARDS:
    - Feature selection: ONLY on training fold
    - Preprocessing: Fit ONLY on training fold
    - Model training: ONLY on training fold
    - Evaluation: ONLY on test fold
    """
    print("=" * 80)
    print("🔬 OVARIAN NESTED CV VALIDATION - Preventing Data Leakage")
    print("=" * 80)
    print()
    
    # Load data
    print("📥 Loading Tier-3 cohort data...")
    with open(TIER3_CHECKPOINT, 'r') as f:
        tier3_data = json.load(f)
    
    # Handle different data structures
    if isinstance(tier3_data, dict):
        patients_dict = tier3_data.get("data", {})
        if not patients_dict:
            patients_dict = tier3_data
    else:
        patients_dict = {}
    
    # Convert to list if needed
    if isinstance(patients_dict, dict):
        patients_list = list(patients_dict.values())
    else:
        patients_list = patients_dict
    
    # Extract labels and filter
    patients_filtered = []
    outcome_labels = []
    patient_ids = []
    
    for patient in patients_list:
        # Tier-3 uses "outcome" field: "sensitive", "resistant", "refractory"
        platinum_response = patient.get("outcome") or patient.get("platinum_response")
        if platinum_response not in ["sensitive", "resistant", "refractory"]:
            continue
        
        # Binary label: resistant/refractory = 1, sensitive = 0
        label = 1 if platinum_response in ["resistant", "refractory"] else 0
        
        patient_id = patient.get("patient_id") or patient.get("id")
        if not patient_id:
            continue
        
        # Only include patients with variants
        if not patient.get("variants"):
            continue
        
        patients_filtered.append(patient)
        outcome_labels.append(label)
        patient_ids.append(patient_id)
    
    print(f"   ✅ Loaded {len(patients_filtered)} patients")
    print(f"   Resistant/Refractory: {sum(outcome_labels)}")
    print(f"   Sensitive: {len(outcome_labels) - sum(outcome_labels)}")
    
    if len(patients_filtered) < 20:
        print("   ❌ Insufficient patients for validation")
        return {}
    
    # Build feature matrix
    print()
    print("🔍 Building feature matrix (mean aggregation)...")
    feature_matrix, feature_patient_ids = build_feature_matrix(patients_filtered)
    X = feature_matrix
    y = np.array(outcome_labels)
    
    print(f"   ✅ Feature matrix: {X.shape}")
    print(f"   ✅ Non-zero features: {np.sum(np.any(X != 0, axis=0))}")
    
    # ✅ ISOLATION CHECK: Verify patient ID alignment
    assert len(patient_ids) == len(feature_patient_ids), "PATIENT ID MISMATCH!"
    
    # Nested CV: Outer loop for model evaluation, inner loop for feature selection
    print()
    print("📊 Running nested cross-validation...")
    print("   Outer CV: 5-fold (model evaluation)")
    print("   Inner CV: Feature selection on training fold only")
    print("   CRITICAL: All feature selection/preprocessing on training only!")
    print()
    
    outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_scores = []
    fold_results = []
    feature_selection_log = []
    
    for fold_idx, (train_idx, test_idx) in enumerate(outer_cv.split(X, y)):
        print(f"   Fold {fold_idx + 1}/5...")
        
        # ✅ ISOLATION CHECK: Verify no overlap
        assert len(set(train_idx) & set(test_idx)) == 0, f"TRAIN/TEST OVERLAP in fold {fold_idx + 1}!"
        
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        train_ids = [patient_ids[i] for i in train_idx]
        test_ids = [patient_ids[i] for i in test_idx]
        
        # ✅ ISOLATION CHECK: Verify no patient ID overlap
        assert len(set(train_ids) & set(test_ids)) == 0, f"PATIENT ID OVERLAP in fold {fold_idx + 1}!"
        
        # ✅ Feature selection on training fold ONLY
        # Debug: Check feature statistics before selection
        n_nonzero_features = np.sum(np.any(X_train != 0, axis=0))
        print(f"      📊 Non-zero features in training: {n_nonzero_features}")
        
        selected_features = select_features_nested(
            X_train, y_train,
            p_threshold=P_VALUE_THRESHOLD,
            cohen_d_threshold=COHENS_D_THRESHOLD,
            max_features=MAX_FEATURES
        )
        
        # ✅ LOG: Track selected features
        feature_selection_log.append({
            'outer_fold': fold_idx + 1,
            'selected_features': selected_features,
            'n_features': len(selected_features)
        })
        
        if len(selected_features) == 0:
            print(f"      ⚠️  No features selected (may need relaxed thresholds)")
            print(f"      💡 Trying relaxed thresholds (p=0.1, d=0.3)...")
            
            # Try relaxed thresholds
            selected_features = select_features_nested(
                X_train, y_train,
                p_threshold=0.1,  # Relaxed
                cohen_d_threshold=0.3,  # Relaxed
                max_features=MAX_FEATURES
            )
            
            if len(selected_features) == 0:
                print(f"      ⚠️  Still no features selected, using random performance")
                cv_scores.append(0.5)
                fold_results.append({
                    'fold': fold_idx + 1,
                    'auroc': 0.5,
                    'n_features': 0,
                    'selected_features': []
                })
                continue
            else:
                print(f"      ✅ Selected {len(selected_features)} features with relaxed thresholds")
        
        print(f"      ✅ Selected {len(selected_features)} features")
        
        # Extract selected features
        X_train_selected = X_train[:, selected_features]
        X_test_selected = X_test[:, selected_features]
        
        # ✅ Preprocessing: Fit ONLY on training fold
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train_selected)  # Fit on train
        X_test_scaled = scaler.transform(X_test_selected)  # Apply to test (don't fit!)
        
        # ✅ Model training: ONLY on training fold
        model = LogisticRegression(
            class_weight='balanced',
            max_iter=1000,
            random_state=42
        )
        model.fit(X_train_scaled, y_train)
        
        # ✅ Evaluation: ONLY on test fold
        y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
        
        if len(set(y_test)) >= 2:  # Both classes present
            auroc = roc_auc_score(y_test, y_pred_proba)
            cv_scores.append(auroc)
            print(f"      ✅ AUROC: {auroc:.3f}")
            
            fold_results.append({
                'fold': fold_idx + 1,
                'auroc': float(auroc),
                'n_features': len(selected_features),
                'selected_features': selected_features
            })
        else:
            print(f"      ⚠️  Only one class in test set")
            cv_scores.append(0.5)
            fold_results.append({
                'fold': fold_idx + 1,
                'auroc': 0.5,
                'n_features': len(selected_features),
                'selected_features': selected_features
            })
    
    # Report results
    mean_auroc = np.mean(cv_scores)
    std_auroc = np.std(cv_scores)
    
    print()
    print("=" * 80)
    print("📊 NESTED CV RESULTS")
    print("=" * 80)
    print(f"   Mean AUROC: {mean_auroc:.3f} ± {std_auroc:.3f}")
    print(f"   Fold scores: {[f'{s:.3f}' for s in cv_scores]}")
    print()
    
    # Compare to original 0.783
    original_auroc = 0.783
    performance_drop = original_auroc - mean_auroc
    
    print("📊 COMPARISON TO ORIGINAL")
    print("=" * 80)
    print(f"   Original (5-fold CV): {original_auroc:.3f}")
    print(f"   Nested CV: {mean_auroc:.3f} ± {std_auroc:.3f}")
    print(f"   Performance drop: {performance_drop:.3f} ({performance_drop*100:.1f} pp)")
    
    if performance_drop < 0.05:
        print("   ✅ No significant leakage detected (drop < 5 pp)")
    elif performance_drop < 0.10:
        print("   ⚠️  Possible leakage (drop 5-10 pp)")
    else:
        print("   🚨 Significant leakage detected (drop > 10 pp)")
    
    print()
    
    # Feature stability analysis
    if len(feature_selection_log) > 0:
        all_selected = []
        for log in feature_selection_log:
            all_selected.extend(log['selected_features'])
        
        from collections import Counter
        feature_counts = Counter(all_selected)
        stable_features = [f for f, count in feature_counts.items() if count >= 3]  # Selected in ≥3 folds
        
        print("📊 FEATURE STABILITY")
        print("=" * 80)
        print(f"   Total unique features selected: {len(feature_counts)}")
        print(f"   Stable features (≥3 folds): {len(stable_features)}")
        print(f"   Most frequently selected: {feature_counts.most_common(10)}")
        print()
    
    # Save results
    results = {
        'timestamp': datetime.now().isoformat(),
        'nested_cv_auroc': float(mean_auroc),
        'nested_cv_std': float(std_auroc),
        'fold_scores': [float(s) for s in cv_scores],
        'fold_results': fold_results,
        'feature_selection_log': feature_selection_log,
        'comparison': {
            'original_auroc': original_auroc,
            'nested_cv_auroc': float(mean_auroc),
            'performance_drop': float(performance_drop),
            'leakage_assessment': 'no_leakage' if performance_drop < 0.05 else 'possible_leakage' if performance_drop < 0.10 else 'significant_leakage'
        },
        'parameters': {
            'p_threshold': P_VALUE_THRESHOLD,
            'cohen_d_threshold': COHENS_D_THRESHOLD,
            'max_features': MAX_FEATURES,
            'n_patients': len(patients_filtered),
            'n_resistant': int(sum(outcome_labels)),
            'n_sensitive': int(len(outcome_labels) - sum(outcome_labels))
        }
    }
    
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"💾 Results saved to: {OUTPUT_FILE}")
    
    return results


if __name__ == "__main__":
    try:
        result = nested_cv_validation()
        if result:
            print("\n✅ NESTED CV VALIDATION COMPLETE")
            sys.exit(0)
        else:
            print("\n❌ VALIDATION FAILED")
            sys.exit(1)
    except Exception as e:
        print(f"\n❌ VALIDATION ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
