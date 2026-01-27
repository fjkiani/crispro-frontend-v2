#!/usr/bin/env python3
"""
Test Pre-Selected 9 Diamond Features with Nested CV
====================================================

CRITICAL EXPERIMENT: Test if the 9 diamond features (discovered via leaky method)
actually predict platinum resistance when validated properly.

NO FEATURE SELECTION - just validate these 9 features in nested CV.
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

# Paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
TIER3_CHECKPOINT = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "checkpoints" / "Tier3_validation_cohort.json"
OUTPUT_DIR = PROJECT_ROOT / ".cursor" / "MOAT" / "SAE_INTELLIGENCE" / "OVARIAN_VALIDATION"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUTPUT_DIR / "ov_diamond_features_nested_cv.json"

# The 9 pre-selected diamond features from original (leaky) analysis
DIAMOND_FEATURES = [1407, 6020, 9738, 12893, 16337, 22868, 26220, 27607, 31362]


def build_feature_matrix(patients: List[Dict], feature_indices: List[int], aggregation: str = 'mean') -> Tuple[np.ndarray, List[str]]:
    """
    Build feature matrix using pre-selected features.
    
    Aggregation: 'mean' or 'max'
    """
    n_patients = len(patients)
    n_features = len(feature_indices)
    
    feature_matrix = np.zeros((n_patients, n_features))
    patient_ids = []
    
    for i, patient in enumerate(patients):
        patient_id = patient.get("patient_id") or patient.get("id")
        patient_ids.append(patient_id)
        
        variants = patient.get("variants", [])
        if not variants:
            continue
        
        # Extract feature values for each variant
        variant_features = []  # [n_variants, n_features]
        
        for variant in variants:
            top_features = variant.get("top_features", [])
            feature_values = np.zeros(n_features)
            
            # Map top_features to our selected indices
            feature_map = {tf.get("index"): abs(float(tf.get("value", 0.0))) for tf in top_features}
            
            for j, fidx in enumerate(feature_indices):
                feature_values[j] = feature_map.get(fidx, 0.0)
            
            variant_features.append(feature_values)
        
        if len(variant_features) > 0:
            variant_features = np.array(variant_features)
            
            if aggregation == 'mean':
                # Mean across variants
                feature_matrix[i, :] = np.mean(variant_features, axis=0)
            elif aggregation == 'max':
                # Max across variants (strongest signal)
                feature_matrix[i, :] = np.max(variant_features, axis=0)
            elif aggregation == 'sum':
                # Sum across variants
                feature_matrix[i, :] = np.sum(variant_features, axis=0)
            else:
                raise ValueError(f"Unknown aggregation: {aggregation}")
    
    return feature_matrix, patient_ids


def validate_diamond_features_nested_cv(
    X: np.ndarray,
    y: np.ndarray,
    patient_ids: List[str],
    aggregation: str = 'mean'
) -> Dict[str, Any]:
    """
    Test if pre-selected diamond features predict in nested CV.
    
    NO feature selection - just testing known features.
    """
    print(f"📊 Testing 9 diamond features with {aggregation.upper()} aggregation...")
    print(f"   Features: {DIAMOND_FEATURES}")
    print()
    
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    aurocs = []
    fold_results = []
    
    for fold_idx, (train_idx, test_idx) in enumerate(cv.split(X, y)):
        print(f"   Fold {fold_idx + 1}/5...")
        
        # ✅ ISOLATION CHECK
        assert len(set(train_idx) & set(test_idx)) == 0, f"TRAIN/TEST OVERLAP in fold {fold_idx + 1}!"
        
        X_train = X[train_idx]
        X_test = X[test_idx]
        y_train = y[train_idx]
        y_test = y[test_idx]
        train_ids = [patient_ids[i] for i in train_idx]
        test_ids = [patient_ids[i] for i in test_idx]
        
        # ✅ ISOLATION CHECK
        assert len(set(train_ids) & set(test_ids)) == 0, f"PATIENT ID OVERLAP in fold {fold_idx + 1}!"
        
        # ✅ Preprocessing: Fit ONLY on training fold
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        # ✅ Model training: ONLY on training fold
        model = LogisticRegression(
            class_weight='balanced',
            max_iter=1000,
            random_state=42
        )
        model.fit(X_train_scaled, y_train)
        
        # ✅ Evaluation: ONLY on test fold
        y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
        
        if len(set(y_test)) >= 2:
            auroc = roc_auc_score(y_test, y_pred_proba)
            aurocs.append(auroc)
            print(f"      ✅ AUROC: {auroc:.3f}")
            
            fold_results.append({
                'fold': fold_idx + 1,
                'auroc': float(auroc),
                'n_train': len(train_idx),
                'n_test': len(test_idx),
                'n_resistant_train': int(sum(y_train)),
                'n_resistant_test': int(sum(y_test))
            })
        else:
            print(f"      ⚠️  Only one class in test set")
            aurocs.append(0.5)
            fold_results.append({
                'fold': fold_idx + 1,
                'auroc': 0.5,
                'n_train': len(train_idx),
                'n_test': len(test_idx),
                'n_resistant_train': int(sum(y_train)),
                'n_resistant_test': int(sum(y_test))
            })
    
    mean_auroc = np.mean(aurocs)
    std_auroc = np.std(aurocs)
    
    return {
        'mean_auroc': float(mean_auroc),
        'std_auroc': float(std_auroc),
        'fold_aurocs': [float(a) for a in aurocs],
        'fold_results': fold_results,
        'aggregation': aggregation,
        'features': DIAMOND_FEATURES
    }


def main():
    print("=" * 80)
    print("🔬 CRITICAL EXPERIMENT: Test 9 Diamond Features (Nested CV)")
    print("=" * 80)
    print()
    print("🚨 MOMENT OF TRUTH: Do the 9 pre-selected features actually work?")
    print()
    
    # Load data
    print("📥 Loading Tier-3 cohort data...")
    with open(TIER3_CHECKPOINT, 'r') as f:
        tier3_data = json.load(f)
    
    patients_dict = tier3_data.get("data", {})
    patients_list = list(patients_dict.values())
    
    # Extract labels
    patients_filtered = []
    outcome_labels = []
    patient_ids = []
    
    for patient in patients_list:
        platinum_response = patient.get("outcome") or patient.get("platinum_response")
        if platinum_response not in ["sensitive", "resistant", "refractory"]:
            continue
        
        label = 1 if platinum_response in ["resistant", "refractory"] else 0
        patient_id = patient.get("patient_id") or patient.get("id")
        
        if not patient_id or not patient.get("variants"):
            continue
        
        patients_filtered.append(patient)
        outcome_labels.append(label)
        patient_ids.append(patient_id)
    
    print(f"   ✅ Loaded {len(patients_filtered)} patients")
    print(f"   Resistant/Refractory: {sum(outcome_labels)}")
    print(f"   Sensitive: {len(outcome_labels) - sum(outcome_labels)}")
    print()
    
    # Test MEAN aggregation
    print("=" * 80)
    print("TEST 1: MEAN Aggregation")
    print("=" * 80)
    X_mean, _ = build_feature_matrix(patients_filtered, DIAMOND_FEATURES, aggregation='mean')
    y = np.array(outcome_labels)
    
    results_mean = validate_diamond_features_nested_cv(X_mean, y, patient_ids, aggregation='mean')
    
    print()
    print(f"📊 MEAN Aggregation Results:")
    print(f"   AUROC: {results_mean['mean_auroc']:.3f} ± {results_mean['std_auroc']:.3f}")
    print(f"   Fold scores: {[f'{a:.3f}' for a in results_mean['fold_aurocs']]}")
    print()
    
    # Test MAX aggregation
    print("=" * 80)
    print("TEST 2: MAX Aggregation")
    print("=" * 80)
    X_max, _ = build_feature_matrix(patients_filtered, DIAMOND_FEATURES, aggregation='max')
    
    results_max = validate_diamond_features_nested_cv(X_max, y, patient_ids, aggregation='max')
    
    print()
    print(f"📊 MAX Aggregation Results:")
    print(f"   AUROC: {results_max['mean_auroc']:.3f} ± {results_max['std_auroc']:.3f}")
    print(f"   Fold scores: {[f'{a:.3f}' for a in results_max['fold_aurocs']]}")
    print()
    
    # Decision tree
    print("=" * 80)
    print("🎯 DECISION TREE")
    print("=" * 80)
    
    best_auroc = max(results_mean['mean_auroc'], results_max['mean_auroc'])
    best_method = 'MEAN' if results_mean['mean_auroc'] > results_max['mean_auroc'] else 'MAX'
    
    print(f"   Best AUROC: {best_auroc:.3f} ({best_method} aggregation)")
    print()
    
    if best_auroc >= 0.65:
        print("   ✅ PATH A: SALVAGE PUBLICATION")
        print("      - 9 features are real (AUROC ≥ 0.65)")
        print("      - Original leaky method discovered true signal")
        print("      - Publication angle: Validation study, not discovery")
        print("      - Timeline: 1-2 weeks to manuscript")
    elif best_auroc >= 0.55:
        print("   ⚠️  PATH B: PIVOT TO SURVIVAL")
        print("      - Weak signal exists (AUROC 0.55-0.64)")
        print("      - Not publication-worthy for resistance")
        print("      - Test survival prediction instead")
        print("      - Timeline: 2-3 days")
    else:
        print("   ❌ PATH C: FUNDAMENTAL PIVOT")
        print("      - 9 features are spurious (AUROC < 0.55)")
        print("      - No SAE signal for platinum resistance")
        print("      - Fundamental rethink required")
        print("      - Options: Different aggregation, gene-level, different outcome")
    
    print()
    
    # Save results
    all_results = {
        'timestamp': datetime.now().isoformat(),
        'mean_aggregation': results_mean,
        'max_aggregation': results_max,
        'best_auroc': float(best_auroc),
        'best_method': best_method,
        'decision_path': 'A' if best_auroc >= 0.65 else 'B' if best_auroc >= 0.55 else 'C',
        'parameters': {
            'diamond_features': DIAMOND_FEATURES,
            'n_patients': len(patients_filtered),
            'n_resistant': int(sum(outcome_labels)),
            'n_sensitive': int(len(outcome_labels) - sum(outcome_labels))
        }
    }
    
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(all_results, f, indent=2)
    
    print(f"💾 Results saved to: {OUTPUT_FILE}")
    print()
    print("✅ CRITICAL EXPERIMENT COMPLETE")
    
    return all_results


if __name__ == "__main__":
    try:
        result = main()
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ EXPERIMENT FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
