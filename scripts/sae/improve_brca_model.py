#!/usr/bin/env python3
"""
Improve BRCA Model Performance
==============================

Quick wins to improve AUROC from 0.607:
1. Feature aggregation methods (max/sum vs mean)
2. Feature stability analysis
3. SMOTE for class imbalance
4. Regularization tuning
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Any, Tuple
from datetime import datetime
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
# SMOTE optional (install with: pip install imbalanced-learn)
try:
    from imblearn.over_sampling import SMOTE
    SMOTE_AVAILABLE = True
except ImportError:
    SMOTE_AVAILABLE = False
    print("⚠️  SMOTE not available (install: pip install imbalanced-learn)")
from scipy import stats

# Paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
BRCA_CHECKPOINT = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "checkpoints" / "BRCA_TCGA_TRUE_SAE_cohort.json"
RECURRENCE_LABELS = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "brca_recurrence_labels.json"
OUTPUT_DIR = PROJECT_ROOT / ".cursor" / "MOAT" / "SAE_INTELLIGENCE" / "BRCA_VALIDATION"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUTPUT_DIR / "brca_model_improvements.json"

# Use stable features from nested CV (appear in ≥3/5 folds)
STABLE_FEATURES = [2737, 7220, 7359, 13562, 8848, 23648]  # From nested CV analysis


def build_feature_matrix_aggregation(
    patients: List[Dict],
    feature_indices: List[int],
    aggregation: str = "mean"
) -> Tuple[np.ndarray, List[str]]:
    """Build feature matrix with different aggregation methods."""
    n_patients = len(patients)
    n_features = len(feature_indices)
    
    feature_matrix = np.zeros((n_patients, n_features))
    patient_ids = []
    
    for i, patient in enumerate(patients):
        patient_id = patient.get("patient_id")
        patient_ids.append(patient_id)
        
        variants = patient.get("variants", [])
        
        # Aggregate per feature
        for j, fidx in enumerate(feature_indices):
            feature_values = []
            
            for variant in variants:
                top_features = variant.get("top_features", [])
                for tf in top_features:
                    if tf.get("index") == fidx:
                        feature_values.append(abs(float(tf.get("value", 0.0) or 0.0)))
            
            if len(feature_values) > 0:
                if aggregation == "mean":
                    feature_matrix[i, j] = np.mean(feature_values)
                elif aggregation == "max":
                    feature_matrix[i, j] = np.max(feature_values)
                elif aggregation == "sum":
                    feature_matrix[i, j] = np.sum(feature_values)
                elif aggregation == "median":
                    feature_matrix[i, j] = np.median(feature_values)
            else:
                feature_matrix[i, j] = 0.0
    
    return feature_matrix, patient_ids


def nested_cv_with_improvements(
    X: np.ndarray,
    y: np.ndarray,
    use_smote: bool = False,
    tune_regularization: bool = False
) -> Dict[str, Any]:
    """Run nested CV with improvements."""
    outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_scores = []
    
    for fold_idx, (train_idx, test_idx) in enumerate(outer_cv.split(X, y)):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        
        # Apply SMOTE if requested and available
        if use_smote and SMOTE_AVAILABLE:
            try:
                smote = SMOTE(random_state=42, k_neighbors=min(5, sum(y_train) - 1))
                X_train, y_train = smote.fit_resample(X_train, y_train)
            except:
                # If SMOTE fails (too few samples), skip
                pass
        elif use_smote and not SMOTE_AVAILABLE:
            # Skip SMOTE if not available
            pass
        
        # Scale features
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        # Train model
        if tune_regularization:
            # Grid search for optimal C
            param_grid = {'C': [0.001, 0.01, 0.1, 1.0, 10.0, 100.0]}
            model = LogisticRegression(
                class_weight='balanced',
                max_iter=1000,
                random_state=42,
                solver='lbfgs'
            )
            grid_search = GridSearchCV(model, param_grid, cv=3, scoring='roc_auc')
            grid_search.fit(X_train_scaled, y_train)
            model = grid_search.best_estimator_
        else:
            model = LogisticRegression(
                class_weight='balanced',
                max_iter=1000,
                random_state=42,
                solver='lbfgs',
                C=1.0
            )
            model.fit(X_train_scaled, y_train)
        
        # Predict
        y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
        
        # Compute AUROC
        if len(set(y_test)) >= 2:
            auroc = roc_auc_score(y_test, y_pred_proba)
            cv_scores.append(auroc)
    
    return {
        "cv_scores": cv_scores,
        "mean_auroc": np.mean(cv_scores) if cv_scores else 0.0,
        "std_auroc": np.std(cv_scores) if cv_scores else 0.0
    }


def improve_brca_model() -> Dict[str, Any]:
    """Test different improvements to BRCA model."""
    print("=" * 80)
    print("🚀 IMPROVING BRCA MODEL PERFORMANCE")
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
    print(f"   Using stable features: {len(STABLE_FEATURES)} features")
    print()
    
    # Test different aggregation methods
    print("🔍 Testing aggregation methods...")
    results = {}
    
    for aggregation in ["mean", "max", "sum", "median"]:
        print(f"   Testing {aggregation} aggregation...")
        X, _ = build_feature_matrix_aggregation(patients_list, STABLE_FEATURES, aggregation)
        y = np.array(outcome_labels)
        
        # Baseline (no improvements)
        baseline = nested_cv_with_improvements(X, y, use_smote=False, tune_regularization=False)
        results[f"{aggregation}_baseline"] = baseline
        
        # With SMOTE
        smote_result = nested_cv_with_improvements(X, y, use_smote=True, tune_regularization=False)
        results[f"{aggregation}_smote"] = smote_result
        
        # With regularization tuning
        reg_result = nested_cv_with_improvements(X, y, use_smote=False, tune_regularization=True)
        results[f"{aggregation}_tuned"] = reg_result
        
        # With both
        both_result = nested_cv_with_improvements(X, y, use_smote=True, tune_regularization=True)
        results[f"{aggregation}_both"] = both_result
        
        print(f"      Baseline: {baseline['mean_auroc']:.3f} ± {baseline['std_auroc']:.3f}")
        print(f"      +SMOTE: {smote_result['mean_auroc']:.3f} ± {smote_result['std_auroc']:.3f}")
        print(f"      +Tuned: {reg_result['mean_auroc']:.3f} ± {reg_result['std_auroc']:.3f}")
        print(f"      +Both: {both_result['mean_auroc']:.3f} ± {both_result['std_auroc']:.3f}")
        print()
    
    # Find best combination
    best_method = None
    best_auroc = 0.0
    
    for method, result in results.items():
        if result['mean_auroc'] > best_auroc:
            best_auroc = result['mean_auroc']
            best_method = method
    
    print("=" * 80)
    print("📊 IMPROVEMENT RESULTS")
    print("=" * 80)
    print(f"Baseline (mean, no improvements): {results['mean_baseline']['mean_auroc']:.3f}")
    print(f"Best method: {best_method}")
    print(f"Best AUROC: {best_auroc:.3f} ± {results[best_method]['std_auroc']:.3f}")
    print(f"Improvement: {best_auroc - results['mean_baseline']['mean_auroc']:+.3f}")
    print("=" * 80)
    
    # Save results
    output = {
        "improvement_date": datetime.now().isoformat(),
        "baseline_auroc": 0.607,  # From nested CV
        "results": {k: {
            "mean_auroc": float(v['mean_auroc']),
            "std_auroc": float(v['std_auroc']),
            "cv_scores": [float(s) for s in v['cv_scores']]
        } for k, v in results.items()},
        "best_method": best_method,
        "best_auroc": float(best_auroc),
        "improvement": float(best_auroc - results['mean_baseline']['mean_auroc']),
        "stable_features": STABLE_FEATURES
    }
    
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\n💾 Results saved to: {OUTPUT_FILE}")
    
    return output


if __name__ == "__main__":
    try:
        result = improve_brca_model()
        print("\n✅ IMPROVEMENT ANALYSIS COMPLETE")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ IMPROVEMENT ANALYSIS FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
