#!/usr/bin/env python3
"""
Build BRCA Multi-Pathway Model - Deliverable #7
===============================================

Mission: Build final validation model using BRCA-specific features
Objective: Train composite model and compare to Oncotype DX baseline
Timeline: 4 hours
Status: EXECUTION READY

Based on: SAE_BRCA_NEXT_10_DELIVERABLES.md (Deliverable #7)
Features: 17 BRCA-specific features (from pathway comparison)
"""

import json
import sys
import pickle
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler
from scipy import stats

# Paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
BRCA_CHECKPOINT = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "checkpoints" / "BRCA_TCGA_TRUE_SAE_cohort.json"
RECURRENCE_LABELS = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "brca_recurrence_labels.json"
PATHWAY_COMPARISON = PROJECT_ROOT / ".cursor" / "MOAT" / "SAE_INTELLIGENCE" / "BRCA_VALIDATION" / "pathway_comparison.json"
OUTPUT_DIR = PROJECT_ROOT / ".cursor" / "MOAT" / "SAE_INTELLIGENCE" / "BRCA_VALIDATION"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
MODEL_OUTPUT = OUTPUT_DIR / "brca_sae_model.pkl"
RESULTS_OUTPUT = OUTPUT_DIR / "brca_sae_validation_results.json"

# Oncotype DX baseline (from literature)
ONCOTYPE_DX_BASELINE_AUROC = 0.65  # Typical performance in validation studies


def compute_pathway_score_mean(patient_data: Dict, feature_indices: List[int]) -> float:
    """Compute pathway score using MEAN aggregation."""
    feature_sums = {fidx: 0.0 for fidx in feature_indices}
    
    for variant in patient_data.get("variants", []):
        for tf in variant.get("top_features", []):
            fidx = tf.get("index")
            if fidx in feature_sums:
                feature_sums[fidx] += abs(float(tf.get("value", 0.0) or 0.0))
    
    return sum(feature_sums.values()) / len(feature_indices) if feature_indices else 0.0


def build_feature_matrix(patients: List[Dict], feature_indices: List[int]) -> np.ndarray:
    """Build feature matrix for model training."""
    n_patients = len(patients)
    n_features = len(feature_indices)
    
    feature_matrix = np.zeros((n_patients, n_features))
    
    for i, patient in enumerate(patients):
        # Compute score for each feature individually
        for j, fidx in enumerate(feature_indices):
            feature_sum = 0.0
            variant_count = 0
            
            for variant in patient.get("variants", []):
                variant_count += 1
                for tf in variant.get("top_features", []):
                    if tf.get("index") == fidx:
                        feature_sum += abs(float(tf.get("value", 0.0) or 0.0))
            
            # Mean activation per variant
            if variant_count > 0:
                feature_matrix[i, j] = feature_sum / variant_count
    
    return feature_matrix


def build_multipathway_model() -> Dict[str, Any]:
    """Build BRCA multi-pathway model using BRCA-specific features."""
    print("=" * 80)
    print("🏗️  BUILDING BRCA MULTI-PATHWAY MODEL - Deliverable #7")
    print("=" * 80)
    print()
    
    # Step 1: Load BRCA SAE features
    print("📥 Step 1: Loading BRCA SAE features...")
    if not BRCA_CHECKPOINT.exists():
        print(f"❌ BRCA checkpoint not found: {BRCA_CHECKPOINT}")
        sys.exit(1)
    
    with open(BRCA_CHECKPOINT, 'r') as f:
        brca_sae = json.load(f)
    
    patients_dict = brca_sae.get("data", {})
    print(f"   ✅ Loaded {len(patients_dict)} patients")
    
    # Step 2: Load recurrence labels
    print()
    print("📥 Step 2: Loading recurrence labels...")
    if not RECURRENCE_LABELS.exists():
        print(f"❌ Recurrence labels not found: {RECURRENCE_LABELS}")
        sys.exit(1)
    
    with open(RECURRENCE_LABELS, 'r') as f:
        labels_data = json.load(f)
    
    outcomes = labels_data.get("outcomes", {})
    print(f"   ✅ Loaded {len(outcomes)} outcome records")
    
    # Step 3: Load pathway comparison to get BRCA feature indices
    print()
    print("📥 Step 3: Loading pathway comparison results...")
    if not PATHWAY_COMPARISON.exists():
        print(f"❌ Pathway comparison not found: {PATHWAY_COMPARISON}")
        sys.exit(1)
    
    with open(PATHWAY_COMPARISON, 'r') as f:
        pathway_data = json.load(f)
    
    brca_feature_indices = pathway_data.get("brca_specific", {}).get("features", [])
    if not brca_feature_indices:
        print("❌ No BRCA feature indices found in pathway comparison")
        sys.exit(1)
    
    print(f"   ✅ Using {len(brca_feature_indices)} BRCA-specific features")
    
    # Step 4: Prepare data for model training
    print()
    print("🔍 Step 4: Preparing data for model training...")
    
    patients_list = []
    outcome_labels = []
    patient_ids = []
    
    for patient_id, patient_data in patients_dict.items():
        outcome = outcomes.get(patient_id, {})
        recurrence = outcome.get("recurrence")
        
        if recurrence is None:
            continue
        
        patient_data["patient_id"] = patient_id
        patients_list.append(patient_data)
        outcome_labels.append(1 if recurrence else 0)
        patient_ids.append(patient_id)
    
    print(f"   ✅ Prepared {len(patients_list)} patients with labels")
    print(f"   Recurrence: {sum(outcome_labels)} True, {len(outcome_labels) - sum(outcome_labels)} False")
    
    # Step 5: Build feature matrix
    print()
    print("🔍 Step 5: Building feature matrix...")
    feature_matrix = build_feature_matrix(patients_list, brca_feature_indices)
    print(f"   ✅ Feature matrix: {feature_matrix.shape}")
    
    # Step 6: Train model with cross-validation
    print()
    print("📊 Step 6: Training model with cross-validation...")
    
    X = feature_matrix
    y = np.array(outcome_labels)
    
    # Standardize features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Train logistic regression with balanced class weights
    model = LogisticRegression(
        class_weight='balanced',
        max_iter=1000,
        random_state=42,
        solver='lbfgs'
    )
    
    # 5-fold stratified cross-validation
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_scores = cross_val_score(model, X_scaled, y, cv=cv, scoring='roc_auc')
    
    print(f"   ✅ CV AUROC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
    print(f"   CV scores: {[f'{s:.3f}' for s in cv_scores]}")
    
    # Train on full dataset
    model.fit(X_scaled, y)
    
    # Step 7: Compute predictions and final AUROC
    print()
    print("📊 Step 7: Computing final model performance...")
    
    y_pred_proba = model.predict_proba(X_scaled)[:, 1]
    final_auroc = roc_auc_score(y, y_pred_proba)
    
    # Bootstrap CI
    n = len(y)
    bootstrap_aurocs = []
    for _ in range(1000):
        indices = np.random.choice(n, n, replace=True)
        y_boot = y[indices]
        y_pred_boot = y_pred_proba[indices]
        if len(set(y_boot)) >= 2:
            auroc_boot = roc_auc_score(y_boot, y_pred_boot)
            bootstrap_aurocs.append(auroc_boot)
    
    ci_lower = np.percentile(bootstrap_aurocs, 2.5)
    ci_upper = np.percentile(bootstrap_aurocs, 97.5)
    
    print(f"   ✅ Final AUROC: {final_auroc:.3f} (95% CI: {ci_lower:.3f}-{ci_upper:.3f})")
    
    # Step 8: Compare to Oncotype DX baseline
    print()
    print("📊 Step 8: Comparing to Oncotype DX baseline...")
    
    improvement = final_auroc - ONCOTYPE_DX_BASELINE_AUROC
    improvement_pct = (improvement / ONCOTYPE_DX_BASELINE_AUROC) * 100
    
    print(f"   Oncotype DX baseline: {ONCOTYPE_DX_BASELINE_AUROC:.3f}")
    print(f"   SAE model: {final_auroc:.3f}")
    print(f"   Improvement: {improvement:+.3f} ({improvement_pct:+.1f}%)")
    
    # Step 9: Save model
    print()
    print("💾 Step 9: Saving model...")
    
    model_data = {
        "model": model,
        "scaler": scaler,
        "feature_indices": brca_feature_indices,
        "n_features": len(brca_feature_indices),
        "model_type": "logistic_regression",
        "class_weight": "balanced"
    }
    
    with open(MODEL_OUTPUT, 'wb') as f:
        pickle.dump(model_data, f)
    
    print(f"   ✅ Saved model to: {MODEL_OUTPUT}")
    
    # Step 10: Save validation results
    print()
    print("💾 Step 10: Saving validation results...")
    
    results = {
        "validation_date": datetime.now().isoformat(),
        "model_type": "brca_specific_multipathway",
        "n_features": len(brca_feature_indices),
        "feature_indices": brca_feature_indices,
        "n_patients": len(patients_list),
        "recurrence_distribution": {
            "true": sum(outcome_labels),
            "false": len(outcome_labels) - sum(outcome_labels)
        },
        "performance": {
            "cv_auroc_mean": float(cv_scores.mean()),
            "cv_auroc_std": float(cv_scores.std()),
            "cv_scores": [float(s) for s in cv_scores],
            "final_auroc": float(final_auroc),
            "ci_lower": float(ci_lower),
            "ci_upper": float(ci_upper)
        },
        "comparison": {
            "oncotype_dx_auroc": ONCOTYPE_DX_BASELINE_AUROC,
            "sae_auroc": float(final_auroc),
            "improvement": float(improvement),
            "improvement_pct": float(improvement_pct)
        },
        "provenance": {
            "brca_checkpoint": str(BRCA_CHECKPOINT),
            "recurrence_labels": str(RECURRENCE_LABELS),
            "pathway_comparison": str(PATHWAY_COMPARISON),
            "model_file": str(MODEL_OUTPUT)
        }
    }
    
    with open(RESULTS_OUTPUT, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"   ✅ Saved results to: {RESULTS_OUTPUT}")
    
    # Step 11: Print summary
    print()
    print("=" * 80)
    print("📊 MODEL VALIDATION SUMMARY")
    print("=" * 80)
    print(f"Model Type: BRCA-specific multi-pathway (17 features)")
    print(f"CV AUROC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
    print(f"Final AUROC: {final_auroc:.3f} (95% CI: {ci_lower:.3f}-{ci_upper:.3f})")
    print()
    print("Comparison to Oncotype DX:")
    print(f"  Oncotype DX: {ONCOTYPE_DX_BASELINE_AUROC:.3f}")
    print(f"  SAE Model: {final_auroc:.3f}")
    print(f"  Improvement: {improvement:+.3f} ({improvement_pct:+.1f}%)")
    print()
    if improvement > 0:
        print("✅ SUCCESS: SAE model outperforms Oncotype DX baseline")
    else:
        print("⚠️  WARNING: SAE model does not outperform Oncotype DX baseline")
    print("=" * 80)
    
    return results


if __name__ == "__main__":
    try:
        result = build_multipathway_model()
        print("\n✅ MODEL BUILDING COMPLETE")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ MODEL BUILDING FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
