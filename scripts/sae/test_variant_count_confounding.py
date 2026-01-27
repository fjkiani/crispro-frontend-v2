#!/usr/bin/env python3
"""
Test Variant Count Confounding - Mars's Red Flag Check
======================================================

Mission: Verify if sum aggregation improvement is spurious due to variant count
Objective: Check if variant count alone predicts recurrence
Status: CRITICAL AUDIT
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Any
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from scipy import stats

# Paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
BRCA_CHECKPOINT = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "checkpoints" / "BRCA_TCGA_TRUE_SAE_cohort.json"
RECURRENCE_LABELS = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "brca_recurrence_labels.json"
OUTPUT_DIR = PROJECT_ROOT / ".cursor" / "MOAT" / "SAE_INTELLIGENCE" / "BRCA_VALIDATION"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUTPUT_DIR / "variant_count_confounding_analysis.json"


def test_variant_count_confounding() -> Dict[str, Any]:
    """Test if variant count is confounding sum aggregation results."""
    print("=" * 80)
    print("🚨 MARS'S RED FLAG CHECK: Variant Count Confounding")
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
    variant_counts = []
    sum_features_list = []
    mean_features_list = []
    
    for patient_id, patient_data in patients_dict.items():
        outcome = outcomes.get(patient_id, {})
        recurrence = outcome.get("recurrence")
        
        if recurrence is None:
            continue
        
        variants = patient_data.get("variants", [])
        variant_count = len(variants)
        variant_counts.append(variant_count)
        
        # Compute sum and mean for stable features
        stable_features = [2737, 7220, 7359, 13562, 8848, 23648]
        sum_features = np.zeros(len(stable_features))
        mean_features = np.zeros(len(stable_features))
        
        for variant in variants:
            for j, fidx in enumerate(stable_features):
                for tf in variant.get("top_features", []):
                    if tf.get("index") == fidx:
                        val = abs(float(tf.get("value", 0.0) or 0.0))
                        sum_features[j] += val
                        mean_features[j] += val
        
        # Normalize mean by variant count
        if variant_count > 0:
            mean_features = mean_features / variant_count
        
        # Use mean of all stable features
        sum_score = np.mean(sum_features)
        mean_score = np.mean(mean_features)
        
        sum_features_list.append(sum_score)
        mean_features_list.append(mean_score)
        
        patients_list.append(patient_data)
        outcome_labels.append(1 if recurrence else 0)
    
    variant_counts = np.array(variant_counts)
    sum_features_array = np.array(sum_features_list)
    mean_features_array = np.array(mean_features_list)
    y = np.array(outcome_labels)
    
    print(f"   ✅ Loaded {len(patients_list)} patients")
    print(f"   Recurrence: {sum(y)} True, {len(y) - sum(y)} False")
    print()
    
    # Test 1: Variant count correlation
    print("🔍 Test 1: Variant count vs recurrence correlation...")
    variant_count_corr, variant_count_p = stats.pearsonr(variant_counts, y)
    print(f"   Correlation: r = {variant_count_corr:.3f}, p = {variant_count_p:.4f}")
    
    if abs(variant_count_corr) > 0.3:
        print(f"   🚨 RED FLAG: High correlation (|r| > 0.3)")
        print(f"   ⚠️  Variant count is confounded with recurrence")
    else:
        print(f"   ✅ Low correlation - variant count not a major confounder")
    print()
    
    # Test 2: Variant count alone as predictor
    print("🔍 Test 2: Variant count alone as predictor...")
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    variant_count_scores = []
    
    for train_idx, test_idx in cv.split(variant_counts.reshape(-1, 1), y):
        X_train = variant_counts[train_idx].reshape(-1, 1)
        X_test = variant_counts[test_idx].reshape(-1, 1)
        y_train = y[train_idx]
        y_test = y[test_idx]
        
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        model = LogisticRegression(class_weight='balanced', max_iter=1000, random_state=42)
        model.fit(X_train_scaled, y_train)
        y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
        
        if len(set(y_test)) >= 2:
            auroc = roc_auc_score(y_test, y_pred_proba)
            variant_count_scores.append(auroc)
    
    variant_count_auroc = np.mean(variant_count_scores) if variant_count_scores else 0.0
    variant_count_std = np.std(variant_count_scores) if variant_count_scores else 0.0
    
    print(f"   Variant count AUROC: {variant_count_auroc:.3f} ± {variant_count_std:.3f}")
    
    if variant_count_auroc > 0.70:
        print(f"   🚨 RED FLAG: Variant count alone predicts recurrence (AUROC > 0.70)")
        print(f"   ⚠️  Sum aggregation is likely just learning variant count")
    elif variant_count_auroc > 0.60:
        print(f"   ⚠️  WARNING: Variant count has moderate predictive power")
        print(f"   ⚠️  Sum aggregation may be partially confounded")
    else:
        print(f"   ✅ Variant count has low predictive power - not a major confounder")
    print()
    
    # Test 3: Sum vs Mean with variant count controlled
    print("🔍 Test 3: Sum vs Mean (controlling for variant count)...")
    
    # Model A: Sum features only
    sum_scores = []
    for train_idx, test_idx in cv.split(sum_features_array.reshape(-1, 1), y):
        X_train = sum_features_array[train_idx].reshape(-1, 1)
        X_test = sum_features_array[test_idx].reshape(-1, 1)
        y_train = y[train_idx]
        y_test = y[test_idx]
        
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        model = LogisticRegression(class_weight='balanced', max_iter=1000, random_state=42)
        model.fit(X_train_scaled, y_train)
        y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
        
        if len(set(y_test)) >= 2:
            auroc = roc_auc_score(y_test, y_pred_proba)
            sum_scores.append(auroc)
    
    sum_auroc = np.mean(sum_scores) if sum_scores else 0.0
    
    # Model B: Mean features only
    mean_scores = []
    for train_idx, test_idx in cv.split(mean_features_array.reshape(-1, 1), y):
        X_train = mean_features_array[train_idx].reshape(-1, 1)
        X_test = mean_features_array[test_idx].reshape(-1, 1)
        y_train = y[train_idx]
        y_test = y[test_idx]
        
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        model = LogisticRegression(class_weight='balanced', max_iter=1000, random_state=42)
        model.fit(X_train_scaled, y_train)
        y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
        
        if len(set(y_test)) >= 2:
            auroc = roc_auc_score(y_test, y_pred_proba)
            mean_scores.append(auroc)
    
    mean_auroc = np.mean(mean_scores) if mean_scores else 0.0
    
    # Model C: Sum + Variant count (to see if sum adds value)
    sum_with_count_scores = []
    for train_idx, test_idx in cv.split(np.column_stack([sum_features_array, variant_counts]), y):
        X_train = np.column_stack([sum_features_array[train_idx], variant_counts[train_idx]])
        X_test = np.column_stack([sum_features_array[test_idx], variant_counts[test_idx]])
        y_train = y[train_idx]
        y_test = y[test_idx]
        
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        model = LogisticRegression(class_weight='balanced', max_iter=1000, random_state=42)
        model.fit(X_train_scaled, y_train)
        y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
        
        if len(set(y_test)) >= 2:
            auroc = roc_auc_score(y_test, y_pred_proba)
            sum_with_count_scores.append(auroc)
    
    sum_with_count_auroc = np.mean(sum_with_count_scores) if sum_with_count_scores else 0.0
    
    print(f"   Sum features only: {sum_auroc:.3f}")
    print(f"   Mean features only: {mean_auroc:.3f}")
    print(f"   Variant count only: {variant_count_auroc:.3f}")
    print(f"   Sum + Variant count: {sum_with_count_auroc:.3f}")
    
    if abs(sum_auroc - variant_count_auroc) < 0.05:
        print(f"   🚨 RED FLAG: Sum AUROC ≈ Variant count AUROC")
        print(f"   ⚠️  Sum aggregation is just variant count in disguise")
    elif sum_with_count_auroc - variant_count_auroc < 0.02:
        print(f"   🚨 RED FLAG: Sum adds < 2 pp over variant count alone")
        print(f"   ⚠️  Sum aggregation provides minimal additional signal")
    else:
        print(f"   ✅ Sum adds significant value over variant count")
    print()
    
    # Test 4: Correlation between sum and variant count
    print("🔍 Test 4: Sum features vs variant count correlation...")
    sum_count_corr, sum_count_p = stats.pearsonr(sum_features_array, variant_counts)
    print(f"   Correlation: r = {sum_count_corr:.3f}, p = {sum_count_p:.4f}")
    
    if abs(sum_count_corr) > 0.7:
        print(f"   🚨 RED FLAG: Sum features highly correlated with variant count (|r| > 0.7)")
        print(f"   ⚠️  Sum is essentially a proxy for variant count")
    elif abs(sum_count_corr) > 0.5:
        print(f"   ⚠️  WARNING: Moderate correlation (|r| > 0.5)")
        print(f"   ⚠️  Sum partially confounded with variant count")
    else:
        print(f"   ✅ Low correlation - sum captures independent signal")
    print()
    
    # Summary
    print("=" * 80)
    print("📊 CONFOUNDING ANALYSIS SUMMARY")
    print("=" * 80)
    
    results = {
        "variant_count_correlation": float(variant_count_corr),
        "variant_count_p_value": float(variant_count_p),
        "variant_count_auroc": float(variant_count_auroc),
        "variant_count_std": float(variant_count_std),
        "sum_features_auroc": float(sum_auroc),
        "mean_features_auroc": float(mean_auroc),
        "sum_with_count_auroc": float(sum_with_count_auroc),
        "sum_vs_count_correlation": float(sum_count_corr),
        "sum_vs_count_p_value": float(sum_count_p),
        "confounding_assessment": {}
    }
    
    # Assessment
    if abs(variant_count_corr) > 0.3 or variant_count_auroc > 0.70:
        results["confounding_assessment"]["status"] = "HIGH_CONFOUNDING"
        results["confounding_assessment"]["verdict"] = "Sum aggregation is confounded by variant count"
        results["confounding_assessment"]["recommendation"] = "Use mean aggregation (unbiased)"
        print("🚨 VERDICT: HIGH CONFOUNDING DETECTED")
        print("   Sum aggregation is confounded by variant count")
        print("   Recommendation: Use mean aggregation (unbiased)")
    elif abs(sum_count_corr) > 0.5:
        results["confounding_assessment"]["status"] = "MODERATE_CONFOUNDING"
        results["confounding_assessment"]["verdict"] = "Sum aggregation partially confounded"
        results["confounding_assessment"]["recommendation"] = "Use mean aggregation or control for variant count"
        print("⚠️  VERDICT: MODERATE CONFOUNDING DETECTED")
        print("   Sum aggregation partially confounded by variant count")
        print("   Recommendation: Use mean aggregation or control for variant count")
    else:
        results["confounding_assessment"]["status"] = "LOW_CONFOUNDING"
        results["confounding_assessment"]["verdict"] = "Sum aggregation appears valid"
        results["confounding_assessment"]["recommendation"] = "Sum aggregation may be acceptable"
        print("✅ VERDICT: LOW CONFOUNDING")
        print("   Sum aggregation appears valid")
    
    print()
    print("Key Findings:")
    print(f"  Variant count correlation: r = {variant_count_corr:.3f}")
    print(f"  Variant count AUROC: {variant_count_auroc:.3f}")
    print(f"  Sum vs count correlation: r = {sum_count_corr:.3f}")
    print(f"  Sum AUROC: {sum_auroc:.3f}")
    print(f"  Mean AUROC: {mean_auroc:.3f}")
    print("=" * 80)
    
    # Save results
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n💾 Results saved to: {OUTPUT_FILE}")
    
    return results


if __name__ == "__main__":
    try:
        result = test_variant_count_confounding()
        print("\n✅ CONFOUNDING ANALYSIS COMPLETE")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ ANALYSIS FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
