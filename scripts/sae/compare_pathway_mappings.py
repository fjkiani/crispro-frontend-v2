#!/usr/bin/env python3
"""
Compare Pathway Mappings - Deliverable #6
==========================================

Mission: Compare OV DDR transfer vs BRCA-specific features
Objective: Decision point for best pathway(s) to use
Timeline: 2 hours
Status: EXECUTION READY

Based on: SAE_BRCA_NEXT_10_DELIVERABLES.md (Deliverable #6)
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime
import numpy as np
from sklearn.metrics import roc_auc_score
from scipy import stats

# Paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
BRCA_CHECKPOINT = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "checkpoints" / "BRCA_TCGA_TRUE_SAE_cohort.json"
RECURRENCE_LABELS = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "brca_recurrence_labels.json"
OV_DDR_RESULTS = PROJECT_ROOT / ".cursor" / "MOAT" / "SAE_INTELLIGENCE" / "BRCA_VALIDATION" / "ov_ddr_transfer_results.json"
BRCA_DIAMONDS = PROJECT_ROOT / ".cursor" / "MOAT" / "SAE_INTELLIGENCE" / "BRCA_VALIDATION" / "brca_diamond_features.json"
OUTPUT_DIR = PROJECT_ROOT / ".cursor" / "MOAT" / "SAE_INTELLIGENCE" / "BRCA_VALIDATION"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUTPUT_DIR / "pathway_comparison.json"

# OV DDR diamond features
OV_DDR_FEATURES = [1407, 6020, 9738, 12893, 16337, 22868, 26220, 27607, 31362]


def compute_pathway_score_mean(patient_data: Dict, feature_indices: List[int]) -> float:
    """
    Compute pathway score using MEAN aggregation (matches manuscript).
    
    Method:
    1. Sum all feature values across all variants
    2. Divide by number of features
    3. Missing features = 0.0 (handled automatically)
    """
    feature_sums = {fidx: 0.0 for fidx in feature_indices}
    
    for variant in patient_data.get("variants", []):
        for tf in variant.get("top_features", []):
            fidx = tf.get("index")
            if fidx in feature_sums:
                # Sum absolute values
                feature_sums[fidx] += abs(float(tf.get("value", 0.0) or 0.0))
    
    # Mean of features (missing features contribute 0.0)
    return sum(feature_sums.values()) / len(feature_indices) if feature_indices else 0.0


def compute_auroc(y_true: List[bool], y_scores: List[float]) -> Tuple[float, Dict[str, Any]]:
    """Compute AUROC with bootstrap confidence intervals."""
    if len(set(y_true)) < 2:
        return 0.0, {"error": "Insufficient outcome diversity"}
    
    auroc = roc_auc_score(y_true, y_scores)
    
    # Bootstrap CI (1000 iterations)
    n = len(y_true)
    bootstrap_aurocs = []
    for _ in range(1000):
        indices = np.random.choice(n, n, replace=True)
        y_true_boot = [y_true[i] for i in indices]
        y_scores_boot = [y_scores[i] for i in indices]
        if len(set(y_true_boot)) >= 2:
            auroc_boot = roc_auc_score(y_true_boot, y_scores_boot)
            bootstrap_aurocs.append(auroc_boot)
    
    if bootstrap_aurocs:
        ci_lower = np.percentile(bootstrap_aurocs, 2.5)
        ci_upper = np.percentile(bootstrap_aurocs, 97.5)
    else:
        ci_lower = ci_upper = auroc
    
    return auroc, {
        "auroc": float(auroc),
        "ci_lower": float(ci_lower),
        "ci_upper": float(ci_upper),
        "n_bootstrap": len(bootstrap_aurocs)
    }


def compare_pathway_mappings() -> Dict[str, Any]:
    """Compare OV DDR transfer vs BRCA-specific features."""
    print("=" * 80)
    print("🔍 COMPARING PATHWAY MAPPINGS - Deliverable #6")
    print("=" * 80)
    print()
    
    # Step 1: Load BRCA SAE features
    print("📥 Step 1: Loading BRCA SAE features...")
    if not BRCA_CHECKPOINT.exists():
        print(f"❌ BRCA checkpoint not found: {BRCA_CHECKPOINT}")
        sys.exit(1)
    
    with open(BRCA_CHECKPOINT, 'r') as f:
        brca_sae = json.load(f)
    
    patients = brca_sae.get("data", {})
    print(f"   ✅ Loaded {len(patients)} patients")
    
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
    
    # Step 3: Load OV DDR results (for reference)
    print()
    print("📥 Step 3: Loading OV DDR transfer results...")
    if not OV_DDR_RESULTS.exists():
        print(f"⚠️  OV DDR results not found: {OV_DDR_RESULTS}")
        ov_ddr_auroc = None
    else:
        with open(OV_DDR_RESULTS, 'r') as f:
            ov_ddr_data = json.load(f)
        ov_ddr_auroc = ov_ddr_data.get("result", {}).get("ov_ddr_auroc")
        print(f"   ✅ OV DDR AUROC: {ov_ddr_auroc:.3f}")
    
    # Step 4: Load BRCA diamond features
    print()
    print("📥 Step 4: Loading BRCA diamond features...")
    if not BRCA_DIAMONDS.exists():
        print(f"❌ BRCA diamonds not found: {BRCA_DIAMONDS}")
        sys.exit(1)
    
    with open(BRCA_DIAMONDS, 'r') as f:
        brca_diamonds_data = json.load(f)
    
    # Extract BRCA feature indices
    # Note: pathway_diamonds are empty (no gene info), so use all significant features
    n_significant = brca_diamonds_data.get("n_significant_features", 0)
    print(f"   ✅ Found {n_significant} significant BRCA features")
    
    # We need to extract the actual feature indices from the analysis
    # For now, we'll compute BRCA score using all features that were significant
    # But we need to re-compute which features those are, or load them from a detailed output
    
    # Step 5: Extract BRCA feature indices by re-running correlation (quick)
    print()
    print("🔍 Step 7: Extracting BRCA feature indices...")
    
    # Inline correlation computation (avoid import issues)
    def build_feature_matrix_inline(patients_list):
        n_patients = len(patients_list)
        n_features = 32768
        feature_matrix = np.zeros((n_patients, n_features))
        patient_ids_list = []
        
        for i, patient in enumerate(patients_list):
            patient_id = patient.get("patient_id")
            patient_ids_list.append(patient_id)
            
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
        
        return feature_matrix, patient_ids_list
    
    def compute_pearson_inline(feature_matrix, outcome_vector):
        n_features = feature_matrix.shape[1]
        pearson_r = np.zeros(n_features)
        pearson_p = np.ones(n_features)
        
        for i in range(n_features):
            if i % 5000 == 0:
                print(f"      Computing Pearson for feature {i}/{n_features}...")
            feature_values = feature_matrix[:, i]
            if np.std(feature_values) == 0:
                continue
            r, p = stats.pearsonr(feature_values, outcome_vector)
            pearson_r[i] = r
            pearson_p[i] = p
        
        return pearson_r, pearson_p
    
    def compute_cohen_d_inline(feature_matrix, outcome_labels):
        n_features = feature_matrix.shape[1]
        cohen_d_values = np.zeros(n_features)
        recurrence_mask = np.array(outcome_labels)
        
        for i in range(n_features):
            if i % 5000 == 0:
                print(f"      Computing Cohen's d for feature {i}/{n_features}...")
            feature_values = feature_matrix[:, i]
            recurrence_values = feature_values[recurrence_mask]
            no_recurrence_values = feature_values[~recurrence_mask]
            
            if len(recurrence_values) == 0 or len(no_recurrence_values) == 0:
                continue
            
            mean_diff = np.mean(recurrence_values) - np.mean(no_recurrence_values)
            pooled_std = np.sqrt((np.var(recurrence_values) + np.var(no_recurrence_values)) / 2)
            
            if pooled_std == 0:
                continue
            
            cohen_d_values[i] = abs(mean_diff / pooled_std)
        
        return cohen_d_values
    
    # Convert patients to list format
    patients_list = []
    for patient_id, patient_data in patients.items():
        patient_data["patient_id"] = patient_id
        patients_list.append(patient_data)
    
    # Build feature matrix and outcome vector
    print("   Building feature matrix...")
    feature_matrix, patient_ids_list = build_feature_matrix_inline(patients_list)
    
    outcome_vector_list = []
    outcome_labels_list = []
    valid_indices = []
    
    for i, pid in enumerate(patient_ids_list):
        outcome = outcomes.get(pid, {})
        recurrence = outcome.get("recurrence")
        if recurrence is not None:
            outcome_vector_list.append(1.0 if recurrence else 0.0)
            outcome_labels_list.append(recurrence)
            valid_indices.append(i)
    
    feature_matrix = feature_matrix[valid_indices, :]
    outcome_vector_array = np.array(outcome_vector_list)
    
    # Quick correlation (only for significant features check)
    print("   Computing correlations for BRCA features...")
    pearson_r_brca, pearson_p_brca = compute_pearson_inline(feature_matrix, outcome_vector_array)
    cohen_d_brca = compute_cohen_d_inline(feature_matrix, outcome_labels_list)
    
    # Extract significant feature indices
    brca_feature_indices = []
    for i in range(feature_matrix.shape[1]):
        if pearson_p_brca[i] < 0.05 and cohen_d_brca[i] >= 0.5:
            brca_feature_indices.append(i)
    
    # Sort by Cohen's d
    brca_features_with_stats = [
        (i, cohen_d_brca[i], pearson_r_brca[i], pearson_p_brca[i])
        for i in brca_feature_indices
    ]
    brca_features_with_stats.sort(key=lambda x: x[1], reverse=True)
    brca_feature_indices = [f[0] for f in brca_features_with_stats]
    
    print(f"   ✅ Found {len(brca_feature_indices)} BRCA feature indices")
    
    # Step 6: Compute pathway scores for each patient
    print()
    print("🔍 Step 6: Computing pathway scores...")
    
    # OV DDR scores
    ov_ddr_scores = []
    brca_scores = []
    recurrence_labels = []
    patient_ids = []
    
    for patient_id, patient_data in patients.items():
        outcome = outcomes.get(patient_id, {})
        recurrence = outcome.get("recurrence")
        
        if recurrence is None:
            continue
        
        # Compute OV DDR score
        ov_ddr_score = compute_pathway_score_mean(patient_data, OV_DDR_FEATURES)
        ov_ddr_scores.append(ov_ddr_score)
        
        # Compute BRCA score using extracted indices
        brca_score = compute_pathway_score_mean(patient_data, brca_feature_indices)
        brca_scores.append(brca_score)
        
        recurrence_labels.append(recurrence)
        patient_ids.append(patient_id)
    
    # Step 7: Compute AUROC for both pathways
    print()
    print("📊 Step 7: Computing AUROC for both pathways...")
    
    # OV DDR AUROC
    ov_ddr_auroc_computed, ov_ddr_stats = compute_auroc(recurrence_labels, ov_ddr_scores)
    print(f"   ✅ OV DDR AUROC: {ov_ddr_auroc_computed:.3f} (95% CI: {ov_ddr_stats['ci_lower']:.3f}-{ov_ddr_stats['ci_upper']:.3f})")
    
    # BRCA AUROC
    brca_auroc_computed, brca_stats = compute_auroc(recurrence_labels, brca_scores)
    print(f"   ✅ BRCA AUROC: {brca_auroc_computed:.3f} (95% CI: {brca_stats['ci_lower']:.3f}-{brca_stats['ci_upper']:.3f})")
    
    # Step 8: Build comparison results
    print()
    print("💾 Step 8: Building comparison results...")
    
    # Determine best pathway
    if brca_auroc_computed > ov_ddr_auroc_computed:
        best_single_pathway = "brca_specific"
        recommendation = "use_brca_specific"
    else:
        best_single_pathway = "ov_ddr"
        recommendation = "use_ov_ddr" if ov_ddr_auroc_computed >= 0.60 else "run_brca_specific"
    
    result = {
        "comparison_date": datetime.now().isoformat(),
        "ov_ddr": {
            "auroc": float(ov_ddr_auroc_computed),
            "ci_lower": ov_ddr_stats["ci_lower"],
            "ci_upper": ov_ddr_stats["ci_upper"],
            "features": OV_DDR_FEATURES,
            "n_features": len(OV_DDR_FEATURES),
            "interpretation": "poor" if ov_ddr_auroc_computed < 0.60 else "moderate" if ov_ddr_auroc_computed < 0.70 else "success"
        },
        "brca_specific": {
            "auroc": float(brca_auroc_computed),
            "ci_lower": brca_stats["ci_lower"],
            "ci_upper": brca_stats["ci_upper"],
            "features": brca_feature_indices,
            "n_features": len(brca_feature_indices),
            "interpretation": "poor" if brca_auroc_computed < 0.60 else "moderate" if brca_auroc_computed < 0.70 else "success"
        },
        "best_single_pathway": best_single_pathway,
        "recommendation": recommendation,
        "n_patients": len(recurrence_labels),
        "recurrence_distribution": {
            "true": sum(recurrence_labels),
            "false": len(recurrence_labels) - sum(recurrence_labels)
        }
    }
    
    # Step 9: Save results
    print()
    print("💾 Step 9: Saving results...")
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(result, f, indent=2)
    
    print(f"   ✅ Saved to: {OUTPUT_FILE}")
    
    # Step 10: Print summary
    print()
    print("=" * 80)
    print("📊 PATHWAY COMPARISON SUMMARY")
    print("=" * 80)
    print(f"OV DDR AUROC: {ov_ddr_auroc_computed:.3f}")
    print(f"95% CI: [{ov_ddr_stats['ci_lower']:.3f}, {ov_ddr_stats['ci_upper']:.3f}]")
    print(f"Interpretation: {result['ov_ddr']['interpretation']}")
    print()
    print(f"BRCA-specific AUROC: {brca_auroc_computed:.3f}")
    print(f"95% CI: [{brca_stats['ci_lower']:.3f}, {brca_stats['ci_upper']:.3f}]")
    print(f"Interpretation: {result['brca_specific']['interpretation']}")
    print(f"Features: {len(brca_feature_indices)}")
    print()
    print(f"Best single pathway: {best_single_pathway}")
    print(f"Recommendation: {recommendation}")
    print("=" * 80)
    
    return result


if __name__ == "__main__":
    try:
        result = compare_pathway_mappings()
        print("\n✅ PATHWAY COMPARISON COMPLETE")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ PATHWAY COMPARISON FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
