#!/usr/bin/env python3
"""
Test OV DDR Transfer to BRCA - Deliverable #4
==============================================

Mission: Test if OV DDR diamond features predict BRCA recurrence
Objective: Core hypothesis test - OV DDR → BRCA recurrence
Timeline: 2 hours
Status: EXECUTION READY

Based on: SAE_BRCA_NEXT_10_DELIVERABLES.md (Deliverable #4)
Corrections: Use mean() aggregation, handle missing features (0.0)
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime
import numpy as np
from sklearn.metrics import roc_auc_score, roc_curve
from scipy import stats
import statistics

# Paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
BRCA_CHECKPOINT = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "checkpoints" / "BRCA_TCGA_TRUE_SAE_cohort.json"
RECURRENCE_LABELS = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "brca_recurrence_labels.json"
OUTPUT_DIR = PROJECT_ROOT / ".cursor" / "MOAT" / "SAE_INTELLIGENCE" / "BRCA_VALIDATION"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUTPUT_DIR / "ov_ddr_transfer_results.json"

# OV DDR diamond features
OV_DDR_FEATURES = [1407, 6020, 9738, 12893, 16337, 22868, 26220, 27607, 31362]


def compute_ddr_bin_mean(patient_data: Dict, diamond_indices: List[int]) -> float:
    """
    Compute DDR_bin using MEAN aggregation (matches manuscript).
    
    Method:
    1. Sum all diamond feature values across all variants
    2. Divide by number of diamond features (9)
    3. Missing features = 0.0 (handled automatically)
    """
    feature_sums = {fidx: 0.0 for fidx in diamond_indices}
    
    for variant in patient_data.get("variants", []):
        for tf in variant.get("top_features", []):
            fidx = tf.get("index")
            if fidx in feature_sums:
                # Sum absolute values (as per publication script)
                feature_sums[fidx] += abs(float(tf.get("value", 0.0) or 0.0))
    
    # Mean of 9 diamond features (missing features contribute 0.0)
    return sum(feature_sums.values()) / len(diamond_indices) if diamond_indices else 0.0


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


def test_ov_ddr_transfer() -> Dict[str, Any]:
    """Test if OV DDR diamond features predict BRCA recurrence."""
    print("=" * 80)
    print("🧪 TESTING OV DDR TRANSFER TO BRCA - Deliverable #4")
    print("=" * 80)
    print(f"OV DDR Features: {OV_DDR_FEATURES}")
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
    
    # Step 3: Compute DDR_bin for each patient
    print()
    print("🔍 Step 3: Computing DDR_bin scores...")
    ddr_bins = []
    recurrence_labels = []
    patient_ids = []
    stats = {
        "total_patients": 0,
        "with_recurrence_label": 0,
        "recurrence_true": 0,
        "recurrence_false": 0,
        "ddr_bin_stats": {
            "mean": 0.0,
            "median": 0.0,
            "min": 0.0,
            "max": 0.0,
            "std": 0.0
        }
    }
    
    for patient_id, patient_data in patients.items():
        stats["total_patients"] += 1
        
        # Get recurrence label
        outcome = outcomes.get(patient_id, {})
        recurrence = outcome.get("recurrence")
        
        if recurrence is None:
            continue  # Skip patients without recurrence labels
        
        stats["with_recurrence_label"] += 1
        if recurrence:
            stats["recurrence_true"] += 1
        else:
            stats["recurrence_false"] += 1
        
        # Compute DDR_bin using mean() aggregation
        ddr_bin = compute_ddr_bin_mean(patient_data, OV_DDR_FEATURES)
        
        ddr_bins.append(ddr_bin)
        recurrence_labels.append(recurrence)
        patient_ids.append(patient_id)
    
    # Compute DDR_bin statistics
    if ddr_bins:
        stats["ddr_bin_stats"]["mean"] = statistics.mean(ddr_bins)
        stats["ddr_bin_stats"]["median"] = statistics.median(ddr_bins)
        stats["ddr_bin_stats"]["min"] = min(ddr_bins)
        stats["ddr_bin_stats"]["max"] = max(ddr_bins)
        stats["ddr_bin_stats"]["std"] = statistics.stdev(ddr_bins) if len(ddr_bins) > 1 else 0.0
    
    print(f"   ✅ Computed DDR_bin for {len(ddr_bins)} patients")
    print(f"   DDR_bin stats: mean={stats['ddr_bin_stats']['mean']:.3f}, std={stats['ddr_bin_stats']['std']:.3f}")
    
    # Step 4: Compute AUROC
    print()
    print("📊 Step 4: Computing AUROC...")
    if len(set(recurrence_labels)) < 2:
        print("❌ Insufficient outcome diversity for AUROC")
        result = {
            "error": "Insufficient outcome diversity",
            "n_patients": len(ddr_bins),
            "recurrence_distribution": {
                "true": stats["recurrence_true"],
                "false": stats["recurrence_false"]
            }
        }
    else:
        auroc, auroc_stats = compute_auroc(recurrence_labels, ddr_bins)
        
        # Interpret results
        if auroc >= 0.70:
            interpretation = "success"
            next_action = "use_ov_mapping"
            scenario = "Scenario A: OV DDR transfers to BRCA (AUROC ≥ 0.70)"
        elif auroc >= 0.60:
            interpretation = "moderate"
            next_action = "refine_brca_ddr"
            scenario = "Scenario B: Moderate transfer (0.60-0.69)"
        else:
            interpretation = "poor"
            next_action = "run_brca_specific"
            scenario = "Scenario C: Poor transfer (< 0.60)"
        
        result = {
            "ov_ddr_auroc": float(auroc),
            "auroc_ci_lower": auroc_stats.get("ci_lower"),
            "auroc_ci_upper": auroc_stats.get("ci_upper"),
            "interpretation": interpretation,
            "next_action": next_action,
            "scenario": scenario,
            "n_patients": len(ddr_bins),
            "recurrence_distribution": {
                "true": stats["recurrence_true"],
                "false": stats["recurrence_false"]
            },
            "ddr_bin_stats": stats["ddr_bin_stats"],
            "aggregation_method": "mean",
            "ov_ddr_features": OV_DDR_FEATURES,
            "missing_features_handled": "0.0 (automatic in mean aggregation)"
        }
        
        print(f"   ✅ AUROC: {auroc:.3f} (95% CI: {auroc_stats.get('ci_lower', 0):.3f}-{auroc_stats.get('ci_upper', 0):.3f})")
        print(f"   Interpretation: {interpretation}")
        print(f"   Next action: {next_action}")
    
    # Step 5: Save results
    print()
    print("💾 Step 5: Saving results...")
    output_data = {
        "test_date": datetime.now().isoformat(),
        "test_name": "OV DDR Transfer to BRCA",
        "result": result,
        "stats": stats,
        "provenance": {
            "brca_checkpoint": str(BRCA_CHECKPOINT),
            "recurrence_labels": str(RECURRENCE_LABELS),
            "aggregation_method": "mean",
            "ov_ddr_features": OV_DDR_FEATURES
        }
    }
    
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"   ✅ Saved to: {OUTPUT_FILE}")
    
    # Step 6: Print summary
    print()
    print("=" * 80)
    print("📊 TRANSFER TEST SUMMARY")
    print("=" * 80)
    if "error" not in result:
        print(f"OV DDR AUROC: {result['ov_ddr_auroc']:.3f}")
        print(f"95% CI: [{result['auroc_ci_lower']:.3f}, {result['auroc_ci_upper']:.3f}]")
        print(f"Interpretation: {result['interpretation']}")
        print(f"Scenario: {result['scenario']}")
        print(f"Next action: {result['next_action']}")
        print(f"Patients: {result['n_patients']}")
        print(f"Recurrence: {result['recurrence_distribution']['true']} True, {result['recurrence_distribution']['false']} False")
    else:
        print(f"Error: {result['error']}")
    print("=" * 80)
    
    return output_data


if __name__ == "__main__":
    try:
        result = test_ov_ddr_transfer()
        print("\n✅ TRANSFER TEST COMPLETE")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ TRANSFER TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
