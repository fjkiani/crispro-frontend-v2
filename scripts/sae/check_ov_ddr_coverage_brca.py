#!/usr/bin/env python3
"""
Check OV DDR Feature Coverage in BRCA - Deliverable #3
=======================================================

Mission: Check how many BRCA variants have OV DDR features in top-64
Objective: Determine if we need to re-extract full 32K for OV DDR indices
Timeline: 2 hours
Status: EXECUTION READY
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Any
from collections import defaultdict
import statistics

# Paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
BRCA_CHECKPOINT = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "checkpoints" / "BRCA_TCGA_TRUE_SAE_cohort.json"
OUTPUT_DIR = PROJECT_ROOT / ".cursor" / "MOAT" / "SAE_INTELLIGENCE" / "BRCA_VALIDATION"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUTPUT_DIR / "ov_ddr_coverage_report.json"

# OV DDR diamond features
OV_DDR_FEATURES = [1407, 6020, 9738, 12893, 16337, 22868, 26220, 27607, 31362]


def check_coverage() -> Dict[str, Any]:
    """Check OV DDR feature coverage in BRCA checkpoint."""
    print("=" * 80)
    print("🔍 CHECKING OV DDR FEATURE COVERAGE - Deliverable #3")
    print("=" * 80)
    print(f"OV DDR Features: {OV_DDR_FEATURES}")
    print(f"BRCA Checkpoint: {BRCA_CHECKPOINT}")
    print()
    
    # Load BRCA checkpoint
    if not BRCA_CHECKPOINT.exists():
        print(f"❌ BRCA checkpoint not found: {BRCA_CHECKPOINT}")
        sys.exit(1)
    
    print("📥 Loading BRCA checkpoint...")
    with open(BRCA_CHECKPOINT, 'r') as f:
        data = json.load(f)
    
    patients = data.get("data", {})
    print(f"   ✅ Loaded {len(patients)} patients")
    print()
    
    # Coverage statistics
    coverage_stats = {
        "total_patients": len(patients),
        "patients_with_any_ddr": 0,
        "patients_with_all_ddr": 0,
        "patients_with_5plus_ddr": 0,
        "mean_features_per_patient": [],
        "median_features_per_patient": None,
        "min_features_per_patient": None,
        "max_features_per_patient": None,
        "missing_features": {f: 0 for f in OV_DDR_FEATURES},
        "feature_presence": {f: 0 for f in OV_DDR_FEATURES},
        "coverage_per_patient": [],
        "patients_by_coverage": {
            "0_features": 0,
            "1-2_features": 0,
            "3-4_features": 0,
            "5-6_features": 0,
            "7-8_features": 0,
            "9_features": 0
        }
    }
    
    # Check each patient
    print("🔍 Checking feature coverage per patient...")
    for patient_id, patient_data in patients.items():
        # Collect all feature indices from all variants
        all_variant_features = set()
        for variant in patient_data.get("variants", []):
            for tf in variant.get("top_features", []):
                idx = tf.get("index")
                if idx is not None:
                    all_variant_features.add(int(idx))
        
        # Check which OV DDR features are present
        present_features = [f for f in OV_DDR_FEATURES if f in all_variant_features]
        num_present = len(present_features)
        
        # Update stats
        if num_present > 0:
            coverage_stats["patients_with_any_ddr"] += 1
            coverage_stats["mean_features_per_patient"].append(num_present)
            coverage_stats["coverage_per_patient"].append(num_present / len(OV_DDR_FEATURES))
        
        if num_present == len(OV_DDR_FEATURES):
            coverage_stats["patients_with_all_ddr"] += 1
        
        if num_present >= 5:
            coverage_stats["patients_with_5plus_ddr"] += 1
        
        # Track missing features
        for f in OV_DDR_FEATURES:
            if f in all_variant_features:
                coverage_stats["feature_presence"][f] += 1
            else:
                coverage_stats["missing_features"][f] += 1
        
        # Categorize by coverage
        if num_present == 0:
            coverage_stats["patients_by_coverage"]["0_features"] += 1
        elif num_present <= 2:
            coverage_stats["patients_by_coverage"]["1-2_features"] += 1
        elif num_present <= 4:
            coverage_stats["patients_by_coverage"]["3-4_features"] += 1
        elif num_present <= 6:
            coverage_stats["patients_by_coverage"]["5-6_features"] += 1
        elif num_present <= 8:
            coverage_stats["patients_by_coverage"]["7-8_features"] += 1
        else:
            coverage_stats["patients_by_coverage"]["9_features"] += 1
    
    # Compute statistics
    if coverage_stats["mean_features_per_patient"]:
        features_list = coverage_stats["mean_features_per_patient"]
        coverage_stats["mean_features_per_patient"] = statistics.mean(features_list)
        coverage_stats["median_features_per_patient"] = statistics.median(features_list)
        coverage_stats["min_features_per_patient"] = min(features_list)
        coverage_stats["max_features_per_patient"] = max(features_list)
        coverage_stats["mean_coverage_percentage"] = statistics.mean(coverage_stats["coverage_per_patient"]) * 100
    else:
        coverage_stats["mean_features_per_patient"] = 0.0
        coverage_stats["median_features_per_patient"] = 0.0
        coverage_stats["min_features_per_patient"] = 0
        coverage_stats["max_features_per_patient"] = 0
        coverage_stats["mean_coverage_percentage"] = 0.0
    
    # Decision
    mean_coverage = coverage_stats["mean_coverage_percentage"]
    if mean_coverage >= 50:
        decision = "proceed"
        recommendation = "Coverage is sufficient (≥50%). Proceed with OV DDR transfer test."
    else:
        decision = "re_extract"
        recommendation = f"Coverage is low ({mean_coverage:.1f}% < 50%). Consider re-extracting full 32K for OV DDR indices only."
    
    # Build report
    report = {
        "coverage_stats": coverage_stats,
        "decision": decision,
        "recommendation": recommendation,
        "ov_ddr_features": OV_DDR_FEATURES,
        "checkpoint_file": str(BRCA_CHECKPOINT),
        "extraction_date": datetime.now().isoformat()
    }
    
    # Save report
    print()
    print("💾 Saving coverage report...")
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(report, f, indent=2)
    print(f"   ✅ Saved to: {OUTPUT_FILE}")
    
    # Print summary
    print()
    print("=" * 80)
    print("📊 COVERAGE SUMMARY")
    print("=" * 80)
    print(f"Total patients: {coverage_stats['total_patients']}")
    print(f"Patients with any DDR feature: {coverage_stats['patients_with_any_ddr']} ({coverage_stats['patients_with_any_ddr']/coverage_stats['total_patients']*100:.1f}%)")
    print(f"Patients with all 9 DDR features: {coverage_stats['patients_with_all_ddr']} ({coverage_stats['patients_with_all_ddr']/coverage_stats['total_patients']*100:.1f}%)")
    print(f"Patients with 5+ DDR features: {coverage_stats['patients_with_5plus_ddr']} ({coverage_stats['patients_with_5plus_ddr']/coverage_stats['total_patients']*100:.1f}%)")
    print()
    print(f"Mean DDR features per patient: {coverage_stats['mean_features_per_patient']:.1f} / 9")
    print(f"Mean coverage percentage: {coverage_stats['mean_coverage_percentage']:.1f}%")
    print()
    print("Feature presence across patients:")
    for f in OV_DDR_FEATURES:
        presence_pct = (coverage_stats['feature_presence'][f] / coverage_stats['total_patients']) * 100
        print(f"  Feature {f}: {coverage_stats['feature_presence'][f]}/{coverage_stats['total_patients']} ({presence_pct:.1f}%)")
    print()
    print("Patients by coverage:")
    for category, count in coverage_stats['patients_by_coverage'].items():
        pct = (count / coverage_stats['total_patients']) * 100
        print(f"  {category}: {count} ({pct:.1f}%)")
    print()
    print("=" * 80)
    print(f"🎯 DECISION: {decision.upper()}")
    print("=" * 80)
    print(recommendation)
    print("=" * 80)
    
    return report


if __name__ == "__main__":
    from datetime import datetime
    import datetime as dt_module
    try:
        result = check_coverage()
        print("\n✅ COVERAGE CHECK COMPLETE")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ COVERAGE CHECK FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
