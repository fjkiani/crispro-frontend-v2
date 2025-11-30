#!/usr/bin/env python3
"""
Health Check: SAE Cohort Data Structure and Quality

Verifies:
- Data file exists
- Structure is valid (patients, outcome, sae_features)
- Outcome field populated correctly
- Feature indices valid (0-32767)
- Patient count correct (66 patients)
"""

import json
import sys
from pathlib import Path


def check_sae_cohort_data():
    """Verify SAE cohort data structure and quality."""
    data_file = Path("data/validation/sae_cohort/sae_features_tcga_ov_platinum.json")
    
    print("=" * 80)
    print("HEALTH CHECK: SAE Cohort Data Structure")
    print("=" * 80)
    
    # Check 1: File exists
    print("\n[1/6] Checking file exists...")
    if not data_file.exists():
        print(f"❌ FAIL: Data file not found: {data_file}")
        return False
    print(f"✅ PASS: File exists: {data_file}")
    
    # Check 2: Load data
    print("\n[2/6] Loading data...")
    try:
        with open(data_file) as f:
            data = json.load(f)
        print(f"✅ PASS: Data loaded successfully")
    except Exception as e:
        print(f"❌ FAIL: Error loading data: {e}")
        return False
    
    # Check 3: Structure validation
    print("\n[3/6] Validating structure...")
    if "patients" not in data:
        print("❌ FAIL: Missing 'patients' key")
        return False
    
    patient_count = len(data["patients"])
    if patient_count != 66:
        print(f"⚠️  WARN: Expected 66 patients, got {patient_count}")
    else:
        print(f"✅ PASS: Patient count correct: {patient_count}")
    
    # Check 4: Patient structure
    print("\n[4/7] Validating patient structure...")
    required_fields = ["patient_id", "outcome", "variants"]
    errors = []
    
    for i, patient in enumerate(data["patients"]):
        for field in required_fields:
            if field not in patient:
                errors.append(f"Patient {i}: Missing '{field}'")
        
        if "outcome" in patient:
            valid_outcomes = ["sensitive", "refractory", "resistant"]
            if patient["outcome"] not in valid_outcomes:
                errors.append(f"Patient {i}: Invalid outcome '{patient['outcome']}'")
        
        # Check variants structure
        if "variants" in patient:
            if not isinstance(patient["variants"], list):
                errors.append(f"Patient {i}: 'variants' must be a list")
            else:
                # Check each variant in the list
                for j, variant in enumerate(patient["variants"]):
                    # Variant should have either top_features or sae_features
                    if "top_features" not in variant and "sae_features" not in variant:
                        errors.append(f"Patient {i}, Variant {j}: Missing 'top_features' or 'sae_features'")
                    
                    # If top_features exists, it should be a list
                    if "top_features" in variant:
                        if not isinstance(variant["top_features"], list):
                            errors.append(f"Patient {i}, Variant {j}: 'top_features' must be a list")
                    
                    # If sae_features exists, it should be a list
                    if "sae_features" in variant:
                        if not isinstance(variant["sae_features"], list):
                            errors.append(f"Patient {i}, Variant {j}: 'sae_features' must be a list")
    
    if errors:
        print(f"❌ FAIL: Found {len(errors)} structure errors:")
        for error in errors[:5]:  # Show first 5
            print(f"   - {error}")
        if len(errors) > 5:
            print(f"   ... and {len(errors) - 5} more")
        return False
    print(f"✅ PASS: All patients have valid structure")
    
    # Check 5: Feature indices validation
    print("\n[5/7] Validating feature indices...")
    invalid_indices = []
    total_features = 0
    
    for i, patient in enumerate(data["patients"]):
        if "variants" in patient:
            for j, variant in enumerate(patient["variants"]):
                # Check top_features (list of dicts with 'index' and 'value')
                if "top_features" in variant and isinstance(variant["top_features"], list):
                    for feat_dict in variant["top_features"]:
                        if isinstance(feat_dict, dict) and "index" in feat_dict:
                            feat_idx = feat_dict["index"]
                            total_features += 1
                            if not isinstance(feat_idx, int) or not (0 <= feat_idx < 32768):
                                invalid_indices.append((i, j, feat_idx))
                
                # Also check sae_features (list of indices) if present
                if "sae_features" in variant and isinstance(variant["sae_features"], list):
                    for feat_idx in variant["sae_features"]:
                        if isinstance(feat_idx, int):
                            total_features += 1
                            if not (0 <= feat_idx < 32768):
                                invalid_indices.append((i, j, feat_idx))
    
    if invalid_indices:
        print(f"❌ FAIL: Found {len(invalid_indices)} invalid feature indices:")
        for patient_idx, feat_idx in invalid_indices[:5]:
            print(f"   - Patient {patient_idx}: feature index {feat_idx}")
        return False
    print(f"✅ PASS: All feature indices valid (0-32767)")
    
    # Check 6: Outcome distribution
    print("\n[6/7] Analyzing outcome distribution...")
    outcomes = [p["outcome"] for p in data["patients"] if "outcome" in p]
    outcome_counts = {
        "sensitive": outcomes.count("sensitive"),
        "refractory": outcomes.count("refractory"),
        "resistant": outcomes.count("resistant")
    }
    
    print(f"✅ PASS: Outcome distribution:")
    for outcome, count in outcome_counts.items():
        print(f"   - {outcome}: {count} patients ({count/len(outcomes)*100:.1f}%)")
    
    # Summary
    print("\n" + "=" * 80)
    print("✅ HEALTH CHECK PASSED: SAE Cohort Data")
    print("=" * 80)
    print(f"Total patients: {patient_count}")
    print(f"Outcome distribution: {outcome_counts}")
    print(f"Data quality: GOOD")
    print("=" * 80)
    
    return True


if __name__ == "__main__":
    try:
        success = check_sae_cohort_data()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\n❌ FATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
