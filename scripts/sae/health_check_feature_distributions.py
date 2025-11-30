#!/usr/bin/env python3
"""
🔍 SAE Feature Distribution Health Check

Purpose: Verify SAE feature distributions are reasonable (sparse, varied)
Agent: Zo (Lead Commander)
Date: January 20, 2025
"""

import json
import numpy as np
from pathlib import Path
import sys

def check_feature_distributions():
    """Verify SAE feature distributions are reasonable."""
    data_file = Path("data/validation/sae_cohort/sae_features_tcga_ov_platinum.json")
    
    print("🔍 Checking SAE feature distributions...")
    
    if not data_file.exists():
        print(f"❌ ERROR: Data file not found: {data_file}")
        return False
    
    # Load data
    with open(data_file) as f:
        data = json.load(f)
    
    # Collect all feature activations from variants
    all_features = []
    patient_features = []
    
    for patient in data["patients"]:
        if "variants" not in patient:
            continue
        
        # Collect all feature values from all variants for this patient
        patient_feat_values = []
        
        for variant in patient["variants"]:
            # Check top_features (list of dicts with 'index' and 'value')
            if "top_features" in variant and isinstance(variant["top_features"], list):
                for feat_dict in variant["top_features"]:
                    if isinstance(feat_dict, dict) and "value" in feat_dict:
                        feat_value = feat_dict["value"]
                        all_features.append(feat_value)
                        patient_feat_values.append(feat_value)
        
        if patient_feat_values:
            patient_features.append(np.array(patient_feat_values))
    
    if not all_features:
        print("❌ ERROR: No feature activations found")
        return False
    
    # Check distributions
    features_array = np.array(all_features)
    
    # Check not all zeros
    if np.all(features_array == 0):
        print("❌ ERROR: All features are zero - check extraction")
        return False
    
    # Check not all ones
    if np.all(features_array == 1.0):
        print("❌ ERROR: All features are one - check extraction")
        return False
    
    # Check reasonable range (SAE activations should be sparse, mostly 0)
    zero_fraction = np.mean(features_array == 0)
    if zero_fraction <= 0.8:
        print(f"⚠️  WARNING: Too few zeros ({zero_fraction:.2%}) - SAE should be sparse")
    else:
        print(f"✅ Zero fraction: {zero_fraction:.2%} (sparse)")
    
    # Check variation across patients
    if len(patient_features) > 1:
        correlations = []
        for i in range(len(patient_features) - 1):
            if len(patient_features[i]) == len(patient_features[i+1]):
                corr = np.corrcoef(patient_features[i], patient_features[i+1])[0, 1]
                if not np.isnan(corr):
                    correlations.append(corr)
        
        if correlations:
            mean_corr = np.mean(correlations)
            if mean_corr >= 0.95:
                print(f"⚠️  WARNING: Patients too similar (mean corr: {mean_corr:.3f})")
            else:
                print(f"✅ Patient variation: mean correlation {mean_corr:.3f}")
    
    # Summary statistics
    print(f"✅ Feature distribution summary:")
    print(f"   - Mean activation: {np.mean(features_array):.4f}")
    print(f"   - Std activation: {np.std(features_array):.4f}")
    print(f"   - Min activation: {np.min(features_array):.4f}")
    print(f"   - Max activation: {np.max(features_array):.4f}")
    print(f"   - Non-zero count: {np.sum(features_array > 0)} / {len(features_array)}")
    
    return True

if __name__ == "__main__":
    success = check_feature_distributions()
    sys.exit(0 if success else 1)
