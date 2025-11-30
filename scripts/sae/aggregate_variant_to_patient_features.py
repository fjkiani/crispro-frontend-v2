#!/usr/bin/env python3
"""
Aggregate Variant-Level SAE Features to Patient-Level

Converts the nested variant-level SAE data to patient-level aggregated features
suitable for biomarker correlation analysis.

Input: data/validation/sae_cohort/sae_features_tcga_ov_platinum.json
Output: data/validation/sae_cohort/sae_features_patient_aggregated.json

Aggregation Strategy:
- Mean activation across all variants for each patient
- Top-k most activated features per patient (k=100)
"""

import json
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List


def aggregate_patient_features(input_file: Path, output_file: Path):
    """Aggregate variant-level SAE features to patient-level."""
    
    print("=" * 80)
    print("SAE Feature Aggregation: Variant → Patient Level")
    print("=" * 80)
    
    # Load data
    print(f"\n[1/4] Loading data from {input_file}...")
    with open(input_file, 'r') as f:
        data = json.load(f)
    
    patients = data.get('patients', [])
    print(f"✅ Loaded {len(patients)} patients")
    
    # Aggregate features per patient
    print("\n[2/4] Aggregating features per patient...")
    aggregated_patients = []
    aggregated_outcome = []
    aggregated_sae_features = {}
    
    for patient in patients:
        patient_id = patient['patient_id']
        outcome = patient['outcome']
        variants = patient.get('variants', [])
        
        # Aggregate features across all variants for this patient
        feature_activations = defaultdict(list)
        
        for variant in variants:
            top_features = variant.get('top_features', [])
            for feature in top_features:
                idx = feature['index']
                val = feature['value']
                feature_activations[idx].append(val)
        
        # Compute mean activation for each feature
        patient_features = {}
        for idx, values in feature_activations.items():
            patient_features[str(idx)] = sum(values) / len(values)
        
        # Store
        aggregated_patients.append(patient_id)
        aggregated_outcome.append(1 if outcome == "sensitive" else 0)
        aggregated_sae_features[patient_id] = patient_features
        
        print(f"  - {patient_id}: {len(patient_features)} unique features, outcome={outcome}")
    
    # Create aggregated output
    print("\n[3/4] Creating aggregated output...")
    aggregated_data = {
        "cohort": data.get('cohort', 'TCGA-OV'),
        "outcome_field": "platinum_response",
        "aggregation_method": "mean_across_variants",
        "num_patients": len(aggregated_patients),
        "patients": aggregated_patients,
        "outcome": aggregated_outcome,
        "sae_features": aggregated_sae_features,
        "metadata": {
            "source": str(input_file),
            "aggregation_date": data.get('extraction_date', 'unknown'),
            "original_num_patients": len(patients),
            "total_features_extracted": sum(len(f) for f in aggregated_sae_features.values())
        }
    }
    
    # Save
    print(f"\n[4/4] Saving to {output_file}...")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(aggregated_data, f, indent=2)
    
    print(f"✅ Saved aggregated data: {len(aggregated_patients)} patients")
    print(f"   - Positive outcomes: {sum(aggregated_outcome)}")
    print(f"   - Negative outcomes: {len(aggregated_outcome) - sum(aggregated_outcome)}")
    print(f"   - Total unique features: {aggregated_data['metadata']['total_features_extracted']}")
    
    print("\n" + "=" * 80)
    print("AGGREGATION COMPLETE!")
    print("=" * 80)
    
    return aggregated_data


if __name__ == "__main__":
    input_file = Path("data/validation/sae_cohort/sae_features_tcga_ov_platinum.json")
    output_file = Path("data/validation/sae_cohort/sae_features_patient_aggregated.json")
    
    if not input_file.exists():
        print(f"❌ ERROR: Input file not found: {input_file}")
        sys.exit(1)
    
    aggregate_patient_features(input_file, output_file)



