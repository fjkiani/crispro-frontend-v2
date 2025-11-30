"""
Merge Jr's platinum response data with Zo's genomic mutations.

Input Files:
- Jr's data: data/validation/tcga_ov_platinum_response_labels.json (469 patients with platinum response)
- Zo's data: tools/benchmarks/hrd_tcga_ov_labeled_sample_use_evo.json (200 patients with mutations)

Output:
- Merged data: tools/benchmarks/tcga_ov_platinum_response_with_genomics.json

Author: Zo
Date: January 14, 2025
For: Resistance Prophet validation
"""

import json
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    """Merge Jr's platinum response with Zo's genomic data"""
    
    # Load Jr's platinum response data (469 patients)
    logger.info("Loading Jr's platinum response data...")
    jr_file = Path("data/validation/tcga_ov_platinum_response_labels.json")
    
    with open(jr_file, 'r') as f:
        jr_data = json.load(f)
    
    jr_patients = jr_data['patients']
    logger.info(f"Jr's data: {len(jr_patients)} patients with platinum response")
    
    # Create lookup by patient_id and sample_id
    jr_lookup = {}
    for patient in jr_patients:
        patient_id = patient['tcga_patient_id']
        sample_id = patient.get('tcga_sample_id')
        
        jr_lookup[patient_id] = patient
        if sample_id:
            jr_lookup[sample_id] = patient
    
    # Load Zo's genomic data (200 patients with mutations)
    logger.info("Loading Zo's genomic data...")
    zo_file = Path("data/validation/tcga_ov_full_validation_dataset.json")
    
    with open(zo_file, 'r') as f:
        zo_data_full = json.load(f)
    
    zo_data = zo_data_full.get('patients', [])
    logger.info(f"Zo's data: {len(zo_data)} patients with mutations")
    
    # Merge datasets
    logger.info("Merging datasets (patient_id/sample_id join)...")
    merged_patients = []
    matched_count = 0
    
    for zo_patient in zo_data:
        patient_id = zo_patient.get('patient_id')
        sample_id = zo_patient.get('sample_id')
        
        # Try to find match in Jr's data
        jr_match = None
        if patient_id and patient_id in jr_lookup:
            jr_match = jr_lookup[patient_id]
        elif sample_id and sample_id in jr_lookup:
            jr_match = jr_lookup[sample_id]
        
        if jr_match:
            # Merge
            merged_patient = {
                **zo_patient,  # Genomic data from Zo
                "platinum_response": jr_match['platinum_response'],  # Outcome from Jr
                "raw_response_value": jr_match['raw_response_value'],
                "platinum_source": jr_match['source_field'],
                "tcga_patient_id": jr_match['tcga_patient_id'],
                "tcga_sample_id": jr_match.get('tcga_sample_id')
            }
            merged_patients.append(merged_patient)
            matched_count += 1
    
    logger.info(f"Matched {matched_count} patients (overlap between Jr's 469 and Zo's {len(zo_data)})")
    logger.info(f"Match rate: {matched_count / len(zo_data):.1%}")
    
    # Create binary resistance label (Manager Q2)
    for patient in merged_patients:
        response = patient['platinum_response']
        patient['resistant_binary'] = 0 if response == 'sensitive' else 1
    
    # Count distribution
    resistant_count = sum(p['resistant_binary'] for p in merged_patients)
    sensitive_count = matched_count - resistant_count
    
    logger.info(f"Outcome distribution:")
    logger.info(f"  Sensitive: {sensitive_count} ({sensitive_count/matched_count:.1%})")
    logger.info(f"  Resistant: {resistant_count} ({resistant_count/matched_count:.1%})")
    
    # Save merged dataset
    output_file = Path("tools/benchmarks/tcga_ov_platinum_response_with_genomics.json")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    output_data = {
        "metadata": {
            "source": "Merged dataset (Zo's genomics + Jr's platinum response)",
            "zo_source": str(zo_file),
            "jr_source": str(jr_file),
            "merge_date": str(Path(__file__).stat().st_mtime),
            "total_patients": matched_count,
            "match_method": "patient_id OR sample_id lookup",
            "outcome_distribution": {
                "sensitive": sensitive_count,
                "resistant": resistant_count
            }
        },
        "patients": merged_patients
    }
    
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    logger.info(f"✅ Merged dataset saved: {output_file}")
    logger.info(f"✅ {matched_count} patients ready for Resistance Prophet validation")
    
    return output_file


if __name__ == "__main__":
    output_file = main()
    print(f"\n🎯 NEXT STEP: Run validation script")
    print(f"python scripts/validate_resistance_prophet.py")

