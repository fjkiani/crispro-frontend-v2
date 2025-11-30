#!/usr/bin/env python3
"""
Transform Ovarian Dataset - Phase 1 Day 1 Task 1

Purpose: Transform the 469-patient dataset to the required validation schema.

Input: data/validation/sae_cohort/tcga_ov_platinum_with_mutations.json
Output: data/validation/tcga_ov_469_transformed.json

Transformation:
1. Map treatment response fields (platinum_response, raw_response_value)
2. Standardize mutation format
3. Add required fields for validation
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Any

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

def map_raw_to_recist(raw_value: str) -> str:
    """Map raw response value to RECIST criteria"""
    mapping = {
        "Complete Remission/Response": "CR",
        "Partial Remission/Response": "PR",
        "Progressive Disease": "PD",
        "Stable Disease": "SD"
    }
    return mapping.get(raw_value, "NE")  # NE = Not Evaluable

def map_treatment_response(patient: Dict[str, Any]) -> Dict[str, Any]:
    """Map treatment response to required schema"""
    return {
        "treatment_type": "chemo",
        "drug_name": "carboplatin/paclitaxel",  # Standard TCGA-OV regimen
        "response": map_raw_to_recist(patient.get("raw_response_value", "")),
        "response_category": patient.get("platinum_response", "unknown"),  # sensitive, resistant, refractory
        "pfs_months": patient.get("pfs_months"),  # If available
        "os_months": patient.get("os_months")  # If available
    }

def standardize_mutation(mut: Dict[str, Any]) -> Dict[str, Any]:
    """Standardize mutation format"""
    return {
        "gene": mut.get("gene", ""),
        "hgvs_p": mut.get("hgvs_p") or mut.get("protein_change", ""),
        "hgvs_c": mut.get("hgvs_c") or mut.get("cDNA_change", ""),
        "variant_classification": mut.get("variant_classification", ""),
        "chromosome": mut.get("chromosome") or mut.get("chrom", ""),
        "position": mut.get("position") or mut.get("pos"),
        "ref_allele": mut.get("ref_allele") or mut.get("ref", ""),
        "alt_allele": mut.get("alt_allele") or mut.get("alt", "")
    }

def transform_patient(patient: Dict[str, Any]) -> Dict[str, Any]:
    """Transform a single patient record"""
    # Standardize mutations
    mutations = patient.get("mutations", [])
    standardized_mutations = [standardize_mutation(m) for m in mutations]
    
    # Map treatment response
    treatment_response = map_treatment_response(patient)
    
    # Build transformed patient record
    transformed = {
        "patient_id": patient.get("patient_id") or patient.get("sample_id", ""),
        "cancer_type": "ovarian_hgsoc",
        "tcga_code": "TCGA-OV",
        "mutations": standardized_mutations,
        "treatment_response": treatment_response,
        "biomarkers": {
            "hrd_score": patient.get("hrd_score"),  # Will be filled in later
            "tmb_score": None,  # Will be computed in Task 2
            "msi_status": patient.get("msi_status", "Unknown"),
            "pd_l1_expression": patient.get("pd_l1_expression")
        },
        "metadata": {
            "source": patient.get("source", "tcga_ov_platinum_with_mutations"),
            "transformation_date": "2025-01-28",
            "original_fields": {
                "platinum_response": patient.get("platinum_response"),
                "raw_response_value": patient.get("raw_response_value")
            }
        }
    }
    
    return transformed

def main():
    """Main transformation function"""
    # Input file
    input_file = Path("data/validation/sae_cohort/tcga_ov_platinum_with_mutations.json")
    
    # Output file
    output_file = Path("data/validation/tcga_ov_469_transformed.json")
    
    # Create output directory if needed
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"📊 Transforming ovarian dataset...")
    print(f"   Input: {input_file}")
    print(f"   Output: {output_file}")
    
    # Check if input file exists
    if not input_file.exists():
        print(f"❌ ERROR: Input file not found: {input_file}")
        print(f"   Please ensure the 469-patient dataset exists.")
        print(f"   Expected location: {input_file.absolute()}")
        sys.exit(1)
    
    # Load input data
    print(f"   Loading {input_file}...")
    with open(input_file, 'r') as f:
        data = json.load(f)
    
    # Handle different data structures
    if isinstance(data, list):
        patients = data
    elif isinstance(data, dict) and "patients" in data:
        patients = data["patients"]
    elif isinstance(data, dict) and "data" in data:
        patients = data["data"]
    else:
        print(f"❌ ERROR: Unexpected data structure in {input_file}")
        print(f"   Expected: list of patients or dict with 'patients' or 'data' key")
        sys.exit(1)
    
    print(f"   Found {len(patients)} patients")
    
    # Transform each patient
    transformed_patients = []
    errors = []
    
    for i, patient in enumerate(patients):
        try:
            transformed = transform_patient(patient)
            transformed_patients.append(transformed)
        except Exception as e:
            errors.append({
                "patient_index": i,
                "patient_id": patient.get("patient_id", "unknown"),
                "error": str(e)
            })
            print(f"   ⚠️  Error transforming patient {i}: {e}")
    
    # Report results
    print(f"\n✅ Transformation complete:")
    print(f"   Transformed: {len(transformed_patients)} patients")
    if errors:
        print(f"   Errors: {len(errors)} patients")
        print(f"   Error details saved to: {output_file.parent / 'transformation_errors.json'}")
        
        # Save errors
        with open(output_file.parent / "transformation_errors.json", 'w') as f:
            json.dump(errors, f, indent=2)
    
    # Save transformed data
    print(f"   Saving to {output_file}...")
    with open(output_file, 'w') as f:
        json.dump(transformed_patients, f, indent=2)
    
    print(f"\n✅ Success! Transformed dataset saved to: {output_file}")
    print(f"   Next step: Run compute_biomarkers.py to compute TMB scores")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

