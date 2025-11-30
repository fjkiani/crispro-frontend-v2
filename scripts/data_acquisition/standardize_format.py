#!/usr/bin/env python3
"""
Standardize Format - Phase 1 Day 1 Task 3

Purpose: Standardize the final dataset format for validation.

Input: data/validation/tcga_ov_469_with_tmb.json
Output: data/validation/tcga_ov_469_validated.json

Standardization:
1. Ensure all required fields are present
2. Validate data quality
3. Rename fields to match validation schema
4. Add validation metadata
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Any, Optional

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

def validate_patient(patient: Dict[str, Any]) -> tuple[bool, List[str]]:
    """
    Validate a patient record
    
    Returns:
        (is_valid, errors)
    """
    errors = []
    
    # Required fields
    required_fields = ["patient_id", "cancer_type", "mutations", "treatment_response"]
    for field in required_fields:
        if field not in patient:
            errors.append(f"Missing required field: {field}")
    
    # Validate mutations
    mutations = patient.get("mutations", [])
    if not isinstance(mutations, list):
        errors.append("mutations must be a list")
    elif len(mutations) == 0:
        errors.append("No mutations found")
    
    # Validate treatment response
    treatment_response = patient.get("treatment_response", {})
    if not isinstance(treatment_response, dict):
        errors.append("treatment_response must be a dict")
    else:
        required_response_fields = ["treatment_type", "drug_name", "response", "response_category"]
        for field in required_response_fields:
            if field not in treatment_response:
                errors.append(f"Missing treatment_response field: {field}")
    
    # Validate biomarkers
    biomarkers = patient.get("biomarkers", {})
    if not isinstance(biomarkers, dict):
        errors.append("biomarkers must be a dict")
    
    is_valid = len(errors) == 0
    return is_valid, errors

def standardize_patient(patient: Dict[str, Any]) -> Dict[str, Any]:
    """Standardize a single patient record"""
    standardized = {
        # Core identifiers
        "patient_id": patient.get("patient_id", ""),
        "cancer_type": patient.get("cancer_type", "ovarian_hgsoc"),
        "tcga_code": patient.get("tcga_code", "TCGA-OV"),
        
        # Mutations (already standardized in Task 1)
        "mutations": patient.get("mutations", []),
        
        # Treatment response (already mapped in Task 1)
        "treatment_response": patient.get("treatment_response", {}),
        
        # Biomarkers (TMB computed in Task 2)
        "biomarkers": patient.get("biomarkers", {}),
        
        # Validation metadata
        "validation_metadata": {
            "dataset_version": "1.0",
            "standardization_date": "2025-01-28",
            "validation_ready": True,
            "data_quality": {
                "has_mutations": len(patient.get("mutations", [])) > 0,
                "has_treatment_response": bool(patient.get("treatment_response")),
                "has_tmb": patient.get("biomarkers", {}).get("tmb_score") is not None,
                "has_hrd": patient.get("biomarkers", {}).get("hrd_score") is not None
            }
        }
    }
    
    # Preserve original metadata if present
    if "metadata" in patient:
        standardized["original_metadata"] = patient["metadata"]
    
    return standardized

def generate_validation_report(patients: List[Dict[str, Any]], validation_results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Generate validation report"""
    total_patients = len(patients)
    valid_patients = sum(1 for r in validation_results if r["is_valid"])
    invalid_patients = total_patients - valid_patients
    
    # Data quality metrics
    patients_with_mutations = sum(1 for p in patients if len(p.get("mutations", [])) > 0)
    patients_with_response = sum(1 for p in patients if p.get("treatment_response"))
    patients_with_tmb = sum(1 for p in patients if p.get("biomarkers", {}).get("tmb_score") is not None)
    patients_with_hrd = sum(1 for p in patients if p.get("biomarkers", {}).get("hrd_score") is not None)
    
    # Response distribution
    response_categories = {}
    for p in patients:
        category = p.get("treatment_response", {}).get("response_category", "unknown")
        response_categories[category] = response_categories.get(category, 0) + 1
    
    report = {
        "validation_summary": {
            "total_patients": total_patients,
            "valid_patients": valid_patients,
            "invalid_patients": invalid_patients,
            "validation_rate": round(valid_patients / total_patients, 3) if total_patients > 0 else 0.0
        },
        "data_quality": {
            "patients_with_mutations": patients_with_mutations,
            "patients_with_treatment_response": patients_with_response,
            "patients_with_tmb": patients_with_tmb,
            "patients_with_hrd": patients_with_hrd,
            "coverage": {
                "mutations": round(patients_with_mutations / total_patients, 3) if total_patients > 0 else 0.0,
                "treatment_response": round(patients_with_response / total_patients, 3) if total_patients > 0 else 0.0,
                "tmb": round(patients_with_tmb / total_patients, 3) if total_patients > 0 else 0.0,
                "hrd": round(patients_with_hrd / total_patients, 3) if total_patients > 0 else 0.0
            }
        },
        "response_distribution": response_categories,
        "validation_errors": [
            r for r in validation_results if not r["is_valid"]
        ]
    }
    
    return report

def main():
    """Main standardization function"""
    # Input file
    input_file = Path("data/validation/tcga_ov_469_with_tmb.json")
    
    # Output file
    output_file = Path("data/validation/tcga_ov_469_validated.json")
    
    # Create output directory if needed
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"📋 Standardizing dataset format...")
    print(f"   Input: {input_file}")
    print(f"   Output: {output_file}")
    
    # Check if input file exists
    if not input_file.exists():
        print(f"❌ ERROR: Input file not found: {input_file}")
        print(f"   Please run compute_biomarkers.py first")
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
        sys.exit(1)
    
    print(f"   Found {len(patients)} patients")
    
    # Validate and standardize each patient
    standardized_patients = []
    validation_results = []
    
    for i, patient in enumerate(patients):
        # Validate
        is_valid, errors = validate_patient(patient)
        validation_results.append({
            "patient_index": i,
            "patient_id": patient.get("patient_id", "unknown"),
            "is_valid": is_valid,
            "errors": errors
        })
        
        if not is_valid:
            print(f"   ⚠️  Patient {i} ({patient.get('patient_id', 'unknown')}): {', '.join(errors)}")
        
        # Standardize
        standardized = standardize_patient(patient)
        standardized_patients.append(standardized)
    
    # Generate validation report
    report = generate_validation_report(patients, validation_results)
    
    # Report results
    print(f"\n✅ Standardization complete:")
    print(f"   Valid patients: {report['validation_summary']['valid_patients']}/{report['validation_summary']['total_patients']}")
    print(f"   Validation rate: {report['validation_summary']['validation_rate']:.1%}")
    print(f"\n   Data Quality:")
    print(f"      Mutations: {report['data_quality']['patients_with_mutations']} ({report['data_quality']['coverage']['mutations']:.1%})")
    print(f"      Treatment Response: {report['data_quality']['patients_with_treatment_response']} ({report['data_quality']['coverage']['treatment_response']:.1%})")
    print(f"      TMB: {report['data_quality']['patients_with_tmb']} ({report['data_quality']['coverage']['tmb']:.1%})")
    print(f"      HRD: {report['data_quality']['patients_with_hrd']} ({report['data_quality']['coverage']['hrd']:.1%})")
    print(f"\n   Response Distribution:")
    for category, count in report['response_distribution'].items():
        print(f"      {category}: {count} patients")
    
    # Save standardized data
    print(f"\n   Saving to {output_file}...")
    with open(output_file, 'w') as f:
        json.dump(standardized_patients, f, indent=2)
    
    # Save validation report
    report_file = output_file.parent / f"{output_file.stem}_validation_report.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\n✅ Success! Standardized dataset saved to: {output_file}")
    print(f"   Validation report saved to: {report_file}")
    print(f"\n   Dataset is now ready for pathway validation!")
    print(f"   Next step: Run pathway validation scripts")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

