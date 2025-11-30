#!/usr/bin/env python3
"""
Compute Biomarkers - Phase 1 Day 1 Task 2

Purpose: Compute TMB scores for all 469 ovarian cancer patients.

Input: data/validation/tcga_ov_469_transformed.json
Output: data/validation/tcga_ov_469_with_tmb.json

Formula:
- TMB = num_mutations / 38.0 (standard WES exome size in Mb)
- Thresholds:
  - TMB-Low: <6 mut/Mb
  - TMB-Intermediate: 6-10 mut/Mb
  - TMB-High: ≥10 mut/Mb
"""

import json
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Any

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

# Standard WES exome size (Mb)
EXOME_SIZE_MB = 38.0

# TMB thresholds
TMB_LOW_THRESHOLD = 6.0
TMB_HIGH_THRESHOLD = 10.0

def compute_tmb(patient: Dict[str, Any]) -> Dict[str, Any]:
    """
    Compute TMB from mutation count
    
    Args:
        patient: Patient record with mutations list
        
    Returns:
        Dict with tmb_score, tmb_category, num_mutations
    """
    mutations = patient.get("mutations", [])
    num_mutations = len(mutations)
    
    # Compute TMB
    tmb_score = round(num_mutations / EXOME_SIZE_MB, 2)
    
    # Categorize TMB
    if tmb_score < TMB_LOW_THRESHOLD:
        tmb_category = "TMB-Low"
    elif tmb_score < TMB_HIGH_THRESHOLD:
        tmb_category = "TMB-Intermediate"
    else:
        tmb_category = "TMB-High"
    
    return {
        "tmb_score": tmb_score,
        "tmb_category": tmb_category,
        "num_mutations": num_mutations,
        "exome_size_mb": EXOME_SIZE_MB
    }

def compute_biomarkers_for_patients(patients: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Compute biomarkers for all patients"""
    updated_patients = []
    
    for patient in patients:
        # Compute TMB
        tmb_result = compute_tmb(patient)
        
        # Update patient record
        updated_patient = patient.copy()
        
        # Update biomarkers section
        if "biomarkers" not in updated_patient:
            updated_patient["biomarkers"] = {}
        
        updated_patient["biomarkers"]["tmb_score"] = tmb_result["tmb_score"]
        updated_patient["biomarkers"]["tmb_category"] = tmb_result["tmb_category"]
        updated_patient["biomarkers"]["num_mutations"] = tmb_result["num_mutations"]
        
        # Add computation metadata
        if "metadata" not in updated_patient:
            updated_patient["metadata"] = {}
        
        if "biomarker_computation" not in updated_patient["metadata"]:
            updated_patient["metadata"]["biomarker_computation"] = {}
        
        updated_patient["metadata"]["biomarker_computation"]["tmb_computed"] = True
        updated_patient["metadata"]["biomarker_computation"]["tmb_date"] = "2025-01-28"
        updated_patient["metadata"]["biomarker_computation"]["tmb_formula"] = f"num_mutations / {EXOME_SIZE_MB}"
        
        updated_patients.append(updated_patient)
    
    return updated_patients

def generate_tmb_summary(patients: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Generate summary statistics for TMB"""
    tmb_scores = [p["biomarkers"]["tmb_score"] for p in patients if p["biomarkers"].get("tmb_score")]
    
    if not tmb_scores:
        return {"error": "No TMB scores computed"}
    
    tmb_categories = [p["biomarkers"]["tmb_category"] for p in patients if p["biomarkers"].get("tmb_category")]
    
    summary = {
        "total_patients": len(patients),
        "patients_with_tmb": len(tmb_scores),
        "tmb_statistics": {
            "mean": round(sum(tmb_scores) / len(tmb_scores), 2),
            "median": round(sorted(tmb_scores)[len(tmb_scores) // 2], 2),
            "min": round(min(tmb_scores), 2),
            "max": round(max(tmb_scores), 2),
            "std": round((sum((x - sum(tmb_scores)/len(tmb_scores))**2 for x in tmb_scores) / len(tmb_scores))**0.5, 2)
        },
        "tmb_distribution": {
            "TMB-Low": tmb_categories.count("TMB-Low"),
            "TMB-Intermediate": tmb_categories.count("TMB-Intermediate"),
            "TMB-High": tmb_categories.count("TMB-High")
        }
    }
    
    return summary

def main():
    """Main biomarker computation function"""
    parser = argparse.ArgumentParser(description="Compute TMB scores for ovarian cancer patients")
    parser.add_argument("--tmb", action="store_true", help="Compute TMB scores (default: all biomarkers)")
    parser.add_argument("--input", type=str, default="data/validation/tcga_ov_469_transformed.json",
                       help="Input file path")
    parser.add_argument("--output", type=str, default="data/validation/tcga_ov_469_with_tmb.json",
                       help="Output file path")
    
    args = parser.parse_args()
    
    # Input file
    input_file = Path(args.input)
    
    # Output file
    output_file = Path(args.output)
    
    # Create output directory if needed
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"🧬 Computing biomarkers...")
    print(f"   Input: {input_file}")
    print(f"   Output: {output_file}")
    
    # Check if input file exists
    if not input_file.exists():
        print(f"❌ ERROR: Input file not found: {input_file}")
        print(f"   Please run transform_ovarian_dataset.py first")
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
    
    # Compute biomarkers
    if args.tmb or True:  # Default to TMB computation
        print(f"\n   Computing TMB scores...")
        updated_patients = compute_biomarkers_for_patients(patients)
        
        # Generate summary
        summary = generate_tmb_summary(updated_patients)
        
        print(f"\n   TMB Summary:")
        print(f"      Patients with TMB: {summary['patients_with_tmb']}/{summary['total_patients']}")
        print(f"      Mean TMB: {summary['tmb_statistics']['mean']} mut/Mb")
        print(f"      Median TMB: {summary['tmb_statistics']['median']} mut/Mb")
        print(f"      Range: {summary['tmb_statistics']['min']} - {summary['tmb_statistics']['max']} mut/Mb")
        print(f"\n      Distribution:")
        print(f"         TMB-Low: {summary['tmb_distribution']['TMB-Low']} patients")
        print(f"         TMB-Intermediate: {summary['tmb_distribution']['TMB-Intermediate']} patients")
        print(f"         TMB-High: {summary['tmb_distribution']['TMB-High']} patients")
        
        # Save updated data
        print(f"\n   Saving to {output_file}...")
        with open(output_file, 'w') as f:
            json.dump(updated_patients, f, indent=2)
        
        # Save summary
        summary_file = output_file.parent / f"{output_file.stem}_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        print(f"\n✅ Success! Biomarkers computed and saved to: {output_file}")
        print(f"   Summary saved to: {summary_file}")
        print(f"   Next step: Run standardize_format.py to finalize dataset")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

