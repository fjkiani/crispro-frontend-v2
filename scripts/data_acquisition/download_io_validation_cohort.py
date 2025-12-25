#!/usr/bin/env python3
"""
IO Validation Cohort Data Acquisition
======================================
Extract TCGA Pan-Cancer Atlas patients for IO prediction validation.

This script extracts patients from 4 TCGA studies:
- Endometrial (UCEC) - high MSI rate
- Colorectal (COAD/READ) - high MSI rate
- Ovarian (OV) - priority for Ayesha context
- Melanoma (SKCM) - high TMB

Output: data/validation/io_validation_cohort.json

Author: Agent
Date: January 28, 2025
"""

import json
import sys
from pathlib import Path
from typing import List, Dict, Optional
from collections import Counter, defaultdict
from datetime import datetime

# Import modular utilities
utils_path = Path(__file__).parent / "utils"
sys.path.insert(0, str(utils_path))
from cbioportal_client import CBioportalClient
from tmb_calculator import (
    extract_tmb_from_clinical,
    calculate_tmb_from_mutations,
    normalize_tmb_value,
    TCGA_EXOME_SIZE_MB
)
from mutation_processor import (
    process_mutations,
    find_primary_tumor_sample,
    has_mbd4_mutation,
    has_tp53_mutation,
    has_ddr_mutation,
    has_mmr_mutation
)
from data_quality import (
    extract_msi_from_clinical,
    generate_data_quality_flags,
    deduplicate_patients,
    validate_patient_record
)

# Study configuration
STUDY_TO_CANCER_TYPE = {
    "ucec_tcga_pan_can_atlas_2018": "endometrial",
    "coadread_tcga_pan_can_atlas_2018": "colorectal",
    "ov_tcga_pan_can_atlas_2018": "ovarian",
    "skcm_tcga_pan_can_atlas_2018": "melanoma"
}

STUDIES = list(STUDY_TO_CANCER_TYPE.keys())

# Output paths
OUTPUT_DIR = Path("data/validation")
OUTPUT_FILE = OUTPUT_DIR / "io_validation_cohort.json"


def extract_study_cohort(client: CBioportalClient, study_id: str, cancer_type: str) -> List[Dict]:
    """
    Extract all patients from a single study.
    
    Args:
        client: cBioPortal client instance
        study_id: cBioPortal study ID
        cancer_type: Cancer type name
    
    Returns:
        List of patient records
    """
    print(f"\n{'='*80}")
    print(f"📊 EXTRACTING: {cancer_type.upper()} ({study_id})")
    print('='*80)
    
    # Step 1: Fetch all patients
    print(f"📥 Fetching patients...")
    patients = client.get_study_patients(study_id)
    print(f"   ✅ Found {len(patients)} patients")
    
    if not patients:
        print(f"   ⚠️  No patients found, skipping study")
        return []
    
    # Step 2: Fetch all samples (for mutation extraction)
    print(f"📥 Fetching samples...")
    samples = client.get_study_samples(study_id)
    print(f"   ✅ Found {len(samples)} samples")
    
    # Step 3: Fetch patient-level clinical data (TMB, MSI)
    print(f"📥 Fetching clinical data...")
    patient_ids = [p["patientId"] for p in patients]
    clinical_data = client.get_clinical_data(patient_ids, entity_type="PATIENT")
    print(f"   ✅ Fetched clinical data for {len(clinical_data)} patients")
    
    # Step 4: Get mutation profile ID
    print(f"📥 Finding mutation profile...")
    mutation_profile_id = client.get_mutation_profile_id(study_id)
    if not mutation_profile_id:
        print(f"   ⚠️  No mutation profile found, mutations will be skipped")
    else:
        print(f"   ✅ Using mutation profile: {mutation_profile_id}")
    
    # Step 5: Extract patient records
    print(f"📥 Extracting patient records...")
    patient_records = []
    
    for i, patient in enumerate(patients):
        patient_id = patient.get("patientId")
        if not patient_id:
            continue
        
        # Get clinical data for this patient
        patient_clinical = clinical_data.get(patient_id, {})
        
        # Extract TMB (provided and calculated)
        tmb_provided = extract_tmb_from_clinical(patient_clinical)
        tmb_provided = normalize_tmb_value(tmb_provided)
        
        # Extract MSI
        msi_status = extract_msi_from_clinical(patient_clinical)
        
        # Find primary tumor sample
        primary_sample_id = find_primary_tumor_sample(samples, patient_id)
        
        # Fetch mutations if sample and profile available
        mutations = []
        if primary_sample_id and mutation_profile_id:
            try:
                raw_mutations = client.get_mutations_for_samples(mutation_profile_id, [primary_sample_id])
                mutations = process_mutations(raw_mutations)
            except Exception as e:
                print(f"   ⚠️  Failed to fetch mutations for {patient_id}: {e}")
        
        # Calculate TMB from mutations
        tmb_calculated = None
        if mutations:
            tmb_calculated = calculate_tmb_from_mutations(mutations, TCGA_EXOME_SIZE_MB)
            tmb_calculated = normalize_tmb_value(tmb_calculated)
        
        # Use provided TMB if available, otherwise calculated
        tmb_final = tmb_provided if tmb_provided is not None else tmb_calculated
        
        # Generate flags
        flags = generate_data_quality_flags(
            has_tmb=(tmb_final is not None),
            has_msi=(msi_status is not None),
            has_mutations=(len(mutations) > 0),
            msi_status=msi_status,
            tmb_value=tmb_final,
            has_mbd4=has_mbd4_mutation(mutations) if mutations else False,
            has_tp53=has_tp53_mutation(mutations) if mutations else False,
            has_ddr=has_ddr_mutation(mutations) if mutations else False
        )
        
        # Build patient record
        patient_record = {
            "patient_id": patient_id,
            "cancer_type": cancer_type,
            "source_study": study_id,
            "sample_id": primary_sample_id,  # For traceability
            "tmb_provided": tmb_provided,
            "tmb_calculated": tmb_calculated,
            "tmb_score": tmb_final,  # Final TMB (provided or calculated)
            "msi_status": msi_status,
            "mutations": mutations,
            "io_treatment": None,  # TCGA predates IO
            "io_response": None,
            "data_quality_flags": flags
        }
        
        # Validate record
        patient_record = validate_patient_record(patient_record)
        patient_records.append(patient_record)
        
        if (i + 1) % 50 == 0:
            print(f"   Processed {i + 1}/{len(patients)} patients...")
    
    print(f"   ✅ Extracted {len(patient_records)} patient records")
    return patient_records


def generate_summary(patients: List[Dict]) -> Dict:
    """Generate summary statistics for the cohort."""
    summary = {
        "total_patients": len(patients),
        "extraction_date": datetime.now().isoformat(),
        "by_cancer_type": {},
        "biomarker_coverage": {},
        "mutation_flags": {},
        "msi_distribution": Counter(),
        "tmb_distribution": {
            "tmb_high_count": 0,
            "tmb_low_count": 0,
            "tmb_missing": 0
        }
    }
    
    # By cancer type
    by_type = defaultdict(list)
    for p in patients:
        by_type[p["cancer_type"]].append(p)
    
    for cancer_type, type_patients in by_type.items():
        summary["by_cancer_type"][cancer_type] = {
            "count": len(type_patients),
            "with_tmb": sum(1 for p in type_patients if p.get("tmb_score") is not None),
            "with_msi": sum(1 for p in type_patients if p.get("msi_status") is not None),
            "with_mutations": sum(1 for p in type_patients if p.get("mutations")),
            "mbd4_mutant": sum(1 for p in type_patients if "mbd4_mutant" in p.get("data_quality_flags", [])),
            "tp53_mutant": sum(1 for p in type_patients if "tp53_mutant" in p.get("data_quality_flags", [])),
            "msi_high": sum(1 for p in type_patients if "msi_high" in p.get("data_quality_flags", [])),
            "tmb_high": sum(1 for p in type_patients if "tmb_high" in p.get("data_quality_flags", []))
        }
    
    # Overall biomarker coverage
    summary["biomarker_coverage"] = {
        "has_tmb": sum(1 for p in patients if p.get("tmb_score") is not None),
        "has_msi": sum(1 for p in patients if p.get("msi_status") is not None),
        "has_both": sum(1 for p in patients if p.get("tmb_score") is not None and p.get("msi_status") is not None),
        "has_mutations": sum(1 for p in patients if p.get("mutations"))
    }
    
    # Mutation flags
    summary["mutation_flags"] = {
        "mbd4_mutant": sum(1 for p in patients if "mbd4_mutant" in p.get("data_quality_flags", [])),
        "tp53_mutant": sum(1 for p in patients if "tp53_mutant" in p.get("data_quality_flags", [])),
        "ddr_mutant": sum(1 for p in patients if "ddr_mutant" in p.get("data_quality_flags", []))
    }
    
    # MSI distribution
    for p in patients:
        msi = p.get("msi_status", "Unknown")
        summary["msi_distribution"][msi] += 1
    
    # TMB distribution
    for p in patients:
        tmb = p.get("tmb_score")
        if tmb is None:
            summary["tmb_distribution"]["tmb_missing"] += 1
        elif tmb >= 10.0:
            summary["tmb_distribution"]["tmb_high_count"] += 1
        else:
            summary["tmb_distribution"]["tmb_low_count"] += 1
    
    return summary


def main():
    """Main extraction function."""
    print("\n" + "="*80)
    print("🧬 IO VALIDATION COHORT DATA ACQUISITION")
    print("="*80)
    print(f"Studies: {len(STUDIES)}")
    print(f"Target: 500+ patients across 4 cancer types")
    print("="*80 + "\n")
    
    # Initialize client
    client = CBioportalClient()
    
    # Extract from all studies
    all_patients = []
    for study_id in STUDIES:
        cancer_type = STUDY_TO_CANCER_TYPE[study_id]
        try:
            patients = extract_study_cohort(client, study_id, cancer_type)
            all_patients.extend(patients)
        except Exception as e:
            print(f"   ❌ ERROR extracting {study_id}: {e}")
            continue
    
    # Deduplicate
    print(f"\n{'='*80}")
    print("🔄 DEDUPLICATING PATIENTS")
    print('='*80)
    before_count = len(all_patients)
    all_patients = deduplicate_patients(all_patients)
    after_count = len(all_patients)
    print(f"   Before: {before_count} patients")
    print(f"   After: {after_count} patients")
    print(f"   Removed: {before_count - after_count} duplicates")
    
    # Generate summary
    print(f"\n{'='*80}")
    print("📊 GENERATING SUMMARY")
    print('='*80)
    summary = generate_summary(all_patients)
    
    # Print summary
    print(f"\nTotal Patients: {summary['total_patients']}")
    print(f"\nBy Cancer Type:")
    for cancer_type, stats in summary["by_cancer_type"].items():
        print(f"  {cancer_type}: {stats['count']} patients")
        print(f"    - TMB: {stats['with_tmb']} ({stats['with_tmb']/stats['count']*100:.1f}%)")
        print(f"    - MSI: {stats['with_msi']} ({stats['with_msi']/stats['count']*100:.1f}%)")
        print(f"    - Mutations: {stats['with_mutations']} ({stats['with_mutations']/stats['count']*100:.1f}%)")
        print(f"    - MSI-H: {stats['msi_high']} ({stats['msi_high']/stats['count']*100:.1f}%)")
        print(f"    - TMB-H: {stats['tmb_high']} ({stats['tmb_high']/stats['count']*100:.1f}%)")
        if stats['mbd4_mutant'] > 0:
            print(f"    - MBD4 mutant: {stats['mbd4_mutant']} ⭐")
        if stats['tp53_mutant'] > 0:
            print(f"    - TP53 mutant: {stats['tp53_mutant']}")
    
    print(f"\nBiomarker Coverage:")
    print(f"  TMB: {summary['biomarker_coverage']['has_tmb']} ({summary['biomarker_coverage']['has_tmb']/summary['total_patients']*100:.1f}%)")
    print(f"  MSI: {summary['biomarker_coverage']['has_msi']} ({summary['biomarker_coverage']['has_msi']/summary['total_patients']*100:.1f}%)")
    print(f"  Both: {summary['biomarker_coverage']['has_both']} ({summary['biomarker_coverage']['has_both']/summary['total_patients']*100:.1f}%)")
    
    print(f"\nMSI Distribution:")
    for status, count in summary['msi_distribution'].items():
        print(f"  {status}: {count}")
    
    print(f"\nTMB Distribution:")
    print(f"  TMB-H (≥10): {summary['tmb_distribution']['tmb_high_count']}")
    print(f"  TMB-L (<10): {summary['tmb_distribution']['tmb_low_count']}")
    print(f"  Missing: {summary['tmb_distribution']['tmb_missing']}")
    
    # Save output
    print(f"\n{'='*80}")
    print("💾 SAVING OUTPUT")
    print('='*80)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    output_data = {
        "metadata": {
            "source": "cBioPortal",
            "studies": STUDIES,
            "extraction_date": summary["extraction_date"],
            "total_patients": summary["total_patients"]
        },
        "summary": summary,
        "patients": all_patients
    }
    
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"   ✅ Saved {len(all_patients)} patients to {OUTPUT_FILE}")
    print(f"\n✅ EXTRACTION COMPLETE!")
    print(f"   Output: {OUTPUT_FILE}")
    print(f"   Next: Run validation script")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

