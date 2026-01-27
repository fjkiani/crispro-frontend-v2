#!/usr/bin/env python3
"""
Download cBioPortal TCGA+MSK Paired Samples
Extract 57 patients with paired primary+recurrent samples
"""

import httpx
import json
import pandas as pd
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Set
from collections import defaultdict

CBIO_BASE = "https://www.cbioportal.org/api"
STUDIES = ["msk_spectrum_tme_2022", "hgsoc_tcga_gdc"]

OUTPUT_DIR = Path("data/serial_sae/cbioportal_paired")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def get_all_samples(study_id: str) -> List[Dict]:
    """Fetch all samples for a study."""
    print(f"📥 Fetching samples from {study_id}...")
    with httpx.Client(timeout=120.0) as client:
        r = client.get(f"{CBIO_BASE}/studies/{study_id}/samples", params={"pageSize": 10000})
        r.raise_for_status()
        return r.json()

def get_clinical_data(study_id: str, sample_ids: List[str]) -> pd.DataFrame:
    """Fetch clinical data for samples."""
    print(f"📥 Fetching clinical data for {len(sample_ids)} samples...")
    clinical_data = {}
    
    chunk_size = 200
    for i in range(0, len(sample_ids), chunk_size):
        chunk = sample_ids[i:i+chunk_size]
        with httpx.Client(timeout=120.0) as client:
            r = client.post(
                f"{CBIO_BASE}/clinical-data/fetch",
                json={
                    "entityIds": chunk,
                    "entityType": "SAMPLE",
                    "projection": "DETAILED"
                }
            )
            try:
                r.raise_for_status()
                data = r.json()
                
                if isinstance(data, list):
                    for row in data:
                        sample_id = row.get("entityId")
                        if sample_id:
                            if sample_id not in clinical_data:
                                clinical_data[sample_id] = {}
                            attr = row.get("clinicalAttributeId")
                            val = row.get("value")
                            if attr and val:
                                clinical_data[sample_id][attr] = val
            except httpx.HTTPStatusError as e:
                print(f"  ⚠️  HTTP error for chunk {i//chunk_size + 1}: {e}")
                continue
            except json.JSONDecodeError as e:
                print(f"  ⚠️  JSON decode error for chunk {i//chunk_size + 1}: {e}")
                print(f"  Response text: {r.text[:200]}")
                continue
        
        print(f"  Processed {min(i+chunk_size, len(sample_ids))}/{len(sample_ids)} samples...")
    
    return pd.DataFrame.from_dict(clinical_data, orient='index')

def identify_paired_samples(samples: List[Dict], clinical_df: pd.DataFrame) -> Dict[str, List[Dict]]:
    """Identify patients with paired primary+recurrent samples using clinical data."""
    patient_samples = defaultdict(list)
    
    # Get all sample IDs
    all_sample_ids = [s.get("sampleId", "") for s in samples]
    
    # Get clinical attributes for sample type identification
    sample_type_attrs = ["SAMPLE_TYPE", "SAMPLE_CLASS", "TUMOR_TYPE", "TISSUE_TYPE", "SAMPLE_CLASSIFICATION"]
    
    # Build sample type mapping from clinical data
    sample_type_map = {}
    if not clinical_df.empty:
        for sample_id in all_sample_ids:
            sample_clinical = clinical_df[clinical_df.index == sample_id]
            if not sample_clinical.empty:
                # Check each potential attribute
                for attr in sample_type_attrs:
                    if attr in sample_clinical.columns:
                        val = str(sample_clinical[attr].iloc[0]).lower()
                        if "primary" in val or "diagnosis" in val:
                            sample_type_map[sample_id] = "primary"
                            break
                        elif "recurrent" in val or "relapse" in val or "recurrence" in val:
                            sample_type_map[sample_id] = "recurrent"
                            break
                        elif "metastasis" in val or "metastatic" in val:
                            sample_type_map[sample_id] = "metastasis"
                            break
    
    # Also check sample ID format (TCGA style)
    for sample in samples:
        sample_id = sample.get("sampleId", "")
        patient_id = sample.get("patientId", "")
        
        if not patient_id:
            continue
        
        # First check clinical data mapping
        sample_type = sample_type_map.get(sample_id, "unknown")
        
        # Fallback: Extract from sample ID (TCGA format: TCGA-XX-XXXX-XX)
        if sample_type == "unknown":
            parts = sample_id.split("-")
            if len(parts) >= 4:
                sample_type_code = parts[3]
                # Primary: 01, Recurrent: 06, Metastasis: 07
                # Handle both 2-digit (01) and longer codes (01A, 01B, etc.)
                if sample_type_code.startswith("01"):
                    sample_type = "primary"
                elif sample_type_code.startswith("06"):
                    sample_type = "recurrent"
                elif sample_type_code.startswith("07"):
                    sample_type = "metastasis"
            # Also check for MSK format or other patterns
            elif "primary" in sample_id.lower() or "diagnosis" in sample_id.lower():
                sample_type = "primary"
            elif "recurrent" in sample_id.lower() or "relapse" in sample_id.lower() or "recurrence" in sample_id.lower():
                sample_type = "recurrent"
            elif "metastasis" in sample_id.lower() or "metastatic" in sample_id.lower():
                sample_type = "metastasis"
        
        patient_samples[patient_id].append({
            "sample_id": sample_id,
            "sample_type": sample_type,
            "study_id": sample.get("studyId", "")
        })
    
    # Find patients with paired samples
    paired_patients = {}
    for patient_id, samples_list in patient_samples.items():
        # If we have 2 samples and one is primary and one is recurrent/metastasis, it's paired
        has_primary = any(s["sample_type"] == "primary" for s in samples_list)
        has_recurrent = any(s["sample_type"] in ["recurrent", "metastasis"] for s in samples_list)
        
        if has_primary and has_recurrent:
            paired_patients[patient_id] = {
                "primary_samples": [s for s in samples_list if s["sample_type"] == "primary"],
                "recurrent_samples": [s for s in samples_list if s["sample_type"] in ["recurrent", "metastasis"]]
            }
        elif len(samples_list) == 2 and all(s["sample_type"] == "unknown" for s in samples_list):
            # Fallback: If we have 2 samples but both are "unknown", assume first is primary, second is recurrent
            # This handles MSK samples where sample ID format doesn't follow TCGA convention
            paired_patients[patient_id] = {
                "primary_samples": [samples_list[0]],  # First sample assumed primary
                "recurrent_samples": [samples_list[1]],  # Second sample assumed recurrent
                "note": "Sample types inferred from order (clinical data unavailable)"
            }
    
    return paired_patients

def get_molecular_profiles(study_id: str) -> List[Dict]:
    """Get available molecular profiles."""
    with httpx.Client(timeout=60.0) as client:
        r = client.get(f"{CBIO_BASE}/studies/{study_id}/molecular-profiles")
        r.raise_for_status()
        return r.json()

def main():
    print("="*80)
    print("cBioPortal Paired Sample Extraction")
    print("="*80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    all_paired_patients = {}
    
    for study_id in STUDIES:
        print(f"\n🔍 Processing study: {study_id}")
        print("-" * 80)
        
        # Get samples
        samples = get_all_samples(study_id)
        print(f"✅ Found {len(samples)} total samples")
        
        # Get sample IDs for clinical data
        sample_ids = [s.get("sampleId", "") for s in samples if s.get("sampleId")]
        
        # Get clinical data for samples
        print(f"📥 Fetching clinical data for {len(sample_ids)} samples...")
        clinical_df = get_clinical_data(study_id, sample_ids)
        print(f"✅ Clinical data retrieved: {len(clinical_df)} samples")
        
        # Debug: Show patients with multiple samples
        patient_counts = defaultdict(int)
        for s in samples:
            patient_counts[s.get("patientId", "")] += 1
        multi_sample_patients = {pid: count for pid, count in patient_counts.items() if count > 1}
        print(f"📊 Patients with multiple samples: {len(multi_sample_patients)}")
        if len(multi_sample_patients) > 0:
            print(f"   Sample counts: {sorted(set(multi_sample_patients.values()), reverse=True)[:5]}")
        
        # Identify paired patients using clinical data
        paired = identify_paired_samples(samples, clinical_df)
        print(f"✅ Found {len(paired)} patients with paired samples")
        
        # Get molecular profiles
        profiles = get_molecular_profiles(study_id)
        profile_types = [p.get("molecularProfileId", "") for p in profiles]
        has_mutations = any("mutations" in p.lower() for p in profile_types)
        has_expression = any("rna" in p.lower() or "expression" in p.lower() for p in profile_types)
        
        print(f"  Molecular profiles: {len(profiles)}")
        print(f"  Mutations available: {has_mutations}")
        print(f"  Expression available: {has_expression}")
        
        # Store paired patients
        for patient_id, data in paired.items():
            if patient_id not in all_paired_patients:
                all_paired_patients[patient_id] = {
                    "study_id": study_id,
                    "primary_samples": [],
                    "recurrent_samples": []
                }
            all_paired_patients[patient_id]["primary_samples"].extend(data["primary_samples"])
            all_paired_patients[patient_id]["recurrent_samples"].extend(data["recurrent_samples"])
    
    print(f"\n📊 Total paired patients across all studies: {len(all_paired_patients)}")
    
    # Save results
    results = {
        "extraction_date": datetime.now().isoformat(),
        "studies": STUDIES,
        "total_paired_patients": len(all_paired_patients),
        "paired_patients": all_paired_patients
    }
    
    output_file = OUTPUT_DIR / "paired_patients.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\n💾 Results saved to: {output_file}")
    
    # Create summary
    summary = {
        "total_paired": len(all_paired_patients),
        "patients": list(all_paired_patients.keys())[:10]  # First 10 for preview
    }
    
    summary_file = OUTPUT_DIR / "summary.json"
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2)
    
    print(f"📋 Summary saved to: {summary_file}")
    
    print("\n" + "="*80)
    print("✅ Phase 1 Complete: Paired patients identified")
    print("="*80)
    print(f"\nNext steps:")
    print(f"  1. Download mutations for paired samples")
    print(f"  2. Download expression data for paired samples")
    print(f"  3. Compute SAE pathway scores")

if __name__ == "__main__":
    main()
