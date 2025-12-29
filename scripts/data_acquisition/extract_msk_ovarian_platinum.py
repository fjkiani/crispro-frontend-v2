#!/usr/bin/env python3
"""
Extract MSK Ovarian Cancer Study for Platinum Response Validation

Extracts data from cBioPortal study: ovarian_msk_2025
- Clinical data (including platinum response if available)
- Somatic mutations with coordinates
- Normalizes to required schemas for DDR_bin validation

Deliverables:
1. Raw data in receipts folder
2. Normalized cBioPortal-style JSON
3. TRUE-SAE cohort JSON (if coordinates available)
"""

import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime
import hashlib

# Add project root to path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root))

# Try to use existing cBioPortal client
try:
    from scripts.data_acquisition.utils.cbioportal_client import CBioportalClient
    USE_CLIENT = True
except ImportError:
    USE_CLIENT = False
    print("⚠️  cBioPortal client not found, using direct API calls")

# Output directories
RAW_DIR = project_root / "data" / "external" / "ov_platinum_non_tcga" / "raw"
RAW_DIR.mkdir(parents=True, exist_ok=True)

BENCHMARK_DIR = project_root / "data" / "benchmarks"
BENCHMARK_DIR.mkdir(parents=True, exist_ok=True)

VALIDATION_DIR = project_root / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "validation" / "sae_cohort" / "checkpoints"
VALIDATION_DIR.mkdir(parents=True, exist_ok=True)

# Study configuration
STUDY_ID = "ovarian_msk_2025"
STUDY_NAME = "Serous Ovarian Cancer (MSK, 2025)"
CBIO_API_BASE = "https://www.cbioportal.org/api"

# Rate limiting
API_DELAY = 1.0


def compute_sha256(file_path: Path) -> str:
    """Compute SHA256 checksum of a file."""
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def fetch_study_info(study_id: str) -> Dict:
    """Fetch study metadata from cBioPortal API."""
    import requests
    
    print(f"📥 Fetching study info for {study_id}...")
    
    try:
        response = requests.get(
            f"{CBIO_API_BASE}/studies/{study_id}",
            timeout=30
        )
        
        if response.status_code == 200:
            return response.json()
        else:
            print(f"   ⚠️  API returned {response.status_code}")
            return {}
    except Exception as e:
        print(f"   ❌ Error: {e}")
        return {}


def fetch_patients(study_id: str) -> List[Dict]:
    """Fetch all patients from study."""
    import requests
    
    print(f"📥 Fetching patients...")
    
    try:
        response = requests.get(
            f"{CBIO_API_BASE}/studies/{study_id}/patients",
            timeout=60
        )
        
        if response.status_code == 200:
            patients = response.json()
            print(f"   ✅ Found {len(patients)} patients")
            return patients
        else:
            print(f"   ⚠️  API returned {response.status_code}")
            return []
    except Exception as e:
        print(f"   ❌ Error: {e}")
        return []


def fetch_clinical_data(study_id: str, entity_type: str = "PATIENT") -> List[Dict]:
    """Fetch clinical data for study."""
    import requests
    
    print(f"📥 Fetching clinical data ({entity_type})...")
    
    try:
        response = requests.get(
            f"{CBIO_API_BASE}/studies/{study_id}/clinical-data",
            params={
                "clinicalDataType": entity_type,
                "projection": "DETAILED"
            },
            timeout=120
        )
        
        if response.status_code == 200:
            clinical = response.json()
            print(f"   ✅ Fetched {len(clinical)} clinical data rows")
            return clinical
        else:
            print(f"   ⚠️  API returned {response.status_code}")
            return []
    except Exception as e:
        print(f"   ❌ Error: {e}")
        return []


def fetch_mutations(study_id: str) -> List[Dict]:
    """Fetch mutations for study."""
    import requests
    
    print(f"📥 Finding mutation profile...")
    
    # First, get molecular profiles
    try:
        response = requests.get(
            f"{CBIO_API_BASE}/studies/{study_id}/molecular-profiles",
            timeout=30
        )
        
        if response.status_code != 200:
            print(f"   ⚠️  Could not fetch molecular profiles: {response.status_code}")
            return []
        
        profiles = response.json()
        mutation_profile = None
        
        for profile in profiles:
            if "mutation" in profile.get("molecularAlterationType", "").lower():
                mutation_profile = profile
                break
        
        if not mutation_profile:
            print(f"   ⚠️  No mutation profile found")
            return []
        
        profile_id = mutation_profile["molecularProfileId"]
        print(f"   ✅ Found mutation profile: {profile_id}")
        
        # Get sample list
        response = requests.get(
            f"{CBIO_API_BASE}/studies/{study_id}/sample-lists",
            timeout=30
        )
        
        if response.status_code != 200:
            print(f"   ⚠️  Could not fetch sample lists: {response.status_code}")
            return []
        
        sample_lists = response.json()
        all_samples_list = None
        
        for sl in sample_lists:
            if "all" in sl.get("sampleListId", "").lower():
                all_samples_list = sl
                break
        
        if not all_samples_list:
            print(f"   ⚠️  No 'all' sample list found")
            return []
        
        sample_list_id = all_samples_list["sampleListId"]
        print(f"   ✅ Using sample list: {sample_list_id}")
        
        # Fetch mutations using mutations/fetch endpoint
        print(f"📥 Fetching mutations...")
        time.sleep(API_DELAY)
        
        # Try mutations/fetch endpoint (POST method)
        try:
            # First, get sample IDs from sample list
            response = requests.get(
                f"{CBIO_API_BASE}/sample-lists/{sample_list_id}",
                timeout=30
            )
            
            if response.status_code == 200:
                sample_list_data = response.json()
                sample_ids = sample_list_data.get("sampleIds", [])
                print(f"   ✅ Found {len(sample_ids)} samples in list")
                
                if sample_ids:
                    # Use mutations/fetch endpoint (POST)
                    time.sleep(API_DELAY)
                    response = requests.post(
                        f"{CBIO_API_BASE}/mutations/fetch",
                        json={
                            "molecularProfileId": profile_id,
                            "sampleIds": sample_ids,
                            "projection": "DETAILED"
                        },
                        timeout=300
                    )
                    
                    if response.status_code == 200:
                        mutations = response.json()
                        print(f"   ✅ Fetched {len(mutations)} mutations")
                        return mutations
                    else:
                        print(f"   ⚠️  Mutations/fetch returned {response.status_code}: {response.text[:200]}")
            else:
                print(f"   ⚠️  Could not fetch sample list: {response.status_code}")
        except Exception as e:
            print(f"   ⚠️  Error fetching mutations: {e}")
        
        # Fallback: try GET endpoint with sample IDs
        print(f"   Trying alternative mutation endpoint...")
        time.sleep(API_DELAY)
        
        response = requests.get(
            f"{CBIO_API_BASE}/molecular-profiles/{profile_id}/mutations",
            params={
                "sampleListId": sample_list_id,
                "projection": "DETAILED",
                "pageSize": 100000
            },
            timeout=300
        )
        
        if response.status_code == 200:
            mutations = response.json()
            print(f"   ✅ Fetched {len(mutations)} mutations (alternative endpoint)")
            return mutations
        else:
            print(f"   ⚠️  Alternative endpoint returned {response.status_code}")
            return []
            
    except Exception as e:
        print(f"   ❌ Error: {e}")
        return []


def map_platinum_response(clinical_value: Any) -> Optional[str]:
    """
    Map clinical data value to platinum response label.
    
    Returns: 'sensitive', 'resistant', 'refractory', or None
    """
    if not clinical_value or not isinstance(clinical_value, str):
        return None
    
    value_lower = clinical_value.lower().strip()
    
    # Sensitive patterns
    if any(term in value_lower for term in ["sensitive", "responder", "response", "cr", "pr"]):
        if "refractory" in value_lower or "resistant" in value_lower:
            return "refractory"
        return "sensitive"
    
    # Resistant patterns
    if any(term in value_lower for term in ["resistant", "non-responder", "progression", "pd"]):
        return "resistant"
    
    # Refractory patterns
    if "refractory" in value_lower:
        return "refractory"
    
    return None


def normalize_to_cbioportal_schema(
    patients: List[Dict],
    clinical_data: List[Dict],
    mutations: List[Dict]
) -> Dict:
    """
    Normalize extracted data to cBioPortal-style schema.
    
    Returns study dict with patients, clinical_outcomes, mutations.
    """
    print(f"\n📊 Normalizing data to cBioPortal schema...")
    
    # Organize clinical data by patient
    clinical_by_patient = {}
    for row in clinical_data:
        patient_id = row.get("patientId", "")
        if not patient_id:
            continue
        
        if patient_id not in clinical_by_patient:
            clinical_by_patient[patient_id] = {}
        
        attr_id = row.get("clinicalAttributeId", "")
        value = row.get("value", "")
        clinical_by_patient[patient_id][attr_id] = value
    
    # Organize mutations by patient
    mutations_by_patient = {}
    for mut in mutations:
        patient_id = mut.get("patientId", "")
        if not patient_id:
            continue
        
        if patient_id not in mutations_by_patient:
            mutations_by_patient[patient_id] = []
        
        # Extract gene (handle both dict and string formats)
        gene_obj = mut.get("gene", {})
        gene = gene_obj.get("hugoGeneSymbol", "") if isinstance(gene_obj, dict) else str(gene_obj) if gene_obj else ""
        
        # Extract chromosome (cBioPortal uses 'chr' field)
        chromosome = mut.get("chr") or mut.get("chromosome", "")
        if chromosome:
            chromosome = str(chromosome).replace("chr", "").replace("Chr", "").replace("CHR", "")
        
        mutations_by_patient[patient_id].append({
            "gene": gene,
            "chromosome": chromosome,
            "start_position": mut.get("startPosition"),
            "end_position": mut.get("endPosition"),
            "reference_allele": mut.get("referenceAllele", ""),
            "variant_allele": mut.get("variantAllele", ""),
            "protein_change": mut.get("proteinChange", ""),
            "variant_type": mut.get("mutationType", "")
        })
    
    # Build patient records
    patient_records = []
    platinum_response_count = 0
    
    for patient in patients:
        patient_id = patient.get("patientId", "")
        if not patient_id:
            continue
        
        clinical = clinical_by_patient.get(patient_id, {})
        patient_mutations = mutations_by_patient.get(patient_id, [])
        
        # Extract platinum response
        platinum_response = None
        for key in ["PLATINUM_RESPONSE", "PLATINUM_SENSITIVITY", "CHEMO_RESPONSE", "RESPONSE"]:
            if key in clinical:
                platinum_response = map_platinum_response(clinical[key])
                if platinum_response:
                    platinum_response_count += 1
                    break
        
        # Build clinical outcomes
        clinical_outcomes = {}
        
        # Map common fields
        if "OS_MONTHS" in clinical:
            clinical_outcomes["OS_MONTHS"] = float(clinical["OS_MONTHS"]) if clinical["OS_MONTHS"] else None
        if "OS_STATUS" in clinical:
            clinical_outcomes["OS_STATUS"] = clinical["OS_STATUS"]
        if "PFS_MONTHS" in clinical:
            clinical_outcomes["PFS_MONTHS"] = float(clinical["PFS_MONTHS"]) if clinical["PFS_MONTHS"] else None
        if "PFS_STATUS" in clinical:
            clinical_outcomes["PFS_STATUS"] = clinical["PFS_STATUS"]
        if "AGE" in clinical:
            clinical_outcomes["AGE"] = int(clinical["AGE"]) if clinical["AGE"] else None
        if "STAGE" in clinical or "CLINICAL_STAGE" in clinical:
            clinical_outcomes["CLINICAL_STAGE"] = clinical.get("STAGE") or clinical.get("CLINICAL_STAGE")
        
        # Build mutations list (with coordinates if available)
        # Note: mutations_by_patient already has full mutation dicts with coordinates
        mutations_list = patient_mutations  # Already normalized in mutations_by_patient
        
        patient_records.append({
            "patient_id": patient_id,
            "clinical_outcomes": clinical_outcomes,
            "mutations": mutations_list,
            "platinum_response": platinum_response  # Add to patient record
        })
    
    print(f"   ✅ Normalized {len(patient_records)} patients")
    print(f"   ✅ Found {platinum_response_count} patients with platinum response labels")
    
    return {
        "study_id": f"external_ov_platinum_{STUDY_ID}",
        "patients": patient_records
    }


def save_raw_data(study_info: Dict, patients: List[Dict], clinical: List[Dict], mutations: List[Dict]):
    """Save raw extracted data to receipts folder."""
    print(f"\n💾 Saving raw data to receipts folder...")
    
    # Save study info
    study_file = RAW_DIR / f"{STUDY_ID}_study_info.json"
    with open(study_file, "w") as f:
        json.dump(study_info, f, indent=2)
    print(f"   ✅ Saved: {study_file.name}")
    
    # Save patients
    patients_file = RAW_DIR / f"{STUDY_ID}_patients.json"
    with open(patients_file, "w") as f:
        json.dump(patients, f, indent=2)
    print(f"   ✅ Saved: {patients_file.name}")
    
    # Save clinical data
    clinical_file = RAW_DIR / f"{STUDY_ID}_clinical_data.json"
    with open(clinical_file, "w") as f:
        json.dump(clinical, f, indent=2)
    print(f"   ✅ Saved: {clinical_file.name}")
    
    # Save mutations
    mutations_file = RAW_DIR / f"{STUDY_ID}_mutations.json"
    with open(mutations_file, "w") as f:
        json.dump(mutations, f, indent=2)
    print(f"   ✅ Saved: {mutations_file.name}")
    
    # Compute checksums
    checksums = {}
    for file in [study_file, patients_file, clinical_file, mutations_file]:
        checksums[file.name] = compute_sha256(file)
    
    # Save checksums
    checksum_file = RAW_DIR / f"{STUDY_ID}_checksums.json"
    with open(checksum_file, "w") as f:
        json.dump(checksums, f, indent=2)
    print(f"   ✅ Saved: {checksum_file.name}")
    
    return {
        "files": [f.name for f in [study_file, patients_file, clinical_file, mutations_file]],
        "checksums": checksums
    }


def main():
    """Main extraction workflow."""
    print("="*80)
    print("MSK OVARIAN CANCER STUDY EXTRACTION")
    print("="*80)
    print(f"Study: {STUDY_NAME} ({STUDY_ID})")
    print()
    
    # Step 1: Fetch study info
    study_info = fetch_study_info(STUDY_ID)
    time.sleep(API_DELAY)
    
    # Step 2: Fetch patients
    patients = fetch_patients(STUDY_ID)
    time.sleep(API_DELAY)
    
    if not patients:
        print("❌ No patients found. Exiting.")
        return
    
    # Step 3: Fetch clinical data
    clinical_data = fetch_clinical_data(STUDY_ID, entity_type="PATIENT")
    time.sleep(API_DELAY)
    
    # Step 4: Fetch mutations
    mutations = fetch_mutations(STUDY_ID)
    time.sleep(API_DELAY)
    
    # Step 5: Save raw data
    raw_metadata = save_raw_data(study_info, patients, clinical_data, mutations)
    
    # Step 6: Normalize to cBioPortal schema
    normalized_study = normalize_to_cbioportal_schema(patients, clinical_data, mutations)
    
    # Step 7: Save normalized dataset
    output_file = BENCHMARK_DIR / f"external_ov_platinum_{STUDY_ID}.json"
    with open(output_file, "w") as f:
        json.dump([normalized_study], f, indent=2)
    print(f"\n✅ Saved normalized dataset: {output_file}")
    
    # Step 8: Update README with extraction metadata
    readme_file = RAW_DIR / "README.md"
    readme_content = f"""# External Ovarian Cancer Platinum Response Dataset

**Source:** {STUDY_NAME}  
**Study ID:** {STUDY_ID}  
**Extraction Date:** {datetime.now().isoformat()}

## Source Information

**Citation:** [To be added - check cBioPortal study page]  
**URL:** https://www.cbioportal.org/study/summary?id={STUDY_ID}  
**License:** [To be verified]

## Extracted Data

**Patients:** {len(patients)}  
**Clinical Data Rows:** {len(clinical_data)}  
**Mutations:** {len(mutations)}  
**Platinum Response Labels:** {sum(1 for p in normalized_study['patients'] if p.get('platinum_response'))}

## Files

{chr(10).join(f"- `{f}`" for f in raw_metadata['files'])}

## Checksums (SHA256)

```json
{json.dumps(raw_metadata['checksums'], indent=2)}
```

## Data Quality Notes

- Platinum response mapping: Check clinical data fields for response labels
- Variant coordinates: Verify chromosome/position/ref/alt availability
- Assembly: [To be verified - check mutation data]

## Next Steps

1. Verify platinum response labels are correctly mapped
2. Check variant coordinate availability for TRUE-SAE extraction
3. Run TRUE-SAE extraction if coordinates available
4. Run validator: `validate_ddr_bin_generic_cohort.py`
"""
    
    with open(readme_file, "w") as f:
        f.write(readme_content)
    print(f"✅ Updated README: {readme_file}")
    
    # Summary
    print("\n" + "="*80)
    print("EXTRACTION SUMMARY")
    print("="*80)
    print(f"✅ Patients extracted: {len(patients)}")
    print(f"✅ Clinical data rows: {len(clinical_data)}")
    print(f"✅ Mutations: {len(mutations)}")
    print(f"✅ Platinum response labels: {sum(1 for p in normalized_study['patients'] if p.get('platinum_response'))}")
    print(f"✅ Normalized dataset: {output_file}")
    print()
    print("Next Steps:")
    print("1. Verify platinum response labels in normalized dataset")
    print("2. Check if mutations have coordinates for TRUE-SAE extraction")
    print("3. If coordinates available, run TRUE-SAE extraction")
    print("4. Run validator on extracted cohort")
    print()


if __name__ == "__main__":
    main()

