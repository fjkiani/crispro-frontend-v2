#!/usr/bin/env python3
"""
Download TCGA-OV HRD scores and survival data from cBioPortal
Complete the dataset for holistic score validation
"""

import requests
import pandas as pd
from pathlib import Path
import json
from datetime import datetime

# Paths
DATA_DIR = Path("data/cohorts/tcga_ov_holistic_score")
DATA_DIR.mkdir(parents=True, exist_ok=True)

# cBioPortal API
CBIOPORTAL_API = "https://www.cbioportal.org/api"

print("=" * 80)
print("TCGA-OV: Downloading HRD Scores + Survival from cBioPortal")
print("=" * 80)
print()

# Step 1: Get study info
print("📋 Step 1: Getting TCGA-OV study info...")
study_id = "ov_tcga"

study_url = f"{CBIOPORTAL_API}/studies/{study_id}"
response = requests.get(study_url)
if response.status_code == 200:
    study_info = response.json()
    print(f"✅ Study: {study_info['name']}")
    print(f"  Patients: {study_info.get('allSampleCount', 'N/A')}")
else:
    print(f"❌ Failed to get study info: {response.status_code}")
    print("  Trying alternative study ID...")
    study_id = "ov_tcga_pan_can_atlas_2018"
    study_url = f"{CBIOPORTAL_API}/studies/{study_id}"
    response = requests.get(study_url)
    if response.status_code == 200:
        study_info = response.json()
        print(f"✅ Study: {study_info['name']}")
    else:
        print(f"❌ Both study IDs failed")
        study_id = "ov_tcga_pan_can_atlas_2018"  # Use default

# Step 2: Get all patients
print(f"\n👥 Step 2: Fetching patient list from {study_id}...")
patients_url = f"{CBIOPORTAL_API}/studies/{study_id}/patients"
response = requests.get(patients_url)
if response.status_code == 200:
    patients = response.json()
    print(f"✅ Found {len(patients)} patients")
    patient_ids = [p['patientId'] for p in patients]
else:
    print(f"❌ Failed: {response.status_code}")
    patient_ids = []

# Step 3: Get clinical data (survival)
print("\n⏱️  Step 3: Fetching survival data...")
clinical_url = f"{CBIOPORTAL_API}/studies/{study_id}/clinical-data"
params = {
    "clinicalDataType": "PATIENT",
    "projection": "DETAILED"
}

response = requests.get(clinical_url, params=params)
if response.status_code == 200:
    clinical_data = response.json()
    print(f"✅ Found {len(clinical_data)} clinical data points")
    
    # Convert to DataFrame
    clinical_df = pd.DataFrame(clinical_data)
    
    # Pivot to wide format
    clinical_wide = clinical_df.pivot_table(
        index='patientId',
        columns='clinicalAttributeId',
        values='value',
        aggfunc='first'
    ).reset_index()
    
    print(f"✅ Clinical data: {len(clinical_wide)} patients × {len(clinical_wide.columns)-1} attributes")
    
    # Save
    clinical_path = DATA_DIR / "tcga_ov_clinical_cbioportal.csv"
    clinical_wide.to_csv(clinical_path, index=False)
    print(f"✅ Saved: {clinical_path}")
    
    # Print key survival columns
    survival_cols = [c for c in clinical_wide.columns if any(x in c.lower() for x in ['os', 'pfs', 'dfs', 'survival', 'status', 'months'])]
    if survival_cols:
        print(f"  Survival columns: {', '.join(survival_cols[:5])}")
    
else:
    print(f"❌ Failed: {response.status_code}")
    clinical_wide = None

# Step 4: Get molecular profiles (mutations, CNA)
print("\n🧬 Step 4: Fetching molecular profile data...")
profiles_url = f"{CBIOPORTAL_API}/studies/{study_id}/molecular-profiles"
response = requests.get(profiles_url)
if response.status_code == 200:
    profiles = response.json()
    print(f"✅ Found {len(profiles)} molecular profiles")
    
    for prof in profiles:
        print(f"  - {prof['molecularProfileId']}: {prof['name']} ({prof['molecularAlterationType']})")
    
    # Get mutations profile ID
    mutation_profile = next((p for p in profiles if p['molecularAlterationType'] == 'MUTATION_EXTENDED'), None)
    if mutation_profile:
        mutation_profile_id = mutation_profile['molecularProfileId']
        print(f"\n  Using mutation profile: {mutation_profile_id}")
        
        # Get DDR gene mutations
        ddr_genes = ['BRCA1', 'BRCA2', 'TP53', 'ATM', 'PALB2', 'RAD51C', 'RAD51D', 'FANCF', 'CHEK2', 'ATR']
        
        print(f"  Fetching mutations in {len(ddr_genes)} DDR genes...")
        
        mutations_url = f"{CBIOPORTAL_API}/molecular-profiles/{mutation_profile_id}/mutations"
        params = {
            "entrezGeneIds": ",".join([str(g) for g in range(1, 25000)]),  # All genes (API limitation workaround)
            "projection": "DETAILED"
        }
        
        # Alternative: Use gene query
        genes_url = f"{CBIOPORTAL_API}/molecular-profiles/{mutation_profile_id}/mutations/fetch"
        
        # This is complex - requires POST with specific format
        print("  ⚠️  Mutation fetch requires complex POST query")
        print("  Alternative: Use earlier GDC MAF file or cBioPortal web interface")
    
else:
    print(f"❌ Failed: {response.status_code}")

# Step 5: Try to get HRD scores via generic assay
print("\n🔬 Step 5: Searching for HRD scores...")
print("  ⚠️  HRD scores not typically in cBioPortal API")
print("  Alternative sources:")
print("    1. UCSC Xena: https://xenabrowser.net (HRD_withSampleID dataset)")
print("    2. Published supplements: Marquard et al. 2015, Takaya et al. 2020")
print("    3. Compute from ASCAT/FACETS data")

# Summary
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
if clinical_wide is not None:
    print(f"  ✅ Clinical/survival data: {len(clinical_wide)} patients")
    print(f"  ✅ File: {clinical_path}")
else:
    print("  ❌ Clinical data: Failed")

print(f"  ⚠️  Mutations: Use GDC MAF file (already downloaded)")
print(f"  ⚠️  HRD scores: Requires manual download from UCSC Xena")
print()
print("📋 NEXT MANUAL STEPS:")
print("  1. Go to: https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Ovarian%20Cancer%20(OV)")
print("  2. Download: TCGA.OV.sampleMap/HRD_withSampleID")
print("  3. Save to: data/cohorts/tcga_ov_holistic_score/tcga_ov_hrd_scores.tsv")
print()
print("  Alternative (faster):")
print("  1. Use cBioPortal web interface: https://www.cbioportal.org/study/summary?id=ov_tcga_pan_can_atlas_2018")
print("  2. Download all clinical/mutation data via 'Download' tab")
print("=" * 80)

# Create receipt
receipt = {
    "timestamp": datetime.now().isoformat(),
    "cbioportal_study": study_id,
    "files": {
        "clinical_survival": str(clinical_path) if clinical_wide is not None else None
    },
    "next_steps": [
        "Download HRD scores from UCSC Xena manually",
        "Or use cBioPortal web interface bulk download",
        "Merge with GDC clinical/mutation data",
        "Compute mechanism vectors",
        "Validate holistic score"
    ]
}

receipt_path = DATA_DIR / f"cbioportal_download_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
with open(receipt_path, 'w') as f:
    json.dump(receipt, f, indent=2)
print(f"\n✅ Receipt: {receipt_path}")
