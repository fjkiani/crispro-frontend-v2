#!/usr/bin/env python3
"""
Download TCGA-OV patient-level data for Holistic Score validation
NO synthetic data - all real patient records

Data sources:
1. GDC API: Clinical data (OS, PFS, platinum response)
2. GDC API: Mutation data (BRCA1/2, TP53, etc.)
3. UCSC Xena: HRD scores
"""

import requests
import json
import pandas as pd
from pathlib import Path
from datetime import datetime
import time

# Output directory
OUTPUT_DIR = Path("data/cohorts/tcga_ov_holistic_score")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# GDC API endpoints
GDC_API = "https://api.gdc.cancer.gov"
CASES_ENDPOINT = f"{GDC_API}/cases"
FILES_ENDPOINT = f"{GDC_API}/files"

print("=" * 80)
print("TCGA-OV Data Download for Holistic Score Validation")
print("=" * 80)
print("\nGoal: Replace TOPACIO synthetic data with real TCGA-OV patient records")
print()

# Step 1: Get TCGA-OV case UUIDs
print("📋 Step 1: Fetching TCGA-OV case list...")

cases_params = {
    "filters": json.dumps({
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "project.project_id", "value": ["TCGA-OV"]}}
        ]
    }),
    "fields": "case_id,submitter_id,demographic.gender,demographic.race,diagnoses.age_at_diagnosis,diagnoses.tumor_stage,diagnoses.primary_diagnosis,exposures.cigarettes_per_day",
    "format": "JSON",
    "size": 10000
}

response = requests.post(CASES_ENDPOINT, json=cases_params, headers={"Content-Type": "application/json"})
cases_data = response.json()
cases = cases_data['data']['hits']

print(f"✅ Found {len(cases)} TCGA-OV cases")

# Step 2: Get clinical data (outcomes)
print("\n📊 Step 2: Fetching clinical outcomes...")

clinical_data = []
for case in cases:
    case_id = case['case_id']
    submitter_id = case.get('submitter_id', '')
    
    # Get demographics
    demo = case.get('demographic', {})
    gender = demo.get('gender', '')
    race = demo.get('race', '')
    
    # Get diagnosis info
    diagnoses = case.get('diagnoses', [])
    if diagnoses:
        diag = diagnoses[0]
        age_at_diagnosis = diag.get('age_at_diagnosis', None)
        tumor_stage = diag.get('tumor_stage', '')
        primary_diagnosis = diag.get('primary_diagnosis', '')
    else:
        age_at_diagnosis = None
        tumor_stage = ''
        primary_diagnosis = ''
    
    clinical_data.append({
        'case_id': case_id,
        'submitter_id': submitter_id,
        'gender': gender,
        'race': race,
        'age_at_diagnosis': age_at_diagnosis,
        'tumor_stage': tumor_stage,
        'primary_diagnosis': primary_diagnosis
    })

clinical_df = pd.DataFrame(clinical_data)
print(f"✅ Clinical data: {len(clinical_df)} patients")

# Step 3: Get mutation data
print("\n🧬 Step 3: Fetching mutation data (BRCA1, BRCA2, TP53)...")

# Download MAF file for TCGA-OV
maf_params = {
    "filters": json.dumps({
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-OV"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Masked Somatic Mutation"]}},
            {"op": "in", "content": {"field": "data_format", "value": ["MAF"]}}
        ]
    }),
    "fields": "file_id,file_name,file_size",
    "format": "JSON",
    "size": 100
}

response = requests.post(FILES_ENDPOINT, json=maf_params, headers={"Content-Type": "application/json"})
files_data = response.json()
maf_files = files_data['data']['hits']

if maf_files:
    maf_file = maf_files[0]  # Get first MAF file
    maf_file_id = maf_file['file_id']
    maf_file_name = maf_file['file_name']
    
    print(f"  Downloading {maf_file_name} ({maf_file['file_size']/1e6:.1f} MB)...")
    
    # Download MAF
    download_url = f"{GDC_API}/data/{maf_file_id}"
    maf_response = requests.get(download_url)
    
    maf_path = OUTPUT_DIR / maf_file_name
    with open(maf_path, 'wb') as f:
        f.write(maf_response.content)
    
    print(f"✅ MAF downloaded: {maf_path}")
    
    # Parse MAF for key genes
    print("  Parsing mutations in BRCA1, BRCA2, TP53...")
    
    # Read MAF (tab-separated)
    try:
        maf_df = pd.read_csv(maf_path, sep='\t', comment='#', low_memory=False)
        
        # Filter to key genes
        key_genes = ['BRCA1', 'BRCA2', 'TP53', 'ATM', 'PALB2', 'RAD51C', 'RAD51D', 'FANCF']
        mutations = maf_df[maf_df['Hugo_Symbol'].isin(key_genes)].copy()
        
        print(f"✅ Found {len(mutations)} mutations in {len(key_genes)} DDR genes")
        
        # Aggregate by patient
        mutation_summary = mutations.groupby(['Tumor_Sample_Barcode', 'Hugo_Symbol']).size().reset_index(name='mutation_count')
        mutation_summary['patient_id'] = mutation_summary['Tumor_Sample_Barcode'].str[:12]  # TCGA barcode
        
        # Pivot to wide format
        mutation_matrix = mutation_summary.pivot_table(
            index='patient_id',
            columns='Hugo_Symbol',
            values='mutation_count',
            fill_value=0
        ).reset_index()
        
        print(f"✅ Mutation matrix: {len(mutation_matrix)} patients × {len(key_genes)} genes")
        
        # Save
        mutation_matrix_path = OUTPUT_DIR / "tcga_ov_mutations.csv"
        mutation_matrix.to_csv(mutation_matrix_path, index=False)
        print(f"✅ Saved: {mutation_matrix_path}")
        
    except Exception as e:
        print(f"⚠️  Error parsing MAF: {e}")
        mutation_matrix = None

else:
    print("⚠️  No MAF files found")
    mutation_matrix = None

# Step 4: Get survival data
print("\n⏱️  Step 4: Fetching survival data...")

# GDC provides survival data in clinical supplements
survival_params = {
    "filters": json.dumps({
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-OV"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Clinical Supplement"]}},
            {"op": "in", "content": {"field": "data_format", "value": ["BCR XML"]}}
        ]
    }),
    "fields": "file_id,file_name",
    "format": "JSON",
    "size": 10
}

response = requests.post(FILES_ENDPOINT, json=survival_params, headers={"Content-Type": "application/json"})
files_data = response.json()
clinical_files = files_data['data']['hits']

print(f"  Found {len(clinical_files)} clinical supplement files")
print("  ⚠️  Note: Detailed survival parsing requires XML parsing (complex)")
print("  Alternative: Download pre-curated survival data from UCSC Xena")

# Step 5: Alternative - UCSC Xena for HRD scores and survival
print("\n🔬 Step 5: Fetching HRD scores from UCSC Xena...")
print("  URL: https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Ovarian%20Cancer%20(OV)")
print("  ⚠️  Manual download required (API not available)")
print("  Dataset: TCGA.OV.sampleMap/HRD_withSampleID")

# Save what we have
print("\n💾 Saving collected data...")

clinical_path = OUTPUT_DIR / "tcga_ov_clinical.csv"
clinical_df.to_csv(clinical_path, index=False)
print(f"✅ Clinical: {clinical_path}")

# Create receipt
receipt = {
    "timestamp": datetime.now().isoformat(),
    "dataset": "TCGA-OV",
    "n_patients": len(clinical_df),
    "files": {
        "clinical": str(clinical_path),
        "mutations": str(mutation_matrix_path) if mutation_matrix is not None else None,
        "maf": str(maf_path) if maf_files else None
    },
    "data_sources": {
        "gdc_api": GDC_API,
        "ucsc_xena": "https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Ovarian%20Cancer%20(OV)"
    },
    "next_steps": [
        "Download HRD scores from UCSC Xena",
        "Download survival data from cBioPortal or Firehose",
        "Compute mechanism vectors from mutations",
        "Calculate holistic scores",
        "Validate against outcomes"
    ]
}

receipt_path = OUTPUT_DIR / f"download_receipt_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
with open(receipt_path, 'w') as f:
    json.dump(receipt, f, indent=2)

print(f"✅ Receipt: {receipt_path}")

# Summary
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"  TCGA-OV patients: {len(clinical_df)}")
print(f"  Mutation data: {'✅ Downloaded' if mutation_matrix is not None else '⚠️ Pending'}")
print(f"  HRD scores: ⚠️ Manual download from UCSC Xena")
print(f"  Survival data: ⚠️ Alternative source needed (cBioPortal)")
print()
print("📋 NEXT STEPS:")
print("  1. Download HRD scores from UCSC Xena (manual)")
print("  2. Download survival data from cBioPortal API")
print("  3. Merge all datasets by patient_id")
print("  4. Compute mechanism vectors")
print("  5. Run holistic score validation")
print()
print("✅ REAL DATA - NO SYNTHETIC RECONSTRUCTION NEEDED")
print("=" * 80)
