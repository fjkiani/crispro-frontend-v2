#!/usr/bin/env python3
"""
FIXED: Download TCGA-OV mutations from cBioPortal
Fetches mutations in batches to ensure we get data for all patients
"""

import requests
import pandas as pd
from pathlib import Path
import json
import time

# Paths
DATA_DIR = Path("data/cohorts/tcga_ov_holistic_score")
DATA_DIR.mkdir(parents=True, exist_ok=True)

# cBioPortal API
CBIOPORTAL_API = "https://www.cbioportal.org/api"
STUDY_ID = "ov_tcga"
PROFILE_ID = "ov_tcga_mutations"

print("=" * 80)
print("TCGA-OV: Force-Downloading Mutations from cBioPortal")
print("=" * 80)

# Step 1: Get patient list again
print("📋 Step 1: Fetching sample list...")
samples_url = f"{CBIOPORTAL_API}/studies/{STUDY_ID}/samples"
response = requests.get(samples_url)
samples = response.json()
sample_ids = [s['sampleId'] for s in samples]
print(f"✅ Found {len(sample_ids)} samples")

# Step 2: Define gene list (DDR + Core Pathways)
genes = [
    'BRCA1', 'BRCA2', 'TP53', 'ATM', 'PALB2', 'RAD51C', 'RAD51D', 'FANCF', 'CHEK2', 'ATR', # DDR
    'KRAS', 'NRAS', 'BRAF', 'MAPK1', # MAPK
    'PIK3CA', 'PTEN', 'AKT1', # PI3K
    'CD274' # PD-L1 (IO)
]

print(f"\n🧬 Step 2: Fetching mutations for {len(genes)} genes in batches...")

# Fetch molecular profile data
# We use the mutation fetch endpoint with a filter
mutations_url = f"{CBIOPORTAL_API}/molecular-profiles/{PROFILE_ID}/mutations/fetch"

all_mutations = []
batch_size = 1000 # Entrez IDs or just gene symbols if supported, but typically we need Entrez IDs

# Map symbols to Entrez IDs (hardcoded for speed)
entrez_map = {
    'BRCA1': 672, 'BRCA2': 675, 'TP53': 7157, 'ATM': 472, 'PALB2': 79728, 
    'RAD51C': 5889, 'RAD51D': 5892, 'FANCF': 2188, 'CHEK2': 11200, 'ATR': 545,
    'KRAS': 3845, 'NRAS': 4893, 'BRAF': 673, 'MAPK1': 5594,
    'PIK3CA': 5290, 'PTEN': 5728, 'AKT1': 207, 'CD274': 29126
}
entrez_ids = list(entrez_map.values())

# Payload matches cBioPortal API expectation
payload = {
    "entrezGeneIds": entrez_ids,
    "sampleIds": sample_ids
}

try:
    print(f"  Querying {len(sample_ids)} samples x {len(genes)} genes...")
    response = requests.post(mutations_url, json=payload)
    
    if response.status_code == 200:
        muts = response.json()
        print(f"✅ Received {len(muts)} mutation records")
        all_mutations.extend(muts)
    else:
        print(f"❌ Failed: {response.status_code} - {response.text}")
        
except Exception as e:
    print(f"❌ Error: {e}")

# Process and save
if all_mutations:
    mut_df = pd.DataFrame(all_mutations)
    
    # Keep useful columns
    cols = ['sampleId', 'gene', 'proteinChange', 'mutationType']
    # Adjust for API response format (nested objects usually flat here)
    if 'entrezGeneId' in mut_df.columns:
        # Map back to symbol
        inv_map = {v: k for k, v in entrez_map.items()}
        mut_df['Hugo_Symbol'] = mut_df['entrezGeneId'].map(inv_map)
    
    output_df = mut_df[['sampleId', 'Hugo_Symbol', 'proteinChange', 'mutationType']].copy()
    output_df = output_df.rename(columns={'sampleId': 'patient_barcode'}) # Close enough for join
    
    # Clean barcodes (TCGA-XX-XXXX-01)
    output_df['patient_barcode'] = output_df['patient_barcode'].apply(lambda x: '-'.join(x.split('-')[:3]))
    
    # Pivot to matrix
    pivot = output_df.pivot_table(
        index='patient_barcode', 
        columns='Hugo_Symbol', 
        values='proteinChange', 
        aggfunc='count'
    ).fillna(0)
    pivot = pivot.reset_index()
    
    # Save
    out_path = DATA_DIR / "tcga_ov_mutations.csv"
    pivot.to_csv(out_path, index=False)
    print(f"\n✅ Saved formatted mutations to {out_path}")
    print(f"  Unique patients with mutations: {len(pivot)}")
else:
    print("\n❌ No mutations found.")

