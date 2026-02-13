import pandas as pd
import requests
import io
import os

# URLs for UCSC Xena TCGA-OV data
# 1. Survival Data (Curated clinical data)
SURVIVAL_URL = "https://gdc.xenahubs.net/download/TCGA-OV.survival.tsv.gz"
# 2. Phenotype/Clinical Data (for Platinum Status)
CLINICAL_URL = "https://gdc.xenahubs.net/download/TCGA-OV.GDC_phenotype.tsv.gz"
# 3. Expression Data (HTSeq - FPKM-UQ) - 
# NOTE: This file is ~200-500MB compressed. We stream it to avoid memory issues if possible, 
# or download to temp.
EXPRESSION_URL = "https://gdc.xenahubs.net/download/TCGA-OV.htseq_fpkm-uq.tsv.gz"

DATA_DIR = "oncology-coPilot/oncology-backend-minimal/data/TCGA-OV/xena_mirror"
os.makedirs(DATA_DIR, exist_ok=True)

def download_file(url, target_path):
    if os.path.exists(target_path):
        print(f"✅ Exists: {target_path}")
        return
    
    print(f"⬇️ Downloading {url} to {target_path}...")
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        with open(target_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print("✅ Download complete.")
    except Exception as e:
        print(f"❌ Failed to download {url}: {e}")

# Execute Downloads
download_file(SURVIVAL_URL, f"{DATA_DIR}/survival.tsv.gz")
download_file(CLINICAL_URL, f"{DATA_DIR}/phenotype.tsv.gz")

print("\n⚠️ Large File Check: Expression Data")
# For expression, we might want to check if we can skip if validation is just PFI
print("Skipping full expression download for this step. We will verify clinical first.")
