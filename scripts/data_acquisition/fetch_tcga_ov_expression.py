import requests
import os
import time
import tarfile
import shutil

# Configuration
# cBioPortal PanCancer Atlas 2018 (Stable, widely used)
CBIO_URL = "https://cbioportal-datahub.s3.amazonaws.com/ov_tcga_pan_can_atlas_2018.tar.gz"
DATA_DIR = "oncology-coPilot/oncology-backend-minimal/data/TCGA-OV/cbioportal"
TARGET_FILE = os.path.join(DATA_DIR, "ov_tcga.tar.gz")
MRI_FILE = os.path.join(DATA_DIR, "data_mrna_seq_v2_rsem.txt")

def download_and_extract_tcga():
    os.makedirs(DATA_DIR, exist_ok=True)
    
    # 1. Download
    if not os.path.exists(TARGET_FILE):
        print(f"⬇️ Downloading TCGA-OV from {CBIO_URL}...")
        try:
            with requests.get(CBIO_URL, stream=True) as r:
                r.raise_for_status()
                with open(TARGET_FILE, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
            print("✅ Download complete.")
        except Exception as e:
            print(f"❌ Download Failed: {e}")
            return

    # 2. Extract
    print("📦 Extracting Expression Data...")
    try:
        with tarfile.open(TARGET_FILE, "r:gz") as tar:
            # Look for the expression file
            # content is usually inside a folder like 'ov_tcga_pan_can_atlas_2018/'
            for member in tar.getmembers():
                if "data_mrna_seq_v2_rsem.txt" in member.name:
                    member.name = os.path.basename(member.name) # flattened
                    tar.extract(member, path=DATA_DIR)
                    print(f"✅ Extracted: {member.name}")
                    
        # Verification
        if os.path.exists(MRI_FILE):
             print(f"🏆 SUCCESS: Expression matrix ready at {MRI_FILE}")
        else:
             print("❌ Extraction failed: Target file not found in archive.")

    except Exception as e:
        print(f"❌ Extraction Error: {e}")

if __name__ == "__main__":
    download_and_extract_tcga()
