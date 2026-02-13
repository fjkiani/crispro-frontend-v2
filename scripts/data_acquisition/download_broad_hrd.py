import requests
import zipfile
import io
import os
import pandas as pd
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# GDC PanCan-DDR-2018 Data Resources
# Source: https://gdc.cancer.gov/about-data/publications/PanCan-DDR-2018
# File: TCGA_DDR_Data_Resources.zip
URL_GDC_ZIP = "https://api.gdc.cancer.gov/data/d9fdd899-bffe-4325-b0ce-58cbf9924c0a"

OUTPUT_DIR = "data/external"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "TCGA_DDR_Data_Resources.zip")
EXTRACT_DIR = os.path.join(OUTPUT_DIR, "TCGA_DDR_Data_Resources")

def download_and_extract():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        
    logger.info(f"Initiating extraction from GDC Node: {URL_GDC_ZIP}")
    
    try:
        response = requests.get(URL_GDC_ZIP, stream=True)
        response.raise_for_status()
        
        logger.info("Payload received. Decrypting/Unzipping...")
        with zipfile.ZipFile(io.BytesIO(response.content)) as z:
            z.extractall(EXTRACT_DIR)
            
        logger.info(f"Extraction complete to {EXTRACT_DIR}")
        
        # Locate the HRD Scores file
        # Usually named something like "TCGA.HRD_with_Sampleinfo.txt" or similar
        # We will list files to be sure
        files = os.listdir(EXTRACT_DIR)
        logger.info(f"Manifest: {files}")
        
        return True
        
    except Exception as e:
        logger.error(f"Extraction failed: {e}")
        return False

if __name__ == "__main__":
    print(">>> ZETA PROTOCOL: DATA INGESTION <<<")
    success = download_and_extract()
    if success:
        print("[SUCCESS] Asset Secured.")
    else:
        print("[FAILURE] Connection Terminated.")
