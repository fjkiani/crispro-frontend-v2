import modal
import os
import logging

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# --- App & Volume Definition ---
# Define the volume to store the large data file
volume = modal.Volume.from_name("alphamissense-data", create_if_missing=True)

# Define a single image for our data preparation pipeline. It needs requests, tqdm, pandas, and pyarrow.
dataprep_image = modal.Image.debian_slim(python_version="3.11").pip_install(
    "requests", "tqdm", "pandas", "pyarrow"
)

app = modal.App("alphamissense-dataprep", image=dataprep_image)

# --- Paths and URLs ---
DATA_URL = "https://zenodo.org/records/8360242/files/AlphaMissense_hg38.tsv.gz"
TSV_PATH = "/data/AlphaMissense_hg38.tsv.gz"
PARQUET_PATH = "/data/AlphaMissense_hg38.parquet"


@app.function(
    volumes={"/data": volume}, 
    timeout=7200, 
    memory=128 * 1024  # Request 128GiB of memory
)
def prepare_data():
    """
    A unified function to download and convert the AlphaMissense data.
    """
    # Defer heavy imports to container/runtime to avoid local import errors
    import requests  # type: ignore
    from tqdm import tqdm  # type: ignore
    import pandas as pd  # type: ignore
    # --- Step 1: Download the TSV if it doesn't exist ---
    if not os.path.exists(TSV_PATH):
        logger.info(f"⬇️ TSV file not found. Downloading from {DATA_URL}...")
        try:
            with requests.get(DATA_URL, stream=True) as r:
                r.raise_for_status()
                total_size = int(r.headers.get('content-length', 0))
                
                with open(TSV_PATH, 'wb') as f, tqdm(
                    desc="AlphaMissense_hg38.tsv.gz",
                    total=total_size,
                    unit='iB',
                    unit_scale=True,
                    unit_divisor=1024,
                ) as bar:
                    for chunk in r.iter_content(chunk_size=8192):
                        size = f.write(chunk)
                        bar.update(size)

            logger.info(f"✅ Download complete. Data saved to volume at {TSV_PATH}.")
            volume.commit()
            
        except Exception as e:
            logger.error(f"❌ An error occurred during download: {e}", exc_info=True)
            if os.path.exists(TSV_PATH):
                os.remove(TSV_PATH)
            raise
    else:
        logger.info(f"✅ TSV file already exists at {TSV_PATH}. Skipping download.")

    # --- Step 2: Convert to Parquet if it doesn't exist ---
    if not os.path.exists(PARQUET_PATH):
        logger.info(f"Parquet file not found. Starting conversion of {TSV_PATH}...")
        try:
            df = pd.read_csv(
                TSV_PATH,
                sep='\\t',
                header=0,
                engine='python',
                skiprows=3,
                usecols=['#CHROM', 'POS', 'REF', 'ALT', 'am_pathogenicity', 'am_class']
            )
            logger.info("TSV file loaded successfully for conversion.")

            df.to_parquet(PARQUET_PATH, engine='pyarrow', index=False)
            
            logger.info(f"✅ Successfully converted data and saved to {PARQUET_PATH}.")
            volume.commit()
            
        except Exception as e:
            logger.error(f"❌ An error occurred during conversion: {e}", exc_info=True)
            if os.path.exists(PARQUET_PATH):
                os.remove(PARQUET_PATH)
            raise
    else:
        logger.info(f"✅ Parquet file already exists at {PARQUET_PATH}. Skipping conversion.")
        
    logger.info("Data preparation complete.")
    return "✅ Data preparation successfully completed."


@app.local_entrypoint()
def main():
    print("🚀 Starting AlphaMissense data preparation pipeline (download & convert)...")
    status = prepare_data.remote()
    print(f"✔️ Remote process finished with status: {status}")
    print("✅ The data should now be in the 'alphamissense-data' volume in Parquet format.") 
 