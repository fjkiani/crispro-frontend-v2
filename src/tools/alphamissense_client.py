import pandas as pd
import logging
import os

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# The canonical path for the pre-computed AlphaMissense data file.
# We now point to the much faster Parquet file.
ALPHAMISSENSE_DATA_PATH = os.getenv("ALPHAMISSENSE_DATA", "/data/AlphaMissense_hg38.parquet")

class AlphaMissenseClient:
    def __init__(self, data_path: str = ALPHAMISSENSE_DATA_PATH):
        """
        Initializes the client by loading the AlphaMissense data from a Parquet file.
        """
        self.data_path = data_path
        self._data = None
        self._load_data()

    def _load_data(self):
        """Loads and indexes the AlphaMissense data from a Parquet file."""
        try:
            logger.info(f"Loading AlphaMissense data from Parquet file: {self.data_path}...")
            # Loading from Parquet is significantly faster than from a gzipped TSV.
            self._data = pd.read_parquet(
                self.data_path,
                engine='pyarrow',
                columns=['#CHROM', 'POS', 'REF', 'ALT', 'am_pathogenicity', 'am_class']
            )
            # Create a more convenient index for fast lookups
            self._data.set_index(['#CHROM', 'POS', 'REF', 'ALT'], inplace=True)
            self._data.sort_index(inplace=True)
            logger.info("AlphaMissense data loaded and indexed successfully from Parquet.")
            logger.info(f"DataFrame shape: {self._data.shape}")

        except FileNotFoundError:
            logger.error(f"AlphaMissense Parquet file not found at {self.data_path}. The client will operate in a degraded mode, returning error scores.")
            self._data = None # Ensure data is None if file not found
        except Exception as e:
            logger.error(f"An error occurred while loading AlphaMissense Parquet data: {e}")
            self._data = None

    def get_score(self, genomic_variant_str: str) -> dict:
        """
        Queries the local AlphaMissense data for a pathogenicity score.

        Args:
            genomic_variant_str: The variant in "chr:pos:ref:alt" format (e.g., "chr7:140753336:A:T").

        Returns:
            A dictionary containing the score and classification.
        """
        logger.info(f"AlphaMissenseClient received raw input: '{genomic_variant_str}'")
        if self._data is None:
            return {"score": -999.0, "classification": "error_data_not_loaded", "error": "AlphaMissense data not available."}
        
        try:
            chrom, pos, ref, alt = genomic_variant_str.split(':')
            pos = int(pos)
            
            # Log the exact key being used for the lookup
            lookup_key = (chrom, pos, ref, alt)
            logger.info(f"Attempting to query AlphaMissense index with key: {lookup_key}")
            
            # Query the multi-index DataFrame
            result = self._data.loc[lookup_key]
            
            # The result could be a Series (single match) or DataFrame (multiple matches)
            if isinstance(result, pd.Series):
                score = result['am_pathogenicity']
                am_class = result['am_class']
            elif isinstance(result, pd.DataFrame) and not result.empty:
                # If multiple matches, take the first one
                score = result['am_pathogenicity'].iloc[0]
                am_class = result['am_class'].iloc[0]
            else:
                score = -998.0 # Not found
                am_class = "Not Found"

            # CRITICAL FIX: The client must return raw values, not pandas objects.
            # .item() extracts the single value from a Series/scalar.
            return {"score": float(score), "classification": str(am_class)}
            
        except KeyError:
            logger.warning(f"Variant {genomic_variant_str} not found in AlphaMissense data.")
            return {"score": -998.0, "classification": "not_found", "error": f"Variant {genomic_variant_str} not found."}
        except Exception as e:
            logger.error(f"An error occurred during lookup for {genomic_variant_str}: {e}")
            return {"score": -997.0, "classification": "error_lookup_failed", "error": str(e)}

if __name__ == '__main__':
    # This example assumes the data file is available locally in Parquet format.
    
    client = AlphaMissenseClient(data_path="data/alphamissense/AlphaMissense_hg38.parquet") # Point to a local path for testing
    
    if client._data is not None:
        # Example for BRAF V600E - NOTE: The coordinates in the original file might be GRCh37.
        # This is a known pathogenic variant.
        score_data = client.get_score("chr7:140453136:A:T") 
        print(f"BRAF V600E (pathogenic): {score_data}")
        
        # Example for a known benign variant
        score_data_benign = client.get_score("chr1:909998:C:T")
        print(f"Benign variant: {score_data_benign}")
    else:
        print("Could not run examples because AlphaMissense data file was not found.")
        print("Please run 'tools/download_alphamissense.py' and 'tools/convert_alphamissense_to_parquet.py' first.")
    
    # Clean up the direct debugging section as it's no longer needed. 