
import os
import sys

# --- DOCTRINE: DYNAMIC PATH CORRECTION ---
# This ensures that the script can find our tools library.
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from src.tools.ncbi_client import NCBIClient
from loguru import logger

def cache_canonical_sequences():
    """
    A one-time execution script to fetch and cache known-good canonical sequences.
    This establishes our ground-truth data assets.
    """
    logger.info("--- OPERATION: GROUND TRUTH ---")
    
    # Instantiate the client. The threat_matrix is not needed for this operation.
    client = NCBIClient(threat_matrix={})
    
    # Define our high-value targets.
    targets = {
        "RUNX1": "NM_001754.5"
    }
    
    output_dir = "data/canonical_sequences"
    os.makedirs(output_dir, exist_ok=True)
    
    for gene, accession_id in targets.items():
        logger.info(f"Acquiring ground truth CDS for {gene} ({accession_id})...")
        
        # Use the new, precise method to get only the coding sequence.
        sequence = client.get_cds_from_accession(accession_id)
        
        if sequence:
            filename = f"{gene}_{accession_id}_CDS.fasta"
            filepath = os.path.join(output_dir, filename)
            
            try:
                with open(filepath, "w") as f:
                    f.write(f">{accession_id} | Homo sapiens {gene} CDS (canonical)\n")
                    f.write(sequence)
                logger.success(f"Ground truth established. Saved CDS to {filepath}")
            except Exception as e:
                logger.error(f"Failed to write canonical CDS for {gene}: {e}")
        else:
            logger.error(f"Failed to acquire ground truth CDS for {gene}. Mission aborted for this target.")
            
    logger.info("--- OPERATION: GROUND TRUTH COMPLETE ---")

if __name__ == "__main__":
    cache_canonical_sequences() 
 