import sys
import os
from loguru import logger

# --- DOCTRINAL CORRECTION: Correct PYTHONPATH ---
# Add the project root to the path to allow imports from `src`.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.tools.ncbi_client import NCBIClient

def run_brca1_recon():
    """
    Uses the NCBIClient to perform reconnaissance on the BRCA1 gene.
    Fetches the full genomic DNA sequence and saves it for analysis.
    """
    logger.info("Initiating Operation: BRCA1 Gauntlet - Reconnaissance Phase")
    
    # We don't need a threat matrix for this targeted mission.
    dummy_threat_matrix = {}
    client = NCBIClient(threat_matrix=dummy_threat_matrix)
    
    gene_symbol = "BRCA1"
    
    logger.info(f"Acquiring genomic sequence for high-value target: {gene_symbol}")
    sequence = client.get_dna_sequence_from_gene_symbol(gene_symbol)
    
    if sequence:
        logger.success(f"Successfully retrieved {len(sequence)}bp for {gene_symbol}.")
        # The client automatically caches it, so we just need to confirm.
        # For clarity, let's print a snippet.
        print("\n--- BEGIN BRCA1 GENOMIC SEQUENCE (first 500bp) ---")
        print(sequence[:500])
        print("--- END BRCA1 GENOMIC SEQUENCE ---\n")
        logger.info("Sequence has been cached locally in data/gene_database/BRCA1.fasta")
    else:
        logger.error(f"Mission failed. Could not retrieve sequence for {gene_symbol}.")

if __name__ == "__main__":
    run_brca1_recon() 
import os
from loguru import logger

# --- DOCTRINAL CORRECTION: Correct PYTHONPATH ---
# Add the project root to the path to allow imports from `src`.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.tools.ncbi_client import NCBIClient

def run_brca1_recon():
    """
    Uses the NCBIClient to perform reconnaissance on the BRCA1 gene.
    Fetches the full genomic DNA sequence and saves it for analysis.
    """
    logger.info("Initiating Operation: BRCA1 Gauntlet - Reconnaissance Phase")
    
    # We don't need a threat matrix for this targeted mission.
    dummy_threat_matrix = {}
    client = NCBIClient(threat_matrix=dummy_threat_matrix)
    
    gene_symbol = "BRCA1"
    
    logger.info(f"Acquiring genomic sequence for high-value target: {gene_symbol}")
    sequence = client.get_dna_sequence_from_gene_symbol(gene_symbol)
    
    if sequence:
        logger.success(f"Successfully retrieved {len(sequence)}bp for {gene_symbol}.")
        # The client automatically caches it, so we just need to confirm.
        # For clarity, let's print a snippet.
        print("\n--- BEGIN BRCA1 GENOMIC SEQUENCE (first 500bp) ---")
        print(sequence[:500])
        print("--- END BRCA1 GENOMIC SEQUENCE ---\n")
        logger.info("Sequence has been cached locally in data/gene_database/BRCA1.fasta")
    else:
        logger.error(f"Mission failed. Could not retrieve sequence for {gene_symbol}.")

if __name__ == "__main__":
    run_brca1_recon() 
 
 
 
 
 