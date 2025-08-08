import sys
import os
from loguru import logger
import argparse

# --- Add src to path to allow for direct NCBIClient import ---
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.tools.ncbi_client import NCBIClient

def populate_armory(gene_symbol: str):
    """
    A dedicated tool to populate our local armory with both DNA and protein
    intelligence for a given gene target. This is a deliberate, manual
    operation to ensure our downstream tests are deterministic and offline-capable.
    """
    logger.info(f"--- OPERATION: ARMORY FIRST ---")
    logger.info(f"Populating local cache for high-value target: {gene_symbol}")
    
    client = NCBIClient(threat_matrix={})
    
    logger.info(f"Acquiring DNA sequence for {gene_symbol}...")
    dna_seq = client.get_dna_sequence_from_gene_symbol(gene_symbol)
    if dna_seq:
        logger.success(f"DNA for {gene_symbol} secured and cached.")
    else:
        logger.error(f"Failed to acquire DNA for {gene_symbol}. Mission failure.")
        sys.exit(1)

    logger.info(f"Acquiring PROTEIN sequence for {gene_symbol}...")
    protein_seq = client.get_protein_sequence_from_gene_symbol(gene_symbol)
    if protein_seq:
        logger.success(f"PROTEIN for {gene_symbol} secured and cached.")
    else:
        logger.error(f"Failed to acquire PROTEIN for {gene_symbol}. Mission failure.")
        sys.exit(1)
        
    logger.info(f"--- ARMORY POPULATED FOR {gene_symbol}. WE ARE MISSION READY. ---")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Populate the local gene armory.")
    parser.add_argument("gene_symbol", type=str, help="The gene symbol to cache (e.g., VEGFA).")
    args = parser.parse_args()
    
    populate_armory(args.gene_symbol) 
import os
from loguru import logger
import argparse

# --- Add src to path to allow for direct NCBIClient import ---
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.tools.ncbi_client import NCBIClient

def populate_armory(gene_symbol: str):
    """
    A dedicated tool to populate our local armory with both DNA and protein
    intelligence for a given gene target. This is a deliberate, manual
    operation to ensure our downstream tests are deterministic and offline-capable.
    """
    logger.info(f"--- OPERATION: ARMORY FIRST ---")
    logger.info(f"Populating local cache for high-value target: {gene_symbol}")
    
    client = NCBIClient(threat_matrix={})
    
    logger.info(f"Acquiring DNA sequence for {gene_symbol}...")
    dna_seq = client.get_dna_sequence_from_gene_symbol(gene_symbol)
    if dna_seq:
        logger.success(f"DNA for {gene_symbol} secured and cached.")
    else:
        logger.error(f"Failed to acquire DNA for {gene_symbol}. Mission failure.")
        sys.exit(1)

    logger.info(f"Acquiring PROTEIN sequence for {gene_symbol}...")
    protein_seq = client.get_protein_sequence_from_gene_symbol(gene_symbol)
    if protein_seq:
        logger.success(f"PROTEIN for {gene_symbol} secured and cached.")
    else:
        logger.error(f"Failed to acquire PROTEIN for {gene_symbol}. Mission failure.")
        sys.exit(1)
        
    logger.info(f"--- ARMORY POPULATED FOR {gene_symbol}. WE ARE MISSION READY. ---")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Populate the local gene armory.")
    parser.add_argument("gene_symbol", type=str, help="The gene symbol to cache (e.g., VEGFA).")
    args = parser.parse_args()
    
    populate_armory(args.gene_symbol) 
 
 
 
 
 