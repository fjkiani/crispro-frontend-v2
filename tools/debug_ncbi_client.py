import argparse
import sys
import os

# This is a hack to allow the script to find the src module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.tools.ncbi_client import NCBIClient
from loguru import logger

def main():
    """
    A local debugging script to test the NCBIClient's ability to retrieve
    a DNA sequence for a given gene symbol.
    
    This helps isolate failures in the hunter-analyst-service from the
    rest of the microservice architecture.
    
    Usage:
        venv/bin/python tools/debug_ncbi_client.py --gene VEGFA
    """
    parser = argparse.ArgumentParser(description="Debug NCBI Client Gene Sequence Fetching.")
    parser.add_argument("--gene", type=str, required=True, help="The gene symbol to test (e.g., VEGFA, MMP9, BRAF).")
    args = parser.parse_args()

    gene_symbol = args.gene
    logger.info(f"🔬 Attempting to fetch DNA sequence for '{gene_symbol}' using NCBIClient...")

    try:
        # The threat_matrix is not used for this function, so we pass an empty dict.
        client = NCBIClient(threat_matrix={})
        dna_sequence = client.get_dna_sequence_from_gene_symbol(gene_symbol)

        if dna_sequence:
            logger.success(f"✅ Successfully retrieved DNA sequence for {gene_symbol}.")
            print("\n--- SEQUENCE ---")
            print(dna_sequence[:200] + "..." if len(dna_sequence) > 200 else dna_sequence)
            print("----------------")
        else:
            logger.error(f"❌ Failed to retrieve DNA sequence for {gene_symbol}. The client returned None or an empty string.")

    except Exception as e:
        logger.error(f"💥 An unexpected error occurred during the fetch for {gene_symbol}:", exc_info=True)

if __name__ == "__main__":
    main() 
 
import sys
import os

# This is a hack to allow the script to find the src module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.tools.ncbi_client import NCBIClient
from loguru import logger

def main():
    """
    A local debugging script to test the NCBIClient's ability to retrieve
    a DNA sequence for a given gene symbol.
    
    This helps isolate failures in the hunter-analyst-service from the
    rest of the microservice architecture.
    
    Usage:
        venv/bin/python tools/debug_ncbi_client.py --gene VEGFA
    """
    parser = argparse.ArgumentParser(description="Debug NCBI Client Gene Sequence Fetching.")
    parser.add_argument("--gene", type=str, required=True, help="The gene symbol to test (e.g., VEGFA, MMP9, BRAF).")
    args = parser.parse_args()

    gene_symbol = args.gene
    logger.info(f"🔬 Attempting to fetch DNA sequence for '{gene_symbol}' using NCBIClient...")

    try:
        # The threat_matrix is not used for this function, so we pass an empty dict.
        client = NCBIClient(threat_matrix={})
        dna_sequence = client.get_dna_sequence_from_gene_symbol(gene_symbol)

        if dna_sequence:
            logger.success(f"✅ Successfully retrieved DNA sequence for {gene_symbol}.")
            print("\n--- SEQUENCE ---")
            print(dna_sequence[:200] + "..." if len(dna_sequence) > 200 else dna_sequence)
            print("----------------")
        else:
            logger.error(f"❌ Failed to retrieve DNA sequence for {gene_symbol}. The client returned None or an empty string.")

    except Exception as e:
        logger.error(f"💥 An unexpected error occurred during the fetch for {gene_symbol}:", exc_info=True)

if __name__ == "__main__":
    main() 
 
 
 
 
 
 
 
 
 