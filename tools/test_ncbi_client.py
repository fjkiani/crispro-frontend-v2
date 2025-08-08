import os
import sys
from loguru import logger

# --- DOCTRINE: DYNAMIC PATH CORRECTION ---
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from src.tools.ncbi_client import NCBIClient

def run_isolated_dna_hunt(gene_symbol: str):
    """
    Performs an isolated test of the NCBIClient's DNA sequence retrieval.
    This allows for rapid debugging without service deployment.
    """
    logger.remove() # Remove default logger to have clean output
    logger.add(sys.stderr, level="DEBUG") # Add a new one with DEBUG level

    print(f"--- 🔬 ISOLATED TEST: DNA HUNT FOR '{gene_symbol}' 🔬 ---")
    
    # Initialize the client with an empty threat matrix for this test
    client = NCBIClient(threat_matrix={})
    
    try:
        # Directly call the method we need to debug
        sequence = client.get_dna_sequence_from_gene_symbol(gene_symbol)
        
        if sequence:
            print("\n--- ✅ TEST SUCCEEDED ---")
            print(f"Successfully retrieved DNA sequence for '{gene_symbol}'.")
            print(f"Sequence Length: {len(sequence)}")
            print(f"Sequence (first 50bp): {sequence[:50]}...")
        else:
            print("\n--- ❌ TEST FAILED ---")
            print(f"Failed to retrieve DNA sequence for '{gene_symbol}'.")
            print("Check the DEBUG logs above for detailed error messages from the client.")

    except Exception as e:
        print(f"\n--- 💥 CRITICAL FAILURE 💥 ---")
        print(f"An unexpected error occurred during the test: {e}")
        logger.exception("Test execution failed.")

if __name__ == "__main__":
    target_gene = "BRAF"
    if len(sys.argv) > 1:
        target_gene = sys.argv[1]
    run_isolated_dna_hunt(target_gene) 
import sys
from loguru import logger

# --- DOCTRINE: DYNAMIC PATH CORRECTION ---
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from src.tools.ncbi_client import NCBIClient

def run_isolated_dna_hunt(gene_symbol: str):
    """
    Performs an isolated test of the NCBIClient's DNA sequence retrieval.
    This allows for rapid debugging without service deployment.
    """
    logger.remove() # Remove default logger to have clean output
    logger.add(sys.stderr, level="DEBUG") # Add a new one with DEBUG level

    print(f"--- 🔬 ISOLATED TEST: DNA HUNT FOR '{gene_symbol}' 🔬 ---")
    
    # Initialize the client with an empty threat matrix for this test
    client = NCBIClient(threat_matrix={})
    
    try:
        # Directly call the method we need to debug
        sequence = client.get_dna_sequence_from_gene_symbol(gene_symbol)
        
        if sequence:
            print("\n--- ✅ TEST SUCCEEDED ---")
            print(f"Successfully retrieved DNA sequence for '{gene_symbol}'.")
            print(f"Sequence Length: {len(sequence)}")
            print(f"Sequence (first 50bp): {sequence[:50]}...")
        else:
            print("\n--- ❌ TEST FAILED ---")
            print(f"Failed to retrieve DNA sequence for '{gene_symbol}'.")
            print("Check the DEBUG logs above for detailed error messages from the client.")

    except Exception as e:
        print(f"\n--- 💥 CRITICAL FAILURE 💥 ---")
        print(f"An unexpected error occurred during the test: {e}")
        logger.exception("Test execution failed.")

if __name__ == "__main__":
    target_gene = "BRAF"
    if len(sys.argv) > 1:
        target_gene = sys.argv[1]
    run_isolated_dna_hunt(target_gene) 
 
 
 
 
 