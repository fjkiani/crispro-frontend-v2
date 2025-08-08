import sys
import os

# Ensure the project root is in the Python path to allow for correct module imports
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from tools.ncbi_client import NCBIClient
from loguru import logger

def main():
    """
    Local test script to validate the Subdomain Hunter (NCBIClient).
    """
    logger.info("--- 🎯 SUBDOMAIN HUNTER LOCAL TEST 🎯 ---")
    
    # Target a known enemy gene for this test
    gene_symbol = "MMP9"
    logger.info(f"Deploying Subdomain Hunter against target: {gene_symbol}")

    try:
        # We need to instantiate the client to test its methods
        # The threat_matrix is a legacy parameter for this test's purpose.
        hunter = NCBIClient(threat_matrix={})
        
        # 1. Test a gene we know is in our local DB
        logger.info("--- Testing successful case (BRAF) ---")
        domains = hunter.find_protein_domains(gene_symbol)

        if not domains:
            logger.warning(f"Reconnaissance complete. No vulnerable subdomains found for {gene_symbol}.")
            return

        logger.success(f"✅ Reconnaissance successful! Found {len(domains)} vulnerable subdomains for {gene_symbol}:")
        for i, domain in enumerate(domains):
            print(f"  [{i+1}] {domain['name']} (Coordinates: {domain['start']} - {domain['end']}) - Priority: {domain['priority']} - Justification: {domain['justification']}")

    except Exception as e:
        logger.error(f"❌❌❌ LOCAL TEST FAILED ❌❌❌")
        logger.error(f"A critical error occurred: {e}")

if __name__ == "__main__":