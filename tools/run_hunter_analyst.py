import sys
import os
import json
from loguru import logger

# --- DOCTRINE: DYNAMIC PATH CORRECTION ---
# Add the project root to the Python path to allow importing from 'src'
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from src.tools.ncbi_client import NCBIClient

def execute_hunt_phase(gene_symbol: str) -> dict | None:
    """
    Executes Phase 1 of the Predator Protocol: The Hunt.
    1. Fetches the full DNA sequence for the target gene.
    2. Hunts for vulnerable protein subdomains.
    3. Identifies the highest-priority target domain.
    4. Slices the full DNA sequence to extract the precise target motif.
    5. Returns the motif and its metadata.
    """
    logger.info(f"--- ⚔️ PREDATOR PROTOCOL: HUNT PHASE INITIATED: {gene_symbol} ⚔️ ---")
    
    try:
        hunter = NCBIClient(threat_matrix={})

        # Step 1: Get the full DNA sequence first. This is our raw material.
        logger.info(f"Acquiring full DNA sequence for {gene_symbol}...")
        full_dna_sequence = hunter.get_gene_sequence(gene_symbol)
        if not full_dna_sequence:
            logger.error(f"Could not retrieve DNA sequence for {gene_symbol}. Hunt aborted.")
            return None
        logger.success(f"✅ Full DNA sequence acquired (Length: {len(full_dna_sequence)} bp).")

        # Step 2: Hunt for protein domains.
        logger.info(f"Deploying Subdomain Hunter against target: {gene_symbol}")
        domains = hunter.find_protein_domains(gene_symbol)

        if not domains:
            logger.warning(f"Reconnaissance complete. No vulnerable subdomains found for {gene_symbol}.")
            return None

        # Step 3: Identify the highest-priority target.
        primary_target_domain = next((d for d in domains if d.get('priority') == 'High'), None)
        
        if not primary_target_domain:
            logger.warning("No 'High' priority target found. Unable to acquire target lock.")
            return None

        logger.success(f"🎯 Primary Target Identified: {primary_target_domain['name']} ({primary_target_domain['justification']})")
        logger.info(f"Domain Amino Acid Coordinates: {primary_target_domain['start']} - {primary_target_domain['end']}")

        # Step 4: Slice the DNA to get the precise motif.
        start_aa = primary_target_domain['start']
        end_aa = primary_target_domain['end']
        
        # Convert amino acid coordinates to nucleotide coordinates (1 AA = 3 bases)
        start_dna = (start_aa - 1) * 3
        end_dna = end_aa * 3
        
        target_motif_dna = full_dna_sequence[start_dna:end_dna]
        
        logger.success(f"✅ Target Lock Acquired. DNA Motif sliced from coordinates {start_dna}-{end_dna}.")
        logger.info(f"Motif Length: {len(target_motif_dna)} bp")

        # Step 5: Return the payload for the next phase.
        return {
            "gene_symbol": gene_symbol,
            "domain_name": primary_target_domain['name'],
            "target_motif_dna": target_motif_dna,
            "start_aa": start_aa,
            "end_aa": end_aa,
            "motif_len_bp": len(target_motif_dna)
        }

    except Exception as e:
        logger.error(f"❌❌❌ HUNT PHASE FAILED ❌❌❌")
        logger.error(f"A critical error occurred: {e}")
        return None


if __name__ == "__main__":
    # We can easily change the target here
    target_gene = "MMP9"
    hunt_result = execute_hunt_phase(target_gene)
    
    if hunt_result:
        print("\n" + "#"*60)
        print("## HUNT PHASE COMPLETE: TARGET PAYLOAD ##")
        print("#"*60 + "\n")
        # Pretty print the JSON output
        print(json.dumps(hunt_result, indent=2))
        print("\n" + "#"*60)

    print("\n\n")
    # Example of another target
    target_gene_2 = "VEGFA"
    hunt_result_2 = execute_hunt_phase(target_gene_2)

    if hunt_result_2:
        print("\n" + "#"*60)
        print("## HUNT PHASE COMPLETE: TARGET PAYLOAD ##")
        print("#"*60 + "\n")
        # Pretty print the JSON output
        print(json.dumps(hunt_result_2, indent=2))
        print("\n" + "#"*60) 