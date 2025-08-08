import modal
import os
import sys
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from loguru import logger

# --- DOCTRINE: DYNAMIC PATH CORRECTION ---
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from src.tools.ncbi_client import NCBIClient

# --- App & Image Definition ---
image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("fastapi", "uvicorn", "pydantic", "loguru", "biopython", "requests", "tenacity")
    .add_local_dir(local_path="src", remote_path="/root/src")
)

app = modal.App("hunter-analyst-service-v2-ghost", image=image) # New name for clarity

class HuntRequest(BaseModel):
    gene_symbol: str

@app.cls(timeout=300)
class HunterAnalyst:
    def __init__(self):
        # We can initialize the client here to be reused.
        self.ncbi_client = NCBIClient(threat_matrix={})

    @modal.web_endpoint(method="POST")
    def hunt(self, request: HuntRequest):
        """
        Executes the Hunt phase using the GHOST PROTOCOL.
        It bypasses the failed NCBI domain analysis and uses an internal
        heuristic to identify a target motif.
        """
        gene_symbol = request.gene_symbol
        logger.info(f"--- ⚔️ GHOST PROTOCOL HUNT INITIATED: {gene_symbol} ⚔️ ---")
        
        try:
            # Step 1: Get sequence using our resilient, local-first client
            full_dna_sequence = self.ncbi_client.get_dna_sequence_from_gene_symbol(gene_symbol)
            
            # --- PREDATOR PROTOCOL UPGRADE ---
            # Step 1.5: Acquire the target's protein sequence for downstream lethality assessment.
            target_protein_sequence = self.ncbi_client.get_protein_sequence_from_gene_symbol(gene_symbol)
            
            if not full_dna_sequence:
                raise HTTPException(status_code=404, detail=f"Could not retrieve DNA sequence for '{gene_symbol}' from local database.")

            if not target_protein_sequence:
                # This is now a mission-critical piece of intelligence.
                raise HTTPException(status_code=404, detail=f"Could not retrieve PROTEIN sequence for '{gene_symbol}'. Lethality assessment is impossible.")

            # Step 2: Heuristic Targeting (Bypass NCBI Domain Hunt)
            # We target a 300bp segment from the center of the gene.
            logger.info("Bypassing external domain analysis. Applying internal targeting heuristic.")
            
            if len(full_dna_sequence) < 300:
                target_motif_dna = full_dna_sequence
                start_aa = 1
                end_aa = len(full_dna_sequence) // 3
            else:
                midpoint = len(full_dna_sequence) // 2
                start_idx = midpoint - 150
                # Ensure the start index is a multiple of 3 for clean codon alignment
                start_idx -= start_idx % 3
                end_idx = start_idx + 300
                target_motif_dna = full_dna_sequence[start_idx:end_idx]
                start_aa = (start_idx // 3) + 1
                end_aa = end_idx // 3

            if not target_motif_dna:
                raise HTTPException(status_code=404, detail=f"Could not determine a target motif for '{gene_symbol}' using internal heuristic.")
            
            selected_domain_name = "Heuristic-Targeted Central Motif"
            
            logger.success(f"🎯 Ghost Protocol Target Lock on {gene_symbol}. Domain: {selected_domain_name}.")
            
            return {
                "gene_symbol": gene_symbol,
                "domain_name": selected_domain_name,
                "target_motif_dna": target_motif_dna,
                "target_protein_sequence": target_protein_sequence, # PREDATOR PROTOCOL UPGRADE
                "full_dna_sequence": full_dna_sequence, # DOCTRINE: Return full context for high-fidelity analysis
                "start_aa": start_aa,
                "end_aa": end_aa,
            }

        except HTTPException as http_exc:
            logger.error(f"Hunt for {gene_symbol} failed: {http_exc.detail} 💀")
            raise http_exc
        except Exception as e:
            logger.error(f"A critical error occurred during the hunt for {gene_symbol}: {e} 💥", exc_info=True)
            raise HTTPException(status_code=500, detail=f"Internal server error during hunt: {str(e)}")

# This is a local entrypoint for testing the class method directly.
@app.local_entrypoint()
def main(gene: str = "RUNX1"):
    """
    A local test harness to verify the HunterAnalyst can acquire full-context intel.
    This now functions as a true client, making an HTTP request to the deployed service.
    """
    import json
    import httpx

    HUNTER_URL = "https://crispro--hunter-analyst-service-v2-ghost-hunteranalyst-hunt.modal.run"

    print(f"--- ⚔️ LOCAL CLIENT TEST FOR HUNTER INTEL ⚔️ ---")
    print(f"  Target Service: {HUNTER_URL}")
    print(f"  Target Gene: {gene}")
    print("-------------------------------------------------")
    
    payload = {"gene_symbol": gene}

    try:
        with httpx.Client(timeout=300.0) as client:
            response = client.post(HUNTER_URL, json=payload)
            response.raise_for_status()
            result = response.json()

        print("\n--- ✅ HUNT REPORT ---")
        # We don't print the whole sequence, just confirm its length
        full_dna = result.get("full_dna_sequence", "")
        result_copy = result.copy()
        if len(full_dna) > 100:
            result_copy["full_dna_sequence"] = f"<DNA Sequence of length {len(full_dna)}>"
        
        print(json.dumps(result_copy, indent=2))
        print("\n---")
        print(f"🧬 Full DNA Sequence Length: {len(full_dna)}")
        print("---")

    except httpx.HTTPStatusError as e:
        print(f"\n--- ❌ TEST FAILED: HTTP Error ---")
        print(f"Status Code: {e.response.status_code}")
        print(f"Response: {e.response.text}")
        print("---------------------------------")
    except Exception as e:
        print(f"\n--- ❌ TEST FAILED: Unexpected Error ---")
        print(f"Error: {e}")
        print("---------------------------------------")
def main(gene: str = "BRAF"):
    hunter = HunterAnalyst()
    result = hunter.hunt.remote(HuntRequest(gene_symbol=gene))
    print(result)
def main(gene: str = "BRAF"):
    hunter = HunterAnalyst()
    result = hunter.hunt.remote(HuntRequest(gene_symbol=gene))
    print(result)