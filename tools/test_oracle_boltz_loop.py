import requests
import sys
import os
import json
import uuid
from loguru import logger
from Bio.Seq import Seq

# --- DOCTRINE: DYNAMIC PATH CORRECTION ---
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from src.tools.ncbi_client import NCBIClient

# --- Configuration ---
UNIFIED_ORACLE_URL = "https://crispro--unified-oracle-service-unifiedoracle-invoke.modal.run"
BOLTZ_URL = "https://crispro--boltz-service-fastapi-app.modal.run/v1/predict_structure"
TARGET_GENE_SYMBOL = "BRAF"

def translate_dna_to_protein(dna_sequence: str) -> str | None:
    """Translates DNA to protein, returns None on failure."""
    try:
        if len(dna_sequence) % 3 != 0:
            remainder = len(dna_sequence) % 3
            dna_sequence = dna_sequence[:-remainder]
        
        protein_sequence = str(Seq(dna_sequence).translate(to_stop=True))
        if not protein_sequence:
            logger.warning("Translation resulted in an empty sequence (immediate stop codon).")
            return None
        return protein_sequence
    except Exception as e:
        logger.error(f"DNA to protein translation failed: {e}")
        return None

def run_oracle_boltz_test():
    """
    Performs an isolated, end-to-end test of the Oracle-to-Boltz loop.
    """
    logger.remove()
    logger.add(sys.stderr, level="INFO")

    logger.info("--- ⚔️  OPERATION FORGED STEEL: INITIATED ⚔️ ---")

    try:
        # --- 1. THE HUNT ---
        logger.info(f"🏹 HUNT: Acquiring long-range context from '{TARGET_GENE_SYMBOL}'...")
        ncbi_client = NCBIClient(threat_matrix={})
        context_dna = ncbi_client.get_dna_sequence_from_gene_symbol(TARGET_GENE_SYMBOL)
        if not context_dna:
            raise RuntimeError("Failed to acquire DNA for context.")
        logger.success(f"Context acquired. Length: {len(context_dna)}bp.")

        # --- 2. THE FORGE (Long-Range) ---
        prompt = context_dna[-1024:]
        logger.info(f"🔥 FORGE: Striking Oracle with {len(prompt)}bp prompt...")
        oracle_payload = {
            "action": "generate",
            "params": {"prompt": prompt, "gen_params": {"n_tokens": 500, "temperature": 0.5}}
        }
        response = requests.post(UNIFIED_ORACLE_URL, json=oracle_payload, timeout=300)
        response.raise_for_status()
        oracle_result = response.json()
        generated_dna = oracle_result.get("completion")
        if not generated_dna:
            raise RuntimeError(f"Oracle failed to generate a sequence. Response: {oracle_result}")
        logger.success(f"Oracle forged a {len(generated_dna)}bp sequence.")

        # --- 3. TRANSLATION ---
        candidate_dna = generated_dna[175:325] # Extract 150bp from the middle
        logger.info(f"🧬 TRANSLATE: Converting 150bp candidate DNA to protein...")
        candidate_protein = translate_dna_to_protein(candidate_dna)
        if not candidate_protein:
            raise RuntimeError("Translation failed.")
        logger.success(f"Translation complete. Protein length: {len(candidate_protein)}aa.")

        # --- 4. VALIDATE (Boltz Strike) ---
        job_id = f"forged_steel_test_{uuid.uuid4().hex[:8]}"
        logger.info(f"🛡️  VALIDATE: Striking Boltz with protein payload (Job ID: {job_id})...")
        boltz_payload = {"protein_sequence": candidate_protein, "job_id": job_id}
        response = requests.post(BOLTZ_URL, json=boltz_payload, timeout=1800)
        response.raise_for_status()
        boltz_result = response.json()
        
        logger.info("--- ✅ BOLTZ RESPONDED ---")
        print(json.dumps(boltz_result, indent=2))
        
        # --- 5. VERDICT ---
        status = boltz_result.get("status")
        plddt = boltz_result.get("plddt_score")
        
        if status == "complete" and plddt is not None and plddt > 70:
            logger.success(f"--- ✅ OPERATION SUCCESS: Weapon is structurally sound (pLDDT: {plddt:.2f}) ---")
        else:
            logger.error(f"--- ❌ OPERATION FAILED: Weapon is a 'wet noodle' (pLDDT: {plddt:.2f}) ---")

    except Exception as e:
        logger.error(f"--- 💥 MISSION CRITICAL FAILURE 💥 ---")
        logger.error(f"An unexpected error occurred: {e}", exc_info=True)

if __name__ == "__main__":
    run_oracle_boltz_test() 
 