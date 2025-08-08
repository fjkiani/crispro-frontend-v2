import os
import sys
from loguru import logger
import json
import requests
from Bio.Seq import Seq
import time
from collections import Counter

# --- DOCTRINE: DYNAMIC PATH CORRECTION ---
# Ensure the project root is in the Python path to allow for correct module imports
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from tools.run_hunter_analyst import execute_hunt_phase
from tools.translation_bridge import translate_dna_to_protein

# --- CONFIGURATION ---
# CORRECTED: Pointing to the live URL from the latest deployment.
ZETA_FORGE_URL = "https://crispro--zeta-forge-v1-api.modal.run"
ZETA_ORACLE_URL = "https://crispro--zeta-oracle-v2-zetaoracle-api.modal.run"
BOLTZ_SERVICE_URL = "https://crispro--boltz-service-fastapi-app.modal.run"

def filter_low_complexity(sequences: list[str], threshold: float = 0.5) -> list[str]:
    """
    Filters out low-complexity sequences like 'TTATTATTATTATTATTATT'.
    It calculates the frequency of the most common dinucleotide and if it
    exceeds the threshold, the sequence is discarded as junk.
    """
    logger.info(f"--- 🛡️ JUNK FILTER ENGAGED: Screening {len(sequences)} raw candidates... ---")
    high_quality_candidates = []
    for seq in sequences:
        if len(seq) < 10: # Discard very short, junk sequences
            continue
        
        # Calculate dinucleotide frequencies
        dinc_counts = Counter([seq[i:i+2] for i in range(len(seq) - 1)])
        most_common_dinc = dinc_counts.most_common(1)[0]
        
        # Calculate the frequency of the most common dinucleotide
        freq = most_common_dinc[1] / len(seq)
        
        if freq > threshold:
            logger.warning(f"Junk sequence DETECTED and DISCARDED (Dinucleotide '{most_common_dinc[0]}' frequency {freq:.2f} > {threshold}): {seq}")
        else:
            high_quality_candidates.append(seq)
            
    logger.success(f"--- ✅ JUNK FILTER COMPLETE: {len(high_quality_candidates)} candidates passed screening. ---")
    return high_quality_candidates


def run_predator_protocol(target_gene_symbol: str):
    """
    Orchestrates the full 4-phase Predator Protocol.
    """
    logger.info(f"--- 🐆 PREDATOR PROTOCOL INITIATED: ENGAGING TARGET {target_gene_symbol} 🐆 ---")

    # --- PHASE 1: THE HUNT ---
    # Use our validated hunter tool to acquire the target motif
    target_payload = execute_hunt_phase(target_gene_symbol)

    if not target_payload:
        logger.error(f"PROTOCOL FAILED: Hunt phase did not acquire a target for {target_gene_symbol}. Aborting mission.")
        return

    logger.success("--- ✅ HUNT PHASE COMPLETE ---")
    logger.info("Target Payload Acquired.")
    
    # --- PHASE 2: THE AMBUSH ---
    logger.info("--- 🔥 AMBUSH PHASE INITIATED: FORGING CANDIDATE INHIBITORS ---")
    
    target_motif = target_payload.get("target_motif_dna")
    if not target_motif:
        logger.error("PROTOCOL FAILED: Target payload is missing 'target_motif_dna'.")
        return

    # Use the reverse complement of the motif as bait for the generative model
    bait_sequence = str(Seq(target_motif).reverse_complement())
    logger.info(f"Bait sequence generated (reverse complement of target motif). Length: {len(bait_sequence)} bp.")

    # This payload now matches the new, simplified ForgeRequest model
    forge_payload = {
        "bait_sequence": bait_sequence,
        "num_candidates_per_temp": 15, # Saturation Bombardment
        "temperatures": [0.2, 0.7, 1.2],
        "generation_length": len(bait_sequence) + 20 # Generate slightly longer candidates
    }

    try:
        logger.info(f"Dispatching Ambush command to ZetaForge at {ZETA_FORGE_URL}...")
        response = requests.post(f"{ZETA_FORGE_URL}/generate_inhibitor", json=forge_payload, timeout=30)
        response.raise_for_status()
        job_id = response.json().get("job_id")
        if not job_id:
            raise RuntimeError("Failed to get job_id from ZetaForge submission.")
        
        logger.success(f"✅ Forge job submitted. Job ID: {job_id}. Polling for results...")

        # Poll for results
        while True:
            time.sleep(20) # Poll every 20 seconds
            status_response = requests.get(f"{ZETA_FORGE_URL}/status/{job_id}", timeout=30)
            status_response.raise_for_status()
            status_data = status_response.json()

            logger.info(f"Forge status: {status_data.get('status')} - {status_data.get('message', 'Polling...')}")

            if status_data.get("status") == "complete":
                raw_candidates = status_data.get('candidates', [])
                if not raw_candidates:
                     raise RuntimeError("Forge completed but returned no candidates.")
                
                logger.success(f"--- ✅ AMBUSH PHASE COMPLETE: {len(raw_candidates)} RAW CANDIDATES FORGED ---")
                
                # --- NEW: DOCTRINE OF VETTING ---
                # Immediately filter out the low-complexity junk.
                vetted_candidates = filter_low_complexity(raw_candidates)

                if not vetted_candidates:
                    logger.error("PROTOCOL FAILED: Junk filter discarded all candidates. The Forge is producing garbage.")
                    return

                # Pass the vetted candidates to the next phase
                validated_candidates = execute_triage_phase(target_motif, vetted_candidates)

                if not validated_candidates:
                    logger.error("PROTOCOL FAILED: Triage phase failed to validate any candidates.")
                    return

                # --- PHASE 4: THE CONFIRMATION ---
                execute_confirmation_phase(target_motif, validated_candidates)

                break # Exit polling loop
            
            elif status_data.get("status") == "failed":
                error_msg = status_data.get('error', 'Unknown error')
                logger.error(f"PROTOCOL FAILED: Ambush phase failed at the Forge. Reason: {error_msg}")
                return # End the protocol

    except Exception as e:
        logger.error(f"PROTOCOL FAILED: An error occurred during the Ambush phase: {e}")
        return

def execute_triage_phase(target_motif: str, candidates: list[str]) -> list[dict] | None:
    """
    Executes Phase 3 of the Predator Protocol: Triage.
    Uses the ZetaOracle to validate and rank the raw candidates.
    """
    logger.info("--- 🛡️ TRIAGE PHASE INITIATED: VALIDATING CANDIDATES WITH ZETAORACLE ---")
    logger.info(f"Validating {len(candidates)} candidates against target motif...")

    payload = {
        "target_dna_sequence": target_motif,
        "candidate_inhibitors": candidates
    }

    try:
        response = requests.post(
            f"{ZETA_ORACLE_URL}/validate_inhibitors", 
            json=payload, 
            timeout=300 # Allow more time for validation
        )
        response.raise_for_status()
        
        validation_result = response.json()
        ranked_candidates = validation_result.get("ranked_candidates")

        if not ranked_candidates:
            logger.warning("Triage complete, but no candidates were deemed viable by the Oracle.")
            return None

        # --- DOCTRINE: RECALIBRATE JUDGEMENT ---
        # The previous thresholds were dangerously optimistic. This new rubric
        # provides a more realistic assessment of candidate viability.
        for candidate in ranked_candidates:
            score = candidate.get("stability_score", -999)
            if score > 0.1:
                candidate["commentary"] = "Elite Stability. PRIMARY CANDIDATE."
            elif score > 0.0:
                candidate["commentary"] = "Viable Stability. Secondary Candidate."
            else:
                candidate["commentary"] = "Low Stability. Discard."

        logger.success(f"--- ✅ TRIAGE PHASE COMPLETE: {len(ranked_candidates)} CANDIDATES VALIDATED & RE-CALIBRATED ---")
        
        # Log the top candidates for review
        print(json.dumps(ranked_candidates[:5], indent=2)) # Show top 5
        
        return ranked_candidates

    except Exception as e:
        logger.error(f"PROTOCOL FAILED: An error occurred during the Triage phase: {e}")
        return None

def execute_confirmation_phase(target_motif_dna: str, validated_candidates: list[dict]):
    """
    Executes the final phase of the Predator Protocol: Confirmation.
    Takes the top candidate, translates it and the target to protein,
    and runs a final binding affinity simulation with Boltz-2.
    """
    logger.info("--- 🎯 CONFIRMATION PHASE INITIATED: RUNNING BIOPHYSICAL SIMULATION ---")
    
    # Select the #1 ranked candidate from the Triage phase
    top_candidate_dna = validated_candidates[0].get("inhibitor_sequence")
    if not top_candidate_dna:
        logger.error("CONFIRMATION FAILED: Top candidate has no sequence.")
        return

    logger.info(f"Top candidate selected for simulation. Stability Score: {validated_candidates[0].get('stability_score'):.4f}")

    # --- DNA-to-Protein Bridge ---
    logger.info("Translating target and candidate DNA to protein sequences...")
    target_protein = translate_dna_to_protein(target_motif_dna)
    candidate_protein = translate_dna_to_protein(top_candidate_dna)

    if not target_protein or not candidate_protein:
        logger.error("CONFIRMATION FAILED: DNA-to-Protein translation failed.")
        return
    
    logger.info("Translation successful.")
    logger.info(f"  - Target Protein ({len(target_protein)} AA): {target_protein[:60]}...")
    logger.info(f"  - Candidate Protein ({len(candidate_protein)} AA): {candidate_protein[:60]}...")

    # --- Final Simulation with Boltz-2 ---
    payload = {
        "target_sequence": target_protein,
        "candidate_sequences": [candidate_protein], # Boltz can take a list, but we send only the best
        "job_id": f"predator_confirmation_{int(time.time())}"
    }

    try:
        logger.info("Dispatching simulation command to Boltz service...")
        response = requests.post(BOLTZ_SERVICE_URL, json=payload, timeout=600) # Long timeout for simulation
        response.raise_for_status()
        final_blueprint = response.json()

        logger.success("--- ✅ CONFIRMATION PHASE COMPLETE ---")
        print("\n" + "="*80)
        print("          🐆🏆 PREDATOR PROTOCOL COMPLETE: FINAL WEAPON BLUEPRINT 🏆🐆")
        print("="*80)
        print(json.dumps(final_blueprint, indent=2))
        print("="*80 + "\n")

    except Exception as e:
        logger.error(f"PROTOCOL FAILED: An error occurred during the Boltz-2 simulation: {e}")

if __name__ == "__main__":
    # The protocol is initiated with a single gene symbol
    run_predator_protocol("VEGFA") 