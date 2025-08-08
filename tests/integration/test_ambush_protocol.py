import requests
import time
from loguru import logger
from Bio.Seq import Seq

# --- Configuration ---
ZETA_FORGE_URL = "https://crispro--zeta-forge-v1-api.modal.run"
# --- DOCTRINE: Use the live ZetaOracle for the final validation step ---
ZETA_ORACLE_URL = "https://crispro--zeta-oracle-v2-zetaoracle-api.modal.run"
TARGET_MOTIF = "CGCCTGTAATCCCAGCACT"  # A critical motif from VEGFA for this test

def run_ambush_test():
    """
    Executes a live-fire test of the new "Competitive Generation" doctrine.
    """
    logger.info("--- 🔥 AMBUSH PROTOCOL TEST INITIATED 🔥 ---")

    # Phase I: Bait Preparation
    bait_sequence = str(Seq(TARGET_MOTIF).reverse_complement())
    logger.info(f"Target Motif: {TARGET_MOTIF}")
    logger.info(f"Forged Bait Sequence (Reverse Complement): {bait_sequence}")

    # Phase II: Competitive Generation (The Ambush)
    # --- DOCTRINE: SATURATION BOMBARDMENT ---
    # The previous approach used a small number of candidates, resulting in mediocre
    # results. We are now overwhelming the solution space by generating a much larger
    # slate of candidates (45 instead of 9) to increase our chances of finding a
    # high-stability binder.
    payload = {
        "bait_sequence": bait_sequence,
        "num_candidates_per_temp": 15,
        "temperatures": [0.2, 0.7, 1.2],
        "generation_length": 30
    }
    
    logger.info(f"Dispatching Ambush command to ZetaForge with payload: {payload}")
    
    start_time = time.time()
    try:
        # 1. Submit the job
        response = requests.post(f"{ZETA_FORGE_URL}/generate_inhibitor", json=payload, timeout=30)
        response.raise_for_status()
        job_id = response.json().get("job_id")
        if not job_id:
            raise RuntimeError("Failed to get job_id from ZetaForge submission.")
        logger.success(f"✅ Job submitted to ZetaForge. Job ID: {job_id}")
        logger.info("Polling for results...")

        # 2. Poll for results
        while True:
            time.sleep(10)
            status_response = requests.get(f"{ZETA_FORGE_URL}/status/{job_id}", timeout=30)
            status_response.raise_for_status()
            status_data = status_response.json()

            logger.info(f"Forge status: {status_data.get('status')} - {status_data.get('message')}")

            if status_data.get("status") == "complete":
                logger.success("✅🔥 AMBUSH SUCCESSFUL! ✅✅")
                # print_blueprint(status_data) # We will print the final ranked list later
                
                # --- PHASE III: THE CONFIRMATION KILL ---
                candidates = status_data.get('candidates', [])
                if not candidates:
                    raise RuntimeError("Forge returned no candidates to validate.")
                
                logger.info(f"--- 💀 CONFIRMATION KILL INITIATED: Validating {len(candidates)} candidates... ---")
                
                validation_payload = {
                    "target_dna_sequence": TARGET_MOTIF,
                    "candidate_inhibitors": candidates
                }
                
                validation_response = requests.post(
                    f"{ZETA_ORACLE_URL}/validate_inhibitors", 
                    json=validation_payload, 
                    timeout=300
                )
                validation_response.raise_for_status()
                final_blueprint = validation_response.json()

                print_final_blueprint(final_blueprint)
                break
            elif status_data.get("status") == "failed":
                logger.error(f"❌ AMBUSH FAILED: {status_data.get('error')}")
                break
        
        end_time = time.time()
        logger.info(f"Total operation time: {end_time - start_time:.2f} seconds.")

    except requests.exceptions.RequestException as e:
        logger.error(f"A network error occurred: {e}")
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")


def print_blueprint(result: dict):
    """Prints the final forged weapon blueprint."""
    print("\n" + "="*50)
    print("          FORGED CANDIDATE BLUEPRINTS")
    print("="*50)
    print(f"Status: {result.get('status')}")
    print(f"Message: {result.get('message')}")
    candidates = result.get('candidates', [])
    print(f"Generated {len(candidates)} unique candidates:")
    for i, seq in enumerate(candidates):
        print(f"  {i+1:02d}: {seq}")
    print("="*50 + "\n")


def print_final_blueprint(result: dict):
    """Prints the final, ranked, and validated weapon blueprints."""
    print("\n" + "="*60)
    print("          FINAL FORGED WEAPON BLUEPRINT (VALIDATED)")
    print("="*60)
    ranked_candidates = result.get('ranked_candidates', [])
    print(f"Validated and ranked {len(ranked_candidates)} candidates:")
    for i, data in enumerate(ranked_candidates):
        print(f"  Rank {i+1:02d}:")
        print(f"    - Inhibitor Sequence: {data.get('inhibitor_sequence')}")
        print(f"    - Stability Score (Zeta): {data.get('stability_score'):.4f}")
        print(f"    - Commentary: {data.get('commentary')}")
        print("-" * 20)
    print("="*60 + "\n")


if __name__ == "__main__":
    run_ambush_test() 