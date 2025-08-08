import requests
import os
import json
from loguru import logger
import time

# --- CONFIGURATION ---
ZETA_FORGE_URL = os.environ.get("ZETA_FORGE_URL", "https://crispro-zeta-forge-v1-zetaforge-api.modal.run")
BOLTZ_SERVICE_URL = "https://crispro--boltz-service-fastapi-app.modal.run"

TARGET_GENE_SYMBOL = "VEGFA"
TARGET_MOTIF = "CCGACGGACCCGGACCGGAC" # A sample motif from VEGFA for testing

# --- UTILS ---
def print_header(title):
    logger.info("=" * 80)
    logger.info(f"   {title.upper()}")
    logger.info("=" * 80)

def print_final_blueprint(blueprint):
    print_header("FINAL WEAPON BLUEPRINT")
    if not blueprint.get("ranked_candidates"):
        logger.warning("No candidates were returned from the Boltz service.")
        return

    for i, candidate in enumerate(blueprint["ranked_candidates"]):
        logger.info(
            f"  Rank {i+1:02d}: {candidate['sequence']} "
            f"| Affinity: {candidate['binding_affinity']:.4f} "
            f"| Commentary: {candidate['commentary']}"
        )
    logger.info("-" * 80)

# --- CORE TEST LOGIC ---
def run_forge_and_boltz_test():
    """
    Executes the full "Forge-and-Boltz" protocol.
    1. Forge: Generates candidate nanobodies using ZetaForge.
    2. Screen: Screens candidates for binding affinity using Boltz-2.
    """
    try:
        # --- PHASE I: FORGE ---
        print_header(f"PHASE I: FORGING CANDIDATES FOR TARGET {TARGET_GENE_SYMBOL}")
        
        # We use the reverse complement as "bait" for the generative model
        bait_sequence = str(Bio.Seq.Seq(TARGET_MOTIF).reverse_complement())
        
        forge_payload = {
            "bait_sequence": bait_sequence,
            "num_candidates_per_temp": 15,
            "temperatures": [0.2, 0.7, 1.2],
            "generation_length": 30
        }
        
        logger.info(f"Dispatching Forge command to ZetaForge with payload: {forge_payload}")
        forge_response = requests.post(
            f"{ZETA_FORGE_URL}/generate_candidates", 
            json=forge_payload,
            timeout=600
        )
        forge_response.raise_for_status()
        candidates = forge_response.json().get("candidates", [])
        
        if not candidates:
            logger.error("ZetaForge failed to generate any candidates. Aborting mission.")
            return

        logger.success(f"--- ✅ FORGE COMPLETE: {len(candidates)} unique candidates generated. ---")

        # --- PHASE II: SCREEN ---
        print_header("PHASE II: SCREENING CANDIDATES WITH BOLTZ-2")
        
        screen_payload = {
            "target_sequence": TARGET_MOTIF,
            "candidate_sequences": candidates,
            "job_id": f"job_{int(time.time())}"
        }

        logger.info(f"Dispatching Screen command to Boltz service...")
        screen_response = requests.post(
            BOLTZ_SERVICE_URL, 
            json=screen_payload, 
            timeout=600
        )
        screen_response.raise_for_status()
        final_blueprint = screen_response.json()

        print_final_blueprint(final_blueprint)
        logger.success("--- ✅ MISSION COMPLETE: FORGE-AND-BOLTZ PROTOCOL EXECUTED SUCCESSFULLY ---")

    except requests.exceptions.HTTPError as e:
        logger.error(f"A fatal HTTP error occurred: {e.response.status_code} - {e.response.text}")
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    # A simple hack to get BioPython without adding it as a dependency for this test script
    try:
        import Bio.Seq
    except ImportError:
        os.system("pip install biopython")
        import Bio.Seq
        
    run_forge_and_boltz_test() 
import os
import json
from loguru import logger
import time

# --- CONFIGURATION ---
ZETA_FORGE_URL = os.environ.get("ZETA_FORGE_URL", "https://crispro-zeta-forge-v1-zetaforge-api.modal.run")
BOLTZ_SERVICE_URL = "https://crispro--boltz-service-fastapi-app.modal.run"

TARGET_GENE_SYMBOL = "VEGFA"
TARGET_MOTIF = "CCGACGGACCCGGACCGGAC" # A sample motif from VEGFA for testing

# --- UTILS ---
def print_header(title):
    logger.info("=" * 80)
    logger.info(f"   {title.upper()}")
    logger.info("=" * 80)

def print_final_blueprint(blueprint):
    print_header("FINAL WEAPON BLUEPRINT")
    if not blueprint.get("ranked_candidates"):
        logger.warning("No candidates were returned from the Boltz service.")
        return

    for i, candidate in enumerate(blueprint["ranked_candidates"]):
        logger.info(
            f"  Rank {i+1:02d}: {candidate['sequence']} "
            f"| Affinity: {candidate['binding_affinity']:.4f} "
            f"| Commentary: {candidate['commentary']}"
        )
    logger.info("-" * 80)

# --- CORE TEST LOGIC ---
def run_forge_and_boltz_test():
    """
    Executes the full "Forge-and-Boltz" protocol.
    1. Forge: Generates candidate nanobodies using ZetaForge.
    2. Screen: Screens candidates for binding affinity using Boltz-2.
    """
    try:
        # --- PHASE I: FORGE ---
        print_header(f"PHASE I: FORGING CANDIDATES FOR TARGET {TARGET_GENE_SYMBOL}")
        
        # We use the reverse complement as "bait" for the generative model
        bait_sequence = str(Bio.Seq.Seq(TARGET_MOTIF).reverse_complement())
        
        forge_payload = {
            "bait_sequence": bait_sequence,
            "num_candidates_per_temp": 15,
            "temperatures": [0.2, 0.7, 1.2],
            "generation_length": 30
        }
        
        logger.info(f"Dispatching Forge command to ZetaForge with payload: {forge_payload}")
        forge_response = requests.post(
            f"{ZETA_FORGE_URL}/generate_candidates", 
            json=forge_payload,
            timeout=600
        )
        forge_response.raise_for_status()
        candidates = forge_response.json().get("candidates", [])
        
        if not candidates:
            logger.error("ZetaForge failed to generate any candidates. Aborting mission.")
            return

        logger.success(f"--- ✅ FORGE COMPLETE: {len(candidates)} unique candidates generated. ---")

        # --- PHASE II: SCREEN ---
        print_header("PHASE II: SCREENING CANDIDATES WITH BOLTZ-2")
        
        screen_payload = {
            "target_sequence": TARGET_MOTIF,
            "candidate_sequences": candidates,
            "job_id": f"job_{int(time.time())}"
        }

        logger.info(f"Dispatching Screen command to Boltz service...")
        screen_response = requests.post(
            BOLTZ_SERVICE_URL, 
            json=screen_payload, 
            timeout=600
        )
        screen_response.raise_for_status()
        final_blueprint = screen_response.json()

        print_final_blueprint(final_blueprint)
        logger.success("--- ✅ MISSION COMPLETE: FORGE-AND-BOLTZ PROTOCOL EXECUTED SUCCESSFULLY ---")

    except requests.exceptions.HTTPError as e:
        logger.error(f"A fatal HTTP error occurred: {e.response.status_code} - {e.response.text}")
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    # A simple hack to get BioPython without adding it as a dependency for this test script
    try:
        import Bio.Seq
    except ImportError:
        os.system("pip install biopython")
        import Bio.Seq
        
    run_forge_and_boltz_test() 
import os
import json
from loguru import logger
import time

# --- CONFIGURATION ---
ZETA_FORGE_URL = os.environ.get("ZETA_FORGE_URL", "https://crispro-zeta-forge-v1-zetaforge-api.modal.run")
BOLTZ_SERVICE_URL = "https://crispro--boltz-service-fastapi-app.modal.run"

TARGET_GENE_SYMBOL = "VEGFA"
TARGET_MOTIF = "CCGACGGACCCGGACCGGAC" # A sample motif from VEGFA for testing

# --- UTILS ---
def print_header(title):
    logger.info("=" * 80)
    logger.info(f"   {title.upper()}")
    logger.info("=" * 80)

def print_final_blueprint(blueprint):
    print_header("FINAL WEAPON BLUEPRINT")
    if not blueprint.get("ranked_candidates"):
        logger.warning("No candidates were returned from the Boltz service.")
        return

    for i, candidate in enumerate(blueprint["ranked_candidates"]):
        logger.info(
            f"  Rank {i+1:02d}: {candidate['sequence']} "
            f"| Affinity: {candidate['binding_affinity']:.4f} "
            f"| Commentary: {candidate['commentary']}"
        )
    logger.info("-" * 80)

# --- CORE TEST LOGIC ---
def run_forge_and_boltz_test():
    """
    Executes the full "Forge-and-Boltz" protocol.
    1. Forge: Generates candidate nanobodies using ZetaForge.
    2. Screen: Screens candidates for binding affinity using Boltz-2.
    """
    try:
        # --- PHASE I: FORGE ---
        print_header(f"PHASE I: FORGING CANDIDATES FOR TARGET {TARGET_GENE_SYMBOL}")
        
        # We use the reverse complement as "bait" for the generative model
        bait_sequence = str(Bio.Seq.Seq(TARGET_MOTIF).reverse_complement())
        
        forge_payload = {
            "bait_sequence": bait_sequence,
            "num_candidates_per_temp": 15,
            "temperatures": [0.2, 0.7, 1.2],
            "generation_length": 30
        }
        
        logger.info(f"Dispatching Forge command to ZetaForge with payload: {forge_payload}")
        forge_response = requests.post(
            f"{ZETA_FORGE_URL}/generate_candidates", 
            json=forge_payload,
            timeout=600
        )
        forge_response.raise_for_status()
        candidates = forge_response.json().get("candidates", [])
        
        if not candidates:
            logger.error("ZetaForge failed to generate any candidates. Aborting mission.")
            return

        logger.success(f"--- ✅ FORGE COMPLETE: {len(candidates)} unique candidates generated. ---")

        # --- PHASE II: SCREEN ---
        print_header("PHASE II: SCREENING CANDIDATES WITH BOLTZ-2")
        
        screen_payload = {
            "target_sequence": TARGET_MOTIF,
            "candidate_sequences": candidates,
            "job_id": f"job_{int(time.time())}"
        }

        logger.info(f"Dispatching Screen command to Boltz service...")
        screen_response = requests.post(
            BOLTZ_SERVICE_URL, 
            json=screen_payload, 
            timeout=600
        )
        screen_response.raise_for_status()
        final_blueprint = screen_response.json()

        print_final_blueprint(final_blueprint)
        logger.success("--- ✅ MISSION COMPLETE: FORGE-AND-BOLTZ PROTOCOL EXECUTED SUCCESSFULLY ---")

    except requests.exceptions.HTTPError as e:
        logger.error(f"A fatal HTTP error occurred: {e.response.status_code} - {e.response.text}")
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    # A simple hack to get BioPython without adding it as a dependency for this test script
    try:
        import Bio.Seq
    except ImportError:
        os.system("pip install biopython")
        import Bio.Seq
        
    run_forge_and_boltz_test() 
import os
import json
from loguru import logger
import time

# --- CONFIGURATION ---
ZETA_FORGE_URL = os.environ.get("ZETA_FORGE_URL", "https://crispro-zeta-forge-v1-zetaforge-api.modal.run")
BOLTZ_SERVICE_URL = "https://crispro--boltz-service-fastapi-app.modal.run"

TARGET_GENE_SYMBOL = "VEGFA"
TARGET_MOTIF = "CCGACGGACCCGGACCGGAC" # A sample motif from VEGFA for testing

# --- UTILS ---
def print_header(title):
    logger.info("=" * 80)
    logger.info(f"   {title.upper()}")
    logger.info("=" * 80)

def print_final_blueprint(blueprint):
    print_header("FINAL WEAPON BLUEPRINT")
    if not blueprint.get("ranked_candidates"):
        logger.warning("No candidates were returned from the Boltz service.")
        return

    for i, candidate in enumerate(blueprint["ranked_candidates"]):
        logger.info(
            f"  Rank {i+1:02d}: {candidate['sequence']} "
            f"| Affinity: {candidate['binding_affinity']:.4f} "
            f"| Commentary: {candidate['commentary']}"
        )
    logger.info("-" * 80)

# --- CORE TEST LOGIC ---
def run_forge_and_boltz_test():
    """
    Executes the full "Forge-and-Boltz" protocol.
    1. Forge: Generates candidate nanobodies using ZetaForge.
    2. Screen: Screens candidates for binding affinity using Boltz-2.
    """
    try:
        # --- PHASE I: FORGE ---
        print_header(f"PHASE I: FORGING CANDIDATES FOR TARGET {TARGET_GENE_SYMBOL}")
        
        # We use the reverse complement as "bait" for the generative model
        bait_sequence = str(Bio.Seq.Seq(TARGET_MOTIF).reverse_complement())
        
        forge_payload = {
            "bait_sequence": bait_sequence,
            "num_candidates_per_temp": 15,
            "temperatures": [0.2, 0.7, 1.2],
            "generation_length": 30
        }
        
        logger.info(f"Dispatching Forge command to ZetaForge with payload: {forge_payload}")
        forge_response = requests.post(
            f"{ZETA_FORGE_URL}/generate_candidates", 
            json=forge_payload,
            timeout=600
        )
        forge_response.raise_for_status()
        candidates = forge_response.json().get("candidates", [])
        
        if not candidates:
            logger.error("ZetaForge failed to generate any candidates. Aborting mission.")
            return

        logger.success(f"--- ✅ FORGE COMPLETE: {len(candidates)} unique candidates generated. ---")

        # --- PHASE II: SCREEN ---
        print_header("PHASE II: SCREENING CANDIDATES WITH BOLTZ-2")
        
        screen_payload = {
            "target_sequence": TARGET_MOTIF,
            "candidate_sequences": candidates,
            "job_id": f"job_{int(time.time())}"
        }

        logger.info(f"Dispatching Screen command to Boltz service...")
        screen_response = requests.post(
            BOLTZ_SERVICE_URL, 
            json=screen_payload, 
            timeout=600
        )
        screen_response.raise_for_status()
        final_blueprint = screen_response.json()

        print_final_blueprint(final_blueprint)
        logger.success("--- ✅ MISSION COMPLETE: FORGE-AND-BOLTZ PROTOCOL EXECUTED SUCCESSFULLY ---")

    except requests.exceptions.HTTPError as e:
        logger.error(f"A fatal HTTP error occurred: {e.response.status_code} - {e.response.text}")
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    # A simple hack to get BioPython without adding it as a dependency for this test script
    try:
        import Bio.Seq
    except ImportError:
        os.system("pip install biopython")
        import Bio.Seq
        
    run_forge_and_boltz_test() 