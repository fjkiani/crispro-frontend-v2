import requests
import json
from loguru import logger

# The local server URL for the CommandCenter
# NOTE: When running with `modal serve`, the default port is 8000.
COMMAND_CENTER_URL = "http://127.0.0.1:8000"
DOSSIER_ENDPOINT = f"{COMMAND_CENTER_URL}/workflow/generate_intelligence_dossier"

def run_dossier_generation_test():
    """
    Calls the new, dedicated dossier generation endpoint to test
    the Sniper's primary weapon.
    """
    logger.info("---  Sniper Mission Rehearsal: Dossier Forge ---")
    
    # This payload matches the DossierRequest Pydantic model in the CommandCenter
    payload = {
        "patient_identifier": "SNIPER-TEST-001",
        "gene": "BRAF",
        "mutation_hgvs_p": "p.Val600Glu", # The infamous V600E mutation
        "protein_sequence": "M...", # In a real test, this would be the full sequence
        "locus": "chr7:140753336" # Locus for V600E
    }

    logger.info(f"Requesting dossier for high-value target: {payload['gene']} {payload['mutation_hgvs_p']}")

    try:
        response = requests.post(DOSSIER_ENDPOINT, json=payload, timeout=400)
        response.raise_for_status()
        
        data = response.json()
        
        logger.success("✅ MISSION SUCCESS: Received valid dossier from CommandCenter.")
        
        # --- Displaying the Dossier ---
        logger.info("\n" + "#" * 50)
        logger.info("## COMMANDER'S INTELLIGENCE DOSSIER ##")
        logger.info("#" * 50 + "\n")
        
        logger.info(f"**TARGET:** {data['request']['gene']} {data['request']['mutation_hgvs_p']}")
        logger.info(f"**VERDICT:** {data['verdict']}\n")
        
        logger.info("**ZETA ORACLE ASSESSMENT:**")
        tm = data['threat_matrix']
        logger.info(f"  - Live Zeta Score: {tm['patient_zeta_score']:.4f}")
        logger.info(f"  - Pathogenic Range: {tm['known_pathogenic_min']} to {tm['known_pathogenic_max']}")
        logger.info(f"  - Benign Range: {tm['known_benign_min']} to {tm['known_benign_max']}\n")

        logger.info("**LLM CLINICAL ANALYSIS:**")
        ca = data['clinical_analysis']
        logger.info(f"  - **Abstract:** {ca['abstract']}")
        logger.info(f"  - **Mechanism:** {ca['mechanism']}")
        logger.info(f"  - **Significance:** {ca['significance']}")
        logger.info(f"  - **Therapeutics:** {ca['therapeutics']}\n")

        logger.info("#" * 50)

    except requests.exceptions.HTTPError as http_err:
        logger.error(f"❌ MISSION FAILED: HTTP error occurred: {http_err}")
        logger.error(f"Response Content: {http_err.response.text}")
    except requests.exceptions.RequestException as req_err:
        logger.error(f"❌ MISSION FAILED: A critical request error occurred: {req_err}")
        logger.error("  Is the CommandCenter service running locally? Try: `modal serve src/services/command_center/main.py`")
    except json.JSONDecodeError:
        logger.error("❌ MISSION FAILED: Failed to decode JSON response from the server.")
        print("Raw Response:", response.text)

if __name__ == "__main__":
    run_dossier_generation_test() 
import json
from loguru import logger

# The local server URL for the CommandCenter
# NOTE: When running with `modal serve`, the default port is 8000.
COMMAND_CENTER_URL = "http://127.0.0.1:8000"
DOSSIER_ENDPOINT = f"{COMMAND_CENTER_URL}/workflow/generate_intelligence_dossier"

def run_dossier_generation_test():
    """
    Calls the new, dedicated dossier generation endpoint to test
    the Sniper's primary weapon.
    """
    logger.info("---  Sniper Mission Rehearsal: Dossier Forge ---")
    
    # This payload matches the DossierRequest Pydantic model in the CommandCenter
    payload = {
        "patient_identifier": "SNIPER-TEST-001",
        "gene": "BRAF",
        "mutation_hgvs_p": "p.Val600Glu", # The infamous V600E mutation
        "protein_sequence": "M...", # In a real test, this would be the full sequence
        "locus": "chr7:140753336" # Locus for V600E
    }

    logger.info(f"Requesting dossier for high-value target: {payload['gene']} {payload['mutation_hgvs_p']}")

    try:
        response = requests.post(DOSSIER_ENDPOINT, json=payload, timeout=400)
        response.raise_for_status()
        
        data = response.json()
        
        logger.success("✅ MISSION SUCCESS: Received valid dossier from CommandCenter.")
        
        # --- Displaying the Dossier ---
        logger.info("\n" + "#" * 50)
        logger.info("## COMMANDER'S INTELLIGENCE DOSSIER ##")
        logger.info("#" * 50 + "\n")
        
        logger.info(f"**TARGET:** {data['request']['gene']} {data['request']['mutation_hgvs_p']}")
        logger.info(f"**VERDICT:** {data['verdict']}\n")
        
        logger.info("**ZETA ORACLE ASSESSMENT:**")
        tm = data['threat_matrix']
        logger.info(f"  - Live Zeta Score: {tm['patient_zeta_score']:.4f}")
        logger.info(f"  - Pathogenic Range: {tm['known_pathogenic_min']} to {tm['known_pathogenic_max']}")
        logger.info(f"  - Benign Range: {tm['known_benign_min']} to {tm['known_benign_max']}\n")

        logger.info("**LLM CLINICAL ANALYSIS:**")
        ca = data['clinical_analysis']
        logger.info(f"  - **Abstract:** {ca['abstract']}")
        logger.info(f"  - **Mechanism:** {ca['mechanism']}")
        logger.info(f"  - **Significance:** {ca['significance']}")
        logger.info(f"  - **Therapeutics:** {ca['therapeutics']}\n")

        logger.info("#" * 50)

    except requests.exceptions.HTTPError as http_err:
        logger.error(f"❌ MISSION FAILED: HTTP error occurred: {http_err}")
        logger.error(f"Response Content: {http_err.response.text}")
    except requests.exceptions.RequestException as req_err:
        logger.error(f"❌ MISSION FAILED: A critical request error occurred: {req_err}")
        logger.error("  Is the CommandCenter service running locally? Try: `modal serve src/services/command_center/main.py`")
    except json.JSONDecodeError:
        logger.error("❌ MISSION FAILED: Failed to decode JSON response from the server.")
        print("Raw Response:", response.text)

if __name__ == "__main__":
    run_dossier_generation_test() 
import json
from loguru import logger

# The local server URL for the CommandCenter
# NOTE: When running with `modal serve`, the default port is 8000.
COMMAND_CENTER_URL = "http://127.0.0.1:8000"
DOSSIER_ENDPOINT = f"{COMMAND_CENTER_URL}/workflow/generate_intelligence_dossier"

def run_dossier_generation_test():
    """
    Calls the new, dedicated dossier generation endpoint to test
    the Sniper's primary weapon.
    """
    logger.info("---  Sniper Mission Rehearsal: Dossier Forge ---")
    
    # This payload matches the DossierRequest Pydantic model in the CommandCenter
    payload = {
        "patient_identifier": "SNIPER-TEST-001",
        "gene": "BRAF",
        "mutation_hgvs_p": "p.Val600Glu", # The infamous V600E mutation
        "protein_sequence": "M...", # In a real test, this would be the full sequence
        "locus": "chr7:140753336" # Locus for V600E
    }

    logger.info(f"Requesting dossier for high-value target: {payload['gene']} {payload['mutation_hgvs_p']}")

    try:
        response = requests.post(DOSSIER_ENDPOINT, json=payload, timeout=400)
        response.raise_for_status()
        
        data = response.json()
        
        logger.success("✅ MISSION SUCCESS: Received valid dossier from CommandCenter.")
        
        # --- Displaying the Dossier ---
        logger.info("\n" + "#" * 50)
        logger.info("## COMMANDER'S INTELLIGENCE DOSSIER ##")
        logger.info("#" * 50 + "\n")
        
        logger.info(f"**TARGET:** {data['request']['gene']} {data['request']['mutation_hgvs_p']}")
        logger.info(f"**VERDICT:** {data['verdict']}\n")
        
        logger.info("**ZETA ORACLE ASSESSMENT:**")
        tm = data['threat_matrix']
        logger.info(f"  - Live Zeta Score: {tm['patient_zeta_score']:.4f}")
        logger.info(f"  - Pathogenic Range: {tm['known_pathogenic_min']} to {tm['known_pathogenic_max']}")
        logger.info(f"  - Benign Range: {tm['known_benign_min']} to {tm['known_benign_max']}\n")

        logger.info("**LLM CLINICAL ANALYSIS:**")
        ca = data['clinical_analysis']
        logger.info(f"  - **Abstract:** {ca['abstract']}")
        logger.info(f"  - **Mechanism:** {ca['mechanism']}")
        logger.info(f"  - **Significance:** {ca['significance']}")
        logger.info(f"  - **Therapeutics:** {ca['therapeutics']}\n")

        logger.info("#" * 50)

    except requests.exceptions.HTTPError as http_err:
        logger.error(f"❌ MISSION FAILED: HTTP error occurred: {http_err}")
        logger.error(f"Response Content: {http_err.response.text}")
    except requests.exceptions.RequestException as req_err:
        logger.error(f"❌ MISSION FAILED: A critical request error occurred: {req_err}")
        logger.error("  Is the CommandCenter service running locally? Try: `modal serve src/services/command_center/main.py`")
    except json.JSONDecodeError:
        logger.error("❌ MISSION FAILED: Failed to decode JSON response from the server.")
        print("Raw Response:", response.text)

if __name__ == "__main__":
    run_dossier_generation_test() 
import json
from loguru import logger

# The local server URL for the CommandCenter
# NOTE: When running with `modal serve`, the default port is 8000.
COMMAND_CENTER_URL = "http://127.0.0.1:8000"
DOSSIER_ENDPOINT = f"{COMMAND_CENTER_URL}/workflow/generate_intelligence_dossier"

def run_dossier_generation_test():
    """
    Calls the new, dedicated dossier generation endpoint to test
    the Sniper's primary weapon.
    """
    logger.info("---  Sniper Mission Rehearsal: Dossier Forge ---")
    
    # This payload matches the DossierRequest Pydantic model in the CommandCenter
    payload = {
        "patient_identifier": "SNIPER-TEST-001",
        "gene": "BRAF",
        "mutation_hgvs_p": "p.Val600Glu", # The infamous V600E mutation
        "protein_sequence": "M...", # In a real test, this would be the full sequence
        "locus": "chr7:140753336" # Locus for V600E
    }

    logger.info(f"Requesting dossier for high-value target: {payload['gene']} {payload['mutation_hgvs_p']}")

    try:
        response = requests.post(DOSSIER_ENDPOINT, json=payload, timeout=400)
        response.raise_for_status()
        
        data = response.json()
        
        logger.success("✅ MISSION SUCCESS: Received valid dossier from CommandCenter.")
        
        # --- Displaying the Dossier ---
        logger.info("\n" + "#" * 50)
        logger.info("## COMMANDER'S INTELLIGENCE DOSSIER ##")
        logger.info("#" * 50 + "\n")
        
        logger.info(f"**TARGET:** {data['request']['gene']} {data['request']['mutation_hgvs_p']}")
        logger.info(f"**VERDICT:** {data['verdict']}\n")
        
        logger.info("**ZETA ORACLE ASSESSMENT:**")
        tm = data['threat_matrix']
        logger.info(f"  - Live Zeta Score: {tm['patient_zeta_score']:.4f}")
        logger.info(f"  - Pathogenic Range: {tm['known_pathogenic_min']} to {tm['known_pathogenic_max']}")
        logger.info(f"  - Benign Range: {tm['known_benign_min']} to {tm['known_benign_max']}\n")

        logger.info("**LLM CLINICAL ANALYSIS:**")
        ca = data['clinical_analysis']
        logger.info(f"  - **Abstract:** {ca['abstract']}")
        logger.info(f"  - **Mechanism:** {ca['mechanism']}")
        logger.info(f"  - **Significance:** {ca['significance']}")
        logger.info(f"  - **Therapeutics:** {ca['therapeutics']}\n")

        logger.info("#" * 50)

    except requests.exceptions.HTTPError as http_err:
        logger.error(f"❌ MISSION FAILED: HTTP error occurred: {http_err}")
        logger.error(f"Response Content: {http_err.response.text}")
    except requests.exceptions.RequestException as req_err:
        logger.error(f"❌ MISSION FAILED: A critical request error occurred: {req_err}")
        logger.error("  Is the CommandCenter service running locally? Try: `modal serve src/services/command_center/main.py`")
    except json.JSONDecodeError:
        logger.error("❌ MISSION FAILED: Failed to decode JSON response from the server.")
        print("Raw Response:", response.text)

if __name__ == "__main__":
    run_dossier_generation_test() 