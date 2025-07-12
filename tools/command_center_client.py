"""
CRISPR Assistant - Command Center Client

This module provides a client interface for interacting with the AI General Command Center.
It handles battle plan formulation, patient data management, and guide RNA design.
"""
import os
import httpx
import json
import random

# Import shared schemas from the service
from services.command_center.schemas import PatientDataPacket, BattlePlan, BattlePlanResponse, GuideDesignResponse

# --- CRUCIAL: The CommandCenter URL should be managed via environment variables ---
# This allows us to switch between local testing and live deployment without code changes.
# For local development, you might set this to "http://127.0.0.1:8000"
COMMAND_CENTER_URL = os.environ.get("COMMAND_CENTER_URL", "https://crispro--command-center-commandcenter-api.modal.run")

def initiate_germline_correction(corrected_sequence: str, target_site_context: str, arm_length: int = 500, log_callback=None):
    """
    Calls the CommandCenter's germline correction blueprint workflow.

    Args:
        corrected_sequence (str): The healthy/reference DNA sequence for the target region.
        target_site_context (str): The genomic coordinates for the broader context, e.g., a ~4kb region.
        arm_length (int, optional): The desired length of the homology arms. Defaults to 500.
        log_callback (function, optional): A function to stream log messages to.

    Returns:
        pandas.DataFrame: A DataFrame containing the prioritized list of guide RNA candidates,
                          enriched with placeholder scores for the UI.
    """
    endpoint = f"{COMMAND_CENTER_URL}/workflow/design_germline_correction"
    
    payload = {
        "corrected_sequence": corrected_sequence,
        "target_site_context": target_site_context,
        "arm_length": arm_length,
    }

    if log_callback:
        log_callback("Connecting to CommandCenter...")
        log_callback(f"  Endpoint: {endpoint}")
        payload_summary = payload.copy()
        payload_summary['corrected_sequence'] = f"{payload_summary['corrected_sequence'][:30]}..."
        log_callback(f"  Payload: {json.dumps(payload_summary, indent=2)}")

    try:
        response = requests.post(endpoint, json=payload, timeout=300, verify=False)
        response.raise_for_status()
        blueprint = response.json()

        if log_callback:
            log_callback("✅ Success! Received response from CommandCenter.")
            raw_response_str = json.dumps(blueprint, indent=2)
            log_callback(f"--- Raw CommandCenter Response ---\n{raw_response_str}")
        
        # --- Pass-through: Return the raw scored guides ---
        scored_guides_data = blueprint.get("scored_guides", [])
        if not scored_guides_data:
            if log_callback: log_callback("NOTE: Blueprint received, but no guide candidates could be scored.")
            return pd.DataFrame(), {}

        # Also return the full blueprint for the UI to use
        return pd.DataFrame(scored_guides_data), blueprint

    except requests.exceptions.RequestException as e:
        if log_callback:
            log_callback(f"--- CommandCenter API Error ---")
            log_callback(f"Failed to connect to {endpoint}. Is the service running?")
            log_callback(f"Error details: {e}")
        return pd.DataFrame(), {}
    except Exception as e:
        if log_callback:
            log_callback(f"--- UNEXPECTED ERROR ---")
            log_callback(f"An error occurred while processing the API response: {e}")
        return pd.DataFrame(), {} 

import requests
import json
import os

# To run this client, you must have the CommandCenter service running.
# You can deploy it with `modal deploy services/command_center/main.py`
# The URL will be provided upon deployment.

BASE_URL = "https://crispro--crispr-assistant-command-center-v2-commandcenter-api.modal.run"

def formulate_battle_plan(patient_data):
    """
    Sends a request to the CommandCenter to formulate a battle plan.
    """
    url = f"{BASE_URL}/v2/workflow/formulate_battle_plan"
    headers = {"Content-Type": "application/json"}
    
    try:
        response = requests.post(url, data=json.dumps(patient_data), headers=headers)
        response.raise_for_status()  # Raises an HTTPError for bad responses (4xx or 5xx)
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")
        return None

def get_all_patients():
    """
    Retrieves a list of all patients from the Command Center database.
    """
    url = f"{BASE_URL}/v2/data/patients"
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"An error occurred while fetching patients: {e}")
        return []

def get_latest_battle_plan_for_patient(patient_identifier: str):
    """
    Retrieves the most recent battle plan for a specific patient.
    """
    url = f"{BASE_URL}/v2/data/battle_plan/{patient_identifier}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"An error occurred while fetching battle plan for {patient_identifier}: {e}")
        return None


if __name__ == '__main__':
    # Example Usage
    print("Fetching all patients...")
    patients = get_all_patients()
    if patients:
        print("Available patients:", patients)
        
        # Fetch the latest battle plan for the first patient
        if patients:
            patient_id_to_test = patients[0]['patient_identifier']
            print(f"\nFetching latest battle plan for patient: {patient_id_to_test}")
            battle_plan = get_latest_battle_plan_for_patient(patient_id_to_test)
            if battle_plan:
                print("Battle Plan Details:")
                print(json.dumps(battle_plan, indent=2))
            else:
                print("Could not retrieve battle plan.")

    else:
        print("No patients found. Formulating a new battle plan as an example...")
        # This is a mock patient data packet. 
        # In a real scenario, this would come from a trusted source.
        mock_patient_data = {
            "patient_identifier": "DEMO_PATIENT_001",
            "mutation_hgvs_c": "c.1907G>A",
            "mutation_hgvs_p": "p.Arg636His",
            "gene": "BRCA1",
            "bam_file_path": "/data/synthetic_brca1.bam",  # This path must be accessible within the Modal container
            "sequence_for_perplexity": "ATGC...",  # A sufficiently long sequence for the model
            "protein_sequence": "MDLSAL...",  # The full protein sequence
            "transcript_id": "ENST00000357654"
        }
        
        new_plan = formulate_battle_plan(mock_patient_data)
        
        if new_plan:
            print("\nSuccessfully formulated new battle plan:")
            print(json.dumps(new_plan, indent=2)) 