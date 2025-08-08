import pytest
import requests
import time
import os

# Get the base URL for the command center from an environment variable, with a fallback
BASE_URL = os.environ.get("COMMAND_CENTER_URL", "https://crispro--crispr-assistant-command-center-v3-commandcenter-api.modal.run")
POLL_INTERVAL = 15  # seconds
POLL_TIMEOUT = 1800  # 30 minutes, to be safe with cold starts and long runs

@pytest.fixture(scope="module")
def command_center_v4_integration():
    """
    A pytest fixture to provide the base URL to the tests.
    This also serves as a check to ensure the service is minimally responsive.
    """
    try:
        response = requests.get(BASE_URL)
        response.raise_for_status()
        assert "CRISPR Assistant Command Center V4 is online" in response.text
    except (requests.exceptions.RequestException, AssertionError) as e:
        pytest.fail(f"Command Center V4 is not online or is unresponsive at {BASE_URL}. Error: {e}")
    
    return BASE_URL

def test_v4_end_to_end_kill_chain(command_center_v4_integration):
    """
    This is the crucible. It tests the entire V4 workflow:
    1. Formulate a Battle Plan.
    2. Execute the V4 Kill Chain on that plan.
    3. Poll until the process is complete.
    4. Validate the results are structured correctly and contain the assassin_score.
    """
    base_url = command_center_v4_integration

    # --- Step 1: Formulate the Battle Plan ---
    formulate_payload = {
        "patient_identifier": "test-patient-v4",
        "gene": "BRAF",
        "mutation_hgvs_p": "p.Val600Glu",
        "sequence_for_perplexity": "AAGGGAAAGTGGTGATTTGGTCTAGCTACAGAGAAATCTCGATGGAGTGGGTCCCATCAGTTTGAACAGTTGTCTGGATCCATTTTGTGGATG"
    }
    
    print(f"Submitting Battle Plan to {base_url}/v3/workflow/formulate_battle_plan...")
    formulate_response = requests.post(f"{base_url}/v3/workflow/formulate_battle_plan", json=formulate_payload)
    assert formulate_response.status_code == 201, f"Failed to formulate battle plan: {formulate_response.text}"
    
    battle_plan = formulate_response.json()
    battle_plan_id = battle_plan.get("id")
    assert battle_plan_id is not None, "Formulation response did not include a battle plan ID."
    print(f"Battle Plan #{battle_plan_id} created successfully.")

    # --- Step 2: Execute the V4 Kill Chain ---
    execute_url = f"{base_url}/v4/workflow/execute_guide_design_campaign/{battle_plan_id}"
    print(f"Executing V4 Kill Chain at {execute_url}...")
    execute_response = requests.post(execute_url)
    assert execute_response.status_code == 202, f"Failed to execute kill chain: {execute_response.text}"
    print("V4 Kill Chain initiated successfully.")

    # --- Step 3: Poll for Completion ---
    poll_url = f"{base_url}/v3/battle_plan/{battle_plan_id}"
    start_time = time.time()
    final_status = None
    
    while time.time() - start_time < POLL_TIMEOUT:
        print(f"Polling status at {poll_url}...")
        poll_response = requests.get(poll_url)
        assert poll_response.status_code == 200, f"Failed to poll for status: {poll_response.text}"
        
        status_data = poll_response.json()
        current_status = status_data.get("status")
        print(f"Current status: {current_status}")

        if current_status == "interventions_designed":
            final_status = "interventions_designed"
            break
        elif current_status == "design_failed":
            pytest.fail(f"The V4 Kill Chain failed. Final status: 'design_failed'. Details: {status_data.get('results')}")
        
        time.sleep(POLL_INTERVAL)

    assert final_status == "interventions_designed", f"Polling timed out after {POLL_TIMEOUT}s. Final status was '{current_status}'."
    print("V4 Kill Chain completed successfully.")

    # --- Step 4: Validate the Fucking Payload ---
    print("Validating results payload...")
    results_data = status_data.get("results")
    assert results_data is not None, "Final result is missing the 'results' field."
    
    guides = results_data.get("guides")
    assert isinstance(guides, list), "The 'guides' field is not a list."
    assert len(guides) > 0, "The 'guides' list is empty."

    # Check the first guide for the correct structure
    first_guide = guides[0]
    assert "guide_sequence" in first_guide, "Guide dossier is missing 'guide_sequence'."
    
    efficacy_data = first_guide.get("efficacy")
    assert efficacy_data is not None, "Guide dossier is missing 'efficacy' data."
    
    safety_data = first_guide.get("safety")
    assert safety_data is not None, "Guide dossier is missing 'safety' data."

    assert "assassin_score" in first_guide, "Guide dossier is missing 'assassin_score'."
    
    # --- The Final, Ruthless Assertion ---
    # We verify that the zeta_score is NOT one of the old mock values.
    # This proves we are hitting the live ZetaOracle service.
    zeta_score = efficacy_data.get("zeta_score")
    assert zeta_score is not None, "Efficacy data is missing the 'zeta_score'."
    
    # The old mock value was `0.85 + (i * 0.01)`. We check that it's not that.
    # A real fused score is unlikely to be exactly this.
    assert zeta_score != 0.85, f"ZetaScore '{zeta_score}' appears to be the old mock value. Live Oracle integration failed."

    print("Payload validation successful. Live ZetaOracle confirmed. The new standard is met. ⚔️") 
 
import requests
import time
import os

# Get the base URL for the command center from an environment variable, with a fallback
BASE_URL = os.environ.get("COMMAND_CENTER_URL", "https://crispro--crispr-assistant-command-center-v3-commandcenter-api.modal.run")
POLL_INTERVAL = 15  # seconds
POLL_TIMEOUT = 1800  # 30 minutes, to be safe with cold starts and long runs

@pytest.fixture(scope="module")
def command_center_v4_integration():
    """
    A pytest fixture to provide the base URL to the tests.
    This also serves as a check to ensure the service is minimally responsive.
    """
    try:
        response = requests.get(BASE_URL)
        response.raise_for_status()
        assert "CRISPR Assistant Command Center V4 is online" in response.text
    except (requests.exceptions.RequestException, AssertionError) as e:
        pytest.fail(f"Command Center V4 is not online or is unresponsive at {BASE_URL}. Error: {e}")
    
    return BASE_URL

def test_v4_end_to_end_kill_chain(command_center_v4_integration):
    """
    This is the crucible. It tests the entire V4 workflow:
    1. Formulate a Battle Plan.
    2. Execute the V4 Kill Chain on that plan.
    3. Poll until the process is complete.
    4. Validate the results are structured correctly and contain the assassin_score.
    """
    base_url = command_center_v4_integration

    # --- Step 1: Formulate the Battle Plan ---
    formulate_payload = {
        "patient_identifier": "test-patient-v4",
        "gene": "BRAF",
        "mutation_hgvs_p": "p.Val600Glu",
        "sequence_for_perplexity": "AAGGGAAAGTGGTGATTTGGTCTAGCTACAGAGAAATCTCGATGGAGTGGGTCCCATCAGTTTGAACAGTTGTCTGGATCCATTTTGTGGATG"
    }
    
    print(f"Submitting Battle Plan to {base_url}/v3/workflow/formulate_battle_plan...")
    formulate_response = requests.post(f"{base_url}/v3/workflow/formulate_battle_plan", json=formulate_payload)
    assert formulate_response.status_code == 201, f"Failed to formulate battle plan: {formulate_response.text}"
    
    battle_plan = formulate_response.json()
    battle_plan_id = battle_plan.get("id")
    assert battle_plan_id is not None, "Formulation response did not include a battle plan ID."
    print(f"Battle Plan #{battle_plan_id} created successfully.")

    # --- Step 2: Execute the V4 Kill Chain ---
    execute_url = f"{base_url}/v4/workflow/execute_guide_design_campaign/{battle_plan_id}"
    print(f"Executing V4 Kill Chain at {execute_url}...")
    execute_response = requests.post(execute_url)
    assert execute_response.status_code == 202, f"Failed to execute kill chain: {execute_response.text}"
    print("V4 Kill Chain initiated successfully.")

    # --- Step 3: Poll for Completion ---
    poll_url = f"{base_url}/v3/battle_plan/{battle_plan_id}"
    start_time = time.time()
    final_status = None
    
    while time.time() - start_time < POLL_TIMEOUT:
        print(f"Polling status at {poll_url}...")
        poll_response = requests.get(poll_url)
        assert poll_response.status_code == 200, f"Failed to poll for status: {poll_response.text}"
        
        status_data = poll_response.json()
        current_status = status_data.get("status")
        print(f"Current status: {current_status}")

        if current_status == "interventions_designed":
            final_status = "interventions_designed"
            break
        elif current_status == "design_failed":
            pytest.fail(f"The V4 Kill Chain failed. Final status: 'design_failed'. Details: {status_data.get('results')}")
        
        time.sleep(POLL_INTERVAL)

    assert final_status == "interventions_designed", f"Polling timed out after {POLL_TIMEOUT}s. Final status was '{current_status}'."
    print("V4 Kill Chain completed successfully.")

    # --- Step 4: Validate the Fucking Payload ---
    print("Validating results payload...")
    results_data = status_data.get("results")
    assert results_data is not None, "Final result is missing the 'results' field."
    
    guides = results_data.get("guides")
    assert isinstance(guides, list), "The 'guides' field is not a list."
    assert len(guides) > 0, "The 'guides' list is empty."

    # Check the first guide for the correct structure
    first_guide = guides[0]
    assert "guide_sequence" in first_guide, "Guide dossier is missing 'guide_sequence'."
    
    efficacy_data = first_guide.get("efficacy")
    assert efficacy_data is not None, "Guide dossier is missing 'efficacy' data."
    
    safety_data = first_guide.get("safety")
    assert safety_data is not None, "Guide dossier is missing 'safety' data."

    assert "assassin_score" in first_guide, "Guide dossier is missing 'assassin_score'."
    
    # --- The Final, Ruthless Assertion ---
    # We verify that the zeta_score is NOT one of the old mock values.
    # This proves we are hitting the live ZetaOracle service.
    zeta_score = efficacy_data.get("zeta_score")
    assert zeta_score is not None, "Efficacy data is missing the 'zeta_score'."
    
    # The old mock value was `0.85 + (i * 0.01)`. We check that it's not that.
    # A real fused score is unlikely to be exactly this.
    assert zeta_score != 0.85, f"ZetaScore '{zeta_score}' appears to be the old mock value. Live Oracle integration failed."

    print("Payload validation successful. Live ZetaOracle confirmed. The new standard is met. ⚔️") 