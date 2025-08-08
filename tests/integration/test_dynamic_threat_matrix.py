import requests
import json
from loguru import logger

# This should be the URL of your successfully deployed command center service
COMMAND_CENTER_URL = "https://crispro--crispr-assistant-command-center-v3-commandcente-70576f.modal.run"

def run_integration_test(gene_symbol: str):
    """
    Calls the live command center to test the full hunter-analyst workflow.
    """
    logger.info(f"--- 🎯 DYNAMIC THREAT MATRIX INTEGRATION TEST 🎯 ---")
    logger.info(f"Connecting to Command Center at: {COMMAND_CENTER_URL}")
    logger.info(f"Deploying Hunter-Analyst against target: {gene_symbol}")

    endpoint_url = f"{COMMAND_CENTER_URL}/workflow/run_hunter_analyst"
    payload = {"gene_symbol": gene_symbol}

    try:
        response = requests.post(endpoint_url, json=payload, timeout=300)
        response.raise_for_status()
        
        data = response.json()

        logger.success("✅ Integration Test Successful! Received response from Command Center.")
        
        print("\n" + "="*50)
        print("          COMMANDER'S STRATEGIC BRIEFING")
        print("="*50)
        print(data.get("briefing", "No briefing content received."))
        print("-"*50)

        primary_target = data.get("primary_target")
        if primary_target:
            logger.info("Primary Target Details:")
            print(json.dumps(primary_target, indent=2))
        else:
            logger.warning("No primary target was identified in the response.")

        print("="*50 + "\n")

    except requests.exceptions.HTTPError as http_err:
        logger.error(f"HTTP error occurred: {http_err}")
        logger.error(f"Response Body: {http_err.response.text}")
    except Exception as err:
        logger.error(f"An error occurred: {err}")

if __name__ == "__main__":
    # We test against VEGFA, which previously failed with the hard-coded intel.
    # A successful run here proves the dynamic system is working.
    run_integration_test("VEGFA") 