import requests
import json
import os

# --- Configuration ---
# Get the CommandCenter URL from an environment variable if it exists, otherwise use the last known good URL.
COMMAND_CENTER_URL = os.getenv("COMMAND_CENTER_ENDPOINT", "https://crispro--command-center-v2-commandcenter-api.modal.run")
ASSESS_THREAT_ENDPOINT = f"{COMMAND_CENTER_URL}/workflow/assess_threat"
HEADERS = {"Content-Type": "application/json"}

def run_test(gene_symbol: str, protein_change: str):
    """Helper function to run a single threat assessment test and print the results."""
    print(f"--- 🧪 TESTING: {gene_symbol} {protein_change} 🧪 ---")
    
    payload = {
        "gene_symbol": gene_symbol,
        "protein_change": protein_change
    }
    
    try:
        response = requests.post(ASSESS_THREAT_ENDPOINT, json=payload, headers=HEADERS)
        response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)
        
        response_data = response.json()
        print("✅ SUCCESS: Received response from CommandCenter.")
        print(json.dumps(response_data, indent=2))
        
        # --- Automated Verification ---
        assessment = response_data.get("assessment", {})
        source = assessment.get("assessment_source")
        
        is_truncating = 'fs' in protein_change or '*' in protein_change or 'Ter' in protein_change
        
        if is_truncating:
            if source == "Truncation Sieve":
                print("\\n✅ VERIFICATION PASSED: Correctly identified by Truncation Sieve.")
            else:
                print(f"\\n❌ VERIFICATION FAILED: Expected 'Truncation Sieve', but got '{source}'.")
        else:
            if source == "Adjudicator":
                print("\\n✅ VERIFICATION PASSED: Correctly routed to the Adjudicator.")
            else:
                print(f"\\n❌ VERIFICATION FAILED: Expected 'Adjudicator', but got '{source}'.")
                
    except requests.exceptions.RequestException as e:
        print(f"❌ FAILED: API request error: {e}")
        if e.response:
            print(f"Response Body: {e.response.text}")

if __name__ == "__main__":
    print("Beginning Triumvirate Protocol Validation...")

    # Test Case 1: Known Pathogenic Missense Mutation
    # This tests the full new pipeline: CommandCenter -> Zeta Oracle -> Adjudicator
    run_test(gene_symbol="BRAF", protein_change="p.V600E")
    
    print("\\n" + "="*60 + "\\n")
    
    # Test Case 2: Known Pathogenic Frameshift Mutation
    # This tests the fast-path "Truncation Sieve" branch of our Triumvirate Protocol.
    run_test(gene_symbol="ASXL1", protein_change="p.Gly646fs")

    print("\\nValidation Script Complete.") 