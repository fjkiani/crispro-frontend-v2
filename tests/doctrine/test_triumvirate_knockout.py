import requests
import json
import os

# --- Configuration ---
# NOTE: This test requires the NEW, DYNAMIC command center to be deployed.
# This URL has been confirmed from the latest successful deployment.
COMMAND_CENTER_URL = "https://crispro--crispr-assistant-command-center-v3-commandcente-70576f.modal.run"
ASSESS_THREAT_ENDPOINT = f"{COMMAND_CENTER_URL}/workflow/assess_threat"
HEADERS = {"Content-Type": "application/json"}

TARGET_GENE = "VEGFA"
# We are simulating a frameshift mutation at the second amino acid.
# This is a catastrophic, truncating event.
FRAMESHIFT_MUTATION = "p.Ala2fs"

def validate_truncation_sieve():
    """
    This test validates the Triumvirate Protocol's first line of defense.
    It sends a known catastrophic mutation to the CommandCenter and verifies
    that the "Truncation Sieve" correctly identifies it, bypassing the flawed Oracle.
    """
    print(f"--- 💥 TRIUMVIRATE DOCTRINE VALIDATION: THE TRUNCATION SIEVE 💥 ---")
    print(f"  TARGET: {TARGET_GENE}")
    print(f"  WEAPON: Catastrophic Frameshift ({FRAMESHIFT_MUTATION})")
    print(f"  COMMAND CENTER: {COMMAND_CENTER_URL}")
    
    payload = {
        "gene_symbol": TARGET_GENE,
        "protein_change": FRAMESHIFT_MUTATION
    }
    
    try:
        response = requests.post(ASSESS_THREAT_ENDPOINT, json=payload, headers=HEADERS, timeout=300)
        response.raise_for_status()
        
        response_data = response.json()
        print("\n--- ✅ SUCCESS: Received Response from CommandCenter ---")
        print(json.dumps(response_data, indent=2))
        
        # --- Automated Verification ---
        assessment = response_data.get("assessment", {})
        source = assessment.get("assessment_source")
        is_pathogenic = assessment.get("is_pathogenic")

        print("\n--- 🧐 Automated Verification ---")
        if source == "Truncation Sieve" and is_pathogenic is True:
            print("  [✅] VERIFICATION PASSED: Correctly identified as Pathogenic by the Truncation Sieve.")
        else:
            print(f"  [❌] VERIFICATION FAILED:")
            print(f"     -> Expected assessment source: 'Truncation Sieve', but got '{source}'.")
            print(f"     -> Expected is_pathogenic: True, but got '{is_pathogenic}'.")
            
    except requests.exceptions.RequestException as e:
        print(f"\n--- ❌ FAILED: API Request Error ---")
        print(f"  ERROR: {e}")
        if e.response:
            print(f"  RESPONSE BODY: {e.response.text}")

if __name__ == "__main__":
    validate_truncation_sieve() 