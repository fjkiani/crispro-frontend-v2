import os
import sys
from dotenv import load_dotenv

# --- Load Environment Variables ---
# This ensures that the most recent .env file is loaded
# especially in environments where the shell might cache variables.
load_dotenv()

# Ensure the script can find the 'tools' directory
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from tools.llm_api import query_llm

def test_generative_endpoint():
    """
    Tests the generative Evo2 endpoint via the llm_api module.
    """
    print("--- Testing Evo2 Generative Endpoint ---")
    
    # A sample 20-nucleotide sequence, similar to a real guide RNA
    sample_sequence = "AATTCTTACCATCCACAAAA"
    prompt = f"Complete the following DNA sequence: {sample_sequence}"
    
    print(f"Sending prompt to local_llm provider: '{prompt}'")
    
    # --- Explicitly pass the endpoint URL ---
    # This bypasses any environment caching issues in subprocesses.
    endpoint = os.environ.get("EVO2_GENERATIVE_ENDPOINT")
    if not endpoint:
        print("\n⚠️ Test Failed: EVO2_GENERATIVE_ENDPOINT is not set in your .env file.")
        return

    # Call the query_llm function, specifically targeting the local provider
    generated_sequence = query_llm(prompt, provider="local_llm", endpoint_url=endpoint)
    
    print("\n--- Response from Endpoint ---")
    print(generated_sequence)
    print("---------------------------------")
    
    if "Error" in generated_sequence or "not configured" in generated_sequence:
        print("\n⚠️ Test Failed. Please check the error message above.")
        print("Ensure your .env file has the EVO2_GENERATIVE_ENDPOINT set correctly.")
    else:
        print("\n✅ Test finished. Review the output to ensure it is a valid DNA sequence.")

if __name__ == "__main__":
    test_generative_endpoint()
