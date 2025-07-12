import os
import requests
import json
from pathlib import Path
from dotenv import load_dotenv
import argparse

# --- Configuration ---
load_dotenv()
COMMAND_CENTER_URL = os.getenv("COMMAND_CENTER_URL")
if not COMMAND_CENTER_URL:
    raise ValueError("COMMAND_CENTER_URL not set in .env file.")

WORKFLOW_ENDPOINT = f"{COMMAND_CENTER_URL}/workflow/run_gene_essentiality_prediction"
OUTPUT_DIR = Path("results/tissue_profiles")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

class StressSimulator:
    """
    Generates a 'stressed' gene essentiality profile for a tissue, simulating the
    impact of a specific pathogenic mutation (the 'seed').
    """
    def __init__(self):
        if not COMMAND_CENTER_URL:
            raise ValueError("COMMAND_CENTER_URL not found. Please set it in your .env file.")
        print(f"🛰️  Stress Simulator initialized. Targeting CommandCenter at: {COMMAND_CENTER_URL}")

    def generate_stressed_profile(self, tissue_name: str, mutation_id: str):
        """
        Sends a request to the CommandCenter to generate a stressed essentiality profile.
        """
        print(f"🔬 Generating STRESSED essentiality profile for '{tissue_name}' with mutation '{mutation_id}'...")

        payload = {
            "tissue_name": tissue_name,
            "stress_factor_mutation": mutation_id,
            "description": f"Stressed gene essentiality profile for {tissue_name} under the influence of {mutation_id}."
        }
        
        headers = {"Content-Type": "application/json"}

        try:
            print(f"   - Sending request to: {WORKFLOW_ENDPOINT}")
            response = requests.post(WORKFLOW_ENDPOINT, headers=headers, json=payload, timeout=1200)
            response.raise_for_status()
            
            profile_data = response.json()
            
            output_filename = f"{tissue_name.lower()}_stressed_{mutation_id.replace(':', '_')}_profile.json"
            output_path = OUTPUT_DIR / output_filename
            
            with open(output_path, "w") as f:
                json.dump(profile_data, f, indent=4)
                
            print(f"   - ✅ Success! Stressed profile saved to: {output_path}")
            print(f"   - Task ID: {profile_data.get('task_id')}")

        except requests.exceptions.HTTPError as http_err:
            print(f"   - ❌ HTTP error occurred: {http_err}")
            print(f"   - Response Body: {response.text}")
        except requests.exceptions.RequestException as req_err:
            print(f"   - ❌ Request error occurred: {req_err}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a stressed tissue profile by simulating a mutation's impact.")
    parser.add_argument("tissue", help="The target tissue to analyze (e.g., Spleen).")
    parser.add_argument("mutation", help="The mutation ID to use as the stress factor (e.g., ASXL1_p.Gly646fs).")
    args = parser.parse_args()
    
    simulator = StressSimulator()
    simulator.generate_stressed_profile(args.tissue, args.mutation) 