import os
import requests
import json
from pathlib import Path
from dotenv import load_dotenv

# --- Configuration ---
load_dotenv()
COMMAND_CENTER_URL = os.getenv("COMMAND_CENTER_URL")
if not COMMAND_CENTER_URL:
    raise ValueError("COMMAND_CENTER_URL not set in .env file.")

# Ensure the URL is for the correct workflow
WORKFLOW_ENDPOINT = f"{COMMAND_CENTER_URL}/workflow/run_gene_essentiality_prediction"
OUTPUT_DIR = Path("results/tissue_profiles")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

class TissueProfiler:
    """
    Characterizes a given tissue by generating a gene essentiality profile.
    This tool invokes the CommandCenter's gene essentiality prediction workflow,
    effectively defining the 'soil' for our Seed & Soil analysis.
    """
    def __init__(self):
        if not COMMAND_CENTER_URL:
            raise ValueError("COMMAND_CENTER_URL not found. Please set it in your .env file.")
        print(f"🛰️  Tissue Profiler initialized. Targeting CommandCenter at: {COMMAND_CENTER_URL}")

    def generate_profile(self, tissue_name: str):
        """
        Sends a request to the CommandCenter to generate an essentiality profile for the specified tissue.
        """
        print(f"🔬 Generating essentiality profile for tissue: '{tissue_name}'...")

        payload = {
            "tissue_name": tissue_name,
            "description": f"Baseline gene essentiality profile for {tissue_name}, to be used as 'soil' in Seed & Soil analysis."
        }
        
        headers = {"Content-Type": "application/json"}

        try:
            print(f"   - Sending request to: {WORKFLOW_ENDPOINT}")
            response = requests.post(WORKFLOW_ENDPOINT, headers=headers, json=payload, timeout=1200) # Long timeout for complex models
            response.raise_for_status()
            
            profile_data = response.json()
            
            # Save the results to a file
            output_filename = f"{tissue_name.lower().replace(' ', '_')}_profile.json"
            output_path = OUTPUT_DIR / output_filename
            
            with open(output_path, "w") as f:
                json.dump(profile_data, f, indent=4)
                
            print(f"   - ✅ Success! Profile saved to: {output_path}")
            print(f"   - Task ID: {profile_data.get('task_id')}")

        except requests.exceptions.HTTPError as http_err:
            print(f"   - ❌ HTTP error occurred: {http_err}")
            print(f"   - Response Body: {response.text}")
        except requests.exceptions.RequestException as req_err:
            print(f"   - ❌ Request error occurred: {req_err}")
        except json.JSONDecodeError:
            print(f"   - ❌ Failed to decode JSON from response. Raw response:")
            print(response.text)

if __name__ == "__main__":
    profiler = TissueProfiler()
    
    # Based on our literature analysis, the spleen is the first high-priority soil.
    target_tissue = "Spleen"
    
    profiler.generate_profile(target_tissue) 