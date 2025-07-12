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

WORKFLOW_ENDPOINT = f"{COMMAND_CENTER_URL}/workflow/generate_threat_report"
OUTPUT_DIR = Path("results/threat_reports")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PROFILES_DIR = Path("results/tissue_profiles")

class ThreatReportGenerator:
    """
    Generates a Threat Report by comparing a baseline and stressed tissue profile
    to identify synthetic lethalities.
    """
    def __init__(self):
        if not COMMAND_CENTER_URL:
            raise ValueError("COMMAND_CENTER_URL not found.")
        print(f"🛰️  Threat Report Generator initialized. Targeting CommandCenter at: {COMMAND_CENTER_URL}")

    def generate_report(self, baseline_profile_path: Path, stressed_profile_path: Path, vulnerability_threshold: float):
        """
        Loads profiles, sends them to the CommandCenter, and saves the resulting report.
        """
        print(f"📄 Loading baseline profile from: {baseline_profile_path}")
        with open(baseline_profile_path, 'r') as f:
            baseline_data = json.load(f)

        print(f"📄 Loading stressed profile from: {stressed_profile_path}")
        with open(stressed_profile_path, 'r') as f:
            stressed_data = json.load(f)
            
        environmental_factor = stressed_data.get("stress_factor", "Unknown Mutation")

        print(f"🛡️  Generating Threat Report for '{baseline_data.get('tissue_name')}' under stress from '{environmental_factor}'...")

        payload = {
            "normal_essentiality_results": baseline_data.get("essentiality_profile", []),
            "stressed_essentiality_results": stressed_data.get("essentiality_profile", []),
            "environmental_factor": environmental_factor,
            "vulnerability_threshold": vulnerability_threshold
        }
        
        headers = {"Content-Type": "application/json"}

        try:
            print(f"   - Sending request to: {WORKFLOW_ENDPOINT}")
            response = requests.post(WORKFLOW_ENDPOINT, headers=headers, json=payload, timeout=300)
            response.raise_for_status()
            
            report_data = response.json()
            
            output_filename = f"threat_report_{baseline_data.get('tissue_name', 'unknown_tissue')}_{environmental_factor.replace(':', '_')}.json"
            output_path = OUTPUT_DIR / output_filename
            
            with open(output_path, "w") as f:
                json.dump(report_data, f, indent=4)
                
            print(f"   - ✅ Success! Threat Report saved to: {output_path}")
            
            # Print a summary of the report to the console
            print("\n--- ⚔️ THREAT REPORT SUMMARY ⚔️ ---")
            print(f"  - Environmental Factor: {report_data.get('environmental_factor')}")
            print(f"  - Synthetic Lethalities Identified: {report_data.get('synthetic_lethalities_identified')}")
            top_hits = report_data.get('synthetic_lethalities', [])
            if top_hits:
                print("  --- Top 5 Vulnerabilities ---")
                for hit in top_hits[:5]:
                    print(f"    - Gene: {hit['gene']}, Vulnerability Increase: {hit['vulnerability_increase']}")
            print("-----------------------------------\n")


        except requests.exceptions.HTTPError as http_err:
            print(f"   - ❌ HTTP error occurred: {http_err}")
            print(f"   - Response Body: {response.text}")
        except requests.exceptions.RequestException as req_err:
            print(f"   - ❌ Request error occurred: {req_err}")

if __name__ == "__main__":
    # For this campaign, the inputs are known.
    baseline_file = PROFILES_DIR / "spleen_profile.json"
    stressed_file = PROFILES_DIR / "spleen_stressed_ASXL1_p.Gly646fs_profile.json"
    
    if not baseline_file.exists() or not stressed_file.exists():
        print(f"❌ Error: Missing profile files. Ensure both '{baseline_file}' and '{stressed_file}' exist.")
    else:
        generator = ThreatReportGenerator()
        generator.generate_report(baseline_file, stressed_file, vulnerability_threshold=0.5) 