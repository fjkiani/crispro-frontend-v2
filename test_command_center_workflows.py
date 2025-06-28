import requests
import json
import os
from dotenv import load_dotenv
import urllib3
import time

# --- Configuration ---
load_dotenv()

# Suppress the InsecureRequestWarning from urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# Fetches the command center URL from .env file, or uses a default.
# You can add a .env file to the project root with the line:
# COMMAND_CENTER_URL="https://your-modal-app-url.modal.run"
COMMAND_CENTER_URL = os.getenv("COMMAND_CENTER_URL")

if not COMMAND_CENTER_URL:
    print("⚠️ WARNING: COMMAND_CENTER_URL not found in .env file.")
    print("Using default URL: https://crispro--command-center-commandcenter-api.modal.run")
    COMMAND_CENTER_URL = "https://crispro--command-center-commandcenter-api.modal.run"


def run_test(name, endpoint, payload):
    """Helper function to run a test against a Command Center endpoint."""
    print(f"\n{'='*20} RUNNING TEST: {name} {'='*20}")
    url = f"{COMMAND_CENTER_URL}{endpoint}"
    print(f"  - Endpoint: {url}")
    print(f"  - Payload: {json.dumps(payload, indent=2)}")
    
    try:
        # Using a long timeout to accomodate generative models
        # verify=False is used to bypass local SSL certificate issues.
        response = requests.post(url, json=payload, timeout=1000, verify=False)
        response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)
        
        print("\n✅ SUCCESS: Request successful.")
        print("--- RESPONSE ---")
        print(json.dumps(response.json(), indent=2))
        print("----------------")
        
    except requests.exceptions.RequestException as e:
        print(f"\n❌ FAILED: {name}")
        print(f"  - Error: {e}")
        if e.response is not None:
            print(f"  - Status Code: {e.response.status_code}")
            try:
                print(f"  - Response Body: {e.response.json()}")
            except json.JSONDecodeError:
                print(f"  - Response Body: {e.response.text}")
    print(f"{'='*58}")

def main():
    """Main function to run all workflow tests."""
    print("🚀 Starting CrisPRO Command Center Integration Tests 🚀")

    # Add a delay to allow services to warm up after deployment.
    # This can help prevent race conditions where the test runs before the app is ready.
    print("--- Waiting 15 seconds for remote services to initialize... ---")
    time.sleep(15)
    
    # --- Test 1: Digital Twin Workflow ---
    run_test(
        name="Digital Twin",
        endpoint="/workflow/digital_twin",
        payload={
            "wild_type_dna": "ATGCGTAC",
            "mutated_dna": "ATGCGTTC"
        }
    )

    # --- Test 2: Second Hit Simulation ---
    run_test(
        name="Second Hit Simulation",
        endpoint="/workflow/second_hit_simulation",
        payload={
            "wild_type_dna": "GATTACAGATTACAGATTACAGATTACAGATTACA", # A longer sequence for generation
            "num_mutations": 2
        }
    )

    # --- Test 3: Gene Essentiality Prediction ---
    run_test(
        name="Gene Essentiality",
        endpoint="/workflow/predict_gene_essentiality",
        payload={
            "genomic_context": "AAABBBCCCDDDEEE",
            "target_gene": "CCC"
        }
    )
    
    # --- Test 4: Germline Correction Blueprint ---
    # NOTE: Using a longer context to test the improved gRNA generation
    run_test(
        name="Germline Correction",
        endpoint="/workflow/design_germline_correction",
        payload={
            "corrected_sequence": "TTT",
            "target_site_context": "AAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAA",
            "arm_length": 50
        }
    )

    # --- Test 5: Clone Assassination ---
    # This gene is essential and should trigger gRNA design
    run_test(
        name="Clone Assassination (Essential Target)",
        endpoint="/workflow/clone_assassination",
        payload={
            "genomic_context": "AAAGGG" + "T"*50 + "CCCTTT",
            "target_gene": "T"*50 # Deleting this long string will have a high impact
        }
    )

    # --- Test 6: Nanobody Inhibitor Design ---
    run_test(
        name="Nanobody Inhibitor Design",
        endpoint="/workflow/design_nanobody",
        payload={
            "target_protein_name": "RUNX1",
            "target_domain": "RUNT"
        }
    )

    # --- Test 7: Drug Repurposing Screen (Placeholder) ---
    run_test(
        name="Drug Repurposing Screen",
        endpoint="/workflow/run_drug_repurposing",
        payload={
            "mutant_protein_structure": "PLACEHOLDER_PDB_FOR_ATGCGTTC"
        }
    )
    
    # --- Test 8: Off-Target Analysis ---
    # This is a live test against the BlastService using the E. coli reference genome.
    run_test(
        name="Off-Target Analysis",
        endpoint="/workflow/off_target_analysis",
        payload={
            "guide_sequence": "GATTACAGATTACAGATTAC" # A 20-bp guide, matches test 2
        }
    )

    # --- Test 9: Threat Report Generation ---
    run_test(
        name="Threat Report Generation",
        endpoint="/workflow/generate_threat_report",
        payload={
            "environmental_factor": "TNF-alpha",
            "vulnerability_threshold": 0.5,
            "normal_essentiality_results": [
                {"target_gene": "GENE_A", "essentiality_score": -0.2},
                {"target_gene": "GENE_B", "essentiality_score": -4.0},
            ],
            "stressed_essentiality_results": [
                # Under stress, GENE_A becomes much more essential (synthetic lethal)
                {"target_gene": "GENE_A", "essentiality_score": -8.5},
                # GENE_B's essentiality is unchanged
                {"target_gene": "GENE_B", "essentiality_score": -4.1},
            ]
        }
    )

    # --- Test 10: Full Patient Assessment Workflow ---
    run_test(
        name="Full Patient Assessment",
        endpoint="/workflow/full_patient_assessment",
        payload={
            "patient_id": "TEST_PATIENT_001",
            "genes_of_interest": ["RUNX1", "FLT3"],
            "germline_target": {
                "gene": "RUNX1",
                "variant": "chr21:36250941 C>T",
                "reference_sequence": "A" * 50 + "C" + "G" * 50,
                "mutated_sequence": "A" * 50 + "T" + "G" * 50
            },
            "somatic_targets": [
                {
                    "id": "FLT3-ITD-1",
                    "ref_seq": "TTT",
                    "mut_seq": "TTTTTT"
                }
            ]
        }
    )
    
    print("\n🏁 All tests complete. 🏁")


if __name__ == "__main__":
    main() 