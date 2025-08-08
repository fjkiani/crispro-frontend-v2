import os
import textwrap
import json

def main():
    """
    Generates the cURL command for the REFORGED Boltz Gauntlet Manual Validation Protocol.

    This script serves as the executable test case for the doctrine outlined in:
    .cursor/rules -> OPERATION: GAUNTLET REFORGED
    """
    
    # --- Test Case: BRAF-Fragment-100 ---
    protein_sequence = "MAALSGGGGGAEPGWVCLHGTLLPEPGGASKICGSGSFPHGLLVFWKVKGLLKMLADRVENTFPIPRVIRKQLHVQMLQYLEKHKIQSWIHVWMRDYHPFDRDQDHYAKGFVRYSVKNT"
    target_name = "manual_gauntlet_reforged_validation"
    # NOTE: This URL must point to the NEWLY DEPLOYED reforged service endpoint
    boltz_service_url = "https://crispro--boltz-service-reforged-api.modal.run/v1/predict_structure"

    # Construct the JSON payload
    payload = {
        "protein_sequence": protein_sequence,
        "target_name": target_name
    }
    json_payload_str = json.dumps(payload)

    # Construct the cURL command
    curl_command = f"""
curl -X POST "{boltz_service_url}" \\
-H "Content-Type: application/json" \\
-d '{json_payload_str}'
"""

    # --- Print the Instructions ---
    print("--- ⚔️  GAUNTLET REFORGED: MANUAL VALIDATION ⚔️  ---")
    print("\\nThis script provides the command to manually test our REFORGED 'boltz-service'.")
    print("It validates that our new '/v1/predict_structure' endpoint is functional.")
    print("\\n------------------------------------------------------------------")
    print("  ✅ 1. DEPLOY THE REFORGED SERVICE")
    print("------------------------------------------------------------------")
    print("  First, ensure you have deployed the modified 'src/services/boltz/main.py'")
    print("  > modal deploy src/services/boltz/main.py")
    print("\\n------------------------------------------------------------------")
    print("  ✅ 2. EXECUTE THIS COMMAND IN YOUR TERMINAL")
    print("------------------------------------------------------------------")
    print(curl_command.strip())
    print("\\n------------------------------------------------------------------")
    print("  ✅ 3. VERIFY THE RESULTS")
    print("------------------------------------------------------------------")
    print("  - The command should return a JSON object, e.g., {'status': 'complete', 'plddt': 85.4}")
    print("  - A 'plddt' score >= 70.0 is considered a SUCCESS.")
    print("------------------------------------------------------------------\\n")

if __name__ == "__main__":
    main() 
import textwrap
import json

def main():
    """
    Generates the cURL command for the REFORGED Boltz Gauntlet Manual Validation Protocol.

    This script serves as the executable test case for the doctrine outlined in:
    .cursor/rules -> OPERATION: GAUNTLET REFORGED
    """
    
    # --- Test Case: BRAF-Fragment-100 ---
    protein_sequence = "MAALSGGGGGAEPGWVCLHGTLLPEPGGASKICGSGSFPHGLLVFWKVKGLLKMLADRVENTFPIPRVIRKQLHVQMLQYLEKHKIQSWIHVWMRDYHPFDRDQDHYAKGFVRYSVKNT"
    target_name = "manual_gauntlet_reforged_validation"
    # NOTE: This URL must point to the NEWLY DEPLOYED reforged service endpoint
    boltz_service_url = "https://crispro--boltz-service-reforged-api.modal.run/v1/predict_structure"

    # Construct the JSON payload
    payload = {
        "protein_sequence": protein_sequence,
        "target_name": target_name
    }
    json_payload_str = json.dumps(payload)

    # Construct the cURL command
    curl_command = f"""
curl -X POST "{boltz_service_url}" \\
-H "Content-Type: application/json" \\
-d '{json_payload_str}'
"""

    # --- Print the Instructions ---
    print("--- ⚔️  GAUNTLET REFORGED: MANUAL VALIDATION ⚔️  ---")
    print("\\nThis script provides the command to manually test our REFORGED 'boltz-service'.")
    print("It validates that our new '/v1/predict_structure' endpoint is functional.")
    print("\\n------------------------------------------------------------------")
    print("  ✅ 1. DEPLOY THE REFORGED SERVICE")
    print("------------------------------------------------------------------")
    print("  First, ensure you have deployed the modified 'src/services/boltz/main.py'")
    print("  > modal deploy src/services/boltz/main.py")
    print("\\n------------------------------------------------------------------")
    print("  ✅ 2. EXECUTE THIS COMMAND IN YOUR TERMINAL")
    print("------------------------------------------------------------------")
    print(curl_command.strip())
    print("\\n------------------------------------------------------------------")
    print("  ✅ 3. VERIFY THE RESULTS")
    print("------------------------------------------------------------------")
    print("  - The command should return a JSON object, e.g., {'status': 'complete', 'plddt': 85.4}")
    print("  - A 'plddt' score >= 70.0 is considered a SUCCESS.")
    print("------------------------------------------------------------------\\n")

if __name__ == "__main__":
    main() 
 
 
 
 
 