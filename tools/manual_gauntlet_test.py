import textwrap
import json

def main():
    """
    Generates the cURL command for the Reforged Boltz Gauntlet Manual Validation Protocol.

    This script serves as the executable test case for the "Gauntlet Reforged" operation.
    """
    
    # --- Test Case: A known stable protein fragment ---
    protein_sequence = "MAALSGGGGGAEPGWVCLHGTLLPEPGGASKICGSGSFPHGLLVFWKVKGLLKMLADRVENTFPIPRVIRKQLHVQMLQYLEKHKIQSWIHVWMRDYHPFDRDQDHYAKGFVRYSVKNT"
    job_id = "manual-gauntlet-validation-1"
    
    # This URL was provided by the Commander after successful deployment.
    boltz_service_url = "https://crispro--boltz-service-fastapi-app.modal.run/v1/predict_structure"

    # Construct the JSON payload
    payload = {
        "protein_sequence": protein_sequence,
        "job_id": job_id
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
    print("  ✅ 1. THE SERVICE IS ALREADY DEPLOYED")
    print("------------------------------------------------------------------")
    print(f"  URL: {boltz_service_url}")
    print("\\n------------------------------------------------------------------")
    print("  ✅ 2. EXECUTE THIS COMMAND IN YOUR TERMINAL")
    print("------------------------------------------------------------------")
    print(curl_command.strip())
    print("\\n------------------------------------------------------------------")
    print("  ✅ 3. INTERPRET THE RESPONSE")
    print("------------------------------------------------------------------")
    print("  - The command will immediately return a JSON object with a 'job_id' and 'pending' status.")
    print("    e.g., {'job_id': 'fu-...', 'status': 'pending', ...}")
    print("  - You must then check the Modal UI logs for the 'boltz-service-gauntlet' app to see the final pLDDT score.")
    print("  - A 'plddt_score' >= 70.0 is considered a SUCCESS.")
    print("------------------------------------------------------------------\\n")

if __name__ == "__main__":
    main() 
import json

def main():
    """
    Generates the cURL command for the Reforged Boltz Gauntlet Manual Validation Protocol.

    This script serves as the executable test case for the "Gauntlet Reforged" operation.
    """
    
    # --- Test Case: A known stable protein fragment ---
    protein_sequence = "MAALSGGGGGAEPGWVCLHGTLLPEPGGASKICGSGSFPHGLLVFWKVKGLLKMLADRVENTFPIPRVIRKQLHVQMLQYLEKHKIQSWIHVWMRDYHPFDRDQDHYAKGFVRYSVKNT"
    job_id = "manual-gauntlet-validation-1"
    
    # This URL was provided by the Commander after successful deployment.
    boltz_service_url = "https://crispro--boltz-service-fastapi-app.modal.run/v1/predict_structure"

    # Construct the JSON payload
    payload = {
        "protein_sequence": protein_sequence,
        "job_id": job_id
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
    print("  ✅ 1. THE SERVICE IS ALREADY DEPLOYED")
    print("------------------------------------------------------------------")
    print(f"  URL: {boltz_service_url}")
    print("\\n------------------------------------------------------------------")
    print("  ✅ 2. EXECUTE THIS COMMAND IN YOUR TERMINAL")
    print("------------------------------------------------------------------")
    print(curl_command.strip())
    print("\\n------------------------------------------------------------------")
    print("  ✅ 3. INTERPRET THE RESPONSE")
    print("------------------------------------------------------------------")
    print("  - The command will immediately return a JSON object with a 'job_id' and 'pending' status.")
    print("    e.g., {'job_id': 'fu-...', 'status': 'pending', ...}")
    print("  - You must then check the Modal UI logs for the 'boltz-service-gauntlet' app to see the final pLDDT score.")
    print("  - A 'plddt_score' >= 70.0 is considered a SUCCESS.")
    print("------------------------------------------------------------------\\n")

if __name__ == "__main__":
    main() 
 
 
 
 
 