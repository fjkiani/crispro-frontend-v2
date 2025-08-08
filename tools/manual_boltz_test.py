import os
import textwrap

def main():
    """
    Generates the cURL command for the Boltz Gauntlet Manual Validation Protocol.

    This script serves as the executable test case for the doctrine outlined in:
    .cursor/rules/boltz_gauntlet_manual_test_plan.mdc
    """
    
    # --- Test Case: BRAF-Fragment-100 ---
    protein_sequence = "MAALSGGGGGAEPGWVCLHGTLLPEPGGASKICGSGSFPHGLLVFWKVKGLLKMLADRVENTFPIPRVIRKQLHVQMLQYLEKHKIQSWIHVWMRDYHPFDRDQDHYAKGFVRYSVKNT"
    target_name = "manual_gauntlet_validation_braf_frag_100"
    boltz_service_url = "https://crispro--boltz-service-fastapi-app.modal.run"

    # Construct the JSON payload
    json_payload = f'{{ "protein_sequence": "{protein_sequence}", "target_name": "{target_name}" }}'

    # Construct the cURL command
    curl_command = f"""
curl -X POST "{boltz_service_url}" \\
-H "Content-Type: application/json" \\
-d '{json_payload}'
"""

    # --- Print the Instructions ---
    print("--- ⚔️  THE GAUNTLET: MANUAL VALIDATION PROTOCOL ⚔️  ---")
    print("\\nThis script provides the command to manually test our 'boltz-service'.")
    print("It validates that our structural prediction service is online and functional.")
    print("\\n------------------------------------------------------------------")
    print("  ✅ 1. TEST CASE DATA")
    print("------------------------------------------------------------------")
    print(f"  - Protein Sequence (Fragment):\\n    {textwrap.fill(protein_sequence, 60, initial_indent='    ', subsequent_indent='    ')}")
    print(f"  - Target Name:\\n    {target_name}")
    print(f"  - Service URL:\\n    {boltz_service_url}")
    print("\\n------------------------------------------------------------------")
    print("  ✅ 2. EXECUTE THIS COMMAND IN YOUR TERMINAL")
    print("------------------------------------------------------------------")
    print(curl_command.strip())
    print("\\n------------------------------------------------------------------")
    print("  ✅ 3. VERIFY THE RESULTS")
    print("------------------------------------------------------------------")
    print("  - The command should return a JSON object.")
    print("  - Look for the 'plddt' key.")
    print("  - A score >= 70.0 is considered a SUCCESS.")
    print("------------------------------------------------------------------\\n")

if __name__ == "__main__":
    main() 
import textwrap

def main():
    """
    Generates the cURL command for the Boltz Gauntlet Manual Validation Protocol.

    This script serves as the executable test case for the doctrine outlined in:
    .cursor/rules/boltz_gauntlet_manual_test_plan.mdc
    """
    
    # --- Test Case: BRAF-Fragment-100 ---
    protein_sequence = "MAALSGGGGGAEPGWVCLHGTLLPEPGGASKICGSGSFPHGLLVFWKVKGLLKMLADRVENTFPIPRVIRKQLHVQMLQYLEKHKIQSWIHVWMRDYHPFDRDQDHYAKGFVRYSVKNT"
    target_name = "manual_gauntlet_validation_braf_frag_100"
    boltz_service_url = "https://crispro--boltz-service-fastapi-app.modal.run"

    # Construct the JSON payload
    json_payload = f'{{ "protein_sequence": "{protein_sequence}", "target_name": "{target_name}" }}'

    # Construct the cURL command
    curl_command = f"""
curl -X POST "{boltz_service_url}" \\
-H "Content-Type: application/json" \\
-d '{json_payload}'
"""

    # --- Print the Instructions ---
    print("--- ⚔️  THE GAUNTLET: MANUAL VALIDATION PROTOCOL ⚔️  ---")
    print("\\nThis script provides the command to manually test our 'boltz-service'.")
    print("It validates that our structural prediction service is online and functional.")
    print("\\n------------------------------------------------------------------")
    print("  ✅ 1. TEST CASE DATA")
    print("------------------------------------------------------------------")
    print(f"  - Protein Sequence (Fragment):\\n    {textwrap.fill(protein_sequence, 60, initial_indent='    ', subsequent_indent='    ')}")
    print(f"  - Target Name:\\n    {target_name}")
    print(f"  - Service URL:\\n    {boltz_service_url}")
    print("\\n------------------------------------------------------------------")
    print("  ✅ 2. EXECUTE THIS COMMAND IN YOUR TERMINAL")
    print("------------------------------------------------------------------")
    print(curl_command.strip())
    print("\\n------------------------------------------------------------------")
    print("  ✅ 3. VERIFY THE RESULTS")
    print("------------------------------------------------------------------")
    print("  - The command should return a JSON object.")
    print("  - Look for the 'plddt' key.")
    print("  - A score >= 70.0 is considered a SUCCESS.")
    print("------------------------------------------------------------------\\n")

if __name__ == "__main__":
    main() 
 
 
 
 
 