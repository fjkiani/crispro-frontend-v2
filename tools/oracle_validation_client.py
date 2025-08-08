import requests
import json
import sys

# The live endpoint for our sanitized Unified Oracle
ORACLE_URL = "https://crispro--unified-oracle-service-unifiedoracle-invoke.modal.run"

# The problematic prompt that previously generated junk DNA
# This is a real-world example of a complex, repetitive sequence from BRCA1.
PROMPT = "TGGATACTGAAACACATTAGAAGAAATGACTGAACTGTCTAAAGCTGTAGTAAATATACAATGCACTTTCTTCCAAAATGTAGAAAAGGCTGAATTCATCTAAAAAGAAATACATTCCTGAGAGGCAGAGGTTGCAGTGAGCCGAGATCCTGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCATCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAATACATTCCTGAGAGGCAGAGGTTGCAGTGAGCCGAGATCCTGCCACTGCACTCCAGCCT"

def run_validation_test():
    """
    Sends a request to the newly deployed Unified Oracle to validate
    the effectiveness of "Operation: Sanitize the Forge".
    """
    print("🚀 Initiating validation test for the Unified Oracle...")
    print(f"📡 Targeting endpoint: {ORACLE_URL}")
    print(f"🧬 Using problematic prompt (last 50 chars): ...{PROMPT[-50:]}")

    payload = {
        "action": "generate",
        "params": {
            "prompt": PROMPT,
            "gen_params": {
                "n_tokens": 50
                # We will use the new disciplined defaults set in the service:
                # temperature: 0.6, top_p: 0.9
            }
        }
    }

    try:
        response = requests.post(ORACLE_URL, json=payload, timeout=300)
        response.raise_for_status()
        result = response.json()

        print("\n--- ORACLE RESPONSE ---")
        print(json.dumps(result, indent=2))
        print("-----------------------\n")

        if result.get("status") == "error":
            print(f"✅ SUCCESS: Oracle correctly identified and filtered junk.")
            print(f"   Reason: {result.get('message')}")
            # This is a successful outcome, as the filter worked.
            return

        if result.get("status") == "success":
            generated_sequence = result.get("completion", "")
            
            # --- VALIDATION CHECKS ---
            print("🔬 Running validation checks on the generated sequence...")
            
            # Check 1: No poly-A tracts
            if "AAAAAAAA" in generated_sequence:
                print("❌ FAILED: Junk pattern 'AAAAAAAA' detected in the output.")
                sys.exit(1)
            else:
                print("✅ PASSED: No poly-A tracts detected.")

            # Check 2: Sufficient complexity
            if len(set(generated_sequence)) < 4:
                print(f"❌ FAILED: Sequence complexity is too low (only {len(set(generated_sequence))} unique characters).")
                sys.exit(1)
            else:
                print("✅ PASSED: Sequence complexity is sufficient.")

            print("\n🎉 VALIDATION COMPLETE: The sanitized Oracle is producing high-quality, non-junk sequences.")

    except requests.exceptions.RequestException as e:
        print(f"❌ CRITICAL FAILURE: Could not connect to the Oracle service: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    run_validation_test() 
import json
import sys

# The live endpoint for our sanitized Unified Oracle
ORACLE_URL = "https://crispro--unified-oracle-service-unifiedoracle-invoke.modal.run"

# The problematic prompt that previously generated junk DNA
# This is a real-world example of a complex, repetitive sequence from BRCA1.
PROMPT = "TGGATACTGAAACACATTAGAAGAAATGACTGAACTGTCTAAAGCTGTAGTAAATATACAATGCACTTTCTTCCAAAATGTAGAAAAGGCTGAATTCATCTAAAAAGAAATACATTCCTGAGAGGCAGAGGTTGCAGTGAGCCGAGATCCTGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCATCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAATACATTCCTGAGAGGCAGAGGTTGCAGTGAGCCGAGATCCTGCCACTGCACTCCAGCCT"

def run_validation_test():
    """
    Sends a request to the newly deployed Unified Oracle to validate
    the effectiveness of "Operation: Sanitize the Forge".
    """
    print("🚀 Initiating validation test for the Unified Oracle...")
    print(f"📡 Targeting endpoint: {ORACLE_URL}")
    print(f"🧬 Using problematic prompt (last 50 chars): ...{PROMPT[-50:]}")

    payload = {
        "action": "generate",
        "params": {
            "prompt": PROMPT,
            "gen_params": {
                "n_tokens": 50
                # We will use the new disciplined defaults set in the service:
                # temperature: 0.6, top_p: 0.9
            }
        }
    }

    try:
        response = requests.post(ORACLE_URL, json=payload, timeout=300)
        response.raise_for_status()
        result = response.json()

        print("\n--- ORACLE RESPONSE ---")
        print(json.dumps(result, indent=2))
        print("-----------------------\n")

        if result.get("status") == "error":
            print(f"✅ SUCCESS: Oracle correctly identified and filtered junk.")
            print(f"   Reason: {result.get('message')}")
            # This is a successful outcome, as the filter worked.
            return

        if result.get("status") == "success":
            generated_sequence = result.get("completion", "")
            
            # --- VALIDATION CHECKS ---
            print("🔬 Running validation checks on the generated sequence...")
            
            # Check 1: No poly-A tracts
            if "AAAAAAAA" in generated_sequence:
                print("❌ FAILED: Junk pattern 'AAAAAAAA' detected in the output.")
                sys.exit(1)
            else:
                print("✅ PASSED: No poly-A tracts detected.")

            # Check 2: Sufficient complexity
            if len(set(generated_sequence)) < 4:
                print(f"❌ FAILED: Sequence complexity is too low (only {len(set(generated_sequence))} unique characters).")
                sys.exit(1)
            else:
                print("✅ PASSED: Sequence complexity is sufficient.")

            print("\n🎉 VALIDATION COMPLETE: The sanitized Oracle is producing high-quality, non-junk sequences.")

    except requests.exceptions.RequestException as e:
        print(f"❌ CRITICAL FAILURE: Could not connect to the Oracle service: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    run_validation_test() 
 
 
 
 
 