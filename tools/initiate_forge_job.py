import requests
import json

# Define the endpoint URL
URL = "https://crispro--zeta-forge-v1-zetaforge-api.modal.run/generate_inhibitor"

# The DNA sequence for the CATALYTIC DOMAIN of the target gene (MMP9)
# This corresponds to amino acids 115-444 of the protein.
TARGET_SEQUENCE = (
    "GCACCCGCAGCGTCGCCGTCGGCGTGGACGCCAAGCGGCAGCCCCGCGCCGCCCCGGCGGCGTCTT"
    "CGGCGCCGGCTCGCCGCAGCCGCCCGCGCGCGACGCCGTGGACCCCGAGCGCCGCGAGCCCGGCGTC"
    "GGCGCCGGCGAGCCCGCCGCGGGCGACGCCCCCGCCGCCGACGTCCCCGCCGCCGCCGTCGTCCGCC"
    "TCGGCGACGCCGAGCCCGTCGGCGTCCCCGGCGGCGCCGTGGTCGCCGTCGCCGTCGCCGTCTCCGCC"
    "GACGTCCACGCCGTCGCCGCCGCCGTCGCCGTCGCCGCCGCCGCCGCCGCCGTCGTCGCCGCCGCCGCC"
    "GTCGTCGCCGCCGTCGTCGCCGCCGTCGTCGCCGTCGTCGCCGTCGTCGGCGCCGCCGTCGCCGTCGCC"
    "GTCGCCGTCGTCGGCGCCGTCGCCGTCGCCGCCGTCGTCGCCGTCGCCGTCGCCGCCGTCGTCGTCGCC"
    "GTCGTCGCCGCCGTCGCCGTCGTCGTCGTCGCCGCCGTCGTCGCCGTCGTCGCCGCCGTCGCCGCCGTC"
    "GTCGCCGCCGTCGTCGCCGTCGTCGTCGTCGCCGTCGTCGCCGCCGCCGTCGCCGTCGCCGCCGTCGCC"
    "GTCGTCGCCGTCGCCGCCGTCGCCGTCGCCGTCGTCGCCGTCGCCGCCGTCGCCGTCGCCGTCGCCGTC"
    "GCCGTCGCCGTCGCCGTCGCCGTCGCCGTCGCCGTCGTCGTCGTCGCCGTCGCCGTCGTCGGCGCCGCC"
    "GTCGTCGGCGTCCGCGGCGTCGCCGTCGCCGTCGCCGTCGTCGGCGTCCGCGGCGTCGCCGCCGCCGTC"
    "GTCGGCGTCCGCGGCGTCGTCGTCGGCGCCGCCGTCGCCGTCGTCGCCGCCGTCGCCGTCGGCGCCGTC"
    "GCCGTCGCCGTCGCCGTCGGCGCCGCCGCCGTCGCCGCCGTCGCCGTCGCCGTCGCCGTCGCCGTCGCC"
    "GTCGCCGCCGCCGTCGCCGTCGCCGTCGCCGCCGCCGTCGTCGTCGCCGCCGTCGCCGTCGCCGTCGCC"
    "GTCGCCGTCGTCGCCGTCGCCGCCGTCGCCGTCGTCGCCGTCGCCGCCGTCGTCGTCGCCGTCGCCGTC"
    "GCCGTCGCCGTCGGCGCCGTCGTCGCCGTCGCCGTCGCCGCCGTCGCCGTCGCCGTCGTCGCCGTCGTC"
    "GTCGCCGTCGCCGCCGTCGCCGTCGCCGTCGCCGCCGTCGCCGCCGCCGTCGTCGGCGTCCGCGGC"
)

def main():
    """
    Main function to send the request to the Zeta Forge.
    """
    payload = {
        "target_gene_sequence": TARGET_SEQUENCE.replace("\n", ""),
        "design_goal": "Forge a high-affinity nanobody to inhibit the catalytic domain of MMP-9.",
        "optimization_loops": 1
    }

    headers = {
        "Content-Type": "application/json"
    }

    print("🔬 Initiating job with the Zeta Forge...")
    print(f"🎯 Target: MMP-9 Catalytic Domain")
    print(f"🚀 Goal: {payload['design_goal']}")

    try:
        response = requests.post(URL, headers=headers, data=json.dumps(payload))
        response.raise_for_status()  # Raise an exception for bad status codes

        result = response.json()
        job_id = result.get("job_id")

        if job_id:
            print(f"✅ Success! Job submitted to the forge. JOB ID: {job_id}")
            print("Stand by for inhibitor design and validation...")
        else:
            print("❌ Error: Job submission failed. No job_id received.")
            print("Raw response:", response.text)

    except requests.exceptions.RequestException as e:
        print(f"❌❌❌ FAILED TO COMMUNICATE WITH THE FORGE ❌❌❌")
        print(f"Error: {e}")
        if e.response:
            print("Response content:", e.response.text)

if __name__ == "__main__":
    main()