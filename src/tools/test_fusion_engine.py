import requests
import json
from Bio import SeqIO

def get_protein_sequence(fasta_file):
    """
    Retrieves the protein sequence from a FASTA file.
    """
    for record in SeqIO.parse(fasta_file, "fasta"):
        return str(record.seq)
    return None

def test_fusion_engine_braf():
    """
    Tests the fusion-engine with a known pathogenic BRAF variant.
    """
    # --- Configuration ---
    base_url = "https://crispro--fusion-engine-v1-fusionengine-api.modal.run"
    endpoint = "/score_variants"
    url = f"{base_url}{endpoint}"
    
    fasta_file = "data/reference/BRAF.fasta"
    
    # --- Get BRAF Protein Sequence ---
    import os
    if not os.path.exists(fasta_file):
        print("BRAF sequence not found, downloading...")
        os.makedirs(os.path.dirname(fasta_file), exist_ok=True)
        # Using a known UniProt URL for BRAF
        import urllib.request
        uniprot_url = "https://www.uniprot.org/uniprot/P15056.fasta"
        urllib.request.urlretrieve(uniprot_url, fasta_file)
        print("Download complete.")

    protein_sequence = get_protein_sequence(fasta_file)
    if not protein_sequence:
        print("Error: Could not retrieve BRAF protein sequence.")
        return

    # --- Construct Payload ---
    variant_data = {
        "variants": [
            {
                "variant_id": "BRAF_test_variant_1",
                "hgvs": "c.XXX", # Placeholder
                "alphamissense_variant_str": "chr7:140734637:G:C"
            }
        ],
        "protein_sequence": protein_sequence
    }

    # --- Send Request ---
    headers = {'Content-Type': 'application/json'}
    print("Sending request to fusion-engine...")
    response = requests.post(url, headers=headers, data=json.dumps(variant_data), verify=False)

    # --- Print Response ---
    if response.status_code == 200:
        print("Request successful!")
        print(json.dumps(response.json(), indent=2))
    else:
        print(f"Request failed with status code: {response.status_code}")
        print("Response:")
        print(response.text)

if __name__ == "__main__":
    test_fusion_engine_braf() 