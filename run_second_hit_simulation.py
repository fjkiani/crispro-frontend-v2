import os
import requests
import json
from dotenv import load_dotenv
from Bio import SeqIO

# --- Configuration ---
# Load environment variables from .env file
load_dotenv()

# Command Center Service URL (fallback to a default if not set)
COMMAND_CENTER_URL = os.getenv("COMMAND_CENTER_URL", "https://crispro--command-center-commandcenter-api.modal.run")
ENDPOINT_URL = f"{COMMAND_CENTER_URL}/workflow/second_hit_simulation"

# Target Gene Configuration
TARGET_GENE_NAME = "ASXL1"
TARGET_CHROM = "chr20"
TARGET_START = 30943441
TARGET_END = 31024429
REFERENCE_GENOME_PATH = "data/reference/hg19.fa"
NUM_SIMULATIONS = 50 # The number of 'second hit' mutations to simulate

def get_reference_sequence(fasta_file, chrom, start, end):
    """
    Extracts a specific sequence from a FASTA file.
    NOTE: This assumes a FASTA index file (.fai) exists.
    """
    print(f"🧬 Extracting reference sequence for {TARGET_GENE_NAME} from {fasta_file}...")
    try:
        with open(fasta_file, "r") as handle:
            fasta = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
            if chrom not in fasta:
                print(f"  ❌ Error: Chromosome '{chrom}' not found in the reference genome.")
                return None
            
            # SeqIO is 1-based, but slicing is 0-based.
            # VCF/Browser coordinates are 1-based.
            sequence = str(fasta[chrom].seq[start-1:end])
            print(f"  ✅ Successfully extracted {len(sequence)}bp sequence for {TARGET_GENE_NAME}.")
            return sequence
    except FileNotFoundError:
        print(f"  ❌ Error: Reference genome file not found at '{fasta_file}'.")
        return None
    except Exception as e:
        print(f"  ❌ An unexpected error occurred while reading the fasta file: {e}")
        return None

def run_simulation(sequence: str, num_mutations: int):
    """
-    Sends the gene sequence to the CommandCenter to run the simulation.
-    """
    print(f"\\n--- 🔥 Initiating Second Hit Simulation for {TARGET_GENE_NAME} 🔥 ---")
    print(f"  - Simulating {num_mutations} potential evolutionary pathways...")
    
    headers = {'Content-Type': 'application/json'}
    payload = {
        "wild_type_dna": sequence,
        "num_mutations": num_mutations
    }
    
    try:
        response = requests.post(ENDPOINT_URL, headers=headers, json=payload, timeout=900) # 15 min timeout
        response.raise_for_status()  # Raises an HTTPError for bad responses (4xx or 5xx)
        print("  - ✅ Simulation request successful. Awaiting results from the Zeta Oracle...")
        return response.json()
    except requests.exceptions.HTTPError as http_err:
        print(f"  ❌ HTTP Error: {http_err}")
        print(f"  - Response Body: {response.text}")
    except requests.exceptions.RequestException as req_err:
        print(f"  ❌ Request Error: {req_err}")
    except json.JSONDecodeError:
        print(f"  ❌ Error: Failed to decode JSON response. Raw response: {response.text}")
    
    return None

def display_results(results: list):
    """
    Displays the simulation results in a structured table.
    """
    if not results:
        print("\\n--- ⚠️ Simulation returned no results. ---")
        return

    print("\\n--- 🗺️ Probabilistic Map of Malignant Evolution: Top Second Hits for ASXL1 ---")
    print("-" * 80)
    print(f"{'Rank':<5} {'Zeta Score':<15} {'Interpretation':<25} {'Mutated Sequence Snippet':<35}")
    print("-" * 80)

    for i, result in enumerate(results, 1):
        score = result.get('zeta_score', 'N/A')
        interpretation = result.get('interpretation', 'N/A')
        sequence = result.get('mutated_sequence', '')
        
        # Format score to 4 decimal places if it's a number
        score_str = f"{score:.4f}" if isinstance(score, (int, float)) else str(score)

        print(f"{i:<5} {score_str:<15} {interpretation:<25} {sequence[:30]}...")

    print("-" * 80)
    print("ANALYSIS: This table ranks potential 'second hit' mutations in ASXL1 by their predicted pathogenicity.")
    print("The mutations at the top represent the most likely and dangerous evolutionary next steps for this malignancy.")


if __name__ == "__main__":
    print("Starting Second Hit Simulation script...")
    
    # 1. Get the wild-type sequence for our target gene
    wt_sequence = get_reference_sequence(REFERENCE_GENOME_PATH, TARGET_CHROM, TARGET_START, TARGET_END)
    
    if wt_sequence:
        # 2. Run the simulation via the Command Center API
        simulation_results = run_simulation(wt_sequence, NUM_SIMULATIONS)
        
        # 3. Display the results
        if simulation_results:
            display_results(simulation_results)

    print("\\n--- 📜 Script Finished ---") 