import requests
import json
import time

# --- Configuration ---
COMMAND_CENTER_URL = "https://crispro--command-center-commandcenter-api.modal.run/workflow/second_hit_simulation"
TARGET_GENES = {
    "ASXL1": "ATGGATTACAGAGGGCCGCGCGGAGGAGGAAGAGGAAGCGGAGCCTCGTCGGCCGAGCAGGAGGAAGAGGAAGAGGAGGAAGAAGAGCGGCCTCGCAGCGGCGGAGCAGGAAGAAGAAGAGGAAGAAGAGGAAGAGGAAGAGCGGAGCCGCGGCGGAGCAGGAAGAAGAAGAGGAAGAAGAGGAAGAGGAAGAGCGGAGCCG",
    "TET2": "ATGCAGCCGGAGGAGCCGCCGCCCCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGC",
    "DNMT3A": "ATGGACGAGCCGGGCGCGGAGGAGGAGGAAGAGGAAGCGGAGCCTCGTCGGCCGAGCAGGAGGAAGAGGAAGAGGAGGAAGAAGAGCGGCCTCGCAGCGGCGGAGCAGGAAGAAGAAGAGGAAGAAGAGGAAGAGGAAGAGCGGAGCCGCGGCGGAGCAGGAAGAAGAAGAGGAAGAAGAGGAAGAGGAAGAGCGGAGCCG"
}
NUM_MUTATIONS_PER_GENE = 5 # Reduced for initial testing
OUTPUT_FILE = "probabilistic_map.json"

# --- Analysis Functions ---
def run_simulation_for_gene(gene_name, wild_type_dna, num_mutations):
    """Calls the Command Center to run the second-hit simulation for a single gene."""
    print(f"--- 💥 Simulating evolution for {gene_name} ({num_mutations} mutations) 💥 ---")
    payload = {
        "wild_type_dna": wild_type_dna,
        "num_mutations": num_mutations
    }
    
    try:
        response = requests.post(COMMAND_CENTER_URL, json=payload, timeout=1800) # Long timeout for large simulations
        response.raise_for_status()
        print(f"--- ✅ Simulation for {gene_name} complete. ---")
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"--- ❌ ERROR: Simulation failed for {gene_name} ---")
        print(f"  Error: {e}")
        return None

def main():
    """Main orchestration script."""
    start_time = time.time()
    full_map = {}

    print("--- 🔥 BEGINNING CLONAL EVOLUTION ANALYSIS 🔥 ---")
    
    for gene_name, wt_dna in TARGET_GENES.items():
        results = run_simulation_for_gene(gene_name, wt_dna, NUM_MUTATIONS_PER_GENE)
        if results:
            full_map[gene_name] = results
            
    print("--- 📊 Aggregating and saving probabilistic map... ---")
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(full_map, f, indent=2)
        
    end_time = time.time()
    print(f"--- 🎉 ANALYSIS COMPLETE! ---")
    print(f"  - Total time: {end_time - start_time:.2f} seconds")
    print(f"  - Probabilistic map saved to: {OUTPUT_FILE}")
    print("---------------------------------")

if __name__ == "__main__":
    main()
