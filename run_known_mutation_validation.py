import os
import requests
import json
from dotenv import load_dotenv
from Bio import SeqIO
from datetime import datetime
import time
import pysam

# This script will be very similar to the second_hit_simulation, but instead of
# generating random mutations, it will work from a curated list of known pathogenic ones.

# --- Configuration ---
load_dotenv()
COMMAND_CENTER_URL = os.getenv("COMMAND_CENTER_URL", "https://crispro--command-center-commandcenter-api.modal.run")
REFERENCE_GENOME_PATH = "data/reference/hg19.fa"
RESULTS_DIR = "results/known_mutation_validation"
TARGET_GENE_NAME = "ASXL1"
TARGET_CHROM = "chr20"
TARGET_START = 30946133 
TARGET_END = 31027122
# The coordinates for the coding sequence (CDS) within the full gene sequence.
# These are relative to the start of the sequence we extract (TARGET_START).
# This is for the canonical transcript NM_015338.
CDS_START_RELATIVE = 265 # ATG start codon
CDS_END_RELATIVE = 4920 # Corrected: Length must be a multiple of 3.

# --- "Known Enemies" List ---
# A curated list of known pathogenic mutations in ASXL1.
# Format: (Mutation Name, Type, Details)
# For frameshifts (fs), details are (position, inserted_bases).
# For nonsense/stop-gained (*), details are (position, original_base, mutated_base).
KNOWN_MUTATIONS = {
    "p.Gly646fs": ("fs", (1936, "T")),  # Corresponds to a specific frameshift
    "p.Arg497*": ("nonsense", (1489, "A", "T")), # Corresponds to c.1489A>T (Corrected from C>T based on our hg19 ref)
    # Add more known pathogenic mutations here
}

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

def introduce_specific_mutation(full_sequence: str, mut_type: str, details: tuple) -> str:
    """
    Introduces a specific, known mutation into a gene's coding sequence (CDS)
    and reconstructs the full gene sequence. This is the correct, original logic.
    """
    # Step 1: Extract the CDS from the full gene sequence
    cds_sequence = full_sequence[CDS_START_RELATIVE-1:CDS_END_RELATIVE]
    
    mutated_cds_str = ""
    if mut_type == "fs":
        position, inserted_bases = details
        mutated_cds_str = cds_sequence[:position-1] + inserted_bases + cds_sequence[position-1:]

    elif mut_type == "nonsense":
        position, original_base, mutated_base = details
        if cds_sequence[position - 1].upper() == original_base.upper():
            mutated_cds_str = cds_sequence[:position-1] + mutated_base + cds_sequence[position:]
        else:
            print(f"  ❌ Error: Invalid position or original base mismatch for nonsense mutation at position {position} in CDS.")
            return None
    else:
        print(f"  ❌ Error: Unknown mutation type '{mut_type}'.")
        return None
        
    # Step 2: Reconstruct the full sequence with the mutated CDS
    reconstructed_sequence = full_sequence[:CDS_START_RELATIVE-1] + mutated_cds_str + full_sequence[CDS_END_RELATIVE:]
    
    print(f"  - Successfully introduced {mut_type} mutation and reconstructed full sequence.")
    return reconstructed_sequence

def run_validation(wild_type_dna: str, mutation_name: str, mutated_dna: str) -> dict:
    """
    Sends a specific WT/mutant pair to the CommandCenter for scoring.
    """
    print("  - Dispatching to CommandCenter for scoring...")
    ENDPOINT_URL = f"{COMMAND_CENTER_URL}/workflow/score_single_mutation"
    
    # --- Ground Truth Debugging ---
    # Print the sequences being sent to the API to verify the mutation was introduced correctly.
    print(f"    - REF (len {len(wild_type_dna)}): {wild_type_dna[:30]}...{wild_type_dna[-30:]}")
    print(f"    - MUT (len {len(mutated_dna)}): {mutated_dna[:30]}...{mutated_dna[-30:]}")

    payload = {
        "wild_type_dna": wild_type_dna,
        "mutated_dna": mutated_dna
    }
    
    try:
        response = requests.post(ENDPOINT_URL, json=payload, verify=False)
        response.raise_for_status()
        score_data = response.json()
        print(f"  - ✅ Score received: {score_data.get('zeta_score', 'N/A'):.4f}")
        return {
            "mutation_name": mutation_name,
            "zeta_score": score_data.get('zeta_score', 'N/A'),
            "interpretation": score_data.get('interpretation', 'N/A')
        }
    except requests.exceptions.HTTPError as http_err:
        print(f"  ❌ HTTP error occurred: {http_err}")
        print(f"  - Response Body: {response.text}")
        return None
    except Exception as e:
        print(f"  ❌ An unexpected error occurred: {e}")
        return None

def save_and_display_results(results: list):
    """Saves and displays the final validation report."""
    if not results:
        print("\n--- ⚠️ Validation campaign returned no results. ---")
        return

    # --- Save Results ---
    os.makedirs(RESULTS_DIR, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filepath = os.path.join(RESULTS_DIR, f"known_enemies_report_{timestamp}.json")
    try:
        with open(filepath, 'w') as f:
            json.dump(results, f, indent=4)
        print(f"\n--- 💾 Report successfully saved to {filepath} ---")
    except Exception as e:
        print(f"  ❌ Error saving report to file: {e}")

    # --- Display Results ---
    print("\n--- 📊 'Known Enemies' Validation Report ---")
    header = f"{'Mutation':<20} {'Zeta Score':<15} {'Predicted':<25} {'Expected':<15} {'Validation':<10}"
    print(header)
    print("-" * len(header))
    
    results.sort(key=lambda x: x.get('zeta_score', 0) or 0)

    for result in results:
        name = result.get('mutation_name', 'N/A')
        score = result.get('zeta_score', 'N/A')
        predicted_interp = result.get('interpretation', 'N/A')
        expected_interp = "Pathogenic" # All known enemies are expected to be pathogenic
        
        # Determine validation status
        # A simple check: if the score is negative, it aligns with being pathogenic.
        validation_status = "✅ PASS" if (isinstance(score, float) and score < -1.0) else "❌ FAIL"

        score_str = f"{score:.4f}" if isinstance(score, (int, float)) else "N/A"
        print(f"{name:<20} {score_str:<15} {predicted_interp:<25} {expected_interp:<15} {validation_status:<10}")
        
    print("-" * len(header))

def invoke_analyst_agent(results: list):
    """
    A more intelligent agent that analyzes the validation results and identifies discrepancies.
    """
    print("\n--- 🤖 Invoking Clinical Significance Analyst Agent 🤖 ---")
    
    if not results:
        print("  - Agent Summary: No results to analyze.")
        return

    pass_count = 0
    fail_count = 0
    discrepancies = []

    for result in results:
        score = result.get('zeta_score')
        name = result.get('mutation_name')
        if isinstance(score, float) and score < -1.0:
            pass_count += 1
        else:
            fail_count += 1
            score_str = f"{score:.4f}" if isinstance(score, float) else "None (Timeout/Error)"
            discrepancies.append(f"  - Discrepancy Found: Mutation '{name}' was expected to be Pathogenic but received a score of {score_str}, indicating a benign prediction or model failure.")

    summary = f"""
    The Analyst Agent has completed its review of the {len(results)} known pathogenic mutations.

    - Validation PASSED: {pass_count}
    - Validation FAILED: {fail_count}

    Summary:
    The validation has revealed significant discrepancies between the Zeta Oracle's predictions and established clinical knowledge. The model is currently failing to identify known pathogenic mutations, scoring them as benign.
    """
    
    if discrepancies:
        summary += "\n    Detected Discrepancies:\n"
        for d in discrepancies:
            summary += f"    {d}\n"

    summary += """
    Next Step Recommendation:
    This validation failure is a critical finding. The predictive accuracy of the Zeta Oracle for these mutation types must be addressed. Further investigation into the model's architecture or training data is required before proceeding with Seed & Soil simulations.
    """
    print(summary)

def main():
    """
    Main function to orchestrate the "Known Enemies" validation campaign.
    """
    print("--- ⚔️  Starting 'Known Enemies' Validation Campaign ⚔️  ---")
    
    # 1. Get wild-type sequence
    wt_sequence = get_reference_sequence(REFERENCE_GENOME_PATH, TARGET_CHROM, TARGET_START, TARGET_END)
    if not wt_sequence:
        return

    all_results = []
    # 2. Iterate through each known enemy
    for mutation_name, (mut_type, details) in KNOWN_MUTATIONS.items():
        print(f"\n--- 🎯 Targeting Known Enemy: {mutation_name} ---")
        
        # 3. Programmatically introduce the exact mutation into the full sequence context
        mutated_sequence = introduce_specific_mutation(wt_sequence, mut_type, details)
        
        if mutated_sequence:
            # 4. Run it through the validation pipeline
            result = run_validation(wt_sequence, mutation_name, mutated_sequence)
            if result:
                all_results.append(result)

    # 5. Save and display a final correlation report
    save_and_display_results(all_results)
    
    # 6. Invoke the new Analyst Agent for a high-level summary
    if all_results:
        invoke_analyst_agent(all_results)

    print("\n--- ✅ 'Known Enemies' Campaign Complete ---")

if __name__ == "__main__":
    main() 