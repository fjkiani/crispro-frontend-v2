import httpx
import json
import os
import sys

# --- Path Correction ---
# Ensure the project root is in the Python path to allow for absolute imports
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from src.tools.pysam_client import PysamClient

# --- Configuration ---
# This is the stable, live Oracle that the CommandCenter is also using.
ORACLE_URL = "https://crispro--zeta-oracle-v3-log-fix-zetaoracle-api.modal.run/invoke"

# From run_patient_assessment.py and VCF analysis
REFERENCE_GENOME_PATH = "data/gene_database/reference/hg19.fa"
RUNX1_CHROM = "chr21"
RUNX1_GENE_REGION_START = 36160098
RUNX1_GENE_REGION_END = 36421599

# Known pathogenic missense mutation from ClinVar
MUTATION_CHROM = "chr21"
MUTATION_POS_GENOMIC = 36250941 # 1-based
REF_ALLELE = "C"
ALT_ALLELE = "T"

def main():
    """
    Executes a single, precise VEP test against the live ZetaOracle
    to verify its scoring logic.
    """
    print("--- 🚀 OPERATION: FINAL BLOW ---")
    print("Objective: Verify ZetaOracle's delta_score calculation.")
    print(f"Targeting Oracle at: {ORACLE_URL}")
    print(f"Known pathogenic variant: {MUTATION_CHROM}:{MUTATION_POS_GENOMIC} {REF_ALLELE}>{ALT_ALLELE}")

    # 1. Acquire reference sequence
    try:
        pysam_client = PysamClient(REFERENCE_GENOME_PATH)
        full_wild_type_dna = pysam_client.fetch_sequence(RUNX1_CHROM, RUNX1_GENE_REGION_START, RUNX1_GENE_REGION_END)
        print(f"\n[1] Reference sequence acquired (Length: {len(full_wild_type_dna)}bp).")
    except Exception as e:
        print(f"\n[1] FAILED: Could not retrieve reference sequence. Error: {e}")
        return

    # 2. Create mutated sequence
    mutation_pos_in_full_gene = MUTATION_POS_GENOMIC - RUNX1_GENE_REGION_START - 1 # 0-indexed
    mutated_dna_list = list(full_wild_type_dna)
    
    # Sanity check reference allele
    if mutated_dna_list[mutation_pos_in_full_gene] != REF_ALLELE:
        print(f"\n[2] FAILED: Reference allele mismatch at position {MUTATION_POS_GENOMIC}.")
        print(f"Expected '{REF_ALLELE}', but found '{mutated_dna_list[mutation_pos_in_full_gene]}'.")
        return

    mutated_dna_list[mutation_pos_in_full_gene] = ALT_ALLELE
    full_mutated_dna = "".join(mutated_dna_list)
    print("\n[2] Mutated sequence forged.")

    # 3. Apply 8k Windowing Protocol
    WINDOW_SIZE = 8192
    half_window = WINDOW_SIZE // 2
    
    slice_start = max(0, mutation_pos_in_full_gene - half_window)
    slice_end = min(len(full_wild_type_dna), mutation_pos_in_full_gene + half_window)
    
    ref_slice = full_wild_type_dna[slice_start:slice_end]
    alt_slice = full_mutated_dna[slice_start:slice_end]
    print(f"\n[3] 8kb windowing protocol applied. Slice length: {len(ref_slice)}bp.")

    # 4. Assemble and Fire Payload
    payload = {
        "action": "score",
        "params": {
            "reference_sequence": ref_slice,
            "alternate_sequence": alt_slice
        }
    }
    print("\n[4] Payload assembled. Firing sniper shot at ZetaOracle...")

    try:
        with httpx.Client(timeout=600.0) as client:
            response = client.post(ORACLE_URL, json=payload)
            response.raise_for_status()
            result = response.json()

        print("\n--- ✅ SNIPER SHOT HIT ---")
        print("--- ORACLE RESPONSE ---")
        print(json.dumps(result, indent=2))
        print("-----------------------")

        delta_score = result.get("delta_score", 0)
        ref_ll = result.get("reference_likelihood", "N/A")
        alt_ll = result.get("alternate_likelihood", "N/A")

        print(f"\n--- ANALYSIS ---")
        print(f"  Reference Log-Likelihood: {ref_ll}")
        print(f"  Alternate Log-Likelihood: {alt_ll}")
        print(f"  Zeta Score (Δ): {delta_score}")

        if isinstance(delta_score, float) and delta_score < -1.0:
            print("\n--- ✅ VICTORY: Oracle returned a significant negative score. The weapon is calibrated. ---")
        else:
            print("\n--- ❌ DEFEAT: Oracle returned a non-significant score. The logic is still flawed. ---")

    except httpx.HTTPStatusError as e:
        print(f"\n--- ❌ ORACLE CALL FAILED ---")
        print(f"  Status Code: {e.response.status_code}")
        print(f"  Response: {e.response.text}")
    except Exception as e:
        print(f"\n--- ❌ TEST FAILED ---")
        print(f"  Error: {e}")

if __name__ == "__main__":
    main() 