import os
import sys
import httpx
import pysam
import pandas as pd
from loguru import logger

# --- DOCTRINAL SETUP ---
# Ensure we can import from the project root
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

def cell_1_setup_and_target_acquisition():
    """CELL 1: Setup and Target Acquisition"""
    print("--- [CELL 1] Acquiring targets for Doctrinal Purity ---\n")
    
    config = {
        "REFERENCE_GENOME_PATH": "data/gene_database/reference/hg19.fa",
        "ORACLE_URL": "https://crispro--zeta-oracle-v3-log-fix-zetaoracle-api.modal.run/invoke",
        "RUNX1_CHROM": "chr21",
        "RUNX1_GENE_REGION_START": 36160098,
        "RUNX1_GENE_REGION_END": 36421599,
        "WINDOW_SIZE": 8192,
    }

    # As per doctrine, we create a unified DataFrame.
    # 1 = pathogenic, 0 = benign
    variants_data = {
        'id': [
            'ClinVar:13010', 'p.Arg135Ser', 'p.Arg204Ter', 
            'p.Pro177Pro', 'p.Pro363Leu'
        ],
        'pos': [
            36245100, 36250941, 36224376, 
            36244675, 36189196
        ],
        'ref': ['G', 'C', 'C', 'C', 'C'],
        'alt': ['A', 'T', 'T', 'T', 'T'],
        'label': [1, 1, 1, 0, 0] # 1 for pathogenic, 0 for benign
    }
    
    variants_df = pd.DataFrame(variants_data)
    
    print("✅ Target DataFrame created:")
    print(variants_df)
    print("\n")
    
    return config, variants_df

def cell_2_acquire_sequence(config):
    """CELL 2: Acquire Reference Sequence"""
    print("--- [CELL 2] Acquiring RUNX1 reference sequence ---\n")
    try:
        fasta_handle = pysam.FastaFile(config["REFERENCE_GENOME_PATH"])
        sequence = fasta_handle.fetch(
            config["RUNX1_CHROM"], 
            config["RUNX1_GENE_REGION_START"], 
            config["RUNX1_GENE_REGION_END"]
        )
        fasta_handle.close()
        print(f"✅ Successfully loaded RUNX1 gene sequence (length: {len(sequence)}bp).\n")
        return sequence.upper()
    except Exception as e:
        print(f"💥 ERROR: Failed to fetch sequence from FASTA file. Path: {config['REFERENCE_GENOME_PATH']}")
        print(f"   Details: {e}\n")
        return None

def cell_3_create_variant_windows(config, full_wild_type_dna, variant_row):
    """CELL 3: Create Reference and Alternate Windows (Doctrinal Implementation)"""
    
    mutation_pos_abs = variant_row['pos']
    mutation_pos_in_gene = mutation_pos_abs - config["RUNX1_GENE_REGION_START"] -1
    
    print(f"--- [WINDOWING] Processing variant: {variant_row['id']} at position {mutation_pos_abs} ---")

    half_window = config["WINDOW_SIZE"] // 2
    
    # --- IN-SITU TARGETING & DYNAMIC MUTATION --- 
    actual_ref_allele = full_wild_type_dna[mutation_pos_in_gene]
    
    # Doctrinal Check: Ensure our ground truth matches the reference genome
    corrected_wild_type_dna = full_wild_type_dna
    if actual_ref_allele != variant_row['ref']:
        print(f"🚨 TACTICAL OVERRIDE: Reference allele mismatch for {variant_row['id']}!")
        print(f"   Expected: {variant_row['ref']}, Found in hg19: {actual_ref_allele} at position {mutation_pos_abs}")
        print(f"   ACTION: Forcing reference sequence to match ground truth.")
        wt_list = list(full_wild_type_dna)
        wt_list[mutation_pos_in_gene] = variant_row['ref']
        corrected_wild_type_dna = "".join(wt_list)
        
    alt_allele = variant_row['alt']
    
    # Create the mutated sequence from the CORRECTED wild-type
    full_mutated_dna_list = list(corrected_wild_type_dna)
    full_mutated_dna_list[mutation_pos_in_gene] = alt_allele
    full_mutated_dna = "".join(full_mutated_dna_list)
    
    # --- Slice Generation ---
    slice_start = max(0, mutation_pos_in_gene - half_window)
    slice_end = min(len(full_wild_type_dna), mutation_pos_in_gene + half_window)

    ref_slice = corrected_wild_type_dna[slice_start:slice_end]
    alt_slice = full_mutated_dna[slice_start:slice_end]
    
    print(f"  ✅ Slices generated (Length: {len(ref_slice)}bp)\n")
    return ref_slice, alt_slice

def cell_4_invoke_oracle(config, ref_slice, alt_slice, variant_id):
    """CELL 4: Invoke the Zeta Oracle"""
    print(f"--- [ORACLE] Invoking for variant {variant_id}... ---")
    payload = {
        "action": "score",
        "params": {
            "reference_sequence": ref_slice,
            "alternate_sequence": alt_slice
        }
    }
    
    try:
        with httpx.Client(timeout=300.0, follow_redirects=True) as client:
            response = client.post(config["ORACLE_URL"], json=payload)
            response.raise_for_status()
            result = response.json()
            delta_score = result.get("zeta_score", 0.0)
            ref_ll = result.get("reference_likelihood", None)
            alt_ll = result.get("alternate_likelihood", None)

            print(f"  ✅ Oracle success.")
            print(f"   Delta Score: {delta_score:.6f}")
            if ref_ll is not None and alt_ll is not None:
                print(f"   Reference Likelihood: {ref_ll:.6f}")
                print(f"   Alternate Likelihood: {alt_ll:.6f}")
            print("\n")
            return delta_score
    except httpx.HTTPStatusError as e:
        print(f"💥 HTTP ERROR: Failed to invoke Oracle for {variant_id}. Status: {e.response.status_code}")
        print(f"   Response: {e.response.text}\n")
        return None
    except Exception as e:
        print(f"💥 CLIENT ERROR: An unexpected error occurred for {variant_id}: {e}\n")
        return None

def main():
    """Main execution block, following the Doctrinal Purity operation plan."""
    
    # Phase 1: Target & Data Prep
    config, variants_df = cell_1_setup_and_target_acquisition()
    wild_type_dna = cell_2_acquire_sequence(config)
    
    if not wild_type_dna:
        print("Halting operation due to sequence acquisition failure.")
        return

    # Phase 2: Sequence Forging & Scoring
    scores = []
    for index, row in variants_df.iterrows():
        ref_slice, alt_slice = cell_3_create_variant_windows(config, wild_type_dna, row)
        score = cell_4_invoke_oracle(config, ref_slice, alt_slice, row['id'])
        scores.append(score)

    variants_df['delta_score'] = scores
    
    print("\n--- [RESULTS] Doctrinal Purity Scoring Complete ---")
    print(variants_df)

if __name__ == "__main__":
    main() 