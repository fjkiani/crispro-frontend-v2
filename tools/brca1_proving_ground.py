import requests
import json
import sys
import os
from Bio import SeqIO
import gzip

# --- Configuration ---
ORACLE_URL = "https://crispro--unified-oracle-service-unifiedoracle-invoke.modal.run"
WINDOW_SIZE = 8192

# --- Target Data (from brca1_target_acquisition.py) ---
TARGET_POS_HG19 = 41197726
REF_ALLELE = 'A'
PATHOGENIC_ALLELE = 'G'
BENIGN_ALLELE = 'T'

def get_brca1_sequence():
    """Reads the reference genome sequence for chromosome 17 (where BRCA1 is located)."""
    # This path is relative to the root of the repository, as seen in the notebook.
    fasta_path = os.path.join('scripts', 'evo2', 'evo2', 'notebooks', 'brca1', 'GRCh37.p13_chr17.fna.gz')
    with gzip.open(fasta_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return str(record.seq)
    raise FileNotFoundError("BRCA1 reference sequence not found.")

def construct_prompt(full_sequence: str, pos: int, variant_allele: str) -> str:
    """
    Constructs a gapped prompt for the Oracle.
    It creates the full sequence with the pathogenic variant, then removes it.
    """
    p = pos - 1 # Convert to 0-indexed
    
    # 1. Create the full sequence with the pathogenic variant
    ref_seq_start = max(0, p - WINDOW_SIZE//2)
    ref_seq_end = min(len(full_sequence), p + WINDOW_SIZE//2)
    context_seq = full_sequence[ref_seq_start:ref_seq_end]
    
    snv_pos_in_context = min(WINDOW_SIZE//2, p)
    
    # Ensure the reference allele is correct
    assert context_seq[snv_pos_in_context] == REF_ALLELE, "Reference allele mismatch!"
    
    pathogenic_sequence = context_seq[:snv_pos_in_context] + variant_allele + context_seq[snv_pos_in_context+1:]
    
    # 2. Create the gapped prompt by removing the pathogenic allele
    prompt = pathogenic_sequence[:snv_pos_in_context]
    
    return prompt, snv_pos_in_context

def run_brca1_proving_ground_test():
    """
    Executes the 'Genetic Repair' test to validate the Oracle's generative intelligence.
    """
    print("🚀 Initiating Operation: BRCA1 Proving Ground...")

    # --- Step 1: Acquire sequence data and construct prompt ---
    print("\n--- STEP 1: Constructing 'Genetic Repair' Prompt ---")
    chr17_sequence = get_brca1_sequence()
    prompt, snv_pos = construct_prompt(chr17_sequence, TARGET_POS_HG19, PATHOGENIC_ALLELE)
    
    print(f"   - Prompt created from pathogenic context around position {TARGET_POS_HG19}")
    print(f"   - Prompt ending: ...{prompt[-30:]}")

    # --- Step 2: Task the Oracle with filling the gap ---
    print("\n--- STEP 2: Commanding Oracle to Generate Repair ---")
    generate_payload = {
        "action": "generate",
        "params": {"prompt": prompt, "gen_params": {"n_tokens": 1}} # Generate the single missing base
    }

    try:
        response = requests.post(ORACLE_URL, json=generate_payload, timeout=300)
        response.raise_for_status()
        gen_result = response.json()

        if gen_result.get("status") != "success":
            print(f"❌ FAILED: Generation did not succeed. Reason: {gen_result.get('message')}")
            sys.exit(1)
        
        generated_allele = gen_result.get("completion", "")
        print(f"✅ Oracle generated allele: '{generated_allele}'")

    except requests.exceptions.RequestException as e:
        print(f"❌ CRITICAL FAILURE: Could not connect to the Oracle service: {e}")
        sys.exit(1)

    # --- Step 3: Final Verdict ---
    print("\n--- STEP 3: FINAL VERDICT ---")
    if generated_allele == BENIGN_ALLELE:
        print(f"🎉 MISSION ACCOMPLISHED! Oracle correctly generated the benign allele ('{BENIGN_ALLELE}').")
    elif generated_allele == REF_ALLELE:
        print(f"✅ SUCCESS (Partial): Oracle generated the reference allele ('{REF_ALLELE}'). This is an acceptable, non-pathogenic outcome.")
    elif generated_allele == PATHOGENIC_ALLELE:
        print(f"❌ FAILED: Oracle regenerated the pathogenic allele ('{PATHOGENIC_ALLELE}').")
        sys.exit(1)
    else:
        print(f"⚠️ UNKNOWN OUTCOME: Oracle generated a novel allele ('{generated_allele}'). Further analysis would be required.")
        sys.exit(1)


if __name__ == "__main__":
    run_brca1_proving_ground_test() 
import json
import sys
import os
from Bio import SeqIO
import gzip

# --- Configuration ---
ORACLE_URL = "https://crispro--unified-oracle-service-unifiedoracle-invoke.modal.run"
WINDOW_SIZE = 8192

# --- Target Data (from brca1_target_acquisition.py) ---
TARGET_POS_HG19 = 41197726
REF_ALLELE = 'A'
PATHOGENIC_ALLELE = 'G'
BENIGN_ALLELE = 'T'

def get_brca1_sequence():
    """Reads the reference genome sequence for chromosome 17 (where BRCA1 is located)."""
    # This path is relative to the root of the repository, as seen in the notebook.
    fasta_path = os.path.join('scripts', 'evo2', 'evo2', 'notebooks', 'brca1', 'GRCh37.p13_chr17.fna.gz')
    with gzip.open(fasta_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return str(record.seq)
    raise FileNotFoundError("BRCA1 reference sequence not found.")

def construct_prompt(full_sequence: str, pos: int, variant_allele: str) -> str:
    """
    Constructs a gapped prompt for the Oracle.
    It creates the full sequence with the pathogenic variant, then removes it.
    """
    p = pos - 1 # Convert to 0-indexed
    
    # 1. Create the full sequence with the pathogenic variant
    ref_seq_start = max(0, p - WINDOW_SIZE//2)
    ref_seq_end = min(len(full_sequence), p + WINDOW_SIZE//2)
    context_seq = full_sequence[ref_seq_start:ref_seq_end]
    
    snv_pos_in_context = min(WINDOW_SIZE//2, p)
    
    # Ensure the reference allele is correct
    assert context_seq[snv_pos_in_context] == REF_ALLELE, "Reference allele mismatch!"
    
    pathogenic_sequence = context_seq[:snv_pos_in_context] + variant_allele + context_seq[snv_pos_in_context+1:]
    
    # 2. Create the gapped prompt by removing the pathogenic allele
    prompt = pathogenic_sequence[:snv_pos_in_context]
    
    return prompt, snv_pos_in_context

def run_brca1_proving_ground_test():
    """
    Executes the 'Genetic Repair' test to validate the Oracle's generative intelligence.
    """
    print("🚀 Initiating Operation: BRCA1 Proving Ground...")

    # --- Step 1: Acquire sequence data and construct prompt ---
    print("\n--- STEP 1: Constructing 'Genetic Repair' Prompt ---")
    chr17_sequence = get_brca1_sequence()
    prompt, snv_pos = construct_prompt(chr17_sequence, TARGET_POS_HG19, PATHOGENIC_ALLELE)
    
    print(f"   - Prompt created from pathogenic context around position {TARGET_POS_HG19}")
    print(f"   - Prompt ending: ...{prompt[-30:]}")

    # --- Step 2: Task the Oracle with filling the gap ---
    print("\n--- STEP 2: Commanding Oracle to Generate Repair ---")
    generate_payload = {
        "action": "generate",
        "params": {"prompt": prompt, "gen_params": {"n_tokens": 1}} # Generate the single missing base
    }

    try:
        response = requests.post(ORACLE_URL, json=generate_payload, timeout=300)
        response.raise_for_status()
        gen_result = response.json()

        if gen_result.get("status") != "success":
            print(f"❌ FAILED: Generation did not succeed. Reason: {gen_result.get('message')}")
            sys.exit(1)
        
        generated_allele = gen_result.get("completion", "")
        print(f"✅ Oracle generated allele: '{generated_allele}'")

    except requests.exceptions.RequestException as e:
        print(f"❌ CRITICAL FAILURE: Could not connect to the Oracle service: {e}")
        sys.exit(1)

    # --- Step 3: Final Verdict ---
    print("\n--- STEP 3: FINAL VERDICT ---")
    if generated_allele == BENIGN_ALLELE:
        print(f"🎉 MISSION ACCOMPLISHED! Oracle correctly generated the benign allele ('{BENIGN_ALLELE}').")
    elif generated_allele == REF_ALLELE:
        print(f"✅ SUCCESS (Partial): Oracle generated the reference allele ('{REF_ALLELE}'). This is an acceptable, non-pathogenic outcome.")
    elif generated_allele == PATHOGENIC_ALLELE:
        print(f"❌ FAILED: Oracle regenerated the pathogenic allele ('{PATHOGENIC_ALLELE}').")
        sys.exit(1)
    else:
        print(f"⚠️ UNKNOWN OUTCOME: Oracle generated a novel allele ('{generated_allele}'). Further analysis would be required.")
        sys.exit(1)


if __name__ == "__main__":
    run_brca1_proving_ground_test() 
 
 
 
 
 