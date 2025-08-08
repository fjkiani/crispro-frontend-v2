import requests
import json
import sys
import os
from Bio import SeqIO
import gzip

# --- Configuration ---
ORACLE_URL = "https://crispro--unified-oracle-service-unifiedoracle-invoke.modal.run"
CONTEXT_SIZE = 1000 # The number of base pairs to provide as a prompt before the variant site.
GENERATION_LENGTH = 1050 # How many tokens to generate (must be > CONTEXT_SIZE to include the variant site)


# --- Target Data ---
TARGET_POS_HG19 = 41197726
REF_ALLELE = 'A'
PATHOGENIC_ALLELE = 'G'
BENIGN_ALLELE = 'T'

def get_brca1_sequence():
    """Reads the reference genome sequence for chromosome 17."""
    fasta_path = os.path.join('scripts', 'evo2', 'evo2', 'notebooks', 'brca1', 'GRCh37.p13_chr17.fna.gz')
    with gzip.open(fasta_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return str(record.seq)
    raise FileNotFoundError("BRCA1 reference sequence not found.")

def construct_completion_prompt(full_sequence: str, pos: int) -> str:
    """Constructs a truncation prompt for the gene completion task."""
    p = pos - 1 # 0-indexed
    
    prompt_start = max(0, p - CONTEXT_SIZE)
    prompt_end = p
    
    prompt = full_sequence[prompt_start:prompt_end]
    return prompt

def run_brca1_gene_completion_test():
    """
    Executes the 'Gene Completion' test, a long-form generative task.
    """
    print("🚀 Initiating Operation: BRCA1 Gene Completion...")

    # --- Step 1: Acquire sequence data and construct prompt ---
    print("\n--- STEP 1: Constructing Truncated Gene Prompt ---")
    chr17_sequence = get_brca1_sequence()
    prompt = construct_completion_prompt(chr17_sequence, TARGET_POS_HG19)
    
    print(f"   - Prompt created by truncating BRCA1 sequence {CONTEXT_SIZE}bp before position {TARGET_POS_HG19}")
    print(f"   - Prompt ending: ...{prompt[-30:]}")

    # --- Step 2: Task the Oracle with completing the gene ---
    print(f"\n--- STEP 2: Commanding Oracle to Generate {GENERATION_LENGTH} tokens ---")
    generate_payload = {
        "action": "generate",
        "params": {"prompt": prompt, "gen_params": {"n_tokens": GENERATION_LENGTH}}
    }

    try:
        response = requests.post(ORACLE_URL, json=generate_payload, timeout=600) # Increased timeout for long generation
        response.raise_for_status()
        gen_result = response.json()

        if gen_result.get("status") != "success":
            print(f"❌ FAILED: Generation did not succeed. Reason: {gen_result.get('message')}")
            sys.exit(1)
        
        generated_completion = gen_result.get("completion", "")
        print(f"✅ Oracle generated completion (first 30bp): {generated_completion[:30]}...")

    except requests.exceptions.RequestException as e:
        print(f"❌ CRITICAL FAILURE: Could not connect to the Oracle service: {e}")
        sys.exit(1)

    # --- Step 3: Final Verdict ---
    print("\n--- STEP 3: FINAL VERDICT ---")
    
    # The variant position is `CONTEXT_SIZE` bases into the generated sequence.
    # We must account for 0-indexing.
    variant_position_in_generation = CONTEXT_SIZE
    if len(generated_completion) <= variant_position_in_generation:
        print(f"❌ FAILED: Generated sequence is too short ({len(generated_completion)}bp) to contain the variant.")
        sys.exit(1)

    generated_allele = generated_completion[variant_position_in_generation]
    
    print(f"   - Allele at target position {TARGET_POS_HG19} is: '{generated_allele}'")

    if generated_allele == BENIGN_ALLELE:
        print(f"🎉 MISSION ACCOMPLISHED! Oracle correctly generated the benign allele ('{BENIGN_ALLELE}').")
    elif generated_allele == REF_ALLELE:
        print(f"✅ SUCCESS (Partial): Oracle generated the reference allele ('{REF_ALLELE}'). This is an acceptable, non-pathogenic outcome.")
    elif generated_allele == PATHOGENIC_ALLELE:
        print(f"❌ FAILED: Oracle generated the pathogenic allele ('{PATHOGENIC_ALLELE}').")
        sys.exit(1)
    else:
        print(f"⚠️ UNKNOWN OUTCOME: Oracle generated a novel allele ('{generated_allele}').")
        sys.exit(1)

if __name__ == "__main__":
    run_brca1_gene_completion_test() 
import json
import sys
import os
from Bio import SeqIO
import gzip

# --- Configuration ---
ORACLE_URL = "https://crispro--unified-oracle-service-unifiedoracle-invoke.modal.run"
CONTEXT_SIZE = 1000 # The number of base pairs to provide as a prompt before the variant site.
GENERATION_LENGTH = 1050 # How many tokens to generate (must be > CONTEXT_SIZE to include the variant site)


# --- Target Data ---
TARGET_POS_HG19 = 41197726
REF_ALLELE = 'A'
PATHOGENIC_ALLELE = 'G'
BENIGN_ALLELE = 'T'

def get_brca1_sequence():
    """Reads the reference genome sequence for chromosome 17."""
    fasta_path = os.path.join('scripts', 'evo2', 'evo2', 'notebooks', 'brca1', 'GRCh37.p13_chr17.fna.gz')
    with gzip.open(fasta_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return str(record.seq)
    raise FileNotFoundError("BRCA1 reference sequence not found.")

def construct_completion_prompt(full_sequence: str, pos: int) -> str:
    """Constructs a truncation prompt for the gene completion task."""
    p = pos - 1 # 0-indexed
    
    prompt_start = max(0, p - CONTEXT_SIZE)
    prompt_end = p
    
    prompt = full_sequence[prompt_start:prompt_end]
    return prompt

def run_brca1_gene_completion_test():
    """
    Executes the 'Gene Completion' test, a long-form generative task.
    """
    print("🚀 Initiating Operation: BRCA1 Gene Completion...")

    # --- Step 1: Acquire sequence data and construct prompt ---
    print("\n--- STEP 1: Constructing Truncated Gene Prompt ---")
    chr17_sequence = get_brca1_sequence()
    prompt = construct_completion_prompt(chr17_sequence, TARGET_POS_HG19)
    
    print(f"   - Prompt created by truncating BRCA1 sequence {CONTEXT_SIZE}bp before position {TARGET_POS_HG19}")
    print(f"   - Prompt ending: ...{prompt[-30:]}")

    # --- Step 2: Task the Oracle with completing the gene ---
    print(f"\n--- STEP 2: Commanding Oracle to Generate {GENERATION_LENGTH} tokens ---")
    generate_payload = {
        "action": "generate",
        "params": {"prompt": prompt, "gen_params": {"n_tokens": GENERATION_LENGTH}}
    }

    try:
        response = requests.post(ORACLE_URL, json=generate_payload, timeout=600) # Increased timeout for long generation
        response.raise_for_status()
        gen_result = response.json()

        if gen_result.get("status") != "success":
            print(f"❌ FAILED: Generation did not succeed. Reason: {gen_result.get('message')}")
            sys.exit(1)
        
        generated_completion = gen_result.get("completion", "")
        print(f"✅ Oracle generated completion (first 30bp): {generated_completion[:30]}...")

    except requests.exceptions.RequestException as e:
        print(f"❌ CRITICAL FAILURE: Could not connect to the Oracle service: {e}")
        sys.exit(1)

    # --- Step 3: Final Verdict ---
    print("\n--- STEP 3: FINAL VERDICT ---")
    
    # The variant position is `CONTEXT_SIZE` bases into the generated sequence.
    # We must account for 0-indexing.
    variant_position_in_generation = CONTEXT_SIZE
    if len(generated_completion) <= variant_position_in_generation:
        print(f"❌ FAILED: Generated sequence is too short ({len(generated_completion)}bp) to contain the variant.")
        sys.exit(1)

    generated_allele = generated_completion[variant_position_in_generation]
    
    print(f"   - Allele at target position {TARGET_POS_HG19} is: '{generated_allele}'")

    if generated_allele == BENIGN_ALLELE:
        print(f"🎉 MISSION ACCOMPLISHED! Oracle correctly generated the benign allele ('{BENIGN_ALLELE}').")
    elif generated_allele == REF_ALLELE:
        print(f"✅ SUCCESS (Partial): Oracle generated the reference allele ('{REF_ALLELE}'). This is an acceptable, non-pathogenic outcome.")
    elif generated_allele == PATHOGENIC_ALLELE:
        print(f"❌ FAILED: Oracle generated the pathogenic allele ('{PATHOGENIC_ALLELE}').")
        sys.exit(1)
    else:
        print(f"⚠️ UNKNOWN OUTCOME: Oracle generated a novel allele ('{generated_allele}').")
        sys.exit(1)

if __name__ == "__main__":
    run_brca1_gene_completion_test() 
 
 
 
 
 