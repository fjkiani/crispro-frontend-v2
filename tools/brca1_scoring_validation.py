import requests
import json
import sys
import os
from Bio import SeqIO
import gzip

# --- Configuration ---
ORACLE_URL = "https://crispro--unified-oracle-service-unifiedoracle-invoke.modal.run"
WINDOW_SIZE = 8192

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

def construct_variant_sequences(full_sequence: str, pos: int):
    """Constructs the full reference and both variant sequences."""
    p = pos - 1 # 0-indexed
    ref_seq_start = max(0, p - WINDOW_SIZE//2)
    ref_seq_end = min(len(full_sequence), p + WINDOW_SIZE//2)
    context_seq = full_sequence[ref_seq_start:ref_seq_end]
    snv_pos_in_context = min(WINDOW_SIZE//2, p)
    
    assert context_seq[snv_pos_in_context] == REF_ALLELE, "Reference allele mismatch!"
    
    pathogenic_seq = context_seq[:snv_pos_in_context] + PATHOGENIC_ALLELE + context_seq[snv_pos_in_context+1:]
    benign_seq = context_seq[:snv_pos_in_context] + BENIGN_ALLELE + context_seq[snv_pos_in_context+1:]
    
    return context_seq, pathogenic_seq, benign_seq

def score_variant(reference_sequence: str, alternate_sequence: str):
    """Scores a single variant against the reference."""
    score_payload = {
        "action": "score",
        "params": {
            "reference_sequence": reference_sequence,
            "alternate_sequence": alternate_sequence
        }
    }
    response = requests.post(ORACLE_URL, json=score_payload, timeout=300)
    response.raise_for_status()
    return response.json()

def run_brca1_scoring_validation():
    """
    Executes a scoring-based validation of the Oracle's understanding of BRCA1 variants.
    """
    print("🚀 Initiating Operation: BRCA1 Scoring Validation...")

    # --- Step 1: Construct variant sequences ---
    print("\n--- STEP 1: Constructing Variant Sequences ---")
    chr17_sequence = get_brca1_sequence()
    ref_seq, pathogenic_seq, benign_seq = construct_variant_sequences(chr17_sequence, TARGET_POS_HG19)
    print("   - Sequences for pathogenic and benign variants constructed.")

    # --- Step 2: Score both variants against the reference ---
    print("\n--- STEP 2: Scoring Variants ---")
    
    try:
        print("   - Scoring pathogenic variant...")
        pathogenic_result = score_variant(ref_seq, pathogenic_seq)
        pathogenic_score = pathogenic_result.get("delta_score")
        print(f"     - Pathogenic Score: {pathogenic_score:.4f}")

        print("   - Scoring benign variant...")
        benign_result = score_variant(ref_seq, benign_seq)
        benign_score = benign_result.get("delta_score")
        print(f"     - Benign Score: {benign_score:.4f}")

    except requests.exceptions.RequestException as e:
        print(f"❌ CRITICAL FAILURE: Could not connect to the Oracle service: {e}")
        sys.exit(1)

    # --- Step 3: Final Verdict ---
    print("\n--- STEP 3: FINAL VERDICT ---")
    if benign_score is None or pathogenic_score is None:
        print("❌ FAILED: Could not retrieve scores for one or both variants.")
        sys.exit(1)
        
    if benign_score > pathogenic_score:
        print(f"🎉 MISSION ACCOMPLISHED! Oracle correctly scored the benign variant ({benign_score:.4f}) as less deleterious than the pathogenic variant ({pathogenic_score:.4f}).")
    else:
        print(f"❌ FAILED: Oracle incorrectly scored the pathogenic variant as less deleterious.")
        sys.exit(1)

if __name__ == "__main__":
    run_brca1_scoring_validation() 
import json
import sys
import os
from Bio import SeqIO
import gzip

# --- Configuration ---
ORACLE_URL = "https://crispro--unified-oracle-service-unifiedoracle-invoke.modal.run"
WINDOW_SIZE = 8192

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

def construct_variant_sequences(full_sequence: str, pos: int):
    """Constructs the full reference and both variant sequences."""
    p = pos - 1 # 0-indexed
    ref_seq_start = max(0, p - WINDOW_SIZE//2)
    ref_seq_end = min(len(full_sequence), p + WINDOW_SIZE//2)
    context_seq = full_sequence[ref_seq_start:ref_seq_end]
    snv_pos_in_context = min(WINDOW_SIZE//2, p)
    
    assert context_seq[snv_pos_in_context] == REF_ALLELE, "Reference allele mismatch!"
    
    pathogenic_seq = context_seq[:snv_pos_in_context] + PATHOGENIC_ALLELE + context_seq[snv_pos_in_context+1:]
    benign_seq = context_seq[:snv_pos_in_context] + BENIGN_ALLELE + context_seq[snv_pos_in_context+1:]
    
    return context_seq, pathogenic_seq, benign_seq

def score_variant(reference_sequence: str, alternate_sequence: str):
    """Scores a single variant against the reference."""
    score_payload = {
        "action": "score",
        "params": {
            "reference_sequence": reference_sequence,
            "alternate_sequence": alternate_sequence
        }
    }
    response = requests.post(ORACLE_URL, json=score_payload, timeout=300)
    response.raise_for_status()
    return response.json()

def run_brca1_scoring_validation():
    """
    Executes a scoring-based validation of the Oracle's understanding of BRCA1 variants.
    """
    print("🚀 Initiating Operation: BRCA1 Scoring Validation...")

    # --- Step 1: Construct variant sequences ---
    print("\n--- STEP 1: Constructing Variant Sequences ---")
    chr17_sequence = get_brca1_sequence()
    ref_seq, pathogenic_seq, benign_seq = construct_variant_sequences(chr17_sequence, TARGET_POS_HG19)
    print("   - Sequences for pathogenic and benign variants constructed.")

    # --- Step 2: Score both variants against the reference ---
    print("\n--- STEP 2: Scoring Variants ---")
    
    try:
        print("   - Scoring pathogenic variant...")
        pathogenic_result = score_variant(ref_seq, pathogenic_seq)
        pathogenic_score = pathogenic_result.get("delta_score")
        print(f"     - Pathogenic Score: {pathogenic_score:.4f}")

        print("   - Scoring benign variant...")
        benign_result = score_variant(ref_seq, benign_seq)
        benign_score = benign_result.get("delta_score")
        print(f"     - Benign Score: {benign_score:.4f}")

    except requests.exceptions.RequestException as e:
        print(f"❌ CRITICAL FAILURE: Could not connect to the Oracle service: {e}")
        sys.exit(1)

    # --- Step 3: Final Verdict ---
    print("\n--- STEP 3: FINAL VERDICT ---")
    if benign_score is None or pathogenic_score is None:
        print("❌ FAILED: Could not retrieve scores for one or both variants.")
        sys.exit(1)
        
    if benign_score > pathogenic_score:
        print(f"🎉 MISSION ACCOMPLISHED! Oracle correctly scored the benign variant ({benign_score:.4f}) as less deleterious than the pathogenic variant ({pathogenic_score:.4f}).")
    else:
        print(f"❌ FAILED: Oracle incorrectly scored the pathogenic variant as less deleterious.")
        sys.exit(1)

if __name__ == "__main__":
    run_brca1_scoring_validation() 
 
 
 
 
 