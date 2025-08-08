import torch
import esm
import logging
import re
import sys
import os

# Ensure tools directory is in the path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tools.esm_client import ESMClient, AA_CODES

# --- Configuration ---
MODEL_NAME = "esm2_t33_650M_UR50D"
WINDOW_SIZE = 1024
PROTEIN_FASTA_PATH = "data/reference/BRCA1.fasta"
VARIANT_HGVS = "p.Cys61Gly"

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def read_fasta_sequence(filepath):
    """Reads a FASTA file and returns the clean protein sequence."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    return "".join(lines[1:]).replace("\n", "").replace("\r", "").strip()

def main():
    """
    Standalone script to debug ESM-2 scoring with a sliding window.
    """
    logger.info("--- 🔬 ESM Client Debugger Initiated 🔬 ---")

    # 1. Load Model
    logger.info(f"Loading ESM model: {MODEL_NAME}...")
    try:
        model, alphabet = esm.pretrained.load_model_and_alphabet(MODEL_NAME)
        model.eval()
        if torch.cuda.is_available():
            model = model.cuda()
            logger.info("✅ Model loaded to GPU.")
        else:
            logger.info("✅ Model loaded to CPU.")
        batch_converter = alphabet.get_batch_converter()
    except Exception as e:
        logger.error(f"💥 Failed to load ESM model: {e}", exc_info=True)
        return

    # 2. Load Sequence
    logger.info(f"Loading protein sequence from: {PROTEIN_FASTA_PATH}")
    protein_sequence = read_fasta_sequence(PROTEIN_FASTA_PATH)
    logger.info(f"✅ Sequence loaded, length: {len(protein_sequence)}")

    # 3. Parse Variant
    logger.info(f"Parsing variant: {VARIANT_HGVS}")
    variant_no_p = VARIANT_HGVS[2:] if VARIANT_HGVS.startswith("p.") else VARIANT_HGVS
    match = re.match(r"([A-Z][a-z]{2}|[A-Z])(\d+)([A-Z][a-z]{2}|[A-Z])", variant_no_p)
    ref_aa_code, pos_str, alt_aa_code = match.groups()
    pos = int(pos_str) - 1
    ref_aa = AA_CODES.get(ref_aa_code, ref_aa_code)
    alt_aa = AA_CODES.get(alt_aa_code, alt_aa_code)
    logger.info(f"✅ Parsed variant: Ref={ref_aa}, Pos={pos+1}, Alt={alt_aa}")

    # 4. Implement Sliding Window
    logger.info(f"Applying sliding window of size {WINDOW_SIZE}...")
    half_window = WINDOW_SIZE // 2
    start_index = max(0, pos - half_window)
    end_index = min(len(protein_sequence), pos + half_window)
    if end_index - start_index < WINDOW_SIZE:
        if start_index == 0:
            end_index = min(len(protein_sequence), WINDOW_SIZE)
        else:
            start_index = max(0, len(protein_sequence) - WINDOW_SIZE)
    
    windowed_sequence = protein_sequence[start_index:end_index]
    window_pos = pos - start_index
    logger.info(f"✅ Window created: {start_index}-{end_index}, new position: {window_pos+1}")

    # 5. Create Mutated Sequence
    mutated_seq_list = list(windowed_sequence)
    mutated_seq_list[window_pos] = alt_aa
    mutated_window = "".join(mutated_seq_list)
    
    # 6. Batch and Score
    logger.info("Scoring wild-type and mutated windows...")
    data = [
        ("wild_type_window", windowed_sequence),
        ("mutated_window", mutated_window)
    ]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    if torch.cuda.is_available():
        batch_tokens = batch_tokens.cuda()

    try:
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[model.num_layers], return_contacts=False)
        
        log_probs = torch.log_softmax(results["logits"], dim=-1)
        
        # WT Score
        wt_tokens = batch_tokens[0, 1 : len(windowed_sequence) + 1]
        wt_log_probs = torch.gather(log_probs[0, 1 : len(windowed_sequence) + 1], 1, wt_tokens.unsqueeze(1)).squeeze(1)
        wt_ll = wt_log_probs.sum().item()

        # Mutant Score
        mut_tokens = batch_tokens[1, 1 : len(mutated_window) + 1]
        mut_log_probs = torch.gather(log_probs[1, 1 : len(mutated_window) + 1], 1, mut_tokens.unsqueeze(1)).squeeze(1)
        mut_ll = mut_log_probs.sum().item()

        delta_ll = mut_ll - wt_ll

        logger.info("\n--- ✅ SUCCESS! ---")
        logger.info(f"Wild-Type Window Log-Likelihood: {wt_ll:.4f}")
        logger.info(f"Mutated Window Log-Likelihood:   {mut_ll:.4f}")
        logger.info(f"Delta Log-Likelihood (ESM Score): {delta_ll:.4f}")

    except torch.cuda.OutOfMemoryError as e:
        logger.error("💥 FAILED: CUDA Out of Memory, even with the sliding window.")
        logger.error(e)
    except Exception as e:
        logger.error(f"💥 An unexpected error occurred during scoring: {e}", exc_info=True)


if __name__ == "__main__":
    main() 