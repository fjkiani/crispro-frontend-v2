import requests
import json
import sys

ORACLE_URL = "https://crispro--unified-oracle-service-unifiedoracle-invoke.modal.run"
# --- DOCTRINE: GAP-FILL PROMPTING ---
# The original prompt's poly-A tract is a pathological attractor. We remove it entirely.
PROMPT = "TGGATACTGAAACACATTAGAAGAAATGACTGAACTGTCTAAAGCTGTAGTAAATATACAATGCACTTTCTTCCAAAATGTAGAAAAGGCTGAATTCATCTAAAAAGAAATACATTCCTGAGAGGCAGAGGTTGCAGTGAGCCGAGATCCTGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCATCTCA"
# The original reference sequence we are trying to replace.
ORIGINAL_REFERENCE_SEQUENCE = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
# We must also update the reference for scoring to be the full, original prompt.
FULL_REFERENCE_SEQUENCE_FOR_SCORING = PROMPT + ORIGINAL_REFERENCE_SEQUENCE

NUM_CANDIDATES_TO_EVALUATE = 15

def generate_candidate():
    """Generates a single diverse candidate from the Oracle."""
    generate_payload = {
        "action": "generate",
        # We ask the model to generate enough tokens to fill the gap.
        "params": {"prompt": PROMPT, "gen_params": {"n_tokens": len(ORIGINAL_REFERENCE_SEQUENCE)}}
    }
    response = requests.post(ORACLE_URL, json=generate_payload, timeout=300)
    response.raise_for_status()
    return response.json()

def score_candidate(generated_sequence: str):
    """Scores a single candidate sequence against the prompt."""
    # The alternate sequence is the prompt plus the generated gap-fill.
    full_alternate_sequence = PROMPT + generated_sequence
    score_payload = {
        "action": "score",
        "params": {
            "reference_sequence": FULL_REFERENCE_SEQUENCE_FOR_SCORING,
            "alternate_sequence": full_alternate_sequence
        }
    }
    response = requests.post(ORACLE_URL, json=score_payload, timeout=300)
    response.raise_for_status()
    return response.json()

def run_intelligent_selection_protocol():
    """
    Implements the full generate-score-select loop to find a high-quality candidate.
    """
    print("🚀 Initiating Intelligent Selection Protocol...")
    
    candidates = []
    
    # --- Step 1: Generate a diverse pool of candidates ---
    print(f"\n--- STEP 1: Generating {NUM_CANDIDATES_TO_EVALUATE} candidates ---")
    for i in range(NUM_CANDIDATES_TO_EVALUATE):
        print(f"   - Generating candidate {i+1}/{NUM_CANDIDATES_TO_EVALUATE}...")
        gen_result = generate_candidate()
        if gen_result.get("status") == "success":
            candidates.append(gen_result.get("completion"))
        else:
            print(f"     - Generation failed or filtered: {gen_result.get('message')}")

    if not candidates:
        print("❌ FAILED: No valid candidates were generated.")
        sys.exit(1)
        
    print(f"✅ Generated {len(candidates)} valid candidates.")

    # --- Step 2: Score all candidates and find the best one ---
    print("\n--- STEP 2: Scoring and selecting the best candidate ---")
    best_candidate = None
    best_score = -float("inf")

    for i, candidate_seq in enumerate(candidates):
        print(f"   - Scoring candidate {i+1}/{len(candidates)}: '{candidate_seq}'")
        score_result = score_candidate(candidate_seq)
        if score_result.get("status") == "success":
            delta_score = score_result.get("delta_score")
            print(f"     - Score: {delta_score:.4f}")
            if delta_score > best_score:
                best_score = delta_score
                best_candidate = candidate_seq
        else:
            print(f"     - Scoring failed: {score_result.get('message')}")
    
    if best_candidate is None:
        print("❌ FAILED: All candidates failed scoring.")
        sys.exit(1)

    print(f"\n✅ Best candidate found: '{best_candidate}' with score: {best_score:.4f}")

    # --- Step 3: Final Verdict ---
    print("\n--- STEP 3: FINAL VERDICT ---")
    if best_score >= 0:
        print(f"🎉 MISSION ACCOMPLISHED! High-quality sequence found with a positive delta_score.")
    else:
        print(f"❌ FAILED: No candidate achieved a positive delta_score. Best score was {best_score:.4f}.")
        sys.exit(1)

if __name__ == "__main__":
    run_intelligent_selection_protocol() 
 
 
 
 
 
 