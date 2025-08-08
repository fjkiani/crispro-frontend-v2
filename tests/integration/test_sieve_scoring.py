import modal
import pytest
from loguru import logger

# --- Test Data ---
# We will use a real-world example from the BRAF gene, focusing on the V600E mutation.
# This provides a clinically relevant and well-understood test case.

# Reference sequence from BRAF (part of exon 15)
REFERENCE_DNA = "ATCATCCACAGAGACAGAGCA" 

# A "good" candidate: This sequence contains the canonical BRAF V600E mutation (GTG -> GAG).
# This is a known, functional, and highly significant oncogenic mutation. 
# The Oracle should recognize this as a biologically plausible change.
GOOD_DNA_V600E = "ATCATCCACAGAGACAGAGAA" 

# A "bad" candidate: This sequence has a single nucleotide deletion, causing a frameshift.
# This is a catastrophic mutation that should be heavily penalized by the Oracle.
BAD_DNA_FRAMESHIFT = "ATCATCCACAGAGACAGAGC" 

# --- Test Case 2: Sieve Scoring Sanity Check ---

def get_command_center():
    """Helper function to look up the deployed CommandCenter service."""
    try:
        # NOTE: The app name must EXACTLY match the name in `src/services/command_center/main.py`
        # Using the modern `from_name` method with both app and class name.
        return modal.Cls.from_name("command-center-v8-override-fix", "CommandCenter")
    except Exception as e:
        logger.error(f"Failed to lookup Modal class: {e}")
        pytest.fail("Could not connect to the deployed CommandCenter Modal service.")

@pytest.mark.integration
def test_sieve_distinguishes_good_from_bad():
    """
    Tests if the Sieve's scoring logic is sane.
    1. Calls the `test_sieve_sanity_check` diagnostic method with good and bad DNA.
    2. Asserts that the score for the good DNA is higher than the score for the bad DNA.
    """
    logger.info("--- 🧪 Running Test Case 2: Sieve Scoring Sanity Check 🧪 ---")
    
    command_center = get_command_center()
    
    logger.info("Requesting sanity check from the Sieve...")
    scores = command_center.test_sieve_sanity_check.remote(
        good_dna=GOOD_DNA_V600E,
        bad_dna=BAD_DNA_FRAMESHIFT,
        reference_dna=REFERENCE_DNA
    )
    
    assert scores, "Sieve failed to return any scores."
    logger.success(f"Sieve returned scores: {scores}")
    
    good_score = scores.get("good_score")
    bad_score = scores.get("bad_score")
    
    assert good_score is not None, "Response from Sieve missing 'good_score'."
    assert bad_score is not None, "Response from Sieve missing 'bad_score'."
    
    # --- THE ASSERTION ---
    assert good_score > bad_score, f"Sieve scoring is insane! Good score ({good_score}) should be greater than bad score ({bad_score})."