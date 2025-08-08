import modal
import pytest
from loguru import logger
from Bio.Seq import Seq

# --- Test Data ---
# A high-quality, long-range bait sequence from the BRCA1 gene.
# This provides a realistic and complex prompt for the Oracle.
BRCA1_BAIT_SEQUENCE = "TGGATACCCCGGATTCTGGTGGTGCAAGATGCTGTGCCCCGAGAGCCGTTTGTGTGTGGAGGAGAAGGCTTTTATTCTGTCACAGGAGGTGCTGCTCCAGAGCAGCAGGAGGGAGGCAGGAAAGGAGGGGGCTGGATGGAGACAGAGGGCTCAGGTGAAGAGGGGCACAGGGTGGCACAGGCTGTGAGCCCTCCCGGACCCCGGCCTTTGCTCAGCTTGTCCTGCAGGATGGGGAGAGGGGAGGCATAGGGACAGGTAGAGGAGGCAGGAGGGGCTCAGGTGCCTCCAGGCAGGAGAGAGGGAGAAAGAGAAGGTGGAGGCTGTTAGGGGCCAGGCCCAGGAGGAAGAGACAGAGGAAGGAGGGGCAGGCAGCTGGCTGTGGGGGAGGCTCCCAGCAGCAGGAGGGTGAGGCTCAGCAGGTGGTCCAGGCCTTAGGAGGAGGGGTGGGGACTGTAGGAGGAGGCAGCTCTGGAGGCAGAACAGGCTGAGGGTGGTGAGGGGTGGACAGGCTGAGGGTGCAGGCTGAGGCGGGGAGAGAGGCTGAGGGTGCAGGCGGTGGCTGGGGCAGAGGGGAGACCTGAGGCTCAGGAGGTGAGGGTGGAGAGAGAGGGCCAGGCAGGCAGCTCTGGTGTGGAGGAGAGGGCTTCAGGTGAGGCCTCAGGTGAGGGAGAAGAAGAAGGTGGTGGCTGTC"

# --- Test Case 1: Forge Generation Quality ---

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
def test_forge_generates_valid_protein():
    """
    Tests if the Forge's raw output is biologically valid.
    1. Calls the `test_generate_raw` diagnostic method.
    2. Translates the resulting DNA.
    3. Asserts that no premature stop codon ('*') is present.
    """
    logger.info("--- 🧪 Running Test Case 1: Forge Generation Quality 🧪 ---")
    
    command_center = get_command_center()
    
    logger.info("Requesting raw DNA generation from the Forge...")
    generated_dna = command_center.test_generate_raw.remote(bait_sequence=BRCA1_BAIT_SEQUENCE)
    
    assert generated_dna, "Forge failed to generate any DNA sequence."
    logger.success(f"Forge produced DNA: {generated_dna[:50]}...")
    
    logger.info("Translating generated DNA to protein...")
    try:
        # Ensure the sequence is a multiple of 3 for translation
        remainder = len(generated_dna) % 3
        if remainder > 0:
            generated_dna = generated_dna[:-remainder]
            
        protein_sequence = str(Seq(generated_dna).translate())
        logger.info(f"Translated Protein: {protein_sequence}")
        
        # --- THE ASSERTION ---
        assert "*" not in protein_sequence, f"Generated DNA contains a premature stop codon!"
        logger.success("✅ PASSED: Generated DNA translates to a valid protein without premature stops.")
        
    except Exception as e: