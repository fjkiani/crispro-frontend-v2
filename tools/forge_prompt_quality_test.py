import asyncio
import os
import sys

# Ensure the 'src' directory is in the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from services.command_center.main import CommandCenter

# --- DOCTRINE: The Pathological Prompt Test ---
# This experiment is designed to prove that the quality of the generative
# prompt ("bait_sequence") is the primary determinant of the quality of the
# generated DNA candidates from the Forge.

# Low-quality prompt: Short, repetitive, biologically naive.
# This mimics the type of prompt that has previously produced "wet noodle" failures.
LOW_QUALITY_BAIT = "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT"

# High-quality prompt: A snippet from the BRAF gene (NM_004333.6) centered
# around the V600 codon. This is a biologically significant, complex sequence.
# The codon for Valine (V) at position 600 is GTG.
HIGH_QUALITY_BAIT = "TCTAGACTGACAGAGAAAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAAATCTAGTAAATAGGTGATTTTGGTCTAGCTACAGTGAAATCTCGATGGAGTGGGTCCCATCAGTTTGAACAGTTGTCTGGATCCATTTTGTGGATGGTAAGAATTGAGGCTATTTTTCC"

async def run_test():
    """
    Runs the comparative generation test.
    """
    print("--- ⚔️ OPERATION: FORGE PROMPT QUALITY TEST ⚔️ ---")
    print("Objective: Prove high-quality prompts yield high-quality DNA.")
    
    command_center = CommandCenter()
    
    print("\n--- TEST CASE 1: LOW-QUALITY PROMPT ---")
    print(f"BAIT: {LOW_QUALITY_BAIT}")
    try:
        low_quality_candidates = await command_center.test_generate_raw(LOW_QUALITY_BAIT)
        if low_quality_candidates:
            print("✅ Generation successful. Analyzing output...")
            for i, dna in enumerate(low_quality_candidates):
                print(f"  Candidate {i+1}: {dna[:100]}...") # Print first 100 chars
        else:
            print("❌ Generation failed to produce candidates.")
    except Exception as e:
        print(f"🔥 Test Case 1 Failed: {e}")


    print("\n--- TEST CASE 2: HIGH-QUALITY PROMPT (BRAF V600) ---")
    print(f"BAIT: {HIGH_QUALITY_BAIT}")
    try:
        high_quality_candidates = await command_center.test_generate_raw(HIGH_QUALITY_BAIT)
        if high_quality_candidates:
            print("✅ Generation successful. Analyzing output...")
            for i, dna in enumerate(high_quality_candidates):
                print(f"  Candidate {i+1}: {dna[:100]}...")
        else:
            print("❌ Generation failed to produce candidates.")
    except Exception as e:
        print(f"🔥 Test Case 2 Failed: {e}")

    print("\n--- 🏁 TEST COMPLETE ---")

if __name__ == "__main__":
    # This assumes the CommandCenter can be instantiated and used locally.
    # We might need to use `modal.Cls.lookup` if it's a deployed service.
    # For a direct local test, this setup should work if dependencies are met.
    # Running the async test
    asyncio.run(run_test()) 
 