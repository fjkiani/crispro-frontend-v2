"""
OPERATION: DIRECT ORACLE STRIKE
OBJECTIVE: To get a definitive delta_score for a known RUNX1 variant
           by bypassing all intermediate services (CommandCenter, Frontend)
           and querying the ZetaOracle directly. This script leverages
           our existing RUNX1DataLoader to construct the necessary sequences.
"""
import asyncio
import os
import sys
import httpx
from loguru import logger

# --- Path Correction ---
# Ensure the 'src' directory is in the Python path for service imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from tools.runx1_data_loader import RUNX1DataLoader

# --- Target Definition ---
VARIANT_ID = "RUNX1_p.Arg135fs"
# This is the live, deployed Unified Oracle.
ORACLE_URL = "https://crispro--unified-oracle-service-unifiedoracle-invoke.modal.run"

def apply_variant(reference_sequence: str, variant_data: dict) -> str:
    """
    Applies a simple SNV or indel to a reference sequence.
    NOTE: This is a simplified function for this specific test case.
    """
    pos = variant_data['position']  # 1-based position from VCF
    ref_allele = variant_data['reference']
    alt_allele = variant_data['alternate']
    
    # Convert to 0-based index for Python strings
    start_index = pos - 1
    
    # Sanity check
    if reference_sequence[start_index : start_index + len(ref_allele)] != ref_allele:
        raise ValueError(
            f"Reference allele mismatch at position {pos}. "
            f"Expected {ref_allele}, found {reference_sequence[start_index : start_index + len(ref_allele)]}"
        )
        
    # Apply the variant
    mutated_sequence = (
        reference_sequence[:start_index] +
        alt_allele +
        reference_sequence[start_index + len(ref_allele):]
    )
    return mutated_sequence

async def run_strike():
    """
    Executes the direct strike on the ZetaOracle.
    """
    logger.info("--- ⚔️  OPERATION: DIRECT ORACLE STRIKE ⚔️  ---")
    
    # --- 1. Acquire Target Intel ---
    logger.info(f"Acquiring intel for variant: {VARIANT_ID} using RUNX1DataLoader...")
    try:
        loader = RUNX1DataLoader()
        variant_data = loader.get_variant_by_id(VARIANT_ID)
        if not variant_data:
            logger.error(f"Intel Failure: Could not find data for variant '{VARIANT_ID}'. Aborting.")
            return

        reference_dna = loader.load_runx1_sequence()
        if not reference_dna:
            logger.error("Intel Failure: Could not load RUNX1 reference sequence. Aborting.")
            return
            
        logger.success("✅ Intel acquired.")
    except Exception as e:
        logger.error(f"Intel Failure: Error during data loading. {e}", exc_info=True)
        return

    # --- 2. Construct Payloads ---
    logger.info("Constructing reference and alternate DNA sequences for scoring...")
    try:
        alternate_dna = apply_variant(reference_dna, variant_data)
        logger.success("✅ Payloads constructed.")
        # logger.debug(f"Reference DNA snippet: ...{reference_dna[variant_data['position']-1-10:variant_data['position']-1+10]}...")
        # logger.debug(f"Alternate DNA snippet: ...{alternate_dna[variant_data['position']-1-10:variant_data['position']-1+10]}...")
    except ValueError as e:
        logger.error(f"Payload Failure: Could not apply variant. {e}")
        return

    # --- 3. Execute Strike ---
    logger.info(f"Striking ZetaOracle at {ORACLE_URL}...")
    
    payload = {
        "action": "score",
        "params": {
            "reference_sequence": reference_dna,
            "alternate_sequence": alternate_dna
        }
    }

    async with httpx.AsyncClient(timeout=300.0) as client:
        try:
            response = await client.post(ORACLE_URL, json=payload)
            response.raise_for_status()
            result = response.json()
            delta_score = result.get("delta_score")
            
            if delta_score is None:
                logger.error(f"Strike Ineffective: Oracle response did not contain 'delta_score'. Full response: {result}")
            else:
                logger.success("✅ STRIKE CONFIRMED. Oracle responded.")

        except httpx.HTTPStatusError as e:
            logger.error(f"Strike Failed: Oracle responded with error {e.response.status_code}. Response: {e.response.text}")
            return
        except httpx.RequestError as e:
            logger.error(f"Strike Failed: Could not connect to Oracle. Error: {e}")
            return

    # --- 4. Deliver Final Intel ---
    logger.info("--- 🏁 FINAL INTELLIGENCE REPORT 🏁 ---")
    print(f"\nVariant ID    : {VARIANT_ID}")
    print(f"Delta Score   : {delta_score}\n")
    logger.info("--- MISSION COMPLETE ---")


if __name__ == "__main__":
    asyncio.run(run_strike()) 
 