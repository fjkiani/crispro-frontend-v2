import pytest
import time
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from tools.runx1_integration_plan import design_runx1_intervention_v3

# --- Mock Data ---
MOCK_TARGET_VARIANT = {
    "gene": "RUNX1",
    "protein_change": "p.Arg135fs",
    "variant_type": "germline",
    "id": "RUNX1_p.Arg135fs",
    "position": 36207648,
    "consequence": "frameshift_variant"
}
MOCK_PATIENT_ID = "LIVE_FIRE_TEST_01"

# --- LIVE-FIRE Integration Test ---

def test_design_intervention_live_fire():
    """
    LIVE-FIRE integration test for the full design_runx1_intervention workflow.
    This test makes REAL requests to the deployed Modal services.
    It will fail if the services are not running or if there are issues
    with the end-to-end communication.
    """
    print("--- ⚔️ INITIATING LIVE-FIRE TEST ⚔️ ---")
    print(f"Patient: {MOCK_PATIENT_ID}, Variant: {MOCK_TARGET_VARIANT['id']}")

    # --- Execute ---
    # We call the new V3 entry point function
    result = design_runx1_intervention_v3(MOCK_TARGET_VARIANT, MOCK_PATIENT_ID)

    # --- Assert ---
    print("--- Final Result from Live Fire Test ---")
    print(result)
    
    # Assert final result structure
    assert result is not None
    assert "error" not in result
    
    guide_design = result.get("guide_design", {})
    assert guide_design is not None
    
    top_guide = guide_design.get("top_guide", {})
    assert top_guide is not None
    assert "guide_sequence" in top_guide
    assert len(top_guide["guide_sequence"]) > 0

    all_guides = guide_design.get("all_guides", [])
    assert len(all_guides) > 0
    print(f"✅ Live-fire test complete. Successfully designed {len(all_guides)} guides.")
    print(f"Top guide: {top_guide['guide_sequence']}") 