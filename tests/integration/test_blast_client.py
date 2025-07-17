import pytest
import asyncio
import sys
from pathlib import Path

# Add project root to the Python path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tools.command_center_client import CommandCenterClient

@pytest.mark.asyncio
async def test_direct_blast_call():
    """
    This test directly invokes the BLAST client to isolate the connection issue.
    """
    client = CommandCenterClient()
    # A known guide sequence that should have some off-targets
    test_guides = ["ATTGCTGCAGTCAGCTCGAT"] 

    print("--- Directly testing check_off_targets_with_blast ---")
    results = await client.check_off_targets_with_blast(test_guides)
    print(f"Results: {results}")

    assert results is not None
    assert len(results) == 1
    assert results[0]["off_target_count"] != 999, "Direct BLAST call failed and returned a dummy value." 