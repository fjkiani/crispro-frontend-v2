
import asyncio
import sys
import os

# Add the BioMed-MCP directory to python path to find modules if needed, 
# though we are running via uv/python in the dir usually.
# But here we are using the mcp client to talk to the server script.

from fastmcp.client import Client
from fastmcp.client.transports import StdioTransport

async def audit_mbd4():
    print("=== Auditing MBD4 in Ovarian Cancer via BioMed-MCP ===")
    
    # Path to the BioMed-MCP server module
    # We assume we are running this from the root or we need to point to the right place.
    # The available tool implies 'biomed_agents' is a module.
    # We will try to run the server via uv run or python -m
    
    server_cwd = "/Users/fahadkiani/Desktop/development/crispr-assistant-main/scripts/data_acquisition/mcp_servers/BioMed-MCP"
    
    # Construct transport
    transport = StdioTransport(
        command="uv",
        args=["run", "python", "-m", "biomed_agents"],
        cwd=server_cwd
    )
    
    client = Client(transport)
    
    async with client:
        print("Connected to BioMed-MCP.")
        
        # 1. Literature Search
        query = "MBD4 loss Ovarian Cancer PARP inhibitor"
        print(f"\nSearching Literature for: '{query}'...")
        try:
            lit_result = await client.call_tool(
                "biomedical_literature_search",
                {
                    "query": query,
                    "max_papers": 5,
                    "synthesize_findings": True
                }
            )
            print(f"\n[Literature Result]\n{lit_result.data}")
        except Exception as e:
            print(f"Literature search failed: {e}")

        # 2. Clinical Trials
        condition = "Ovarian Cancer"
        print(f"\nSearching Clinical Trials for Condition: '{condition}' with keyword 'MBD4'...")
        # Note: The tool 'clinical_trials_research' might not have a keyword filter in the high level args shown in test_client.py
        # But let's try to pass it or just search for the condition and we'll see if we can filter.
        # Actually, let's search for "MBD4" as the condition/term if possible, or rely on the lit search.
        # The test_client.py showed 'condition', 'study_phase'.
        
        try:
            # We'll try to use a more specific condition string if it allows free text
            trials_result = await client.call_tool(
                "clinical_trials_research",
                {
                    "condition": "Ovarian Cancer MBD4", # Trying to combine
                    "max_studies": 5
                }
            )
            print(f"\n[Clinical Trials (MBD4) Result]\n{trials_result.data}")
        except Exception as e:
             # Fallback to just Ovarian Cancer and we scan manually/visually
             print(f"Specific trial search failed ({e}), trying generic Ovarian Cancer check...")


if __name__ == "__main__":
    asyncio.run(audit_mbd4())
