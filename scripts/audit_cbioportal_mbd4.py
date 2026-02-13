import asyncio
import sys
import os
import json
from fastmcp.client import Client
from fastmcp.client.transports import StdioTransport

async def audit_cbioportal_mbd4():
    print("=== Auditing MBD4 in Ovarian Cancer via cBioPortal-MCP ===")
    
    server_cwd = "/Users/fahadkiani/Desktop/development/crispr-assistant-main/tools/cbioportal-mcp"
    
    # Transport for cbioportal-mcp
    transport = StdioTransport(
        command="uv",
        args=["run", "cbioportal-mcp"], # As per README
        cwd=server_cwd
    )
    
    client = Client(transport)
    
    async with client:
        print("Connected to cBioPortal-MCP.")
        
        # 1. Search for Ovarian Cancer Studies (Verification)
        print("\nSearching for Ovarian Cancer Studies...")
        studies_result = await client.call_tool(
            "search_studies",
            {"keyword": "Ovarian", "limit": 50}
        )
        
        try:
            # FastMCP returns a TextContent object in a list
            studies_text = studies_result.content[0].text
            studies_data = json.loads(studies_text)
            print(f"Found studies count: {len(studies_data.get('studies', []))}")
        except Exception as e:
            print(f"Failed to parse studies result: {e}")
            # print(studies_result) # Debug if needed

        # 2. Check MBD4 Mutations in specific high-value studies
        # 'hgsoc_tcga_gdc': High-Grade Serous Ovarian Cancer (TCGA)
        # 'ovarian_msk_2025': Serous Ovarian Cancer (MSK 2025)
        # 'hgsoc_tcga_pan_can_atlas_2018': TCGA PanCancer Atlas (often most comprehensive)
        target_studies = ["hgsoc_tcga_gdc", "ovarian_msk_2025", "ov_tcga_pan_can_atlas_2018"]
        
        for study in target_studies:
            print(f"\nChecking MBD4 in {study}...")
            # Sample List ID is typically study_id + "_all"
            sample_list_id = f"{study}_all"
            
            try:
                # Check study existence first
                # details_res = await client.call_tool("get_study_details", {"study_id": study})
                # ... skipping detail check to just try querying ...

                print(f"  Querying mutations for MBD4 in {sample_list_id}...")
                mutations_res = await client.call_tool(
                    "get_mutations_in_gene",
                    {
                        "gene_id": "MBD4",
                        "study_id": study,
                        "sample_list_id": sample_list_id
                    }
                )
                
                try:
                    mutations_text = mutations_res.content[0].text
                    mutations = json.loads(mutations_text)
                except Exception as e:
                    print(f"  Failed to parse mutations response: {e}")
                    # print(mutations_res)
                    continue
                
                if isinstance(mutations, dict) and "error" in mutations:
                     print(f"  Error querying API: {mutations['error']}")
                elif isinstance(mutations, dict) and "mutations" in mutations:
                     m_list = mutations["mutations"]
                     print(f"  Found {len(m_list)} mutations in MBD4 in {study}.")
                     if len(m_list) > 0:
                         # Print specific mutations to see if they are pathogenic
                         for m in m_list:
                             print(f"    - {m.get('proteinChange')} ({m.get('mutationType')}) Status: {m.get('mutationStatus')}")
                else:
                     print(f"  Unexpected response structure. Keys: {mutations.keys() if isinstance(mutations, dict) else 'Not a dict'}")
                     
            except Exception as e:
                print(f"Failed to process {study}: {e}")

if __name__ == "__main__":
    asyncio.run(audit_cbioportal_mbd4())
