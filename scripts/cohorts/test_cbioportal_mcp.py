#!/usr/bin/env python3
"""
Test cBioPortal MCP for clinical data with outcomes
Focus: Find studies with treatment outcomes for PGx validation
"""

import asyncio
import sys
sys.path.insert(0, 'tools/cbioportal-mcp')

from cbioportal_mcp.api_client import APIClient
from cbioportal_mcp.endpoints.studies import StudiesEndpoints
from cbioportal_mcp.endpoints.molecular_profiles import MolecularProfilesEndpoints
from cbioportal_mcp.endpoints.genes import GenesEndpoints

async def main():
    print("=" * 60)
    print("cBioPortal MCP: Exploring Clinical Data with Outcomes")
    print("=" * 60)
    
    # Initialize API client
    client = APIClient()
    studies = StudiesEndpoints(client)
    molecular = MolecularProfilesEndpoints(client)
    genes = GenesEndpoints(client)
    
    # Search for colorectal cancer studies (fluoropyrimidine treatment)
    print("\n1. Searching for colorectal cancer studies...")
    search_result = await studies.search_studies("colorectal", limit=20)
    
    if "error" not in search_result:
        colorectal_studies = search_result.get("studies", [])
        print(f"   Found {len(colorectal_studies)} colorectal studies")
        
        for study in colorectal_studies[:5]:
            study_id = study.get("studyId", "")
            name = study.get("name", "")[:50]
            print(f"   - {study_id}: {name}")
    
    # Get clinical data from a large colorectal study
    target_study = "coadread_tcga_pan_can_atlas_2018"
    print(f"\n2. Getting clinical data from: {target_study}")
    
    clinical_result = await molecular.get_clinical_data(
        study_id=target_study,
        limit=0  # Get all
    )
    
    if "error" not in clinical_result:
        patients = clinical_result.get("clinical_data_by_patient", {})
        print(f"   Patients found: {len(patients)}")
        
        # Check what clinical attributes are available
        all_attributes = set()
        for patient_id, attrs in patients.items():
            all_attributes.update(attrs.keys())
        
        print(f"\n3. Clinical attributes available ({len(all_attributes)} total):")
        
        # Look for outcome-related attributes
        outcome_attrs = [a for a in all_attributes if any(x in a.lower() for x in 
            ['survival', 'os', 'pfs', 'dfs', 'response', 'status', 'vital', 'death', 'outcome', 'treatment'])]
        
        print(f"   Outcome-related attributes ({len(outcome_attrs)}):")
        for attr in sorted(outcome_attrs):
            print(f"     - {attr}")
        
        # Look for PGx-related attributes
        pgx_attrs = [a for a in all_attributes if any(x in a.lower() for x in 
            ['dpyd', 'tpmt', 'ugt1a1', 'cyp', 'pgx', 'pharmacogen', 'genotype', 'metabolizer'])]
        
        if pgx_attrs:
            print(f"\n   PGx-related attributes ({len(pgx_attrs)}):")
            for attr in sorted(pgx_attrs):
                print(f"     - {attr}")
        else:
            print(f"\n   No explicit PGx attributes found")
        
        # Show sample patient data
        print(f"\n4. Sample patient data (first 3 patients):")
        for i, (patient_id, attrs) in enumerate(list(patients.items())[:3]):
            print(f"\n   Patient: {patient_id}")
            for key in ['OS_STATUS', 'OS_MONTHS', 'DFS_STATUS', 'DFS_MONTHS', 'TREATMENT_TYPE', 'SAMPLE_TYPE']:
                if key in attrs:
                    print(f"     {key}: {attrs[key]}")
    else:
        print(f"   Error: {clinical_result}")
    
    # Check for mutations in PGx genes
    print(f"\n5. Checking for DPYD mutations in colorectal studies...")
    dpyd_mutations = await genes.get_mutations_in_gene(
        study_id=target_study,
        gene_id="1806",  # DPYD Entrez ID
        limit=50
    )
    
    if "error" not in dpyd_mutations:
        mutations = dpyd_mutations.get("mutations", [])
        print(f"   DPYD mutations found: {len(mutations)}")
        if mutations:
            for mut in mutations[:5]:
                print(f"     - {mut.get('proteinChange', 'N/A')} ({mut.get('mutationType', 'N/A')})")
    else:
        print(f"   Error: {dpyd_mutations}")
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print("cBioPortal MCP provides:")
    print("  ✅ Access to clinical data with outcomes (OS, DFS, etc.)")
    print("  ✅ Mutation data for PGx genes")
    print("  ⚠️  Note: TCGA data is somatic, not germline PGx variants")
    print("\nFor true PGx validation, we need:")
    print("  - Germline variant data (not somatic)")
    print("  - Treatment-specific outcomes")
    print("  - Partner datasets with PGx screening")

if __name__ == "__main__":
    asyncio.run(main())
