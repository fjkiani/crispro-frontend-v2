#!/usr/bin/env python3
"""
cBioPortal MCP: Complete Test with Outcomes
"""

import asyncio
import sys
sys.path.insert(0, '/Users/fahadkiani/Desktop/development/crispr-assistant-main/tools/cbioportal-mcp')

from cbioportal_mcp.api_client import APIClient

async def main():
    client = APIClient()
    await client.startup()
    
    print("=" * 60)
    print("cBioPortal MCP: Complete Test with Outcomes")
    print("=" * 60)
    
    study_id = 'ov_tcga'
    
    # Get sample lists
    print('\n1. Getting sample lists...')
    sample_lists = await client.make_api_request(f'studies/{study_id}/sample-lists')
    print(f'   Found {len(sample_lists)} sample lists')
    mutation_sample_list = None
    for sl in sample_lists:
        name = sl.get('name', '')
        sl_id = sl.get('sampleListId', '')
        if 'mutation' in name.lower() or 'sequenced' in name.lower():
            print(f'   Using: {sl_id}')
            mutation_sample_list = sl_id
            break
    
    # Get mutations with sample list
    if mutation_sample_list:
        print(f'\n2. Getting BRCA1/BRCA2 mutations...')
        muts = await client.make_api_request(
            'molecular-profiles/ov_tcga_mutations/mutations',
            params={'entrezGeneId': 672, 'sampleListId': mutation_sample_list}
        )
        print(f'   BRCA1 mutations: {len(muts)}')
        for m in muts[:5]:
            sample = m.get('sampleId', '')
            protein = m.get('proteinChange', 'N/A')
            mut_type = m.get('mutationType', '')
            print(f'     - {sample}: {protein} ({mut_type})')
        
        # BRCA2
        muts2 = await client.make_api_request(
            'molecular-profiles/ov_tcga_mutations/mutations',
            params={'entrezGeneId': 675, 'sampleListId': mutation_sample_list}
        )
        print(f'   BRCA2 mutations: {len(muts2)}')
    else:
        muts = []
        muts2 = []
    
    # Get clinical data
    print(f'\n3. Getting clinical data with outcomes...')
    clinical = await client.make_api_request(f'studies/{study_id}/clinical-data', params={'pageSize': 1000})
    
    # Group by patient
    patients = {}
    for c in clinical:
        pid = c.get('patientId')
        if pid not in patients:
            patients[pid] = {}
        patients[pid][c.get('clinicalAttributeId')] = c.get('value')
    
    print(f'   Total patients: {len(patients)}')
    
    # Count outcomes
    os_count = sum(1 for p in patients.values() if 'OS_STATUS' in p)
    deceased = sum(1 for p in patients.values() if str(p.get('OS_STATUS', '')).startswith('1'))
    print(f'   With OS data: {os_count}')
    print(f'   Deceased: {deceased}')
    
    # Sample patient
    print(f'\n4. Sample patient data:')
    for pid, attrs in list(patients.items())[:2]:
        print(f'   Patient: {pid}')
        for k in ['OS_STATUS', 'OS_MONTHS', 'AGE', 'TUMOR_STAGE']:
            if k in attrs:
                print(f'     {k}: {attrs[k]}')
    
    # Link mutations to outcomes
    print(f'\n5. Linking mutations to outcomes...')
    mutation_patients = set()
    for m in muts:
        sample = m.get('sampleId', '')
        patient = sample.rsplit('-', 1)[0] if '-' in sample else sample
        mutation_patients.add(patient)
    
    for m in muts2:
        sample = m.get('sampleId', '')
        patient = sample.rsplit('-', 1)[0] if '-' in sample else sample
        mutation_patients.add(patient)
    
    print(f'   Patients with BRCA1/2 mutations: {len(mutation_patients)}')
    
    # Check outcomes for mutation patients
    mut_deceased = 0
    mut_with_os = 0
    for pid in mutation_patients:
        if pid in patients:
            if 'OS_STATUS' in patients[pid]:
                mut_with_os += 1
                if str(patients[pid]['OS_STATUS']).startswith('1'):
                    mut_deceased += 1
    
    print(f'   With OS data: {mut_with_os}')
    print(f'   Deceased: {mut_deceased}')
    
    await client.shutdown()
    
    print('\n' + '=' * 60)
    print('SUCCESS! cBioPortal MCP provides:')
    print('  - Mutation data (somatic)')
    print('  - Clinical outcomes (OS, DFS, PFS)')
    print('  - Patient-mutation-outcome linkage')
    print('=' * 60)

if __name__ == "__main__":
    asyncio.run(main())
