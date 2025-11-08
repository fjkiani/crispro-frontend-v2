#!/usr/bin/env python3
"""
Test 2: Ensembl Gene Sequence Fetching

Goal: Verify we can reliably fetch gene sequences
"""

import httpx
import asyncio

async def test_ensembl_fetch():
    """Test fetching gene sequences from Ensembl."""
    
    print("🧪 TEST 2: Ensembl Gene Sequence Fetching\n")
    
    test_genes = ["VDR", "TP53", "NFKB1"]
    
    async with httpx.AsyncClient(timeout=30.0) as client:
        for gene in test_genes:
            print(f"📥 Fetching {gene}...")
            
            try:
                # Step 1: Lookup
                lookup_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}?content-type=application/json;expand=0"
                lookup_resp = await client.get(lookup_url)
                
                if lookup_resp.status_code != 200:
                    print(f"   ❌ Lookup failed: {lookup_resp.status_code}")
                    continue
                
                gene_data = lookup_resp.json()
                gene_id = gene_data.get('id')
                
                if not gene_id:
                    print(f"   ❌ No gene ID found")
                    continue
                
                print(f"   ✅ Gene ID: {gene_id}")
                
                # Step 2: Get sequence
                seq_url = f"https://rest.ensembl.org/sequence/id/{gene_id}?content-type=application/json"
                seq_resp = await client.get(seq_url)
                
                if seq_resp.status_code != 200:
                    print(f"   ❌ Sequence fetch failed: {seq_resp.status_code}")
                    continue
                
                seq_data = seq_resp.json()
                sequence = seq_data.get('seq', '')
                
                if not sequence:
                    print(f"   ❌ No sequence returned")
                    continue
                
                print(f"   ✅ Sequence: {len(sequence)}bp")
                print(f"   ✅ First 50bp: {sequence[:50]}")
                print()
                
            except Exception as e:
                print(f"   ❌ Error: {e}")
                print()
    
    print("✅ TEST 2 COMPLETE: Ensembl fetching tested")

if __name__ == "__main__":
    asyncio.run(test_ensembl_fetch())

