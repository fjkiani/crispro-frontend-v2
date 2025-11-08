#!/usr/bin/env python3
"""
Test 3: Baseline vs Intervention - Core Hypothesis Test

CRITICAL TEST: Can Evo2 differentiate baseline (disease) vs intervention (compound) states?

This test determines if our approach is viable.
"""

import httpx
import asyncio
import json
import os

API_BASE = os.getenv("API_BASE", "http://127.0.0.1:8000")

async def fetch_gene_sequence(gene_symbol: str) -> str:
    """Fetch gene sequence from Ensembl REST API."""
    
    async with httpx.AsyncClient(timeout=30.0) as client:
        # Step 1: Lookup gene ID
        lookup_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}?content-type=application/json;expand=0"
        
        try:
            lookup_resp = await client.get(lookup_url)
            if lookup_resp.status_code == 200:
                gene_data = lookup_resp.json()
                gene_id = gene_data.get('id')
                
                if gene_id:
                    # Step 2: Get sequence
                    seq_url = f"https://rest.ensembl.org/sequence/id/{gene_id}?content-type=application/json"
                    seq_resp = await client.get(seq_url)
                    
                    if seq_resp.status_code == 200:
                        seq_data = seq_resp.json()
                        return seq_data.get('seq', '')
        except Exception as e:
            print(f"⚠️  Error fetching {gene_symbol}: {e}")
    
    return ""

async def test_baseline_vs_intervention():
    """Test if Evo2 can differentiate baseline vs intervention using same sequence."""
    
    print("🧪 TEST 3: Baseline vs Intervention - Core Hypothesis\n")
    
    # Test compound: Vitamin D
    compound = "Vitamin D"
    target_gene = "VDR"
    
    print(f"📋 Test Setup:")
    print(f"   Compound: {compound}")
    print(f"   Target Gene: {target_gene}\n")
    
    # Step 1: Fetch VDR sequence
    print(f"📥 Fetching {target_gene} sequence from Ensembl...")
    vdr_sequence = await fetch_gene_sequence(target_gene)
    
    if not vdr_sequence:
        print(f"❌ FAILED: Could not fetch {target_gene} sequence")
        return False
    
    print(f"✅ Got sequence: {len(vdr_sequence)}bp")
    print(f"   First 50bp: {vdr_sequence[:50]}...\n")
    
    # Use first 1000bp for testing (Evo2 may have length limits)
    test_sequence = vdr_sequence[:1000]
    
    async with httpx.AsyncClient(timeout=120.0) as client:
        # Step 2: Score baseline (disease-active state)
        print(f"📤 Scoring baseline (disease-active state)...")
        baseline_resp = await client.post(
            f"{API_BASE}/api/evo/score",
            json={
                "sequence": test_sequence,
                "model_id": "evo2_1b"
            }
        )
        
        if baseline_resp.status_code != 200:
            print(f"❌ Baseline scoring failed: {baseline_resp.status_code}")
            print(f"   {baseline_resp.text}")
            return False
        
        baseline_data = baseline_resp.json()
        baseline_score = baseline_data.get('score', baseline_data.get('likelihood', 0))
        
        print(f"✅ Baseline score: {baseline_score}")
        
        # Step 3: Score intervention (same sequence!)
        # ⚠️ CRITICAL: If Evo2 only scores sequences (not context), these will be identical
        print(f"\n📤 Scoring intervention ({compound} active state)...")
        print(f"   ⚠️  Using SAME sequence (this is the test!)")
        
        intervention_resp = await client.post(
            f"{API_BASE}/api/evo/score",
            json={
                "sequence": test_sequence,  # SAME sequence
                "model_id": "evo2_1b"
            }
        )
        
        if intervention_resp.status_code != 200:
            print(f"❌ Intervention scoring failed: {intervention_resp.status_code}")
            print(f"   {intervention_resp.text}")
            return False
        
        intervention_data = intervention_resp.json()
        intervention_score = intervention_data.get('score', intervention_data.get('likelihood', 0))
        
        print(f"✅ Intervention score: {intervention_score}")
        
        # Step 4: Compute delta
        delta = abs(baseline_score - intervention_score)
        
        print(f"\n📊 RESULTS:")
        print(f"   Baseline: {baseline_score}")
        print(f"   Intervention: {intervention_score}")
        print(f"   Delta: {delta}\n")
        
        # Analysis
        if delta == 0:
            print("❌ CRITICAL FINDING: Delta = 0")
            print("   ⚠️  Same sequence → Same score")
            print("   ⚠️  Evo2 cannot differentiate context using sequence-only approach")
            print("\n📝 CONCLUSION:")
            print("   ❌ This approach DOES NOT WORK")
            print("   ✅ Need alternative: Test 6 (variant scoring) or Test 7 (generation)")
            return False
        elif delta < 0.01:
            print("⚠️  WARNING: Delta very small (< 0.01)")
            print("   May be noise, not meaningful biological difference")
            print("\n📝 CONCLUSION:")
            print("   ⚠️  Approach is questionable - delta too small")
            return False
        else:
            print("✅ Delta is meaningful (> 0.01)")
            print("   ✅ Evo2 CAN differentiate states (somehow)")
            print("\n📝 CONCLUSION:")
            print("   ✅ Approach MIGHT WORK (but need to understand HOW)")
            print("   🤔 Question: Why did delta > 0 if sequence was identical?")
            print("      - Non-deterministic scoring?")
            print("      - Evo2 has internal context?")
            print("      - Need to investigate further")
            return True

if __name__ == "__main__":
    success = asyncio.run(test_baseline_vs_intervention())
    
    if success:
        print("\n✅ TEST 3: APPROACH VIABLE (with caveats)")
    else:
        print("\n❌ TEST 3: APPROACH NOT VIABLE")
        print("\n⚠️  NEXT STEPS:")
        print("   1. Run Test 6: Variant scoring alternative")
        print("   2. Run Test 7: Generation alternative")
        print("   3. Reconsider approach if both fail")

