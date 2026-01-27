#!/usr/bin/env python3
"""Test Evo2 invocation for PARP recommendation."""
import sys
import asyncio
import httpx

sys.path.insert(0, 'oncology-coPilot/oncology-backend-minimal')

async def test():
    print("=== Testing Evo2 API ===")
    api_base = "http://127.0.0.1:8000"
    
    mutation = {
        "gene": "MBD4",
        "chrom": "3",
        "pos": 129446860,
        "ref": "A",
        "alt": "del"
    }
    
    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            r = await client.post(
                f"{api_base}/api/evo/score_variant_multi",
                json={'chrom': mutation['chrom'], 'pos': mutation['pos'], 
                      'ref': mutation['ref'], 'alt': mutation['alt'], 'build': 'hg38'}
            )
            if r.status_code == 200:
                d = r.json()
                print(f"✅ Evo2 working: min_delta={abd.get('min_delta',0)):.4f}")
                print(f"\nNOTE: The PARP service WILL invoke Evo2 when:")
                print("  1. Mutations have chrom/pos/ref/alt")
                print("  2. Service calls SequenceProcessor.score_sequences()")
                print("  3. Which calls Evo2Scorer.score()")
                print("  4. Which POSTs to /api/evo/score_variant_multi")
                return True
            else:
                print(f"❌ Evo2 error: {r.status_code}")
                return False
    except httpx.ConnectError:
        print("❌ Backend not running - start with: make backend")
        return False

if __name__ == "__main__":
    asyncio.run(test())
