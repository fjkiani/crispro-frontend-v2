#!/usr/bin/env python3
"""
Test 1: Evo2 API Format Discovery

Goal: Understand what /api/evo/score actually returns
"""

import httpx
import asyncio
import json
import os

API_BASE = os.getenv("API_BASE", "http://127.0.0.1:8000")

async def test_evo2_api_format():
    """Test basic Evo2 score endpoint to understand response format."""
    
    print("🧪 TEST 1: Evo2 API Format Discovery\n")
    
    async with httpx.AsyncClient(timeout=60.0) as client:
        # Test with a simple DNA sequence
        test_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        
        print(f"📤 Request: POST {API_BASE}/api/evo/score")
        print(f"   Sequence: {test_sequence[:20]}...")
        print(f"   Model: evo2_1b\n")
        
        try:
            response = await client.post(
                f"{API_BASE}/api/evo/score",
                json={
                    "sequence": test_sequence,
                    "model_id": "evo2_1b"
                }
            )
            
            print(f"📥 Response Status: {response.status_code}")
            
            if response.status_code == 200:
                data = response.json()
                print(f"📥 Response Body:")
                print(json.dumps(data, indent=2))
                
                # Analyze response structure
                print("\n🔍 Analysis:")
                print(f"   Keys: {list(data.keys())}")
                
                # Look for score/likelihood fields
                if "score" in data:
                    score = data["score"]
                    print(f"   ✅ Found 'score': {score} (type: {type(score).__name__})")
                    print(f"   📊 Score range: Need to determine range")
                elif "likelihood" in data:
                    likelihood = data["likelihood"]
                    print(f"   ✅ Found 'likelihood': {likelihood} (type: {type(likelihood).__name__})")
                    print(f"   📊 Likelihood range: Need to determine range")
                else:
                    print(f"   ⚠️  No 'score' or 'likelihood' field found!")
                    print(f"   📋 Available fields: {list(data.keys())}")
                
                return data
            else:
                print(f"❌ Error: {response.status_code}")
                print(f"   {response.text}")
                return None
                
        except Exception as e:
            print(f"❌ Exception: {e}")
            return None

if __name__ == "__main__":
    result = asyncio.run(test_evo2_api_format())
    
    if result:
        print("\n✅ TEST 1 PASSED: Got response from Evo2 API")
        print("\n📝 Next Steps:")
        print("   1. Note the response format")
        print("   2. Determine score range (0-1? -inf to inf?)")
        print("   3. Proceed to Test 2")
    else:
        print("\n❌ TEST 1 FAILED: Could not get response from Evo2 API")
        print("\n⚠️  BLOCKING: Need to fix Evo2 API connection before proceeding")

