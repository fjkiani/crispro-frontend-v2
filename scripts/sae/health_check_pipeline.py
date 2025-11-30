#!/usr/bin/env python3
"""
🔍 SAE Pipeline Health Check

Purpose: Test end-to-end Evo2 → SAE extraction pipeline
Agent: Zo (Lead Commander)
Date: January 20, 2025
"""

import asyncio
import httpx
import json
import sys
from pathlib import Path

async def test_evo2_to_sae_pipeline():
    """Test end-to-end Evo2 → SAE extraction pipeline."""
    
    print("🔍 Testing Evo2 → SAE extraction pipeline...")
    
    # Test variant: MBD4 c.1239delA (from MBD4+TP53 plan)
    test_variant = {
        "assembly": "GRCh38",
        "chrom": "3",
        "pos": 129430456,
        "ref": "A",
        "alt": "",
        "window": 8192,
        "model_id": "evo2_7b"
    }
    
    try:
        async with httpx.AsyncClient(timeout=60.0) as client:
            # Step 1: Get Evo2 activations
            print("  Step 1: Getting Evo2 activations...")
            evo_response = await client.post(
                "http://localhost:8000/api/evo/score_variant_with_activations",
                json=test_variant
            )
            
            if evo_response.status_code != 200:
                print(f"❌ ERROR: Evo2 failed: {evo_response.status_code}")
                print(f"   Response: {evo_response.text[:200]}")
                return False
            
            evo_data = evo_response.json()
            
            if "activations" not in evo_data:
                print(f"❌ ERROR: Missing activations in Evo2 response")
                return False
            
            activation_dim = len(evo_data["activations"])
            if activation_dim != 4096:
                print(f"⚠️  WARNING: Expected 4096-dim, got {activation_dim}")
            else:
                print(f"✅ Evo2 activations: {activation_dim}-dim")
            
            # Step 2: Extract SAE features
            print("  Step 2: Extracting SAE features...")
            sae_response = await client.post(
                "http://localhost:8000/api/sae/extract_features",
                json={
                    "activations": evo_data["activations"],
                    "model_id": "evo2_7b"
                }
            )
            
            if sae_response.status_code != 200:
                print(f"❌ ERROR: SAE failed: {sae_response.status_code}")
                print(f"   Response: {sae_response.text[:200]}")
                return False
            
            sae_data = sae_response.json()
            
            if "features" not in sae_data:
                print(f"❌ ERROR: Missing features in SAE response")
                return False
            
            feature_dim = len(sae_data["features"])
            if feature_dim != 32768:
                print(f"⚠️  WARNING: Expected 32768 features, got {feature_dim}")
            else:
                print(f"✅ SAE features: {feature_dim}-dim")
            
            # Check sparsity (top K=64 active)
            active_features = [f for f in sae_data["features"] if f > 0]
            active_count = len(active_features)
            
            if active_count > 64:
                print(f"⚠️  WARNING: Expected ≤64 active features, got {active_count}")
            else:
                print(f"✅ Active features: {active_count} (sparse)")
            
            # Check provenance
            if "provenance" in sae_data:
                provenance = sae_data["provenance"]
                print(f"✅ Provenance: {provenance.get('model', 'unknown')}, weights: {provenance.get('weights', 'unknown')}")
            
            print("✅ Evo2 → SAE pipeline working")
            return True
            
    except httpx.TimeoutException:
        print("❌ ERROR: Request timed out (60s)")
        return False
    except Exception as e:
        print(f"❌ ERROR: Pipeline test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = asyncio.run(test_evo2_to_sae_pipeline())
    sys.exit(0 if success else 1)

