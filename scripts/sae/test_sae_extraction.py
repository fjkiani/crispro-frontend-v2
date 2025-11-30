#!/usr/bin/env python3
"""
🧪 SAE Extraction Test Script
==============================

Quick test script to verify SAE extraction with evo2_7b and trained weights.

Usage:
    python3 scripts/sae/test_sae_extraction.py

Environment Variables:
    - BACKEND_URL: Backend URL (default: http://localhost:8000)
    - ENABLE_TRUE_SAE: Must be 1 (set in backend .env)
    - ENABLE_EVO2_SAE: Must be 1 (set in backend .env)
"""

import asyncio
import httpx
import json
import os
import sys
from pathlib import Path
from loguru import logger

# Configure logger
logger.remove()
logger.add(sys.stderr, level="INFO")

# Configuration
BACKEND_URL = os.getenv("BACKEND_URL", "http://localhost:8000")
SAE_ENDPOINT = f"{BACKEND_URL}/api/sae/extract_features"

# Test variants (known oncogenic drivers)
TEST_VARIANTS = [
    {
        "name": "BRCA1 (DDR pathway)",
        "chrom": "17",
        "pos": 43044295,
        "ref": "T",
        "alt": "G",
        "assembly": "GRCh38",
        "expected_pathway": "DDR"
    },
    {
        "name": "KRAS G12D (MAPK pathway)",
        "chrom": "12",
        "pos": 25398284,
        "ref": "C",
        "alt": "A",
        "assembly": "GRCh38",
        "expected_pathway": "MAPK"
    },
    {
        "name": "BRAF V600E (MAPK pathway)",
        "chrom": "7",
        "pos": 140753336,
        "ref": "T",
        "alt": "A",
        "assembly": "GRCh38",
        "expected_pathway": "MAPK"
    }
]


async def test_sae_extraction(variant: dict) -> dict:
    """Test SAE extraction for a single variant."""
    payload = {
        "chrom": variant["chrom"],
        "pos": variant["pos"],
        "ref": variant["ref"],
        "alt": variant["alt"],
        "assembly": variant["assembly"],
        "window": 8192,
        "model_id": "evo2_7b"
    }
    
    async with httpx.AsyncClient(timeout=180.0) as client:
        try:
            logger.info(f"Testing: {variant['name']}")
            logger.info(f"  Variant: {variant['chrom']}:{variant['pos']} {variant['ref']}>{variant['alt']}")
            
            response = await client.post(SAE_ENDPOINT, json=payload)
            
            if response.status_code == 200:
                result = response.json()
                prov = result.get("provenance", {})
                
                # Check model and weights
                model = prov.get("model", "")
                d_in = prov.get("d_in", 0)
                top_features = result.get("top_features", [])
                
                logger.info(f"  ✅ Status: 200 OK")
                logger.info(f"  Model: {model}")
                logger.info(f"  d_in: {d_in}")
                logger.info(f"  Top features: {len(top_features)}")
                
                # Verify trained weights
                if "trained weights" in model.lower():
                    logger.info(f"  🎉 TRAINED WEIGHTS LOADED!")
                    weights_status = "✅ Trained"
                elif "random init" in model.lower():
                    logger.warning(f"  ⚠️  Still using random weights")
                    weights_status = "❌ Random"
                else:
                    logger.warning(f"  ⚠️  Unknown weights status")
                    weights_status = "❓ Unknown"
                
                # Verify dimension
                if d_in == 4096:
                    dim_status = "✅ Correct (4096)"
                elif d_in == 1920:
                    dim_status = "❌ Wrong (1920 - old model)"
                else:
                    dim_status = f"❓ Unexpected ({d_in})"
                
                return {
                    "variant": variant["name"],
                    "status": "success",
                    "weights": weights_status,
                    "dimension": dim_status,
                    "d_in": d_in,
                    "model": model,
                    "top_features_count": len(top_features),
                    "provenance": prov
                }
            else:
                error_text = response.text[:200] if response.text else "No error message"
                logger.error(f"  ❌ Status: {response.status_code}")
                logger.error(f"  Error: {error_text}")
                return {
                    "variant": variant["name"],
                    "status": "error",
                    "status_code": response.status_code,
                    "error": error_text
                }
        except httpx.TimeoutException:
            logger.error(f"  ❌ Timeout (180s)")
            return {
                "variant": variant["name"],
                "status": "timeout"
            }
        except Exception as e:
            logger.error(f"  ❌ Exception: {e}")
            return {
                "variant": variant["name"],
                "status": "exception",
                "error": str(e)
            }


async def main():
    """Run all tests."""
    logger.info("=" * 80)
    logger.info("🧪 SAE Extraction Test - evo2_7b with Trained Weights")
    logger.info("=" * 80)
    logger.info(f"Backend URL: {BACKEND_URL}")
    logger.info(f"SAE Endpoint: {SAE_ENDPOINT}")
    logger.info("")
    
    # Check if backend is reachable
    try:
        async with httpx.AsyncClient(timeout=5.0) as client:
            health_response = await client.get(f"{BACKEND_URL}/api/ayesha/complete_care_v2/health")
            if health_response.status_code == 200:
                logger.info("✅ Backend is reachable")
            else:
                logger.warning(f"⚠️  Backend returned status {health_response.status_code}")
    except Exception as e:
        logger.error(f"❌ Backend not reachable: {e}")
        logger.error("   Please start the backend first:")
        logger.error("   cd oncology-coPilot/oncology-backend-minimal")
        logger.error("   uvicorn api.main:app --reload")
        sys.exit(1)
    
    logger.info("")
    
    # Run tests
    results = []
    for variant in TEST_VARIANTS:
        result = await test_sae_extraction(variant)
        results.append(result)
        logger.info("")
    
    # Summary
    logger.info("=" * 80)
    logger.info("📊 TEST SUMMARY")
    logger.info("=" * 80)
    
    success_count = sum(1 for r in results if r.get("status") == "success")
    trained_count = sum(1 for r in results if "✅ Trained" in r.get("weights", ""))
    correct_dim_count = sum(1 for r in results if r.get("d_in") == 4096)
    
    logger.info(f"Total tests: {len(results)}")
    logger.info(f"Successful: {success_count}/{len(results)}")
    logger.info(f"Trained weights: {trained_count}/{success_count}")
    logger.info(f"Correct dimension (4096): {correct_dim_count}/{success_count}")
    logger.info("")
    
    # Detailed results
    for result in results:
        logger.info(f"{result['variant']}:")
        logger.info(f"  Status: {result.get('status', 'unknown')}")
        if result.get("status") == "success":
            logger.info(f"  Weights: {result.get('weights', 'unknown')}")
            logger.info(f"  Dimension: {result.get('dimension', 'unknown')}")
        else:
            logger.info(f"  Error: {result.get('error', 'unknown')}")
    
    # Final verdict
    logger.info("")
    if success_count == len(results) and trained_count == success_count and correct_dim_count == success_count:
        logger.info("✅ ALL TESTS PASSED - Ready for cohort extraction!")
        return 0
    else:
        logger.warning("⚠️  Some tests failed or using wrong configuration")
        logger.warning("   Check backend .env file and restart backend")
        return 1


if __name__ == "__main__":
    exit_code = asyncio.run(main())
    sys.exit(exit_code)

