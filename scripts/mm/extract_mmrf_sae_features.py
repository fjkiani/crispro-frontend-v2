#!/usr/bin/env python3
"""
🧬 MM RESISTANCE PREDICTION - Work Item 1: Extract MMRF TRUE SAE Features
==========================================================================

Mission: Extract TRUE SAE features for MMRF CoMMpass cohort
Agent: Plumber (Implementation)
Date: January 29, 2025
Priority: P0 (Critical)

Objective:
- Load MMRF CoMMpass cohort (from download_mmrf_commpass.py)
- Extract TRUE SAE features using same pipeline as OV
- Use Modal SAE service (same as OV)
- Aggregate to patient-level (mean/max)
- Save to: data/validation/mm_cohort/mmrf_sae_features.json

Success Criteria:
- TRUE SAE features extracted for ≥500 patients
- Patient-level aggregation complete
- Compatible with validation pipeline

Guardrails:
- RUO/validation-only
- Requires ENABLE_EVO2_SAE=1 and ENABLE_TRUE_SAE=1
- Cost control: MAX_PATIENTS, MAX_TOTAL_VARIANTS
"""

import asyncio
import httpx
import json
import sys
import os
from pathlib import Path
from typing import List, Dict, Any, Optional
from datetime import datetime
from collections import defaultdict
import numpy as np
from loguru import logger

# Configure logger
logger.remove()
logger.add(sys.stderr, level="INFO")

# Paths
MMRF_DATA_FILE = Path("data/validation/mm_cohort/mmrf_commpass.json")
OUTPUT_DIR = Path("data/validation/mm_cohort")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUTPUT_DIR / "mmrf_sae_features.json"
CHECKPOINT_FILE = OUTPUT_DIR / "mmrf_sae_features_extraction_checkpoint.json"

# API Configuration
BACKEND_URL = os.getenv("BACKEND_URL", "http://localhost:8000")
SAE_FEATURES_ENDPOINT = f"{BACKEND_URL}/api/sae/extract_features"

# Extraction Configuration (COST CONTROL)
MAX_VARIANTS_PER_PATIENT = 50
BATCH_SIZE = 10
MAX_PATIENTS = int(os.getenv("MAX_PATIENTS", "500"))
MAX_TOTAL_VARIANTS = int(os.getenv("MAX_TOTAL_VARIANTS", "25000"))
REQUEST_TIMEOUT = 180.0
MAX_RETRIES = 3
RETRY_DELAY = 5.0

# Circuit breaker
MAX_ERROR_RATE = 0.30
MIN_CALLS_BEFORE_CHECK = 20

# Safety: require explicit flag
ENABLE_SAE_COHORT_RUN = os.getenv("ENABLE_SAE_COHORT_RUN", "0").lower() in ("true", "1", "yes")

def check_feature_flags() -> bool:
    """Check if required feature flags are enabled."""
    enable_evo2_sae = os.getenv("ENABLE_EVO2_SAE", "0").lower() in ("true", "1", "yes")
    enable_true_sae = os.getenv("ENABLE_TRUE_SAE", "0").lower() in ("true", "1", "yes")
    
    if not enable_evo2_sae:
        logger.error("Environment variable ENABLE_EVO2_SAE is not set to '1' or 'true'.")
    if not enable_true_sae:
        logger.error("Environment variable ENABLE_TRUE_SAE is not set to '1' or 'true'.")
    if not ENABLE_SAE_COHORT_RUN:
        logger.error("Environment variable ENABLE_SAE_COHORT_RUN is not set to '1' or 'true'.")
        logger.info("To run, set ENABLE_SAE_COHORT_RUN=1 ENABLE_EVO2_SAE=1 ENABLE_TRUE_SAE=1 python scripts/mm/extract_mmrf_sae_features.py")
    
    return enable_evo2_sae and enable_true_sae and ENABLE_SAE_COHORT_RUN

def load_mmrf_data() -> Dict[str, Any]:
    """Load MMRF CoMMpass patient data."""
    if not MMRF_DATA_FILE.exists():
        logger.error(f"❌ MMRF data file not found: {MMRF_DATA_FILE}. Please run download_mmrf_commpass.py first.")
        sys.exit(1)
    
    logger.info(f"📥 Loading MMRF data from {MMRF_DATA_FILE}")
    with open(MMRF_DATA_FILE, 'r') as f:
        data = json.load(f)
    
    patients = data.get('patients', [])
    logger.info(f"   ✅ Loaded {len(patients)} patients")
    return data

async def extract_sae_features_for_variant(
    client: httpx.AsyncClient,
    variant: Dict[str, Any],
    model_id: str = "evo2_7b"
) -> Optional[Dict[str, Any]]:
    """Extract SAE features for a single variant."""
    payload = {
        "chrom": variant.get("chrom"),
        "pos": variant.get("pos"),
        "ref": variant.get("ref"),
        "alt": variant.get("alt"),
        "gene": variant.get("gene"),
        "hgvs_p": variant.get("hgvs_p"),
        "model_id": model_id,
        "assembly": "GRCh38"
    }
    
    for attempt in range(MAX_RETRIES):
        try:
            response = await client.post(SAE_FEATURES_ENDPOINT, json=payload, timeout=REQUEST_TIMEOUT)
            response.raise_for_status()
            result = response.json()
            
            if 'features' not in result or not isinstance(result['features'], (list, dict)):
                logger.warning(f"      ⚠️  SAE features missing or malformed for {variant.get('gene')}")
                return None

            return {
                "variant": variant,
                "sae_features": result["features"],
                "provenance": result.get("provenance", {})
            }
        except Exception as e:
            logger.error(f"      ❌ Error for {variant.get('gene')}: {e}")
            if attempt < MAX_RETRIES - 1:
                await asyncio.sleep(RETRY_DELAY)
            else:
                return None
    return None

async def extract_cohort_sae_features():
    """Main extraction pipeline for MMRF CoMMpass cohort."""
    logger.info("=" * 80)
    logger.info("🧬 MMRF CoMMpass SAE FEATURE EXTRACTION")
    logger.info("=" * 80)
    
    if not check_feature_flags():
        logger.error("❌ Feature flags not enabled. Exiting.")
        sys.exit(1)
    
    mmrf_data = load_mmrf_data()
    all_patients = mmrf_data.get("patients", [])

    patients_with_mutations = [p for p in all_patients if p.get("mutations")]
    logger.info(f"📊 Patients with mutations: {len(patients_with_mutations)}")
    logger.info(f"   Processing limit: {MAX_PATIENTS}")

    patients_to_process = patients_with_mutations[:MAX_PATIENTS]
    
    logger.info(f"✅ Starting extraction for {len(patients_to_process)} patients")
    logger.info("⚠️  Full implementation requires Modal SAE service integration")
    logger.info("   See scripts/sae/extract_sae_features_cohort.py for reference")

if __name__ == "__main__":
    asyncio.run(extract_cohort_sae_features())
