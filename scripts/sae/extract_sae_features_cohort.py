#!/usr/bin/env python3
"""
🧬 SAE PHASE 2: COHORT FEATURE EXTRACTION
==========================================

Mission: Extract SAE features for TCGA-OV platinum response cohort
Agent: Zo (Lead Commander)
Date: January 14, 2025
Priority: HIGH (Phase 2, Task 1)

Objective:
- Load TCGA-OV platinum response labels (from extract_platinum_response_labels.py)
- For each patient with labeled variants:
  - Extract Evo2 layer 26 activations via /api/evo/score_variant_with_activations
  - Extract SAE features (32K-dim) via /api/sae/extract_features
  - Aggregate per-patient SAE features (mean, max, top features)
- Save cohort file: data/validation/sae_features_tcga_ov_platinum.json

Success Criteria:
- Process N≈200 fully labeled TCGA-OV patients
- Extract SAE features for ~30-50 variants per patient
- Wall-time: ~10-20 minutes with modest parallelism
- Output structured JSON for correlation analysis

Guardrails:
- RUO/validation-only (no production impact)
- Requires ENABLE_EVO2_SAE=1 and ENABLE_TRUE_SAE=1
- All features tagged with provenance
"""

import asyncio
import httpx
import json
import sys
import time
import os
from pathlib import Path
from typing import List, Dict, Optional, Any
from datetime import datetime
from collections import defaultdict
import numpy as np
from loguru import logger

# Configure logger
logger.remove()
logger.add(sys.stderr, level="INFO")

# Paths
PLATINUM_WITH_MUTATIONS_FILE = Path("data/validation/sae_cohort/tcga_ov_platinum_with_mutations.json")
OUTPUT_DIR = Path("data/validation/sae_cohort")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUTPUT_DIR / "sae_features_tcga_ov_platinum.json"
CHECKPOINT_FILE = OUTPUT_DIR / "sae_features_extraction_checkpoint.json"

# API Configuration
BACKEND_URL = os.getenv("BACKEND_URL", "http://localhost:8000")
EVO_ACTIVATIONS_ENDPOINT = f"{BACKEND_URL}/api/evo/score_variant_with_activations"
SAE_FEATURES_ENDPOINT = f"{BACKEND_URL}/api/sae/extract_features"

# Extraction Configuration (COST CONTROL)
MAX_VARIANTS_PER_PATIENT = 50  # Limit variants per patient for initial validation
BATCH_SIZE = 10  # Process variants in batches
MAX_PATIENTS = int(os.getenv("MAX_PATIENTS", "50"))  # Default to 50; set env var to override
MAX_TOTAL_VARIANTS = int(os.getenv("MAX_TOTAL_VARIANTS", "2500"))  # Hard cap on total variants processed
REQUEST_TIMEOUT = 180.0  # 3 minutes per request
MAX_RETRIES = 3
RETRY_DELAY = 5.0

# Circuit breaker: stop if error rate exceeds threshold
MAX_ERROR_RATE = 0.30  # Stop if >30% of calls fail
MIN_CALLS_BEFORE_CHECK = 20  # Wait for at least 20 calls before checking error rate

# Safety: require explicit flag to run cohort extraction
ENABLE_SAE_COHORT_RUN = os.getenv("ENABLE_SAE_COHORT_RUN", "0").lower() in ("true", "1", "yes")

# Optional: wire in pyBioPortal to fetch real TCGA‑OV mutations when variants
# are not already embedded in the labels file.
PROJECT_ROOT = Path(__file__).resolve().parent.parent  # /scripts
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

PYBIOPORTAL_PATH = PROJECT_ROOT / "oncology-coPilot" / "oncology-backend" / "tests" / "pyBioPortal-master"
if PYBIOPORTAL_PATH.exists() and str(PYBIOPORTAL_PATH) not in sys.path:
    sys.path.insert(0, str(PYBIOPORTAL_PATH))

try:
    from pybioportal import molecular_profiles as mp  # type: ignore
    from pybioportal import sample_lists as sl  # type: ignore
    from pybioportal import mutations as mut  # type: ignore
    # Reuse helper for TCGA patient IDs
    from extract_platinum_response_labels import extract_tcga_patient_id  # type: ignore
except Exception:
    mp = sl = mut = None  # type: ignore


def build_patient_variants_from_pybioportal(
    study_id: str = "ov_tcga_pan_can_atlas_2018",
    max_page_size: int = 10000,
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Build a mapping of TCGA patient_id → list of variant dicts by querying cBioPortal
    via the local pyBioPortal client.

    This is RUO / validation‑only and is only used if pybioportal is available.
    """
    variants_by_patient: Dict[str, List[Dict[str, Any]]] = defaultdict(list)

    if mp is None or sl is None or mut is None:
        logger.warning("pybioportal not available; skipping variant extraction from cBioPortal.")
        return {}

    try:
        logger.info(f"📡 Fetching mutation profile metadata for {study_id} via pyBioPortal...")
        df_prof = mp.get_all_molecular_profiles_in_study(study_id)
        profile_id: Optional[str] = None

        if isinstance(df_prof, np.ndarray):
            # defensive: pybioportal usually returns DataFrame, but guard anyway
            logger.warning("Unexpected numpy array from pybioportal.molecular_profiles; skipping.")
            return {}

        import pandas as pd  # local import to avoid hard dependency at module import

        if isinstance(df_prof, pd.DataFrame) and "molecularProfileId" in df_prof.columns:
            for pid in df_prof["molecularProfileId"].astype(str).tolist():
                if pid.endswith("_mutations"):
                    profile_id = pid
                    break

        if not profile_id:
            logger.error(f"❌ No mutation profile found for {study_id}")
            return {}

        logger.info(f"   ✅ Using mutation profile: {profile_id}")

        df_lists = sl.get_all_sample_lists_in_study(study_id)
        sample_list_id: Optional[str] = None

        if isinstance(df_lists, pd.DataFrame) and "sampleListId" in df_lists.columns:
            sample_list_candidates = df_lists["sampleListId"].astype(str).tolist()
            for sid in sample_list_candidates:
                if sid.endswith("_all"):
                    sample_list_id = sid
                    break
            if sample_list_id is None and sample_list_candidates:
                sample_list_id = sample_list_candidates[0]

        logger.info(f"   ✅ Using sample list: {sample_list_id}")

        df_muts = mut.get_muts_in_mol_prof_by_sample_list_id(
            profile_id,
            sample_list_id,
            projection="DETAILED",
            pageSize=max_page_size,
        )

        if not isinstance(df_muts, pd.DataFrame) or df_muts.empty:
            logger.error("❌ No mutations returned from pyBioPortal")
            return {}

        # Determine gene column
        gene_cols = [c for c in df_muts.columns if "hugo" in c.lower() or c.lower() == "gene"]
        if not gene_cols:
            logger.error("❌ Could not find gene column in mutation dataframe")
            return {}
        gene_col = gene_cols[0]

        logger.info(f"   📊 Retrieved {len(df_muts)} mutation rows; building patient→variants map...")

        for _, row in df_muts.iterrows():
            sample_id = str(row.get("sampleId") or row.get("sample_id") or "")
            if not sample_id:
                continue
            patient_id = extract_tcga_patient_id(sample_id) or sample_id

            chrom = row.get("chromosome") or row.get("chr") or row.get("chrom")
            pos = row.get("startPosition") or row.get("start_position") or row.get("position")
            ref = row.get("referenceAllele") or row.get("ref")
            alt = row.get("tumorSeqAllele2") or row.get("alt")

            if not (chrom and pos and ref and alt):
                continue

            variant: Dict[str, Any] = {
                "gene": row[gene_col],
                "chrom": str(chrom).replace("chr", ""),
                "pos": int(pos),
                "ref": str(ref).upper(),
                "alt": str(alt).upper(),
                "hgvs_p": row.get("proteinChange") or row.get("hgvsp"),
            }
            variants_by_patient[patient_id].append(variant)

        logger.info(f"   ✅ Built variants for {len(variants_by_patient)} TCGA patients from pyBioPortal.")
        return variants_by_patient

    except Exception as e:
        logger.error(f"❌ Error while building patient variants from pyBioPortal: {e}")
        return {}

# Feature Flags Check
def check_feature_flags():
    """Verify that required feature flags are enabled."""
    enable_evo2_sae = os.getenv("ENABLE_EVO2_SAE", "false").lower() in ("true", "1", "yes")
    enable_true_sae = os.getenv("ENABLE_TRUE_SAE", "false").lower() in ("true", "1", "yes")
    
    if not enable_evo2_sae:
        logger.error("❌ ENABLE_EVO2_SAE=0. Set ENABLE_EVO2_SAE=1 to enable Evo2 activations endpoint.")
        return False
    
    if not enable_true_sae:
        logger.error("❌ ENABLE_TRUE_SAE=0. Set ENABLE_TRUE_SAE=1 to enable true SAE features.")
        return False
    
    if not ENABLE_SAE_COHORT_RUN:
        logger.error("❌ ENABLE_SAE_COHORT_RUN=0. This is a safety check to prevent accidental credit burn.")
        logger.error("   Set ENABLE_SAE_COHORT_RUN=1 explicitly when you want to run the cohort extraction.")
        return False
    
    logger.info("✅ Feature flags enabled: ENABLE_EVO2_SAE=1, ENABLE_TRUE_SAE=1, ENABLE_SAE_COHORT_RUN=1")
    logger.info(f"✅ Cost limits: MAX_PATIENTS={MAX_PATIENTS}, MAX_TOTAL_VARIANTS={MAX_TOTAL_VARIANTS}")
    return True

def load_platinum_labels() -> Dict[str, Any]:
    """Load TCGA-OV platinum response labels with mutations."""
    if not PLATINUM_WITH_MUTATIONS_FILE.exists():
        logger.error(f"❌ Platinum labels with mutations file not found: {PLATINUM_WITH_MUTATIONS_FILE}")
        logger.error(f"   Run scripts/sae/extract_patient_mutations_for_cohort.py first")
        sys.exit(1)
    
    logger.info(f"📥 Loading platinum labels with mutations from {PLATINUM_WITH_MUTATIONS_FILE}")
    with open(PLATINUM_WITH_MUTATIONS_FILE, 'r') as f:
        data = json.load(f)
    
    # The new file structure is a direct list of patient dicts (not wrapped in {"patients": [...]}])
    if isinstance(data, list):
        patients = data
        logger.info(f"   ✅ Loaded {len(patients)} patients directly from list")
        return {"patients": patients}  # Wrap for compatibility
    else:
        # Fallback if it's still a dict with 'patients' key
        patients = data.get("patients", [])
        logger.info(f"   ✅ Loaded {len(patients)} patients from dict")
        return data

def load_checkpoint() -> Dict[str, Any]:
    """Load extraction checkpoint if it exists."""
    if CHECKPOINT_FILE.exists():
        logger.info(f"📥 Loading checkpoint from {CHECKPOINT_FILE}")
        with open(CHECKPOINT_FILE, 'r') as f:
            checkpoint = json.load(f)
        logger.info(f"   ✅ Checkpoint loaded: {len(checkpoint.get('completed_patients', []))} patients completed")
        return checkpoint
    
    return {
        "completed_patients": [],
        "failed_patients": [],
        "extraction_start": datetime.now().isoformat(),
        "last_checkpoint": None
    }

def save_checkpoint(checkpoint: Dict[str, Any]):
    """Save extraction checkpoint."""
    checkpoint["last_checkpoint"] = datetime.now().isoformat()
    with open(CHECKPOINT_FILE, 'w') as f:
        json.dump(checkpoint, f, indent=2)

async def extract_sae_features_for_variant(
    client: httpx.AsyncClient,
    variant: Dict[str, Any],
    model_id: str = "evo2_7b"
) -> Optional[Dict[str, Any]]:
    """
    Extract SAE features for a single variant.
    
    Args:
        client: httpx AsyncClient
        variant: Variant dict with chrom, pos, ref, alt
        model_id: Evo2 model to use
    
    Returns:
        {
            "variant": {...},
            "sae_features": [...],  # 32K-dim
            "top_features": [...],  # Top k=64
            "stats": {...},
            "provenance": {...}
        }
    """
    chrom = variant.get("chrom") or variant.get("chromosome")
    pos = variant.get("pos") or variant.get("position")
    ref = variant.get("ref") or variant.get("reference_allele")
    alt = variant.get("alt") or variant.get("alternate_allele")
    
    if not all([chrom, pos, ref, alt]):
        logger.warning(f"   ⚠️  Skipping variant with missing coordinates: {variant}")
        return None
    
    # Call SAE endpoint (which internally calls Evo2 + SAE)
    # NOTE: TCGA-OV mutations from cBioPortal are on GRCh37/hg19.
    # Using GRCh38 here causes Ensembl 400 errors when positions exceed chr length.
    payload = {
        "chrom": str(chrom).replace("chr", ""),
        "pos": int(pos),
        "ref": str(ref).upper(),
        "alt": str(alt).upper(),
        "model_id": model_id,
        "assembly": "GRCh37",
        "window": 8192,
    }
    
    for attempt in range(MAX_RETRIES):
        try:
            response = await client.post(
                SAE_FEATURES_ENDPOINT,
                json=payload,
                headers={"Content-Type": "application/json"},
                timeout=REQUEST_TIMEOUT
            )
            
            if response.status_code == 200:
                result = response.json()
                return {
                    "variant": {
                        "chrom": payload["chrom"],
                        "pos": payload["pos"],
                        "ref": payload["ref"],
                        "alt": payload["alt"],
                        "gene": variant.get("gene"),
                        "hgvs_p": variant.get("hgvs_p")
                    },
                    "sae_features": result.get("features", []),
                    "top_features": result.get("top_features", []),
                    "stats": result.get("stats", {}),
                    "provenance": result.get("provenance", {})
                }
            elif response.status_code == 403:
                logger.error(f"❌ Feature flags not enabled: {response.text}")
                return None
            elif response.status_code == 503:
                logger.error(f"❌ SAE service not configured: {response.text}")
                return None
            else:
                logger.warning(f"   ⚠️  HTTP {response.status_code} for {chrom}:{pos} {ref}>{alt}: {response.text}")
                if attempt < MAX_RETRIES - 1:
                    await asyncio.sleep(RETRY_DELAY)
                    continue
                return None
        
        except httpx.TimeoutException:
            logger.warning(f"   ⏱️  Timeout for {chrom}:{pos} {ref}>{alt} (attempt {attempt + 1}/{MAX_RETRIES})")
            if attempt < MAX_RETRIES - 1:
                await asyncio.sleep(RETRY_DELAY)
                continue
            return None
        
        except Exception as e:
            logger.error(f"   ❌ Error extracting SAE features for {chrom}:{pos} {ref}>{alt}: {e}")
            if attempt < MAX_RETRIES - 1:
                await asyncio.sleep(RETRY_DELAY)
                continue
            return None
    
    return None

def aggregate_patient_sae_features(variant_features: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Aggregate SAE features across variants for a single patient.
    
    Aggregation methods:
    - Mean activation per SAE feature (32K-dim)
    - Max activation per SAE feature (32K-dim)
    - Top-k most activated features across all variants
    - Sparsity statistics
    
    Args:
        variant_features: List of variant-level SAE feature dicts
    
    Returns:
        {
            "num_variants": int,
            "mean_features": List[float],  # 32K-dim
            "max_features": List[float],   # 32K-dim
            "top_features_aggregated": List[{"index": int, "value": float}],
            "sparsity_mean": float,
            "sparsity_std": float
        }
    """
    if not variant_features:
        return {
            "num_variants": 0,
            "mean_features": [],
            "max_features": [],
            "top_features_aggregated": [],
            "sparsity_mean": 0.0,
            "sparsity_std": 0.0
        }
    
    # Extract feature arrays (32K-dim each)
    all_features = []
    all_sparsities = []
    feature_activation_counts = defaultdict(float)
    
    for vf in variant_features:
        features = vf.get("sae_features", [])
        if features:
            all_features.append(features)
            
            # Track sparsity
            stats = vf.get("stats", {})
            sparsity = stats.get("sparsity", 0.0)
            all_sparsities.append(sparsity)
            
            # Track top features
            for top_feat in vf.get("top_features", []):
                idx = top_feat.get("index")
                val = top_feat.get("value", 0.0)
                if idx is not None:
                    feature_activation_counts[idx] += val
    
    # Compute mean and max across variants
    if all_features:
        features_array = np.array(all_features)  # Shape: [num_variants, 32768]
        mean_features = np.mean(features_array, axis=0).tolist()
        max_features = np.max(features_array, axis=0).tolist()
    else:
        mean_features = []
        max_features = []
    
    # Top features aggregated (most activated across all variants)
    top_features_aggregated = [
        {"index": int(idx), "value": float(val)}
        for idx, val in sorted(feature_activation_counts.items(), key=lambda x: x[1], reverse=True)[:64]
    ]
    
    # Sparsity statistics
    sparsity_mean = float(np.mean(all_sparsities)) if all_sparsities else 0.0
    sparsity_std = float(np.std(all_sparsities)) if all_sparsities else 0.0
    
    return {
        "num_variants": len(variant_features),
        "mean_features": mean_features,
        "max_features": max_features,
        "top_features_aggregated": top_features_aggregated,
        "sparsity_mean": sparsity_mean,
        "sparsity_std": sparsity_std
    }

async def extract_sae_features_for_patient(
    client: httpx.AsyncClient,
    patient: Dict[str, Any],
    model_id: str = "evo2_7b"
) -> Optional[Dict[str, Any]]:
    """
    Extract SAE features for all variants of a single patient.
    
    Args:
        client: httpx AsyncClient
        patient: Patient dict with variants and outcome
        model_id: Evo2 model to use
    
    Returns:
        {
            "patient_id": str,
            "outcome": str,  # "sensitive", "resistant", "refractory"
            "variants": [...],  # Variant-level SAE features
            "aggregated_features": {...},  # Patient-level aggregated features
            "provenance": {...}
        }
    """
    patient_id = patient.get("patient_id") or patient.get("tcga_patient_id")
    outcome = patient.get("platinum_response") or patient.get("outcome")
    variants = patient.get("variants", [])
    
    if not patient_id or not outcome:
        logger.warning(f"   ⚠️  Skipping patient with missing ID or outcome: {patient}")
        return None
    
    # Limit variants per patient
    variants_to_process = variants[:MAX_VARIANTS_PER_PATIENT]
    
    logger.info(f"   🧬 Patient {patient_id}: {len(variants_to_process)} variants, outcome={outcome}")
    
    # Extract SAE features for each variant
    variant_features = []
    for i, variant in enumerate(variants_to_process):
        logger.debug(f"      Variant {i+1}/{len(variants_to_process)}: {variant.get('gene')} {variant.get('hgvs_p')}")
        
        sae_result = await extract_sae_features_for_variant(client, variant, model_id)
        if sae_result:
            variant_features.append(sae_result)
    
    if not variant_features:
        logger.warning(f"   ⚠️  No SAE features extracted for patient {patient_id}")
        return None
    
    # Aggregate patient-level features
    aggregated = aggregate_patient_sae_features(variant_features)
    
    logger.info(f"   ✅ Patient {patient_id}: {len(variant_features)}/{len(variants_to_process)} variants successful")
    
    return {
        "patient_id": patient_id,
        "outcome": outcome,
        "variants": variant_features,
        "aggregated_features": aggregated,
        "provenance": {
            "extraction_date": datetime.now().isoformat(),
            "model_id": model_id,
            "num_variants_attempted": len(variants_to_process),
            "num_variants_successful": len(variant_features),
            "sae_version": "v1",
            "phase": "phase2_validation"
        }
    }

async def extract_cohort_sae_features():
    """Main extraction pipeline for TCGA-OV cohort."""
    logger.info("=" * 80)
    logger.info("🧬 SAE PHASE 2: COHORT FEATURE EXTRACTION")
    logger.info("=" * 80)
    
    # Check feature flags
    if not check_feature_flags():
        logger.error("❌ Feature flags not enabled. Exiting.")
        sys.exit(1)
    
    # Load platinum labels
    platinum_data = load_platinum_labels()
    all_patients = platinum_data.get("patients", [])

    # The new patient records have "mutations" key (from extract_patient_mutations_for_cohort.py)
    # Map to "variants" for compatibility with existing SAE extraction logic
    logger.info("🔍 Mapping mutations to variants for SAE extraction...")
    for p in all_patients:
        if p.get("mutations") and not p.get("variants"):
            p["variants"] = p["mutations"]  # Use mutations as variants
    
    has_any_variants = any(bool(p.get("variants")) for p in all_patients)
    if not has_any_variants:
        logger.error("❌ No mutations/variants found in patient records. Cannot proceed.")
        sys.exit(1)
    
    logger.info(f"   ✅ Mapped mutations to variants for {sum(1 for p in all_patients if p.get('variants'))} patients.")
    
    # Filter: only fully labeled patients with variants
    labeled_patients = [
        p for p in all_patients
        if p.get("platinum_response") in ["sensitive", "resistant", "refractory"]
        and p.get("variants")
    ]
    
    logger.info(f"📊 Cohort stats:")
    logger.info(f"   Total patients: {len(all_patients)}")
    logger.info(f"   Fully labeled with variants: {len(labeled_patients)}")
    logger.info(f"   Processing limit: {MAX_PATIENTS}")
    
    # Limit to MAX_PATIENTS for Phase 2
    patients_to_process = labeled_patients[:MAX_PATIENTS]
    
    # Load checkpoint
    checkpoint = load_checkpoint()
    completed_patient_ids = set(checkpoint.get("completed_patients", []))
    failed_patient_ids = set(checkpoint.get("failed_patients", []))
    
    # Filter out already completed patients AND previously failed patients
    remaining_patients = [
        p for p in patients_to_process
        if p.get("patient_id") not in completed_patient_ids
        and p.get("tcga_patient_id") not in completed_patient_ids
        and p.get("patient_id") not in failed_patient_ids
        and p.get("tcga_patient_id") not in failed_patient_ids
    ]
    
    logger.info(f"   Already completed: {len(completed_patient_ids)}")
    logger.info(f"   Previously failed (skipping): {len(failed_patient_ids)}")
    logger.info(f"   Remaining: {len(remaining_patients)}")
    
    if not remaining_patients:
        logger.info("✅ All patients already processed!")
        return
    
    # Extract SAE features for each patient
    cohort_results = []
    failed_patients = checkpoint.get("failed_patients", [])
    
    # Circuit breaker tracking
    total_variants_processed = 0
    total_variants_failed = 0
    
    async with httpx.AsyncClient() as client:
        for i, patient in enumerate(remaining_patients):
            patient_id = patient.get("patient_id") or patient.get("tcga_patient_id")
            logger.info(f"\n[{i+1}/{len(remaining_patients)}] Processing patient: {patient_id}")
            
            # Check total variant cap
            if total_variants_processed >= MAX_TOTAL_VARIANTS:
                logger.warning(f"⚠️  Reached MAX_TOTAL_VARIANTS limit ({MAX_TOTAL_VARIANTS}). Stopping extraction.")
                break
            
            try:
                patient_result = await extract_sae_features_for_patient(client, patient, model_id="evo2_7b")
                
                if patient_result:
                    cohort_results.append(patient_result)
                    checkpoint["completed_patients"].append(patient_id)
                    
                    # Track variant counts for circuit breaker
                    num_variants = patient_result.get("num_variants", 0)
                    total_variants_processed += num_variants
                    
                    # Save checkpoint every 10 patients
                    if (i + 1) % 10 == 0:
                        save_checkpoint(checkpoint)
                        logger.info(f"   💾 Checkpoint saved: {len(checkpoint['completed_patients'])} patients completed")
                        logger.info(f"   📊 Variants processed: {total_variants_processed}/{MAX_TOTAL_VARIANTS}")
                else:
                    failed_patients.append(patient_id)
                    total_variants_failed += patient.get("variants", [])[:MAX_VARIANTS_PER_PATIENT].__len__()
                    logger.warning(f"   ⚠️  Failed to extract SAE features for {patient_id}")
                
                # Circuit breaker: check error rate
                if total_variants_processed + total_variants_failed >= MIN_CALLS_BEFORE_CHECK:
                    error_rate = total_variants_failed / (total_variants_processed + total_variants_failed)
                    if error_rate > MAX_ERROR_RATE:
                        logger.error(f"🚨 Circuit breaker triggered! Error rate: {error_rate:.1%} (>{MAX_ERROR_RATE:.1%})")
                        logger.error(f"   Variants processed: {total_variants_processed}, failed: {total_variants_failed}")
                        logger.error("   Stopping extraction to prevent credit burn. Check configuration and logs.")
                        break
            
            except Exception as e:
                logger.error(f"   ❌ Error processing patient {patient_id}: {e}")
                failed_patients.append(patient_id)
                total_variants_failed += len(patient.get("variants", [])[:MAX_VARIANTS_PER_PATIENT])
                continue
    
    # Final checkpoint
    checkpoint["failed_patients"] = failed_patients
    save_checkpoint(checkpoint)
    
    # Save cohort results
    logger.info(f"\n💾 Saving cohort results to {OUTPUT_FILE}")
    
    output_data = {
        "cohort": "TCGA-OV",
        "outcome": "platinum_response",
        "extraction_date": datetime.now().isoformat(),
        "num_patients": len(cohort_results),
        "num_failed": len(failed_patients),
        "patients": cohort_results,
        "provenance": {
            "phase": "phase2_validation",
            "model_id": "evo2_7b",  # Default to evo2_7b for trained weights
            "sae_version": "v1",
            "max_variants_per_patient": MAX_VARIANTS_PER_PATIENT,
            "feature_flags": {
                "ENABLE_EVO2_SAE": True,
                "ENABLE_TRUE_SAE": True
            }
        }
    }
    
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    logger.info("=" * 80)
    logger.info("✅ COHORT EXTRACTION COMPLETE")
    logger.info(f"   Patients processed: {len(cohort_results)}")
    logger.info(f"   Patients failed: {len(failed_patients)}")
    logger.info(f"   Output: {OUTPUT_FILE}")
    logger.info("=" * 80)

if __name__ == "__main__":
    asyncio.run(extract_cohort_sae_features())

