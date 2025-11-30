#!/usr/bin/env python3
"""
⚔️ MISSION: EXTRACT SOMATIC MUTATIONS FOR TCGA-OV PLATINUM COHORT
==================================================================

Agent: Zo (Lead Commander)
Date: November 21, 2025
Priority: CRITICAL - BLOCKING REAL SAE COHORT
Timeline: 1-2 hours

Objective: Fetch somatic mutations for all 469 labeled TCGA-OV patients
using pyBioPortal and merge them into the platinum labels file.

Inputs:
- data/validation/tcga_ov_platinum_response_labels.json (469 patients with labels)
- pyBioPortal API (via oncology-backend/tests/pyBioPortal-master)

Outputs:
- data/validation/sae_cohort/tcga_ov_platinum_with_mutations.json
  (469 patients with labels + somatic mutations ready for SAE extraction)

Success Criteria:
- Extract mutations for ≥90% of labeled patients (≥422/469)
- Each patient has chrom, pos, ref, alt, gene, hgvs_p for each mutation
- Output schema compatible with extract_sae_features_cohort.py
- Full provenance tracking
"""

import sys
from pathlib import Path
import pandas as pd
import json
import time
from typing import Dict, List, Any, Optional
from datetime import datetime
from loguru import logger

# Add pyBioPortal to path
pybioportal_parent = Path("/Users/fahadkiani/Desktop/development/crispr-assistant-main/oncology-coPilot/oncology-backend/tests/pyBioPortal-master")
if str(pybioportal_parent) not in sys.path:
    sys.path.insert(0, str(pybioportal_parent))

from pybioportal import molecular_profiles as mp
from pybioportal import mutations as mut

# --- Configuration ---
PLATINUM_LABELS_FILE = Path("data/validation/tcga_ov_platinum_response_labels.json")
OUTPUT_DIR = Path("data/validation/sae_cohort")
OUTPUT_FILE = OUTPUT_DIR / "tcga_ov_platinum_with_mutations.json"
CHECKPOINT_FILE = OUTPUT_DIR / "mutations_extraction_checkpoint.json"

# TCGA study IDs for ovarian cancer
STUDY_ID = "ov_tcga_pan_can_atlas_2018"
MUTATION_PROFILE_SUFFIX = "_mutations"

# Extraction parameters
MAX_RETRIES = 3
REQUEST_DELAY_SECONDS = 0.5  # Rate limit
BATCH_SIZE = 100  # Process samples in batches to avoid timeouts

# --- Logging Setup ---
logger.remove()
logger.add(sys.stderr, level="INFO")
logger.add(OUTPUT_DIR / "mutation_extraction_log.log", rotation="10 MB", level="DEBUG")

# --- Helper Functions ---

def load_platinum_labels() -> Dict[str, Any]:
    """Load TCGA-OV platinum response labels."""
    if not PLATINUM_LABELS_FILE.exists():
        logger.error(f"❌ Platinum labels file not found: {PLATINUM_LABELS_FILE}")
        sys.exit(1)
    
    logger.info(f"📥 Loading platinum labels from {PLATINUM_LABELS_FILE}")
    with open(PLATINUM_LABELS_FILE, 'r') as f:
        data = json.load(f)
    
    patients = data.get('patients', [])
    logger.info(f"   ✅ Loaded {len(patients)} patients")
    return data


def get_mutation_profile_id(study_id: str) -> Optional[str]:
    """Get the mutation molecular profile ID for a study."""
    try:
        logger.info(f"🔍 Looking up mutation profile for {study_id}...")
        df_prof = mp.get_all_molecular_profiles_in_study(study_id)
        
        if isinstance(df_prof, pd.DataFrame) and "molecularProfileId" in df_prof.columns:
            candidates = df_prof["molecularProfileId"].astype(str).tolist()
            for pid in candidates:
                if pid.endswith(MUTATION_PROFILE_SUFFIX):
                    logger.info(f"   ✅ Found profile: {pid}")
                    return pid
        
        logger.error(f"❌ No mutation profile found for {study_id}")
        return None
        
    except Exception as e:
        logger.error(f"❌ Error getting mutation profile: {e}")
        return None


def fetch_mutations_for_samples(
    profile_id: str,
    sample_ids: List[str],
    max_retries: int = MAX_RETRIES
) -> pd.DataFrame:
    """
    Fetch mutations for specific sample IDs.
    Uses pyBioPortal mutations.fetch_muts_in_mol_prof with sample_ids filter.
    """
    logger.info(f"📊 Fetching mutations for {len(sample_ids)} samples from {profile_id}...")
    
    for attempt in range(max_retries):
        try:
            # Use fetch_muts_in_mol_prof with sample_ids filter
            df_muts = mut.fetch_muts_in_mol_prof(
                molecular_profile_id=profile_id,
                sample_ids=sample_ids,
                projection="DETAILED",
                pageSize=10000000,  # Get all mutations
                pageNumber=0,
                direction="ASC",
                sortBy=None
            )
            
            if not isinstance(df_muts, pd.DataFrame) or df_muts.empty:
                logger.warning(f"⚠️  No mutations found (attempt {attempt + 1}/{max_retries})")
                if attempt < max_retries - 1:
                    time.sleep(2 ** attempt)
                    continue
                return pd.DataFrame()
            
            logger.info(f"✅ Fetched {len(df_muts)} mutations across {df_muts['sampleId'].nunique()} samples")
            return df_muts
            
        except Exception as e:
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt
                logger.warning(f"⚠️  Error (attempt {attempt + 1}/{max_retries}): {e}")
                logger.info(f"   Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                logger.error(f"❌ Failed to fetch mutations after {max_retries} attempts: {e}")
                return pd.DataFrame()
    
    return pd.DataFrame()


def normalize_mutation_to_schema(row: pd.Series) -> Dict[str, Any]:
    """
    Normalize a pyBioPortal mutation row to our standard schema.
    
    Expected output schema:
    {
        "gene": str,
        "chrom": str,
        "pos": int,
        "ref": str,
        "alt": str,
        "hgvs_p": str (optional),
        "variant_type": str,
        "assembly": "GRCh37" (TCGA Pan-Can uses GRCh37)
    }
    """
    # pyBioPortal actual column names (verified via DETAILED projection):
    # - gene_hugoGeneSymbol
    # - chr
    # - startPosition
    # - referenceAllele
    # - variantAllele
    # - proteinChange
    # - variantType / mutationType
    # - ncbiBuild (e.g., "GRCh37")
    
    gene = row.get('gene_hugoGeneSymbol', 'UNKNOWN')
    chrom = str(row.get('chr', '')).replace('chr', '')
    pos = int(row.get('startPosition', 0))
    ref = row.get('referenceAllele', '')
    alt = row.get('variantAllele', '')
    hgvs_p = row.get('proteinChange', '')
    variant_type = row.get('variantType', row.get('mutationType', 'SNP'))
    assembly = row.get('ncbiBuild', 'GRCh37')  # Use actual build from data
    
    return {
        "gene": gene,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "hgvs_p": hgvs_p,
        "variant_type": variant_type,
        "assembly": assembly
    }


def extract_mutations_for_patient(
    patient_id: str,
    sample_id: str,
    all_mutations_df: pd.DataFrame
) -> List[Dict[str, Any]]:
    """Extract and normalize mutations for a single patient."""
    # Filter mutations for this sample
    patient_muts = all_mutations_df[all_mutations_df['sampleId'] == sample_id]
    
    if patient_muts.empty:
        return []
    
    # Normalize each mutation to our schema
    normalized_muts = []
    for _, row in patient_muts.iterrows():
        try:
            mut_data = normalize_mutation_to_schema(row)
            # Basic validation
            if mut_data['gene'] != 'UNKNOWN' and mut_data['chrom'] and mut_data['pos'] > 0:
                normalized_muts.append(mut_data)
        except Exception as e:
            logger.debug(f"  Skipping mutation for {patient_id} due to normalization error: {e}")
    
    return normalized_muts


def load_checkpoint() -> List[Dict[str, Any]]:
    """Load extraction checkpoint if it exists."""
    if CHECKPOINT_FILE.exists():
        logger.info(f"🔄 Loading checkpoint from {CHECKPOINT_FILE}")
        with open(CHECKPOINT_FILE, 'r') as f:
            return json.load(f)
    return []


def save_checkpoint(patients_with_mutations: List[Dict[str, Any]]):
    """Save checkpoint after processing patients."""
    with open(CHECKPOINT_FILE, 'w') as f:
        json.dump(patients_with_mutations, f, indent=2)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    logger.info("=" * 80)
    logger.info("⚔️ TCGA-OV PLATINUM COHORT MUTATION EXTRACTION")
    logger.info("=" * 80)
    logger.info(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    # 1. Load platinum labels
    platinum_data = load_platinum_labels()
    patients = platinum_data.get('patients', [])
    
    if not patients:
        logger.error("❌ No patients found in platinum labels file")
        sys.exit(1)
    
    # 2. Get mutation profile ID
    profile_id = get_mutation_profile_id(STUDY_ID)
    if not profile_id:
        logger.error(f"❌ Cannot proceed without mutation profile ID for {STUDY_ID}")
        sys.exit(1)
    
    # 3. Extract sample IDs
    sample_ids = [p['tcga_sample_id'] for p in patients]
    logger.info(f"📋 Target: {len(sample_ids)} samples")
    
    # 4. Fetch all mutations for these samples
    logger.info(f"🔬 Fetching mutations from cBioPortal...")
    all_mutations_df = fetch_mutations_for_samples(profile_id, sample_ids)
    
    if all_mutations_df.empty:
        logger.error("❌ No mutations fetched. Cannot proceed.")
        sys.exit(1)
    
    logger.info(f"✅ Total mutations fetched: {len(all_mutations_df)}")
    logger.info(f"   Covering {all_mutations_df['sampleId'].nunique()} unique samples")
    
    # 5. Load checkpoint if exists
    patients_with_mutations = load_checkpoint()
    processed_sample_ids = set(p['tcga_sample_id'] for p in patients_with_mutations)
    logger.info(f"Loaded checkpoint with {len(processed_sample_ids)} already processed patients.")
    
    # 6. Process each patient and extract their mutations
    logger.info(f"\n{'='*80}")
    logger.info("📊 Processing mutations per patient...")
    logger.info(f"{'='*80}\n")
    
    for i, patient in enumerate(patients):
        patient_id = patient['tcga_patient_id']
        sample_id = patient['tcga_sample_id']
        
        if sample_id in processed_sample_ids:
            logger.debug(f"Skipping {patient_id} (already processed)")
            continue
        
        if (i + 1) % 50 == 0:
            logger.info(f"   Processing patient {i+1}/{len(patients)}...")
        
        # Extract mutations for this patient
        mutations = extract_mutations_for_patient(patient_id, sample_id, all_mutations_df)
        
        # Create patient record with mutations
        patient_with_mutations = {
            "tcga_patient_id": patient_id,
            "tcga_sample_id": sample_id,
            "platinum_response": patient['platinum_response'],
            "raw_response_value": patient.get('raw_response_value'),
            "source_field": patient.get('source_field'),
            "mutations": mutations,
            "num_mutations": len(mutations),
            "provenance": {
                "mutation_extraction_script": "scripts/sae/extract_patient_mutations_for_cohort.py",
                "timestamp": datetime.now().isoformat(),
                "study_id": STUDY_ID,
                "mutation_profile_id": profile_id,
                "assembly": "GRCh38"
            }
        }
        
        patients_with_mutations.append(patient_with_mutations)
        
        # Save checkpoint every 50 patients
        if (i + 1) % 50 == 0:
            save_checkpoint(patients_with_mutations)
            logger.info(f"   💾 Checkpoint saved ({len(patients_with_mutations)} patients)")
    
    # 7. Save final results
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(patients_with_mutations, f, indent=2)
    
    # 8. Compute and log statistics
    patients_with_muts = [p for p in patients_with_mutations if p['num_mutations'] > 0]
    total_mutations = sum(p['num_mutations'] for p in patients_with_mutations)
    
    logger.info(f"\n{'='*80}")
    logger.info("✅ MUTATION EXTRACTION COMPLETE")
    logger.info(f"{'='*80}")
    logger.info(f"Output: {OUTPUT_FILE}")
    logger.info(f"Patients processed: {len(patients_with_mutations)}/{len(patients)}")
    logger.info(f"Patients with mutations: {len(patients_with_muts)} ({len(patients_with_muts)/len(patients)*100:.1f}%)")
    logger.info(f"Total mutations extracted: {total_mutations}")
    logger.info(f"Average mutations per patient: {total_mutations/len(patients_with_muts):.1f}" if patients_with_muts else "N/A")
    logger.info(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    logger.info("📝 Next step: Run scripts/sae/extract_sae_features_cohort.py to extract SAE features")


if __name__ == "__main__":
    main()

