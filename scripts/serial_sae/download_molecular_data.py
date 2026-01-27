#!/usr/bin/env python3
"""
Download Mutations and Expression Data for Paired Samples
Extract molecular data for primary and recurrent samples
"""

import httpx
import json
import pandas as pd
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional
import time

CBIO_BASE = "https://www.cbioportal.org/api"
OUTPUT_DIR = Path("data/serial_sae/cbioportal_paired")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def get_molecular_profile_id(study_id: str, profile_type: str) -> Optional[str]:
    """Get molecular profile ID for a study."""
    with httpx.Client(timeout=60.0) as client:
        r = client.get(f"{CBIO_BASE}/studies/{study_id}/molecular-profiles")
        r.raise_for_status()
        profiles = r.json()
        
        for profile in profiles:
            profile_id = profile.get("molecularProfileId", "").lower()
            if profile_type.lower() in profile_id:
                return profile.get("molecularProfileId")
    
    return None

def download_mutations(study_id: str, sample_ids: List[str], profile_id: str) -> pd.DataFrame:
    """Download mutation data for samples."""
    print(f"📥 Downloading mutations for {len(sample_ids)} samples...")
    
    # cBioPortal mutation API requires POST with sample list
    # Limit to 50 samples per request (API limit)
    all_mutations = []
    for i in range(0, len(sample_ids), 50):
        chunk = sample_ids[i:i+50]
        print(f"  Fetching chunk {i//50 + 1}/{(len(sample_ids)-1)//50 + 1} ({len(chunk)} samples)...")
        
        try:
            with httpx.Client(timeout=300.0) as client:
                r = client.post(
                    f"{CBIO_BASE}/mutations/fetch",
                    json={
                        "molecularProfileId": profile_id,
                        "sampleIds": chunk,
                        "projection": "DETAILED"
                    }
                )
                r.raise_for_status()
                
                if r.text:
                    mutations = r.json()
                    if mutations:
                        all_mutations.extend(mutations)
                else:
                    print(f"    ⚠️  Empty response for chunk")
        except Exception as e:
            print(f"    ⚠️  Error fetching chunk: {e}")
            continue
    
    if not all_mutations:
        print("⚠️  No mutations found")
        return pd.DataFrame()
    
    df = pd.DataFrame(all_mutations)
    print(f"✅ Downloaded {len(df)} mutations")
    return df

def download_expression(study_id: str, sample_ids: List[str], profile_id: str) -> pd.DataFrame:
    """Download expression data for samples."""
    print(f"📥 Downloading expression for {len(sample_ids)} samples...")
    
    # cBioPortal expression API
    with httpx.Client(timeout=300.0) as client:
        r = client.post(
            f"{CBIO_BASE}/molecular-data/fetch",
            json={
                "molecularProfileIds": [profile_id],
                "sampleIds": sample_ids,
                "projection": "DETAILED"
            }
        )
        r.raise_for_status()
        data = r.json()
    
    if not data or "data" not in data:
        return pd.DataFrame()
    
    # Convert to DataFrame
    rows = []
    for row in data["data"]:
        rows.append({
            "sample_id": row.get("sampleId"),
            "gene": row.get("entrezGeneId"),
            "value": row.get("value")
        })
    
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.pivot_table(index="gene", columns="sample_id", values="value")
    
    print(f"✅ Downloaded expression for {len(df.columns)} samples, {len(df)} genes")
    return df

def main():
    print("="*80)
    print("Download Molecular Data for Paired Samples")
    print("="*80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Load paired patients
    paired_file = OUTPUT_DIR / "paired_patients.json"
    if not paired_file.exists():
        print("❌ Error: paired_patients.json not found. Run download_cbioportal_paired.py first.")
        return
    
    with open(paired_file) as f:
        data = json.load(f)
    
    # Extract paired_patients dict
    paired_data = data.get("paired_patients", {})
    print(f"📊 Loaded {len(paired_data)} paired patients")
    
    # Collect all sample IDs by study
    study_samples = {}
    for patient_id, patient_data in paired_data.items():
        study_id = patient_data.get("study_id", "")
        if study_id not in study_samples:
            study_samples[study_id] = {"primary": [], "recurrent": []}
        
        for sample in patient_data.get("primary_samples", []):
            study_samples[study_id]["primary"].append(sample.get("sample_id", ""))
        for sample in patient_data.get("recurrent_samples", []):
            study_samples[study_id]["recurrent"].append(sample.get("sample_id", ""))
    
    # Download data for each study
    for study_id, samples in study_samples.items():
        print(f"\n🔍 Processing study: {study_id}")
        print("-" * 80)
        
        all_sample_ids = samples["primary"] + samples["recurrent"]
        print(f"📊 Total samples: {len(all_sample_ids)} ({len(samples['primary'])} primary + {len(samples['recurrent'])} recurrent)")
        
        # Get mutation profile
        mut_profile = get_molecular_profile_id(study_id, "mutations")
        if mut_profile:
            print(f"✅ Mutation profile: {mut_profile}")
            mutations_df = download_mutations(study_id, all_sample_ids, mut_profile)
            if not mutations_df.empty:
                mut_file = OUTPUT_DIR / f"{study_id}_mutations.csv"
                mutations_df.to_csv(mut_file, index=False)
                print(f"💾 Saved: {mut_file}")
        else:
            print("⚠️  No mutation profile found")
        
        # Get expression profile
        expr_profile = get_molecular_profile_id(study_id, "rna")
        if not expr_profile:
            expr_profile = get_molecular_profile_id(study_id, "expression")
        
        if expr_profile:
            print(f"✅ Expression profile: {expr_profile}")
            expression_df = download_expression(study_id, all_sample_ids, expr_profile)
            if not expression_df.empty:
                expr_file = OUTPUT_DIR / f"{study_id}_expression.csv"
                expression_df.to_csv(expr_file)
                print(f"💾 Saved: {expr_file}")
        else:
            print("⚠️  No expression profile found")
        
        time.sleep(1)  # Rate limiting
    
    print("\n" + "="*80)
    print("✅ Molecular data download complete")
    print("="*80)

if __name__ == "__main__":
    main()
