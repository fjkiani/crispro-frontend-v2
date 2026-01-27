#!/usr/bin/env python3
"""
Download Mutations and Expression Data for Paired Samples
Using cBioPortal MCP API Client
"""

import asyncio
import json
import pandas as pd
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional
import sys

# Add cBioPortal MCP to path
mcp_path = Path(__file__).parent.parent.parent / "tools" / "cbioportal-mcp"
sys.path.insert(0, str(mcp_path))

# Import API client directly (avoiding full server import)
import importlib.util
spec = importlib.util.spec_from_file_location(
    "api_client",
    mcp_path / "cbioportal_mcp" / "api_client.py"
)
api_client_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(api_client_module)
APIClient = api_client_module.APIClient

# Simple config loader
def load_config():
    """Load configuration with defaults."""
    return {
        "server.base_url": "https://www.cbioportal.org/api",
        "server.client_timeout": 480.0
    }

DATA_DIR = Path("data/serial_sae/cbioportal_paired")
OUTPUT_DIR = DATA_DIR
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

async def get_molecular_profiles(api_client: APIClient, study_id: str) -> List[Dict]:
    """Get available molecular profiles for a study."""
    try:
        profiles = await api_client.make_api_request(
            f"studies/{study_id}/molecular-profiles"
        )
        return profiles if isinstance(profiles, list) else []
    except Exception as e:
        print(f"⚠️  Error fetching molecular profiles: {e}")
        return []

async def download_mutations_for_samples(
    api_client: APIClient,
    study_id: str,
    sample_ids: List[str],
    profile_id: str
) -> pd.DataFrame:
    """Download mutation data for specific samples."""
    print(f"📥 Downloading mutations for {len(sample_ids)} samples...")
    
    all_mutations = []
    
    # cBioPortal mutation API requires POST with sample list
    # Process in batches of 50 (API limit)
    batch_size = 50
    for i in range(0, len(sample_ids), batch_size):
        chunk = sample_ids[i:i+batch_size]
        print(f"  Fetching chunk {i//batch_size + 1}/{(len(sample_ids)-1)//batch_size + 1} ({len(chunk)} samples)...")
        
        try:
            mutations = await api_client.make_api_request(
                "mutations/fetch",
                method="POST",
                json_data={
                    "molecularProfileId": profile_id,
                    "sampleIds": chunk,
                    "projection": "DETAILED"
                }
            )
            
            if isinstance(mutations, list):
                all_mutations.extend(mutations)
                print(f"    ✅ Got {len(mutations)} mutations")
            else:
                print(f"    ⚠️  Unexpected response type: {type(mutations)}")
        except Exception as e:
            print(f"    ⚠️  Error fetching chunk: {e}")
            continue
    
    if not all_mutations:
        print("⚠️  No mutations found")
        return pd.DataFrame()
    
    df = pd.DataFrame(all_mutations)
    print(f"✅ Downloaded {len(df)} total mutations")
    return df

async def download_expression_for_samples(
    api_client: APIClient,
    study_id: str,
    sample_ids: List[str],
    profile_id: str
) -> pd.DataFrame:
    """Download expression data for specific samples."""
    print(f"📥 Downloading expression for {len(sample_ids)} samples...")
    
    try:
        data = await api_client.make_api_request(
            "molecular-data/fetch",
            method="POST",
            json_data={
                "molecularProfileIds": [profile_id],
                "sampleIds": sample_ids,
                "projection": "DETAILED"
            }
        )
        
        if not data or not isinstance(data, dict) or "data" not in data:
            print("⚠️  No expression data found")
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
    except Exception as e:
        print(f"⚠️  Error downloading expression: {e}")
        return pd.DataFrame()

async def main():
    print("="*80)
    print("Download Molecular Data for Paired Samples (via cBioPortal MCP)")
    print("="*80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Load configuration
    config = load_config()
    
    # Initialize API client
    api_client = APIClient(
        base_url=config.get("server.base_url", "https://www.cbioportal.org/api"),
        client_timeout=config.get("server.client_timeout", 480.0)
    )
    
    # Start the client
    await api_client.startup()
    
    # Load paired patients
    paired_file = DATA_DIR / "paired_patients.json"
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
        
        # Get molecular profiles
        profiles = await get_molecular_profiles(api_client, study_id)
        print(f"📋 Available molecular profiles: {len(profiles)}")
        
        # Find mutation profile
        mut_profile = None
        for profile in profiles:
            profile_id = profile.get("molecularProfileId", "").lower()
            if "mutation" in profile_id:
                mut_profile = profile.get("molecularProfileId")
                print(f"✅ Mutation profile: {mut_profile}")
                break
        
        if mut_profile:
            mutations_df = await download_mutations_for_samples(
                api_client, study_id, all_sample_ids, mut_profile
            )
            if not mutations_df.empty:
                mut_file = OUTPUT_DIR / f"{study_id}_mutations.csv"
                mutations_df.to_csv(mut_file, index=False)
                print(f"💾 Saved: {mut_file}")
        else:
            print("⚠️  No mutation profile found")
        
        # Find expression profile
        expr_profile = None
        for profile in profiles:
            profile_id = profile.get("molecularProfileId", "").lower()
            if "rna" in profile_id or "expression" in profile_id or "mrna" in profile_id:
                expr_profile = profile.get("molecularProfileId")
                print(f"✅ Expression profile: {expr_profile}")
                break
        
        if expr_profile:
            expression_df = await download_expression_for_samples(
                api_client, study_id, all_sample_ids, expr_profile
            )
            if not expression_df.empty:
                expr_file = OUTPUT_DIR / f"{study_id}_expression.csv"
                expression_df.to_csv(expr_file)
                print(f"💾 Saved: {expr_file}")
        else:
            print("⚠️  No expression profile found")
    
    # Shutdown client
    await api_client.shutdown()
    
    print("\n" + "="*80)
    print("✅ Molecular data download complete")
    print("="*80)

if __name__ == "__main__":
    asyncio.run(main())
