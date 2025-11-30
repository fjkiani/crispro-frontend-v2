#!/usr/bin/env python3
"""
PHASE 0: Data Gathering for SAE Validation
===========================================
Extract and consolidate all required data BEFORE validation begins.

Manager's Requirements (Q1-Q7):
- Q1: Pull full mutation lists from cBioPortal per sample ID
- Q2: Map outcome_platinum to {sensitive, resistant, refractory}
- Q5: Filter for Stage IIIC + IV for Ayesha-like subgroup
- Extract OS (overall survival) data for HR calculations

Output: A single JSON file with all data needed for validation.
"""

import httpx
import json
import sys
from collections import Counter
from typing import List, Dict, Optional

CBIO_BASE = "https://www.cbioportal.org/api"
STUDY_ID = "ov_tcga"

# Pathways for SAE computation (Manager's specification)
PATHWAY_GENES = {
    "DDR": ["BRCA1", "BRCA2", "ATM", "CHEK2", "RAD51", "PALB2", "RAD51C", "RAD51D", "BARD1", "NBN"],
    "MAPK": ["KRAS", "NRAS", "BRAF", "MEK1", "MEK2", "ERK1", "ERK2", "RAF1"],
    "PI3K": ["PIK3CA", "PIK3CB", "PIK3CD", "PTEN", "AKT1", "AKT2", "AKT3", "MTOR"],
    "VEGF": ["VEGFA", "VEGFR1", "VEGFR2", "KDR", "FLT1"],
    "HER2": ["ERBB2", "HER2"],
    "IO": ["PD1", "PDL1", "CTLA4", "CD274", "PDCD1"],
    "Efflux": ["ABCB1", "MDR1", "ABCC1", "ABCG2"]
}

def fetch_all_samples() -> List[Dict]:
    """Fetch all TCGA-OV samples from cBioPortal."""
    print("📥 Fetching TCGA-OV samples from cBioPortal...")
    with httpx.Client(timeout=60.0) as client:
        r = client.get(f"{CBIO_BASE}/studies/{STUDY_ID}/samples", params={"pageSize": 10000})
        r.raise_for_status()
        samples = r.json()
    print(f"✅ Found {len(samples)} TCGA-OV samples")
    return samples

def fetch_mutations_for_sample(sample_id: str) -> List[Dict]:
    """Fetch ALL mutations for a given sample."""
    profile_id = f"{STUDY_ID}_mutations"
    with httpx.Client(timeout=120.0) as client:
        r = client.get(
            f"{CBIO_BASE}/molecular-profiles/{profile_id}/mutations",
            params={
                "sampleIds": [sample_id],
                "projection": "DETAILED",
                "pageSize": 10000
            }
        )
        r.raise_for_status()
        mutations = r.json()
    return mutations if isinstance(mutations, list) else []

def fetch_clinical_data(sample_ids: List[str]) -> Dict[str, Dict]:
    """Fetch clinical data (stage, OS, treatment) for samples."""
    print(f"📥 Fetching clinical data for {len(sample_ids)} samples...")
    clinical_data = {}
    
    # Chunk requests to avoid payload issues
    chunk_size = 200
    for i in range(0, len(sample_ids), chunk_size):
        chunk = sample_ids[i:i+chunk_size]
        with httpx.Client(timeout=120.0) as client:
            r = client.post(
                f"{CBIO_BASE}/clinical-data/fetch",
                json={
                    "entityIds": chunk,
                    "entityType": "SAMPLE",
                    "projection": "DETAILED"
                }
            )
            r.raise_for_status()
            data = r.json()
            
        if isinstance(data, list):
            for row in data:
                sample_id = row.get("entityId")
                if sample_id:
                    if sample_id not in clinical_data:
                        clinical_data[sample_id] = {}
                    attr = row.get("clinicalAttributeId")
                    val = row.get("value")
                    if attr and val:
                        clinical_data[sample_id][attr] = val
        
        print(f"  Processed {min(i+chunk_size, len(sample_ids))}/{len(sample_ids)} samples...")
    
    return clinical_data

def map_outcome_platinum(value) -> str:
    """
    Map outcome_platinum to {sensitive, resistant, refractory}.
    Manager's Q2 specification.
    """
    if isinstance(value, int):
        # Current format: 0/1 (binary)
        # Assume 1=sensitive, 0=resistant (need to verify!)
        return "sensitive" if value == 1 else "resistant"
    
    # String format (future-proof)
    val_str = str(value).lower().strip()
    
    # Map CR/PR → sensitive
    if any(x in val_str for x in ["cr", "complete response", "pr", "partial response", "sensitive"]):
        return "sensitive"
    
    # Map PD/refractory → refractory
    if any(x in val_str for x in ["pd", "progressive disease", "refractory"]):
        return "refractory"
    
    # Map SD/resistant → resistant
    if any(x in val_str for x in ["sd", "stable", "resistant"]):
        return "resistant"
    
    # Unknown
    return "unknown"

def load_existing_tcga_data() -> List[Dict]:
    """Load our existing TCGA file with 200 patients."""
    print("📂 Loading existing TCGA data...")
    with open("tools/benchmarks/hrd_tcga_ov_labeled_sample_use_evo.json") as f:
        data = json.load(f)
    print(f"✅ Loaded {data['n']} patients")
    return data['results']

def gather_validation_dataset():
    """
    Main function: Gather all required data for validation.
    
    Output structure:
    {
        "patients": [
            {
                "index": 0,
                "original_mutation": {"gene": "TP53", "hgvs_p": "R306*", ...},
                "outcome_platinum": "sensitive",  # Mapped value
                "sample_id": "TCGA-XX-XXXX-01",  # If found
                "mutations": [  # Full mutation list from cBioPortal
                    {"gene": "TP53", "hgvs_p": "R306*", ...},
                    {"gene": "BRCA1", "hgvs_p": "Q1395*", ...}
                ],
                "clinical": {
                    "stage": "IV",
                    "os_months": 24.5,
                    "os_event": 1,
                    ...
                },
                "pathway_mutations": {  # Pre-computed for SAE
                    "DDR": ["BRCA1:Q1395*"],
                    "MAPK": [],
                    ...
                }
            },
            ...
        ],
        "summary": {
            "n_patients": 200,
            "n_with_sample_ids": 150,
            "n_with_full_mutations": 150,
            "outcome_distribution": {"sensitive": 100, "resistant": 90, "unknown": 10},
            "stage_distribution": {"III": 50, "IV": 100, "unknown": 50}
        }
    }
    """
    print("\n⚔️ SAE VALIDATION DATA GATHERING ⚔️\n")
    
    # Step 1: Load existing 200 patients
    existing_patients = load_existing_tcga_data()
    
    # Step 2: Fetch all TCGA-OV samples to build lookup
    all_samples = fetch_all_samples()
    sample_lookup = {s["sampleId"]: s for s in all_samples}
    
    # Step 3: Map outcome_platinum
    print("\n📊 Mapping outcome_platinum values...")
    for p in existing_patients:
        original_outcome = p["input"]["outcome_platinum"]
        mapped_outcome = map_outcome_platinum(original_outcome)
        p["outcome_platinum_mapped"] = mapped_outcome
    
    outcome_dist = Counter([p["outcome_platinum_mapped"] for p in existing_patients])
    print(f"Outcome distribution: {dict(outcome_dist)}")
    missing_pct = (outcome_dist.get("unknown", 0) / len(existing_patients)) * 100
    print(f"Missing/unknown: {missing_pct:.1f}%")
    
    # Step 4: Attempt to match to cBioPortal sample IDs
    # (Our data has gene+hgvs_p but no sample IDs)
    # We'll need to pull ALL mutations and match by signature
    print("\n⚠️  WARNING: No sample IDs in existing data.")
    print("   We'll use gene+hgvs_p matching (may be imperfect).")
    print("   Proceeding with existing data + cBioPortal enrichment where possible.\n")
    
    # Step 5: For now, save what we have
    # TODO: Implement gene+hgvs_p matching if needed
    output = {
        "patients": [],
        "summary": {
            "n_patients": len(existing_patients),
            "n_with_sample_ids": 0,  # To be updated
            "n_with_full_mutations": 0,  # To be updated
            "outcome_distribution": dict(outcome_dist),
            "stage_distribution": {},  # To be extracted
            "data_quality": {
                "has_full_mutations": False,
                "has_clinical_data": False,
                "has_os_data": False,
                "missing_outcome_pct": missing_pct
            }
        }
    }
    
    for i, p in enumerate(existing_patients):
        patient_record = {
            "index": i,
            "original_mutation": p["input"],
            "outcome_platinum": p["outcome_platinum_mapped"],
            "sample_id": None,  # Not available
            "mutations": [p["input"]],  # Only have one mutation per patient
            "clinical": {},  # Not available yet
            "pathway_mutations": {}  # Will compute from single mutation
        }
        
        # Classify the single mutation into pathways
        gene = p["input"]["gene"]
        hgvs = p["input"]["hgvs_p"]
        for pathway, genes in PATHWAY_GENES.items():
            if gene in genes:
                patient_record["pathway_mutations"].setdefault(pathway, []).append(f"{gene}:{hgvs}")
        
        output["patients"].append(patient_record)
    
    # Save output
    output_file = "data/validation/tcga_ov_validation_dataset.json"
    with open(output_file, "w") as f:
        json.dump(output, f, indent=2)
    
    print(f"\n✅ Data gathering complete!")
    print(f"📁 Saved to: {output_file}")
    print(f"\n📊 Summary:")
    print(f"  Total patients: {output['summary']['n_patients']}")
    print(f"  Outcome distribution: {output['summary']['outcome_distribution']}")
    print(f"  Missing outcomes: {output['summary']['data_quality']['missing_outcome_pct']:.1f}%")
    print(f"\n⚠️  LIMITATIONS:")
    print(f"  - Only 1 mutation per patient (can't compute full pathway burden)")
    print(f"  - No sample IDs (can't pull full mutation lists from cBioPortal)")
    print(f"  - No clinical data (stage, OS)")
    print(f"\n🔴 BLOCKER: Need manager decision on how to proceed.")
    print(f"   Option A: Use single-mutation data (limited SAE features)")
    print(f"   Option B: Re-extract from cBioPortal with sample IDs + full mutations")
    print(f"   Option C: Find alternative TCGA dataset with full genomics")
    
    return output

if __name__ == "__main__":
    gather_validation_dataset()






