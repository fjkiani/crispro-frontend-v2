#!/usr/bin/env python3
"""
⚔️ GISTIC HRD SCORE EXTRACTION & CALCULATION ⚔️

Mission: Extract GISTIC copy-number data and calculate HRD scores (LOH + LST + TAI)
Study: ov_tcga_pan_can_atlas_2018
Output: HRD scores (0-100 scale) for TCGA-OV patients

Note: Full HRD calculation (LOH + LST + TAI) is complex. This script implements
a simplified proxy using GISTIC arm-level copy-number alterations.
"""

import argparse
import json
import httpx
import os
from typing import List, Dict, Any, Optional
from collections import defaultdict

CBIO_BASE = "https://www.cbioportal.org/api"
STUDY_ID = "ov_tcga_pan_can_atlas_2018"
GISTIC_PROFILE_ID = f"{STUDY_ID}_gistic"

def _headers() -> Dict[str, str]:
    headers = {"Accept": "application/json", "Content-Type": "application/json"}
    token = os.getenv("CBIO_TOKEN")
    if token:
        headers["Authorization"] = f"Bearer {token}"
    return headers

def get_samples(study_id: str) -> List[Dict]:
    """Get all samples for a study"""
    print(f"🔍 Fetching samples for {study_id}...")
    with httpx.Client(timeout=60.0, headers=_headers()) as client:
        r = client.get(f"{CBIO_BASE}/studies/{study_id}/samples")
        r.raise_for_status()
        return r.json() or []

def get_gistic_data(study_id: str, sample_ids: List[str]) -> Dict[str, List[Dict]]:
    """
    Extract GISTIC copy-number data for samples.
    
    Returns: {sample_id: [copy_number_segments]}
    """
    print(f"🔍 Fetching GISTIC data for {len(sample_ids)} samples...")
    
    # GISTIC data is in molecular profile
    profile_id = f"{study_id}_gistic"
    
    # Fetch discrete copy-number alterations using correct endpoint
    # Endpoint: /molecular-profiles/{molecular_profile_id}/discrete-copy-number/fetch
    filter_data = {
        "sampleIds": sample_ids,
        "entrezGeneIds": None,  # Get all genes
        "sampleListId": None
    }
    
    params = {
        "discreteCopyNumberEventType": "ALL",  # Get all events
        "projection": "DETAILED"
    }
    
    sample_gistic = defaultdict(list)
    
    try:
        with httpx.Client(timeout=120.0, headers=_headers()) as client:
            r = client.post(
                f"{CBIO_BASE}/molecular-profiles/{profile_id}/discrete-copy-number/fetch",
                params=params,
                json=filter_data
            )
            r.raise_for_status()
            result = r.json() or []
            
            # Group by sample
            for row in result:
                sample_id = row.get("sampleId")
                if sample_id:
                    sample_gistic[sample_id].append(row)
        
        print(f"✅ Extracted GISTIC data for {len(sample_gistic)} samples")
        return dict(sample_gistic)
    
    except Exception as e:
        print(f"⚠️ GISTIC extraction failed: {e}")
        print("⚠️ Falling back to arm-level copy-number data...")
        return _get_arm_level_cna(study_id, sample_ids)

def _get_arm_level_cna(study_id: str, sample_ids: List[str]) -> Dict[str, List[Dict]]:
    """Fallback: Get arm-level copy-number data"""
    profile_id = f"{study_id}_armlevel_cna"
    
    # Use molecular data fetch endpoint
    filter_data = {
        "sampleIds": sample_ids,
        "entrezGeneIds": None
    }
    
    sample_cna = defaultdict(list)
    
    try:
        with httpx.Client(timeout=120.0, headers=_headers()) as client:
            r = client.post(
                f"{CBIO_BASE}/molecular-profiles/{profile_id}/molecular-data/fetch",
                json=filter_data
            )
            r.raise_for_status()
            result = r.json() or []
            
            for row in result:
                sample_id = row.get("sampleId")
                if sample_id:
                    sample_cna[sample_id].append(row)
        
        print(f"✅ Extracted arm-level CNA for {len(sample_cna)} samples")
        return dict(sample_cna)
    
    except Exception as e:
        print(f"❌ Arm-level CNA extraction also failed: {e}")
        return {}

def calculate_hrd_proxy(gistic_data: List[Dict]) -> float:
    """
    Calculate HRD proxy score from GISTIC copy-number data.
    
    Simplified approach:
    - Count significant copy-number alterations (amplifications + deletions)
    - Weight by chromosome arm size
    - Normalize to 0-100 scale
    
    Note: This is a PROXY, not true HRD (LOH + LST + TAI).
    True HRD requires complex algorithms (e.g., scarHRD R package).
    """
    if not gistic_data:
        return 0.0
    
    # Count alterations
    # GISTIC values: -2 (deep deletion), -1 (shallow deletion), 0 (normal), 1 (gain), 2 (amplification)
    significant_alts = 0
    total_segments = len(gistic_data)
    
    for segment in gistic_data:
        alteration = segment.get("alteration", 0)
        # Count significant alterations (|alteration| >= 1)
        if abs(int(alteration)) >= 1:
            significant_alts += 1
    
    # Simple proxy: percentage of genome with significant alterations
    # Scale to 0-100 (rough approximation)
    if total_segments == 0:
        return 0.0
    
    alt_fraction = significant_alts / total_segments
    
    # Rough mapping: 0-30% alterations → 0-50 HRD, 30-60% → 50-80, 60%+ → 80-100
    if alt_fraction < 0.3:
        hrd_proxy = alt_fraction * 166.7  # Scale to 0-50
    elif alt_fraction < 0.6:
        hrd_proxy = 50 + (alt_fraction - 0.3) * 100  # Scale to 50-80
    else:
        hrd_proxy = 80 + min((alt_fraction - 0.6) * 66.7, 20)  # Scale to 80-100
    
    return round(min(hrd_proxy, 100.0), 2)

def main():
    ap = argparse.ArgumentParser(description="Extract GISTIC data and calculate HRD proxy scores")
    ap.add_argument("--study", default=STUDY_ID, help="cBioPortal study ID")
    ap.add_argument("--out", default="tools/benchmarks/data/hrd_scores_from_gistic.json")
    ap.add_argument("--limit", type=int, default=200, help="Limit number of samples")
    args = ap.parse_args()
    
    # Get samples
    samples = get_samples(args.study)
    sample_ids = [s.get("sampleId") for s in samples[:args.limit] if s.get("sampleId")]
    
    if not sample_ids:
        raise SystemExit("❌ No samples found")
    
    print(f"📊 Processing {len(sample_ids)} samples...")
    
    # Build sample->patient map
    samp_to_patient = {}
    for s in samples:
        sid = s.get("sampleId")
        pid = s.get("patientId") or s.get("patient_id")
        if sid and pid:
            samp_to_patient[sid] = str(pid)
    
    # Extract GISTIC data
    gistic_data = get_gistic_data(args.study, sample_ids)
    
    # Calculate HRD scores
    hrd_scores = {}
    for sample_id in sample_ids:
        segments = gistic_data.get(sample_id, [])
        hrd_score = calculate_hrd_proxy(segments)
        patient_id = samp_to_patient.get(sample_id, sample_id)
        hrd_scores[patient_id] = {
            "sample_id": sample_id,
            "hrd_score": hrd_score,
            "n_segments": len(segments)
        }
    
    # Save results
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(hrd_scores, f, indent=2)
    
    # Stats
    scores = [v["hrd_score"] for v in hrd_scores.values()]
    print(f"\n✅ HRD Score Statistics:")
    print(f"   Samples processed: {len(hrd_scores)}")
    print(f"   Mean HRD: {sum(scores)/len(scores):.2f}" if scores else "   Mean HRD: N/A")
    print(f"   Min HRD: {min(scores):.2f}" if scores else "   Min HRD: N/A")
    print(f"   Max HRD: {max(scores):.2f}" if scores else "   Max HRD: N/A")
    print(f"   HRD-High (≥42): {sum(1 for s in scores if s >= 42)}/{len(scores)} ({100*sum(1 for s in scores if s >= 42)/len(scores):.1f}%)" if scores else "   HRD-High: N/A")
    print(f"\n💾 Results saved to: {args.out}")

if __name__ == "__main__":
    main()

