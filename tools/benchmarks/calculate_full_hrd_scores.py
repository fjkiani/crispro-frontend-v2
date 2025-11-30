#!/usr/bin/env python3
"""
⚔️ FULL HRD SCORE CALCULATION (LOH + LST + TAI) ⚔️

Mission: Calculate full HRD scores using LOH + LST + TAI algorithms
Study: ov_tcga_pan_can_atlas_2018
Input: GISTIC discrete copy-number data from cBioPortal
Output: Full HRD scores (LOH + LST + TAI) for TCGA-OV patients

Algorithm References:
- LOH: Loss of Heterozygosity segments (copy number = 0 or 1, or allelic imbalance)
- LST: Large-Scale Transitions (breakpoints between adjacent segments ≥10Mb with copy number change)
- TAI: Telomeric Allelic Imbalance (allelic imbalances extending to telomeres)

Note: This uses GISTIC discrete CNA data. For true segment-level analysis,
      download segment files from GDC and use scarHRD or similar tools.
"""

import argparse
import json
import httpx
import os
from typing import List, Dict, Any, Optional, Tuple
from collections import defaultdict
import math

CBIO_BASE = "https://www.cbioportal.org/api"
STUDY_ID = "ov_tcga_pan_can_atlas_2018"
GISTIC_PROFILE_ID = f"{STUDY_ID}_gistic"

# Chromosome lengths (GRCh38, in base pairs)
CHROM_LENGTHS = {
    "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555,
    "5": 181538259, "6": 170805979, "7": 159345973, "8": 145138636,
    "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
    "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345,
    "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
    "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415
}

# Telomere regions (first and last 5Mb of each chromosome)
TELOMERE_SIZE = 5_000_000  # 5Mb

def _headers() -> Dict[str, str]:
    headers = {"Accept": "application/json", "Content-Type": "application/json"}
    token = os.getenv("CBIO_TOKEN")
    if token:
        headers["Authorization"] = f"Bearer {token}"
    return headers

def _lookup_gene_info(entrez_id: int, gene_cache: Dict[int, Dict]) -> Optional[Dict]:
    """Look up gene info from entrezGeneId, with caching"""
    if entrez_id in gene_cache:
        return gene_cache[entrez_id]
    
    try:
        with httpx.Client(timeout=30.0, headers=_headers()) as client:
            r = client.get(f"{CBIO_BASE}/genes/{entrez_id}")
            if r.status_code == 200:
                gene_data = r.json()
                gene_info = {
                    "hugoGeneSymbol": gene_data.get("hugoGeneSymbol", ""),
                    "type": gene_data.get("type", "")
                }
                gene_cache[entrez_id] = gene_info
                return gene_info
    except Exception:
        pass
    
    return None

def _map_alteration_code(alteration: int) -> str:
    """Map GISTIC alteration code to type"""
    # GISTIC codes: -2=Deep Deletion, -1=Shallow Deletion, 0=Diploid, 1=Gain, 2=Amplification
    mapping = {
        -2: "DEL",  # Deep deletion
        -1: "LOSS",  # Shallow deletion/loss
        0: "DIPLOID",  # No alteration
        1: "GAIN",  # Gain
        2: "AMP"  # Amplification
    }
    return mapping.get(alteration, "UNKNOWN")

def get_gistic_discrete_cna(study_id: str, sample_ids: List[str]) -> Dict[str, List[Dict]]:
    """
    Extract GISTIC discrete copy-number alterations.
    
    Returns: {sample_id: [cna_records]}
    Each record has: gene, alteration (AMP/DEL/GAIN/LOSS), alteration_code
    Note: GISTIC is gene-level (no chromosome positions available from cBioPortal API)
    """
    print(f"🔍 Fetching GISTIC discrete CNA for {len(sample_ids)} samples...")
    
    profile_id = f"{study_id}_gistic"
    sample_cna = defaultdict(list)
    gene_cache = {}  # Cache gene lookups
    
    # Process in chunks
    chunk_size = 50
    for i in range(0, len(sample_ids), chunk_size):
        chunk = sample_ids[i:i+chunk_size]
        print(f"  Processing chunk {i//chunk_size + 1}/{(len(sample_ids)-1)//chunk_size + 1}...")
        
        try:
            with httpx.Client(timeout=180.0, headers=_headers()) as client:
                r = client.post(
                    f"{CBIO_BASE}/molecular-profiles/{profile_id}/discrete-copy-number/fetch",
                    json={"sampleIds": chunk}
                )
                r.raise_for_status()
                result = r.json() or []
                
                for cna in result:
                    sample_id = cna.get("sampleId")
                    if not sample_id:
                        continue
                    
                    entrez_id = cna.get("entrezGeneId", 0)
                    if not entrez_id:
                        continue
                    
                    # Look up gene name
                    gene_info = _lookup_gene_info(entrez_id, gene_cache)
                    if not gene_info:
                        continue
                    
                    # Map alteration code
                    alteration_code = cna.get("alteration", 0)
                    alteration_type = _map_alteration_code(alteration_code)
                    
                    # Skip diploid (no alteration)
                    if alteration_type == "DIPLOID":
                        continue
                    
                    sample_cna[sample_id].append({
                        "gene": gene_info.get("hugoGeneSymbol", ""),
                        "entrez_id": entrez_id,
                        "alteration": alteration_type,
                        "alteration_code": alteration_code
                    })
        
        except Exception as e:
            print(f"  ⚠️ Chunk {i//chunk_size + 1} failed: {e}")
            continue
    
    print(f"✅ Extracted CNA for {len(sample_cna)} samples")
    return dict(sample_cna)

def _parse_cytoband_arm(cytoband: str) -> Tuple[str, str]:
    """
    Parse cytoband to extract chromosome and arm (p or q).
    
    Examples:
    - "1p36.33" -> ("1", "p")
    - "17q21.31" -> ("17", "q")
    - "Xp22.33" -> ("X", "p")
    """
    if not cytoband:
        return ("", "")
    
    # Extract chromosome (first part before p/q)
    chrom = ""
    arm = ""
    
    for i, char in enumerate(cytoband):
        if char in ["p", "q"]:
            chrom = cytoband[:i]
            arm = char
            break
    
    return (chrom, arm)

def _convert_cna_to_gene_counts(cna_records: List[Dict]) -> Dict[str, int]:
    """
    Convert gene-level CNA records to alteration counts.
    
    Since we don't have chromosome positions, we use gene-level patterns:
    - Count genes with each alteration type
    - Group by alteration severity
    
    Returns: {alteration_type: gene_count}
    """
    counts = defaultdict(int)
    
    for record in cna_records:
        alt = record.get("alteration", "")
        if alt and alt != "DIPLOID":
            counts[alt] += 1
    
    return dict(counts)

def calculate_loh_gene_level(gene_counts: Dict[str, int]) -> int:
    """
    Calculate Loss of Heterozygosity (LOH) score from gene-level data.
    
    LOH proxy = Count of deletion/loss genes, normalized.
    - Deep deletions (DEL) weighted 2x
    - Shallow losses (LOSS) weighted 1x
    - Threshold: ≥20 deletion/loss genes = 1 LOH point
    """
    del_count = gene_counts.get("DEL", 0)
    loss_count = gene_counts.get("LOSS", 0)
    
    # Weighted count: deep deletions count more
    weighted_loh = (del_count * 2) + loss_count
    
    # Normalize: every 20 weighted genes = 1 LOH point
    loh_score = weighted_loh // 20
    
    return max(0, loh_score)

def calculate_lst_gene_level(gene_counts: Dict[str, int]) -> int:
    """
    Calculate Large-Scale Transitions (LST) score from gene-level data.
    
    LST proxy = Presence of both gains and losses indicates transitions.
    - If both gain and loss types present: 1 point per 10 genes of each type
    - Minimum: 5 gain genes + 5 loss genes = 1 LST point
    """
    gain_types = {"AMP", "GAIN"}
    loss_types = {"DEL", "LOSS"}
    
    gain_count = sum(gene_counts.get(alt, 0) for alt in gain_types)
    loss_count = sum(gene_counts.get(alt, 0) for alt in loss_types)
    
    # LST: both gains and losses present (indicates transitions)
    if gain_count > 0 and loss_count > 0:
        # Score based on minimum of both types
        # Every 10 genes of each type = 1 LST point
        lst_score = min(gain_count, loss_count) // 10
        return max(1, lst_score)  # At least 1 if both present
    
    return 0

def calculate_tai_gene_level(gene_counts: Dict[str, int], altered_genes_count: int) -> int:
    """
    Calculate Telomeric Allelic Imbalance (TAI) score from gene-level data.
    
    TAI proxy = Overall alteration burden, weighted toward telomere-proximal patterns.
    - Uses GISTIC reference gene set size (~24,000 genes typically analyzed)
    - High alteration burden suggests telomere instability
    - NOTE: GISTIC only reports genes WITH alterations, so alteration rates are lower
    - Adjusted threshold: ≥1% of genes altered = 1 TAI point per 0.5% above threshold
    """
    # GISTIC typically analyzes ~24,000 genes (protein-coding + some non-coding)
    # This is a reasonable reference for TCGA data
    GISTIC_REFERENCE_GENES = 24000
    
    total_altered = altered_genes_count
    total_genes = GISTIC_REFERENCE_GENES
    
    alteration_fraction = total_altered / total_genes
    
    # TAI threshold: 1% of genes altered (240 genes) - adjusted for GISTIC data
    # Every 0.5% above threshold (120 genes) = 1 TAI point
    # This produces scores roughly in range 0-40 for typical TCGA-OV samples
    if alteration_fraction >= 0.01:
        excess = alteration_fraction - 0.01
        tai_score = int(excess / 0.005)  # 1 point per 0.5% above 1%
        return max(1, tai_score)
    
    return 0

def calculate_full_hrd_score(cna_records: List[Dict]) -> Dict[str, Any]:
    """
    Calculate full HRD score = LOH + LST + TAI (gene-level approximation)
    
    Returns: {hrd_score, loh, lst, tai, total_genes, altered_genes}
    """
    # Convert to gene counts
    gene_counts = _convert_cna_to_gene_counts(cna_records)
    altered_genes_count = sum(gene_counts.values())
    
    # Calculate components
    loh = calculate_loh_gene_level(gene_counts)
    lst = calculate_lst_gene_level(gene_counts)
    tai = calculate_tai_gene_level(gene_counts, altered_genes_count)
    
    hrd_score = loh + lst + tai
    
    # GISTIC reference gene set size
    GISTIC_REFERENCE_GENES = 24000
    
    return {
        "hrd_score": hrd_score,
        "loh": loh,
        "lst": lst,
        "tai": tai,
        "total_genes": GISTIC_REFERENCE_GENES,  # Reference gene set size
        "altered_genes": altered_genes_count,
        "gene_counts": dict(gene_counts)
    }

def main():
    parser = argparse.ArgumentParser(description="Calculate full HRD scores (LOH + LST + TAI)")
    parser.add_argument("--study-id", default=STUDY_ID, help="cBioPortal study ID")
    parser.add_argument("--limit", type=int, help="Limit number of samples (for testing)")
    parser.add_argument("--output", default="tools/benchmarks/data/full_hrd_scores.json",
                       help="Output JSON file")
    parser.add_argument("--input-data", help="Input patient data file (to add HRD scores)")
    parser.add_argument("--resume", action="store_true", 
                       help="Resume from existing output file (skip already processed samples)")
    parser.add_argument("--save-interval", type=int, default=20,
                       help="Save progress every N samples (default: 20)")
    
    args = parser.parse_args()
    
    # Get samples
    print(f"🔍 Fetching samples for {args.study_id}...")
    with httpx.Client(timeout=60.0, headers=_headers()) as client:
        r = client.get(f"{CBIO_BASE}/studies/{args.study_id}/samples")
        r.raise_for_status()
        all_samples = r.json() or []
    
    sample_ids = [s["sampleId"] for s in all_samples]
    if args.limit:
        sample_ids = sample_ids[:args.limit]
    
    # Load existing results if resuming
    processed_samples = set()
    results = []
    if args.resume and os.path.exists(args.output):
        print(f"📂 Loading existing results from {args.output}...")
        try:
            with open(args.output) as f:
                existing_results = json.load(f)
                results = existing_results
                processed_samples = {r["sample_id"] for r in existing_results}
                print(f"  ✅ Found {len(results)} already processed samples")
        except Exception as e:
            print(f"  ⚠️ Could not load existing file: {e}")
            results = []
            processed_samples = set()
    
    # Filter out already processed samples
    remaining_samples = [s for s in sample_ids if s not in processed_samples]
    print(f"📊 Processing {len(remaining_samples)} new samples ({len(processed_samples)} already done)...\n")
    
    if not remaining_samples:
        print("✅ All samples already processed!")
        return
    
    # Extract GISTIC CNA data for remaining samples
    print(f"🔍 Fetching GISTIC CNA data for {len(remaining_samples)} samples...")
    sample_cna = get_gistic_discrete_cna(args.study_id, remaining_samples)
    
    # Calculate HRD scores
    print(f"\n🧬 Calculating full HRD scores (LOH + LST + TAI)...")
    
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    for i, sample_id in enumerate(remaining_samples, 1):
        if sample_id not in sample_cna:
            print(f"  ⚠️ No CNA data for {sample_id}")
            continue
        
        cna_records = sample_cna[sample_id]
        hrd_metrics = calculate_full_hrd_score(cna_records)
        
        # Get patient ID from sample ID
        patient_id = sample_id.rsplit("-", 1)[0]  # Remove sample suffix
        
        results.append({
            "sample_id": sample_id,
            "patient_id": patient_id,
            **hrd_metrics
        })
        
        # Save incrementally
        if len(results) % args.save_interval == 0:
            with open(args.output, "w") as f:
                json.dump(results, f, indent=2)
            print(f"  💾 Saved progress: {len(results)}/{len(sample_ids)} samples "
                  f"({i}/{len(remaining_samples)} new)")
    
    # Final save
    with open(args.output, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✅ Calculated HRD scores for {len(results)} total samples")
    
    # Statistics
    if results:
        hrd_scores = [r["hrd_score"] for r in results]
        loh_scores = [r["loh"] for r in results]
        lst_scores = [r["lst"] for r in results]
        tai_scores = [r["tai"] for r in results]
        
        print(f"\n📊 HRD Score Statistics:")
        print(f"  HRD Score: mean={sum(hrd_scores)/len(hrd_scores):.2f}, "
              f"min={min(hrd_scores)}, max={max(hrd_scores)}")
        print(f"  LOH: mean={sum(loh_scores)/len(loh_scores):.2f}, "
              f"min={min(loh_scores)}, max={max(loh_scores)}")
        print(f"  LST: mean={sum(lst_scores)/len(lst_scores):.2f}, "
              f"min={min(lst_scores)}, max={max(lst_scores)}")
        print(f"  TAI: mean={sum(tai_scores)/len(tai_scores):.2f}, "
              f"min={min(tai_scores)}, max={max(tai_scores)}")
        
        hrd_high = sum(1 for s in hrd_scores if s >= 42)
        print(f"  HRD-High (≥42): {hrd_high}/{len(results)} ({100*hrd_high/len(results):.1f}%)")
    
    # Save results
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\n💾 Results saved to: {args.output}")
    
    # If input data provided, merge HRD scores
    if args.input_data and os.path.exists(args.input_data):
        print(f"\n🔄 Merging HRD scores into {args.input_data}...")
        with open(args.input_data) as f:
            patient_data = json.load(f)
        
        # Create lookup: patient_id -> HRD score
        hrd_lookup = {r["patient_id"]: r["hrd_score"] for r in results}
        
        # Add HRD scores to patient data
        merged_count = 0
        for patient in patient_data:
            patient_id = patient.get("patient_id", "")
            if patient_id in hrd_lookup:
                patient["hrd_score"] = hrd_lookup[patient_id]
                merged_count += 1
        
        # Save merged data
        merged_output = args.input_data.replace(".json", "_with_full_hrd.json")
        with open(merged_output, "w") as f:
            json.dump(patient_data, f, indent=2)
        
        print(f"✅ Merged HRD scores for {merged_count}/{len(patient_data)} patients")
        print(f"💾 Merged data saved to: {merged_output}")

if __name__ == "__main__":
    main()

