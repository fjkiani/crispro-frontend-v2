#!/usr/bin/env python3
"""
Platinum Response Data Hunt - Modular Orchestrator

Coordinates all extraction services to find ≥100 TCGA-OV patients with platinum response labels.
"""

import sys
from pathlib import Path
from typing import Dict
import json
from datetime import datetime
from collections import Counter

# Add services to path
sys.path.insert(0, str(Path(__file__).parent))

# Add parent directory to path for imports
import sys
from pathlib import Path
parent_dir = Path(__file__).parent.parent
sys.path.insert(0, str(parent_dir))

from platinum_hunt.services.gdc_xml_downloader import extract_platinum_response_from_gdc_xml_aggressive
from platinum_hunt.services.tcga_cdr_extractor import extract_platinum_response_from_tcga_cdr
from platinum_hunt.services.pybioportal_treatments_extractor import extract_platinum_response_from_pybioportal
from platinum_hunt.services.broad_firehose_extractor import extract_platinum_response_from_broad_firehose
from extract_platinum_response_labels import (
    get_zo_sample_ids,
    extract_tcga_patient_id,
    search_tcga_ov_original,
    search_pancancer_atlas
)

OUTPUT_DIR = Path("data/validation")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUTPUT_DIR / "tcga_ov_platinum_response_labels.json"

def main():
    print("="*80)
    print("⚔️ PLATINUM RESPONSE DATA HUNT - MODULAR EXECUTION")
    print("="*80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Output: {OUTPUT_FILE}")
    print()
    
    # Get Zo's sample IDs for matching
    zo_sample_ids = get_zo_sample_ids()
    print(f"✅ Target: Match to {len(zo_sample_ids)} Zo's sample IDs")
    
    # Build patient ID → sample ID mapping
    patient_to_samples = {}
    for sample_id in zo_sample_ids:
        patient_id = extract_tcga_patient_id(sample_id)
        if patient_id:
            if patient_id not in patient_to_samples:
                patient_to_samples[patient_id] = []
            patient_to_samples[patient_id].append(sample_id)
    
    print(f"✅ Mapped {len(patient_to_samples)} patient IDs to sample IDs\n")
    
    # Search all sources (priority order)
    all_results = {}
    
    # Source 1: GDC XML files (AGGRESSIVE DOWNLOAD)
    print("🔥 PRIORITY 1: GDC XML Clinical Files (Aggressive Download)")
    import asyncio
    gdc_results = asyncio.run(extract_platinum_response_from_gdc_xml_aggressive())
    all_results.update(gdc_results)
    print(f"   GDC XML: Found {len(gdc_results)} patients")
    print(f"   Total so far: {len(all_results)} patients\n")
    
    # Source 1.5: TCGA-CDR (Pan-Cancer Clinical Data Resource)
    if len(all_results) < 100:
        print("🔥 PRIORITY 1.5: TCGA-CDR (Pan-Cancer Clinical Data Resource)")
        cdr_results = extract_platinum_response_from_tcga_cdr()
        # Merge (GDC takes precedence)
        for patient_id, data in cdr_results.items():
            if patient_id not in all_results:
                all_results[patient_id] = data
        print(f"   TCGA-CDR: Found {len(cdr_results)} patients")
        print(f"   Total so far: {len(all_results)} patients\n")
    
    # Source 2: pyBioPortal Treatments (DEEP EXTRACTION)
    if len(all_results) < 100:
        print("🔥 PRIORITY 2: pyBioPortal Treatments (Deep Extraction)")
        pybioportal_results = extract_platinum_response_from_pybioportal()
        # Merge (GDC takes precedence)
        for patient_id, data in pybioportal_results.items():
            if patient_id not in all_results:
                all_results[patient_id] = data
        print(f"   Total so far: {len(all_results)} patients\n")
    
    # Source 3: Broad Firehose (DEEP EXTRACTION)
    if len(all_results) < 100:
        print("🔥 PRIORITY 3: Broad Firehose (Deep Extraction)")
        broad_results = extract_platinum_response_from_broad_firehose()
        # Merge (earlier sources take precedence)
        for patient_id, data in broad_results.items():
            if patient_id not in all_results:
                all_results[patient_id] = data
        print(f"   Total so far: {len(all_results)} patients\n")
    
    # Source 4: cBioPortal API (SURFACE-LEVEL - fallback)
    if len(all_results) < 100:
        print("📡 FALLBACK: cBioPortal API (Surface-Level)")
        tcga_ov_results = search_tcga_ov_original()
        pancancer_results = search_pancancer_atlas()
        # Merge (earlier sources take precedence)
        for patient_id, data in tcga_ov_results.items():
            if patient_id not in all_results:
                all_results[patient_id] = data
        for patient_id, data in pancancer_results.items():
            if patient_id not in all_results:
                all_results[patient_id] = data
        print(f"   Total so far: {len(all_results)} patients\n")
    
    # Structure output
    print("="*80)
    print("STRUCTURING OUTPUT")
    print("="*80)
    
    patients_list = []
    matched_count = 0
    
    for patient_id, response_data in all_results.items():
        # Find matching sample IDs
        sample_ids = patient_to_samples.get(patient_id, [])
        
        if sample_ids:
            matched_count += 1
            primary_sample_id = sample_ids[0]
        else:
            # Generate sample ID (assume -01 suffix)
            primary_sample_id = f"{patient_id}-01"
        
        # Handle different response key names from different sources
        response_value = response_data.get("platinum_response") or response_data.get("response", "unknown")
        raw_value = response_data.get("raw_response_value") or response_data.get("original_value") or response_data.get("response", "")
        source_field = response_data.get("source_field") or response_data.get("file_id", "unknown")
        
        # Normalize response if needed
        if response_value not in ["sensitive", "resistant", "refractory"]:
            # Import normalize function
            from extract_platinum_response_labels import normalize_response
            response_value = normalize_response(response_value)
        
        patients_list.append({
            "tcga_patient_id": patient_id,
            "tcga_sample_id": primary_sample_id,
            "platinum_response": response_value,
            "raw_response_value": raw_value,
            "source_field": source_field
        })
    
    # Calculate statistics
    response_dist = Counter(p["platinum_response"] for p in patients_list)
    source_dist = Counter(p["source_field"] for p in patients_list)
    
    output = {
        "metadata": {
            "source": "Multi-source modular extraction (GDC XML, pyBioPortal, Broad Firehose, cBioPortal)",
            "extraction_date": datetime.now().isoformat(),
            "n_patients": len(patients_list),
            "n_matched_to_zo_samples": matched_count,
            "match_rate": f"{matched_count / max(len(zo_sample_ids), 1) * 100:.1f}%",
            "response_distribution": dict(response_dist),
            "source_distribution": dict(source_dist),
            "sources_searched": [
                "GDC XML Clinical Files (Deep)",
                "pyBioPortal Treatments (Deep)",
                "Broad Firehose (Deep)",
                "TCGA-OV Original (cBioPortal API)",
                "TCGA PanCancer Atlas (cBioPortal API)"
            ]
        },
        "patients": patients_list
    }
    
    # Save output
    print(f"\n💾 Saving results to {OUTPUT_FILE}...")
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"✅ Saved {len(patients_list)} patients with response labels")
    print(f"✅ Matched {matched_count}/{len(zo_sample_ids)} to Zo's samples ({matched_count/max(len(zo_sample_ids),1)*100:.1f}%)")
    print(f"\n📊 Response Distribution:")
    for response, count in response_dist.most_common():
        print(f"   {response}: {count} ({count/len(patients_list)*100:.1f}%)")
    
    print(f"\n📊 Source Distribution:")
    for source, count in source_dist.most_common():
        print(f"   {source}: {count}")
    
    # Success criteria check
    print("\n" + "="*80)
    print("SUCCESS CRITERIA CHECK")
    print("="*80)
    
    success = True
    if len(patients_list) < 100:
        print(f"   ❌ Minimum N: {len(patients_list)} < 100 (FAIL)")
        success = False
    else:
        print(f"   ✅ Minimum N: {len(patients_list)} ≥ 100 (PASS)")
    
    match_rate = matched_count / max(len(zo_sample_ids), 1)
    if match_rate < 0.80:
        print(f"   ❌ Match Rate: {match_rate*100:.1f}% < 80% (FAIL)")
        success = False
    else:
        print(f"   ✅ Match Rate: {match_rate*100:.1f}% ≥ 80% (PASS)")
    
    if len(response_dist) == 1:
        print(f"   ❌ Response Distribution: Only 1 category (FAIL)")
        success = False
    else:
        print(f"   ✅ Response Distribution: {len(response_dist)} categories (PASS)")
    
    if success:
        print("\n🎉 ALL SUCCESS CRITERIA MET!")
    else:
        print("\n⚠️  SOME SUCCESS CRITERIA NOT MET - Consider additional sources")
    
    print(f"\n📄 Output file: {OUTPUT_FILE}")
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())

