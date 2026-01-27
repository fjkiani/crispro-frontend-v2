#!/usr/bin/env python3
"""
Test Post-Treatment SAE on GSE165897

CRITICAL REQUIREMENT: We need POST-TREATMENT variant data (WES/WGS), not scRNA-seq.

GSE165897 Status:
- ✅ Has paired pre/post-NACT scRNA-seq (expression data)
- ✅ Has resistance labels (PFI < 180 days = resistant)
- ❌ Does NOT have variant data (WES/WGS) for SAE extraction

This script:
1. Checks if variant data exists for GSE165897
2. If not, documents the requirement
3. Suggests alternative datasets with post-treatment variants
"""

import json
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd

# Known GSE165897 patient resistance labels
RESISTANCE_LABELS = {
    "EOC1005": "resistant",  # PFI: 65 days
    "EOC136": "sensitive",   # PFI: 520 days
    "EOC153": "sensitive",   # PFI: 393 days
    "EOC227": "resistant",   # PFI: 230 days
    "EOC3": "resistant",     # PFI: 14 days
    "EOC349": "resistant",   # PFI: 36 days
    "EOC372": "sensitive",   # PFI: 460 days
    "EOC443": "resistant",   # PFI: 177 days
    "EOC540": "resistant",   # PFI: 126 days
    "EOC733": "resistant",   # PFI: 83 days
    "EOC87": "resistant",    # PFI: 30 days
}

def check_variant_data_availability() -> Dict[str, any]:
    """
    Check if variant data (WES/WGS) is available for GSE165897.
    
    Returns:
        Status dict with availability info
    """
    print("=" * 80)
    print("GSE165897 Post-Treatment SAE Test - Data Availability Check")
    print("=" * 80)
    
    # Check known data locations
    data_dir = Path("scripts/data_acquisition/sae")
    ega_info = {
        "dataset_id": "EGAS00001005010",
        "access": "Controlled (requires EGA DAC approval)",
        "url": "https://ega-archive.org/datasets/EGAD00001005010",
        "data_type": "Raw FASTQ files (scRNA-seq, not WES/WGS)"
    }
    
    # Check if we have any variant files
    variant_files = list(data_dir.glob("*GSE165897*variant*"))
    variant_files.extend(list(data_dir.glob("*GSE165897*WES*")))
    variant_files.extend(list(data_dir.glob("*GSE165897*WGS*")))
    variant_files.extend(list(data_dir.glob("*GSE165897*mutation*")))
    
    status = {
        "has_variant_data": len(variant_files) > 0,
        "variant_files_found": [str(f) for f in variant_files],
        "has_scRNA_seq": True,  # We know this exists
        "has_resistance_labels": True,
        "n_patients": len(RESISTANCE_LABELS),
        "n_resistant": sum(1 for v in RESISTANCE_LABELS.values() if v == "resistant"),
        "n_sensitive": sum(1 for v in RESISTANCE_LABELS.values() if v == "sensitive"),
        "ega_info": ega_info,
        "recommendation": None
    }
    
    if not status["has_variant_data"]:
        status["recommendation"] = (
            "GSE165897 does not have variant data (WES/WGS) for SAE extraction. "
            "The dataset contains scRNA-seq (expression) data only. "
            "To test post-treatment SAE, we need:\n"
            "1. Post-treatment WES/WGS variant data from GSE165897 (if available in EGA)\n"
            "2. OR alternative dataset with paired pre/post-treatment variants\n"
            "3. OR extract variants from scRNA-seq (not recommended, low quality)"
        )
    
    return status

def document_requirement():
    """Document the requirement for post-treatment variant data."""
    
    status = check_variant_data_availability()
    
    print("\n📊 Dataset Status:")
    print(f"  Patients: {status['n_patients']}")
    print(f"  Resistant: {status['n_resistant']}")
    print(f"  Sensitive: {status['n_sensitive']}")
    print(f"  Has variant data: {status['has_variant_data']}")
    
    if status['variant_files_found']:
        print(f"\n✅ Variant files found:")
        for f in status['variant_files_found']:
            print(f"  - {f}")
    else:
        print("\n❌ No variant data files found")
    
    print(f"\n📋 EGA Controlled Access:")
    print(f"  Dataset: {status['ega_info']['dataset_id']}")
    print(f"  Access: {status['ega_info']['access']}")
    print(f"  URL: {status['ega_info']['url']}")
    print(f"  Data type: {status['ega_info']['data_type']}")
    
    if status['recommendation']:
        print(f"\n⚠️  RECOMMENDATION:")
        print(status['recommendation'])
    
    # Save status
    output_path = Path("data/serial_sae/gse165897/gse165897_sae_requirement.json")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump(status, f, indent=2)
    
    print(f"\n💾 Status saved to: {output_path}")
    
    return status

if __name__ == "__main__":
    document_requirement()
