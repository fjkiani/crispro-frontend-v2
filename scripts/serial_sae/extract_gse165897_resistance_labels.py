#!/usr/bin/env python3
"""
Extract and create resistance labels for GSE165897 based on PFI data.

PFI < 6 months (< 180 days) = resistant
PFI >= 6 months (>= 180 days) = sensitive
"""

import json
import pandas as pd
from pathlib import Path
from typing import Dict, Optional

# Known PFI values from user's investigation
# Format: patient_id: (pfi_days, source)
KNOWN_PFI = {
    "EOC1005": (126, "user_investigation"),  # < 180 = resistant
    "EOC87": (274, "user_investigation"),    # >= 180 = sensitive
    "EOC136": (210, "user_investigation"),  # >= 180 = sensitive (borderline)
    # Remaining patients need to be extracted from Table S1
    "EOC153": None,
    "EOC227": None,
    "EOC3": None,
    "EOC349": None,
    "EOC372": None,
    "EOC443": None,
    "EOC540": None,
    "EOC733": None,
}

# Median PFI: 4.2 months (127 days) - predominantly resistant cohort
# This suggests most remaining patients are likely resistant

def classify_resistance(pfi_days: Optional[int]) -> Optional[str]:
    """Classify resistance based on PFI."""
    if pfi_days is None:
        return None
    return "resistant" if pfi_days < 180 else "sensitive"

def create_resistance_labels() -> Dict[str, Optional[str]]:
    """Create resistance labels dictionary."""
    labels = {}
    
    for patient_id, pfi_data in KNOWN_PFI.items():
        if pfi_data is None:
            labels[patient_id] = None
        else:
            pfi_days, source = pfi_data
            labels[patient_id] = classify_resistance(pfi_days)
    
    return labels

def main():
    print("="*80)
    print("GSE165897 Resistance Labels Extraction")
    print("="*80)
    
    labels = create_resistance_labels()
    
    # Count classifications
    resistant_count = sum(1 for v in labels.values() if v == "resistant")
    sensitive_count = sum(1 for v in labels.values() if v == "sensitive")
    missing_count = sum(1 for v in labels.values() if v is None)
    
    print(f"\nResistance Classifications:")
    print(f"  Resistant (PFI < 180 days): {resistant_count}")
    print(f"  Sensitive (PFI >= 180 days): {sensitive_count}")
    print(f"  Missing: {missing_count}")
    
    print(f"\nKnown Classifications:")
    for patient_id, label in sorted(labels.items()):
        if label:
            pfi_data = KNOWN_PFI.get(patient_id)
            if pfi_data:
                pfi_days, source = pfi_data
                print(f"  {patient_id}: {label} (PFI: {pfi_days} days, source: {source})")
            else:
                print(f"  {patient_id}: {label}")
        else:
            print(f"  {patient_id}: MISSING (need PFI from Table S1)")
    
    # Save JSON file
    output_path = Path("data/serial_sae/gse165897/resistance_labels.json")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Add metadata
    labels_with_meta = {
        "_comment": "Resistance labels for GSE165897 based on PFI data",
        "_classification": "PFI < 6 months (< 180 days) = resistant, PFI >= 6 months (>= 180 days) = sensitive",
        "_source": "DECIDER cohort (Lahtinen et al., 2023), Table S1 from Science Advances supplementary",
        "_median_pfi": "4.2 months (127 days) - predominantly resistant cohort",
        "_known_values": {
            "EOC1005": {"pfi_days": 126, "classification": "resistant"},
            "EOC87": {"pfi_days": 274, "classification": "sensitive"},
            "EOC136": {"pfi_days": 210, "classification": "sensitive"},
        },
        **labels
    }
    
    with open(output_path, "w") as f:
        json.dump(labels_with_meta, f, indent=2)
    
    print(f"\n💾 Saved resistance labels to: {output_path}")
    
    # Also create CSV version
    csv_path = Path("data/serial_sae/gse165897/resistance_labels.csv")
    csv_data = []
    for patient_id, label in sorted(labels.items()):
        pfi_data = KNOWN_PFI.get(patient_id)
        pfi_days = pfi_data[0] if pfi_data else None
        csv_data.append({
            "patient_id": patient_id,
            "pfi_days": pfi_days,
            "resistance_label": label
        })
    
    csv_df = pd.DataFrame(csv_data)
    csv_df.to_csv(csv_path, index=False)
    print(f"💾 Saved CSV version to: {csv_path}")
    
    if missing_count > 0:
        print(f"\n⚠️  Warning: {missing_count} patients missing PFI data.")
        print("   Need to extract from Table S1 in Science Advances supplementary PDF.")
        print("   Or access GenomeSpy tool: https://csbi.ltdk.helsinki.fi/p/lahtinen_et_al_2022/")
    
    print("\n" + "="*80)
    print("Next: Re-run pathway_kinetics_gse165897.py to generate resistance-stratified analysis")
    print("="*80)

if __name__ == "__main__":
    main()
