#!/usr/bin/env python3
"""
Check overlap between Jr2's platinum response data and Zo's mutation dataset.
"""

import json
from pathlib import Path

def main():
    # Load datasets
    jr_file = Path("data/validation/tcga_ov_platinum_response_labels.json")
    zo_file = Path("data/validation/tcga_ov_full_validation_dataset.json")
    
    if not jr_file.exists():
        print(f"❌ Jr2's file not found: {jr_file}")
        return
    
    if not zo_file.exists():
        print(f"❌ Zo's file not found: {zo_file}")
        return
    
    print("="*80)
    print("🔍 CHECKING OVERLAP BETWEEN DATASETS")
    print("="*80)
    print()
    
    # Load Jr2's data
    with open(jr_file, 'r') as f:
        jr_data = json.load(f)
    
    jr_patients = jr_data.get("patients", [])
    print(f"📊 Jr2's Dataset:")
    print(f"   Total patients: {len(jr_patients)}")
    print(f"   Response distribution: {jr_data.get('metadata', {}).get('response_distribution', {})}")
    print()
    
    # Load Zo's data
    with open(zo_file, 'r') as f:
        zo_data = json.load(f)
    
    zo_patients = zo_data.get("patients", [])
    print(f"📊 Zo's Dataset:")
    print(f"   Total patients: {len(zo_patients)}")
    print(f"   Patients with mutations: {sum(1 for p in zo_patients if p.get('mutations'))}")
    print()
    
    # Extract sample IDs
    jr_sample_ids = {p.get("tcga_sample_id") for p in jr_patients if p.get("tcga_sample_id")}
    zo_sample_ids = {p.get("sample_id") for p in zo_patients if p.get("sample_id")}
    
    print(f"🔑 Sample ID Sets:")
    print(f"   Jr2 sample IDs: {len(jr_sample_ids)}")
    print(f"   Zo sample IDs: {len(zo_sample_ids)}")
    print()
    
    # Find exact matches
    exact_overlap = jr_sample_ids & zo_sample_ids
    print(f"✅ Exact Sample ID Overlap: {len(exact_overlap)} patients")
    print(f"   Match rate: {len(exact_overlap)/len(jr_sample_ids)*100:.1f}%")
    print()
    
    # Try patient ID matching (TCGA-XX-XXXX without -01 suffix)
    jr_patient_ids = set()
    for sample_id in jr_sample_ids:
        if sample_id:
            # Extract patient ID (remove -01, -02, etc.)
            patient_id = sample_id.split('-01')[0] if '-01' in sample_id else sample_id
            jr_patient_ids.add(patient_id)
    
    zo_patient_ids = set()
    for sample_id in zo_sample_ids:
        if sample_id:
            patient_id = sample_id.split('-01')[0] if '-01' in sample_id else sample_id
            zo_patient_ids.add(patient_id)
    
    patient_overlap = jr_patient_ids & zo_patient_ids
    print(f"🔑 Patient ID Sets (without sample suffix):")
    print(f"   Jr2 patient IDs: {len(jr_patient_ids)}")
    print(f"   Zo patient IDs: {len(zo_patient_ids)}")
    print(f"   Patient ID Overlap: {len(patient_overlap)} patients")
    print(f"   Match rate: {len(patient_overlap)/len(jr_patient_ids)*100:.1f}%")
    print()
    
    # Show some examples
    if exact_overlap:
        print(f"📋 Example Exact Matches (first 5):")
        for sample_id in list(exact_overlap)[:5]:
            print(f"   ✅ {sample_id}")
        print()
    
    if patient_overlap and len(patient_overlap) > len(exact_overlap):
        print(f"📋 Example Patient ID Matches (first 5):")
        for patient_id in list(patient_overlap)[:5]:
            # Find corresponding sample IDs
            jr_samples = [s for s in jr_sample_ids if s and s.startswith(patient_id)]
            zo_samples = [s for s in zo_sample_ids if s and s.startswith(patient_id)]
            print(f"   ✅ {patient_id}: Jr2={jr_samples}, Zo={zo_samples}")
        print()
    
    # Recommendation
    print("="*80)
    print("🎯 RECOMMENDATION")
    print("="*80)
    
    if len(exact_overlap) >= 40:
        print(f"✅ SUFFICIENT OVERLAP: {len(exact_overlap)} patients")
        print("   → Proceed with validation using exact sample ID matches")
    elif len(patient_overlap) >= 40:
        print(f"⚠️  PARTIAL OVERLAP: {len(patient_overlap)} patients (patient ID level)")
        print("   → Proceed with validation using patient ID matching")
        print("   → May need to handle multiple samples per patient")
    else:
        print(f"❌ INSUFFICIENT OVERLAP: {len(exact_overlap)} exact, {len(patient_overlap)} patient-level")
        print("   → Need ≥40 patients for statistical validation")
        print("   → RECOMMENDED: Re-run GDC extraction (early exit removed)")
        print("   → Expected: 200-300 patients with response data (if all 597 files processed)")
    
    print()

if __name__ == "__main__":
    main()


