#!/usr/bin/env python3
"""
Extract PTRC-HGSOC Collection Data

PTRC-HGSOC (Proteogenomic Tumor Research Consortium - High-Grade Serous Ovarian Cancer)
Collection from The Cancer Imaging Archive (TCIA)

This script processes downloaded PTRC-HGSOC data:
- MAF files (mutation data)
- Clinical annotations (platinum response labels)
- Normalizes to standard format
"""

import json
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime
import hashlib

# Add project root to path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root))

OUTPUT_DIR = project_root / "data" / "external" / "ov_platinum_non_tcga" / "raw"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

NORMALIZED_OUTPUT = project_root / "data" / "benchmarks" / "external_ov_platinum_ptrc_hgsoc.json"

# Expected data structure
# PTRC-HGSOC typically provides:
# - MAF files (Mutation Annotation Format)
# - Clinical data files (CSV/TSV)
# - Platinum response labels (refractory/sensitive)


def compute_checksum(file_path: Path) -> str:
    """Compute SHA256 checksum of a file."""
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def parse_maf_file(maf_path: Path) -> List[Dict[str, Any]]:
    """
    Parse MAF (Mutation Annotation Format) file.
    
    MAF format typically has columns:
    - Hugo_Symbol (gene)
    - Chromosome
    - Start_Position
    - End_Position
    - Reference_Allele
    - Tumor_Seq_Allele2 (variant allele)
    - Tumor_Sample_Barcode (sample ID)
    - Variant_Classification
    - HGVSp_Short (protein change)
    """
    mutations = []
    
    if not maf_path.exists():
        print(f"   ⚠️  MAF file not found: {maf_path}")
        return mutations
    
    print(f"   📥 Parsing MAF file: {maf_path.name}")
    
    try:
        with open(maf_path, "r", encoding="utf-8") as f:
            # MAF files can have comment lines starting with #
            lines = f.readlines()
            
            # Find header line
            header_idx = None
            for i, line in enumerate(lines):
                if line.startswith("Hugo_Symbol") or line.startswith("Gene"):
                    header_idx = i
                    break
            
            if header_idx is None:
                print(f"   ⚠️  Could not find MAF header")
                return mutations
            
            # Parse header
            header_line = lines[header_idx].strip()
            delimiter = "\t" if "\t" in header_line else ","
            headers = [h.strip() for h in header_line.split(delimiter)]
            
            # Parse data rows
            for line in lines[header_idx + 1:]:
                if line.strip().startswith("#") or not line.strip():
                    continue
                
                values = [v.strip() for v in line.split(delimiter)]
                if len(values) < len(headers):
                    continue
                
                row = dict(zip(headers, values))
                
                # Extract mutation data
                mutation = {
                    "gene": row.get("Hugo_Symbol", row.get("Gene", "")),
                    "chromosome": row.get("Chromosome", "").replace("chr", "").replace("Chr", ""),
                    "start_position": int(row.get("Start_Position", 0)) if row.get("Start_Position") else None,
                    "end_position": int(row.get("End_Position", 0)) if row.get("End_Position") else None,
                    "reference_allele": row.get("Reference_Allele", ""),
                    "variant_allele": row.get("Tumor_Seq_Allele2", row.get("Tumor_Allele", "")),
                    "protein_change": row.get("HGVSp_Short", row.get("Protein_Change", "")),
                    "variant_type": row.get("Variant_Classification", ""),
                    "sample_id": row.get("Tumor_Sample_Barcode", row.get("Sample_ID", "")),
                    "patient_id": row.get("Matched_Norm_Sample_Barcode", row.get("Patient_ID", ""))
                }
                
                # Extract patient ID from sample ID if needed
                if not mutation["patient_id"] and mutation["sample_id"]:
                    # Sample IDs often have format: PATIENT-SAMPLE
                    parts = mutation["sample_id"].split("-")
                    if len(parts) >= 2:
                        mutation["patient_id"] = "-".join(parts[:2])
                
                mutations.append(mutation)
        
        print(f"   ✅ Parsed {len(mutations)} mutations from MAF")
        return mutations
    
    except Exception as e:
        print(f"   ❌ Error parsing MAF file: {e}")
        return mutations


def parse_clinical_file(clinical_path: Path) -> Dict[str, Dict[str, Any]]:
    """
    Parse clinical data file (CSV/TSV).
    
    Expected columns:
    - Patient ID
    - Platinum response (refractory/sensitive/resistant)
    - Other clinical data
    """
    clinical_data = {}
    
    if not clinical_path.exists():
        print(f"   ⚠️  Clinical file not found: {clinical_path}")
        return clinical_data
    
    print(f"   📥 Parsing clinical file: {clinical_path.name}")
    
    try:
        with open(clinical_path, "r", encoding="utf-8") as f:
            # Detect delimiter
            first_line = f.readline()
            delimiter = "\t" if "\t" in first_line else ","
            f.seek(0)
            
            reader = csv.DictReader(f, delimiter=delimiter)
            
            for row in reader:
                # Find patient ID column (various possible names)
                patient_id = None
                for key in ["Patient_ID", "PatientID", "patient_id", "PATIENT_ID", "Sample_ID"]:
                    if key in row and row[key]:
                        patient_id = row[key].strip()
                        break
                
                if not patient_id:
                    continue
                
                # Extract platinum response
                platinum_response = None
                for key in row.keys():
                    key_lower = key.lower()
                    if "platinum" in key_lower or "response" in key_lower or "sensitive" in key_lower or "refractory" in key_lower:
                        value = row[key].strip() if row[key] else ""
                        if value:
                            # Map to standard format
                            value_lower = value.lower()
                            if "refractory" in value_lower or "resistant" in value_lower:
                                platinum_response = "refractory"
                            elif "sensitive" in value_lower or "responder" in value_lower:
                                platinum_response = "sensitive"
                            break
                
                # Store all clinical data
                clinical_data[patient_id] = {
                    "platinum_response": platinum_response,
                    "all_fields": dict(row)
                }
        
        print(f"   ✅ Parsed clinical data for {len(clinical_data)} patients")
        if any(d["platinum_response"] for d in clinical_data.values()):
            response_count = sum(1 for d in clinical_data.values() if d["platinum_response"])
            print(f"   ✅ Found platinum response labels for {response_count} patients")
        
        return clinical_data
    
    except Exception as e:
        print(f"   ❌ Error parsing clinical file: {e}")
        return clinical_data


def normalize_ptrc_hgsoc_data(
    mutations: List[Dict],
    clinical_data: Dict[str, Dict],
    study_id: str = "external_ov_platinum_ptrc_hgsoc"
) -> Dict[str, Any]:
    """
    Normalize PTRC-HGSOC data to standard format.
    """
    print(f"\n📊 Normalizing PTRC-HGSOC data...")
    
    # Organize mutations by patient
    mutations_by_patient = {}
    for mut in mutations:
        patient_id = mut.get("patient_id", "")
        if not patient_id:
            # Try to extract from sample ID
            sample_id = mut.get("sample_id", "")
            if sample_id:
                parts = sample_id.split("-")
                patient_id = "-".join(parts[:2]) if len(parts) >= 2 else sample_id
        
        if not patient_id:
            continue
        
        if patient_id not in mutations_by_patient:
            mutations_by_patient[patient_id] = []
        
        # Clean mutation dict (remove None values)
        clean_mut = {k: v for k, v in mut.items() if v is not None and v != ""}
        if clean_mut:
            mutations_by_patient[patient_id].append(clean_mut)
    
    # Build normalized dataset
    normalized = {
        "study_id": study_id,
        "source": "PTRC-HGSOC (TCIA)",
        "extraction_date": datetime.now().isoformat(),
        "patients": []
    }
    
    # Get all unique patients
    all_patient_ids = set(mutations_by_patient.keys()) | set(clinical_data.keys())
    
    platinum_response_count = 0
    
    for patient_id in sorted(all_patient_ids):
        patient_mutations = mutations_by_patient.get(patient_id, [])
        clinical = clinical_data.get(patient_id, {})
        
        # Extract platinum response
        platinum_response = clinical.get("platinum_response") if isinstance(clinical, dict) else None
        if platinum_response:
            platinum_response_count += 1
        
        # Build patient record
        patient_record = {
            "patient_id": patient_id,
            "mutations": patient_mutations,
            "platinum_response": platinum_response,
            "clinical_outcomes": clinical.get("all_fields", {}) if isinstance(clinical, dict) else {}
        }
        
        normalized["patients"].append(patient_record)
    
    print(f"   ✅ Normalized {len(normalized['patients'])} patients")
    print(f"   ✅ Found platinum response labels for {platinum_response_count} patients")
    print(f"   ✅ Total mutations: {sum(len(p['mutations']) for p in normalized['patients'])}")
    
    return normalized


def main():
    """Main extraction function."""
    print("="*80)
    print("PTRC-HGSOC COLLECTION EXTRACTION")
    print("="*80)
    print()
    print("⚠️  NOTE: This script processes downloaded PTRC-HGSOC data.")
    print("   You must first download the data from TCIA:")
    print("   1. Register at https://www.cancerimagingarchive.net/registration/")
    print("   2. Access https://www.cancerimagingarchive.net/collection/ptrc-hgsoc/")
    print("   3. Download MAF files and clinical data")
    print("   4. Place files in: data/external/ov_platinum_non_tcga/raw/ptrc_hgsoc/")
    print()
    
    # Expected data directory
    data_dir = OUTPUT_DIR / "ptrc_hgsoc"
    
    if not data_dir.exists():
        print(f"❌ Data directory not found: {data_dir}")
        print(f"   Please download PTRC-HGSOC data and place it in this directory.")
        print()
        print("Expected structure:")
        print(f"   {data_dir}/")
        print(f"   {data_dir}/maf/  (MAF files)")
        print(f"   {data_dir}/clinical/  (Clinical data files)")
        return
    
    print(f"📁 Data directory: {data_dir}")
    print()
    
    # Find MAF files
    maf_dir = data_dir / "maf"
    maf_files = list(maf_dir.glob("*.maf")) + list(maf_dir.glob("*.txt")) if maf_dir.exists() else []
    
    if not maf_files:
        print("⚠️  No MAF files found. Searching in data directory...")
        maf_files = list(data_dir.glob("*.maf")) + list(data_dir.glob("*maf*.txt"))
    
    print(f"📥 Found {len(maf_files)} MAF file(s)")
    
    # Parse all MAF files
    all_mutations = []
    for maf_file in maf_files:
        mutations = parse_maf_file(maf_file)
        all_mutations.extend(mutations)
    
    print(f"✅ Total mutations parsed: {len(all_mutations)}")
    
    # Find clinical files
    clinical_dir = data_dir / "clinical"
    clinical_files = list(clinical_dir.glob("*.csv")) + list(clinical_dir.glob("*.tsv")) + list(clinical_dir.glob("*.txt")) if clinical_dir.exists() else []
    
    if not clinical_files:
        print("⚠️  No clinical files found in clinical directory. Searching in data directory...")
        clinical_files = list(data_dir.glob("*.csv")) + list(data_dir.glob("*.tsv"))
    
    print(f"\n📥 Found {len(clinical_files)} clinical file(s)")
    
    # Parse all clinical files
    all_clinical_data = {}
    for clinical_file in clinical_files:
        clinical_data = parse_clinical_file(clinical_file)
        # Merge clinical data (later files override earlier ones)
        for patient_id, data in clinical_data.items():
            all_clinical_data[patient_id] = data
    
    print(f"✅ Total patients with clinical data: {len(all_clinical_data)}")
    
    # Normalize data
    normalized = normalize_ptrc_hgsoc_data(all_mutations, all_clinical_data)
    
    # Save normalized dataset
    print(f"\n💾 Saving normalized dataset...")
    with open(NORMALIZED_OUTPUT, "w") as f:
        json.dump(normalized, f, indent=2)
    
    print(f"✅ Saved: {NORMALIZED_OUTPUT}")
    
    # Summary
    print("\n" + "="*80)
    print("EXTRACTION SUMMARY")
    print("="*80)
    print(f"✅ Patients: {len(normalized['patients'])}")
    print(f"✅ Mutations: {sum(len(p['mutations']) for p in normalized['patients'])}")
    print(f"✅ Platinum response labels: {sum(1 for p in normalized['patients'] if p.get('platinum_response'))}")
    print(f"✅ Normalized dataset: {NORMALIZED_OUTPUT}")
    print()
    
    # Check requirements
    print("Requirements Check:")
    patient_count = len(normalized['patients'])
    platinum_count = sum(1 for p in normalized['patients'] if p.get('platinum_response'))
    mutations_with_coords = sum(
        1 for p in normalized['patients']
        for m in p['mutations']
        if m.get('chromosome') and m.get('start_position')
    )
    
    print(f"   N ≥ 50 patients: {'✅' if patient_count >= 50 else '❌'} ({patient_count})")
    print(f"   Platinum response labels: {'✅' if platinum_count > 0 else '❌'} ({platinum_count})")
    print(f"   Mutations with coordinates: {'✅' if mutations_with_coords > 0 else '❌'} ({mutations_with_coords})")
    print()


if __name__ == "__main__":
    main()

