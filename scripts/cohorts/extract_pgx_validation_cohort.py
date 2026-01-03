#!/usr/bin/env python3
"""
Extract PGx validation cohort from cBioPortal
Following cohort_context_concept.mdc framework

Purpose: Validate PGx prevention claims with real-world data
"""

import requests
import json
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Any

# PGx genes (Entrez IDs)
PGX_GENES = {
    "DPYD": 1806,
    "TPMT": 7172,
    "UGT1A1": 54658,
    "CYP2D6": 1565,
}

# PGx-relevant drugs
PGX_DRUGS = {
    "DPYD": ["fluorouracil", "5-FU", "capecitabine", "tegafur"],
    "TPMT": ["mercaptopurine", "6-MP", "azathioprine", "thioguanine"],
    "UGT1A1": ["irinotecan", "camptosar"],
    "CYP2D6": ["tamoxifen", "codeine", "tramadol"],
}

# Output directory
OUTPUT_DIR = Path("oncology-coPilot/oncology-backend-minimal/data/cohorts")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
RECEIPTS_DIR = OUTPUT_DIR / "receipts"
RECEIPTS_DIR.mkdir(exist_ok=True)


def discover_pgx_studies() -> List[Dict[str, Any]]:
    """
    Step 1: Discovery - Find studies with PGx-relevant data
    """
    print("🔍 Step 1: Discovering PGx-relevant studies...")
    
    studies_url = "https://www.cbioportal.org/api/studies"
    resp = requests.get(studies_url, timeout=60)
    all_studies = resp.json()
    
    # Filter for cancer studies (PGx is most relevant in oncology)
    cancer_studies = [
        s for s in all_studies 
        if s.get("cancerTypeId") and s.get("cancerTypeId") != "mixed"
    ]
    
    print(f"  Found {len(cancer_studies)} cancer studies")
    
    # Save discovery receipt
    discovery_receipt = {
        "timestamp": datetime.now().isoformat(),
        "total_studies": len(all_studies),
        "cancer_studies": len(cancer_studies),
        "pgx_genes": list(PGX_GENES.keys()),
        "method": "cBioPortal API /studies endpoint"
    }
    
    receipt_path = RECEIPTS_DIR / f"pgx_discovery_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(receipt_path, 'w') as f:
        json.dump(discovery_receipt, f, indent=2)
    
    print(f"  ✅ Discovery receipt saved: {receipt_path}")
    
    return cancer_studies


def extract_pgx_mutations(study_id: str, gene_entrez_id: int) -> List[Dict[str, Any]]:
    """
    Step 2: Extraction - Get PGx mutations for a study
    Uses mutations/fetch POST endpoint (per cohort_context_concept.mdc)
    """
    mutation_profile = f"{study_id}_mutations"
    
    mutations_url = f"https://www.cbioportal.org/api/molecular-profiles/{mutation_profile}/mutations/fetch"
    
    body = {
        "sampleListId": f"{study_id}_all",
        "entrezGeneIds": [gene_entrez_id]
    }
    
    try:
        resp = requests.post(mutations_url, json=body, timeout=120)
        if resp.status_code == 200:
            return resp.json()
        elif resp.status_code == 404:
            # Try alternative profile name
            alt_profile = f"{study_id}_sequenced"
            mutations_url = f"https://www.cbioportal.org/api/molecular-profiles/{alt_profile}/mutations/fetch"
            resp = requests.post(mutations_url, json=body, timeout=120)
            if resp.status_code == 200:
                return resp.json()
        return []
    except Exception as e:
        print(f"  ⚠️  Error extracting mutations: {e}")
        return []


def extract_clinical_data(study_id: str) -> Dict[str, Any]:
    """
    Step 3: Extraction - Get clinical data (outcomes, treatment)
    """
    clinical_url = f"https://www.cbioportal.org/api/studies/{study_id}/clinical-data"
    
    # Get patient-level clinical data
    params = {"clinicalDataType": "PATIENT", "projection": "DETAILED"}
    
    try:
        resp = requests.get(clinical_url, params=params, timeout=120)
        if resp.status_code == 200:
            return resp.json()
        return []
    except Exception as e:
        print(f"  ⚠️  Error extracting clinical data: {e}")
        return []


def transform_to_cohort_schema(
    mutations: Dict[str, List[Dict]],
    clinical_data: List[Dict],
    study_id: str
) -> Dict[str, Any]:
    """
    Step 4: Transformation - Map to validation schema
    """
    # Group mutations by patient
    patient_mutations = {}
    for gene, muts in mutations.items():
        for m in muts:
            patient_id = m.get("patientId")
            if patient_id:
                if patient_id not in patient_mutations:
                    patient_mutations[patient_id] = {}
                patient_mutations[patient_id][gene] = {
                    "mutation_type": m.get("mutationType"),
                    "protein_change": m.get("proteinChange"),
                    "variant_classification": m.get("variantClassification")
                }
    
    # Group clinical data by patient
    patient_clinical = {}
    for cd in clinical_data:
        patient_id = cd.get("patientId")
        attr_id = cd.get("clinicalAttributeId")
        value = cd.get("value")
        
        if patient_id:
            if patient_id not in patient_clinical:
                patient_clinical[patient_id] = {}
            patient_clinical[patient_id][attr_id] = value
    
    # Build cohort
    all_patients = set(patient_mutations.keys()) | set(patient_clinical.keys())
    
    patients = []
    for patient_id in all_patients:
        patient = {
            "patient_id": patient_id,
            "pgx_variants": patient_mutations.get(patient_id, {}),
            "clinical": patient_clinical.get(patient_id, {})
        }
        patients.append(patient)
    
    cohort = {
        "cohort_id": f"pgx_validation_{study_id}",
        "source": "cbioportal",
        "study_id": study_id,
        "patients": patients,
        "metadata": {
            "extraction_date": datetime.now().isoformat(),
            "n_patients": len(patients),
            "n_with_pgx_variants": len([p for p in patients if p["pgx_variants"]]),
            "pgx_genes_searched": list(PGX_GENES.keys())
        }
    }
    
    return cohort


def validate_cohort(cohort: Dict[str, Any]) -> Dict[str, Any]:
    """
    Step 5: Validation - Quality checks
    """
    validation = {
        "n_patients": cohort["metadata"]["n_patients"],
        "n_with_pgx_variants": cohort["metadata"]["n_with_pgx_variants"],
        "pgx_coverage": cohort["metadata"]["n_with_pgx_variants"] / cohort["metadata"]["n_patients"] if cohort["metadata"]["n_patients"] > 0 else 0,
        "has_outcomes": False,
        "has_treatment": False,
        "validation_status": "PENDING"
    }
    
    # Check for outcome data
    sample_patient = cohort["patients"][0] if cohort["patients"] else {}
    clinical = sample_patient.get("clinical", {})
    
    outcome_keys = ["OS_MONTHS", "OS_STATUS", "PFS_MONTHS", "PFS_STATUS", "DFS_MONTHS"]
    if any(k in clinical for k in outcome_keys):
        validation["has_outcomes"] = True
    
    # Check for treatment data
    treatment_keys = ["DRUG", "TREATMENT", "THERAPY"]
    if any(k in str(clinical).upper() for k in treatment_keys):
        validation["has_treatment"] = True
    
    # Determine validation status
    if validation["pgx_coverage"] > 0.01 and validation["has_outcomes"]:
        validation["validation_status"] = "VALID"
    elif validation["pgx_coverage"] > 0.01:
        validation["validation_status"] = "PARTIAL"
    else:
        validation["validation_status"] = "INSUFFICIENT"
    
    return validation


def main():
    """
    Main extraction pipeline following cohort_context_concept.mdc pattern
    """
    print("=" * 60)
    print("PGx Validation Cohort Extraction")
    print("Following cohort_context_concept.mdc framework")
    print("=" * 60)
    print()
    
    # Step 1: Discovery
    studies = discover_pgx_studies()
    
    # Focus on TCGA PanCancer Atlas (best biomarker coverage)
    target_study = "pancan_tcga_pan_can_atlas_2018"
    
    print(f"\n📊 Step 2: Extracting PGx mutations from {target_study}...")
    
    # Extract mutations for each PGx gene
    all_mutations = {}
    for gene_name, entrez_id in PGX_GENES.items():
        print(f"  Extracting {gene_name} (Entrez: {entrez_id})...")
        muts = extract_pgx_mutations(target_study, entrez_id)
        all_mutations[gene_name] = muts
        print(f"    Found {len(muts)} mutations")
    
    print(f"\n📋 Step 3: Extracting clinical data...")
    clinical_data = extract_clinical_data(target_study)
    print(f"  Found {len(clinical_data)} clinical data points")
    
    print(f"\n🔄 Step 4: Transforming to cohort schema...")
    cohort = transform_to_cohort_schema(all_mutations, clinical_data, target_study)
    
    print(f"\n✅ Step 5: Validating cohort...")
    validation = validate_cohort(cohort)
    
    print(f"\n📊 Validation Results:")
    print(f"  Patients: {validation['n_patients']}")
    print(f"  With PGx variants: {validation['n_with_pgx_variants']}")
    print(f"  PGx coverage: {validation['pgx_coverage']:.2%}")
    print(f"  Has outcomes: {validation['has_outcomes']}")
    print(f"  Status: {validation['validation_status']}")
    
    # Save cohort
    output_path = OUTPUT_DIR / f"pgx_validation_cohort_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(output_path, 'w') as f:
        json.dump(cohort, f, indent=2)
    
    print(f"\n💾 Cohort saved: {output_path}")
    
    # Save validation receipt
    validation_receipt = {
        "timestamp": datetime.now().isoformat(),
        "cohort_file": str(output_path),
        "validation": validation,
        "extraction_method": "cBioPortal API following cohort_context_concept.mdc"
    }
    
    receipt_path = RECEIPTS_DIR / f"pgx_validation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(receipt_path, 'w') as f:
        json.dump(validation_receipt, f, indent=2)
    
    print(f"📋 Validation receipt saved: {receipt_path}")
    
    return cohort, validation


if __name__ == "__main__":
    try:
        cohort, validation = main()
        sys.exit(0 if validation["validation_status"] != "INSUFFICIENT" else 1)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
