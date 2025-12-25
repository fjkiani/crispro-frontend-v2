#!/usr/bin/env python3
"""
Mutation Processing Utilities
==============================
Modular utilities for processing and extracting mutation data.

Author: Agent
Date: January 28, 2025
"""

from typing import List, Dict, Set, Optional

# DDR pathway genes (for flagging)
DDR_GENES = {
    "BRCA1", "BRCA2", "ATM", "ATR", "CHEK2", "RAD51", "PALB2", 
    "RAD51C", "RAD51D", "BARD1", "NBN", "BRIP1", "MBD4"
}

# dMMR genes for MSI prediction
MMR_GENES = {"MLH1", "MSH2", "MSH6", "PMS2", "EPCAM"}

# Special genes of interest
TP53_GENE = "TP53"
MBD4_GENE = "MBD4"


def extract_mutation_fields(mutation: Dict) -> Dict:
    """
    Extract required mutation fields from cBioPortal mutation object.
    
    Args:
        mutation: Raw mutation dict from cBioPortal
    
    Returns:
        Standardized mutation dict with required fields
    """
    return {
        "gene": mutation.get("gene", ""),
        "hgvs_p": mutation.get("proteinChange", "") or mutation.get("hgvs_p", ""),
        "variant_classification": mutation.get("variantType", "") or mutation.get("variant_classification", ""),
        "chromosome": mutation.get("chromosome", ""),
        "start_position": mutation.get("startPosition", "") or mutation.get("start_position", ""),
        "reference_allele": mutation.get("referenceAllele", "") or mutation.get("reference_allele", ""),
        "variant_allele": mutation.get("variantAllele", "") or mutation.get("variant_allele", ""),
        # Optional fields
        "hgvs_c": mutation.get("hgvs_c", "") or mutation.get("cDNAChange", ""),
        "protein_change": mutation.get("proteinChange", "") or mutation.get("protein_change", ""),
        "tumor_sample_barcode": mutation.get("sampleId", "") or mutation.get("tumor_sample_barcode", "")
    }


def process_mutations(mutations: List[Dict]) -> List[Dict]:
    """
    Process raw mutations into standardized format.
    
    Args:
        mutations: List of raw mutation dicts from cBioPortal
    
    Returns:
        List of standardized mutation dicts
    """
    return [extract_mutation_fields(mut) for mut in mutations if mut.get("gene")]


def get_mutated_genes(mutations: List[Dict]) -> Set[str]:
    """Extract set of mutated gene names from mutation list."""
    return {mut.get("gene", "").upper() for mut in mutations if mut.get("gene")}


def has_ddr_mutation(mutations: List[Dict]) -> bool:
    """Check if patient has any DDR pathway mutation."""
    genes = get_mutated_genes(mutations)
    return bool(genes & DDR_GENES)


def has_mmr_mutation(mutations: List[Dict]) -> bool:
    """Check if patient has any dMMR gene mutation (MSI-related)."""
    genes = get_mutated_genes(mutations)
    return bool(genes & MMR_GENES)


def has_tp53_mutation(mutations: List[Dict]) -> bool:
    """Check if patient has TP53 mutation."""
    genes = get_mutated_genes(mutations)
    return TP53_GENE.upper() in genes


def has_mbd4_mutation(mutations: List[Dict]) -> bool:
    """Check if patient has MBD4 mutation (Ayesha-relevant)."""
    genes = get_mutated_genes(mutations)
    return MBD4_GENE.upper() in genes


def find_primary_tumor_sample(samples: List[Dict], patient_id: str) -> Optional[str]:
    """
    Find primary tumor sample ID for a patient.
    Prefers -01 (primary tumor) over -02 (recurrence) or -06 (metastatic).
    
    Args:
        samples: List of sample dicts from cBioPortal
        patient_id: Patient ID (e.g., "TCGA-XX-XXXX")
    
    Returns:
        Sample ID (e.g., "TCGA-XX-XXXX-01") or None
    """
    patient_samples = [s for s in samples if s.get("patientId") == patient_id]
    
    if not patient_samples:
        return None
    
    # Prefer primary tumor (-01)
    for sample in patient_samples:
        sample_id = sample.get("sampleId", "")
        if sample_id.endswith("-01"):
            return sample_id
    
    # Fallback to any sample
    return patient_samples[0].get("sampleId")

