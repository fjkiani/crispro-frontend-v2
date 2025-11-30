#!/usr/bin/env python3
"""
🗺️ Create SAE Feature→Pathway Mapping
======================================

Creates feature→pathway mapping file from biomarker analysis results.

Strategy: Gene→Pathway Inference
1. Load biomarker analysis results (top significant features)
2. For each feature, find which patients have it activated
3. Extract genes/mutations from those patients
4. Use gene→pathway mapping to infer pathway
5. Build mapping JSON file

Usage:
    python3 scripts/sae/create_feature_pathway_mapping.py \
        --biomarkers data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers.json \
        --cohort data/validation/sae_cohort/sae_features_tcga_ov_platinum.json \
        --output api/resources/sae_feature_mapping.json \
        --top-n 100

Inputs:
    - Biomarker analysis results (top features with statistics)
    - SAE cohort file (patient features and mutations)

Output:
    - sae_feature_mapping.json (ready for use in SAE service)
"""

import json
import argparse
from pathlib import Path
from typing import Dict, List, Any, Set
from collections import defaultdict, Counter
from datetime import datetime
from loguru import logger
import sys

# Gene→Pathway Mappings (from existing codebase)
GENE_TO_PATHWAY = {
    # DDR pathway
    "BRCA1": ["ddr"],
    "BRCA2": ["ddr"],
    "RAD51": ["ddr"],
    "ATM": ["ddr", "tp53"],
    "ATR": ["ddr", "tp53"],
    "CHEK2": ["ddr", "tp53"],
    "PALB2": ["ddr"],
    "FANCA": ["ddr"],
    "FANCD2": ["ddr"],
    "PARP1": ["ddr"],
    "XRCC1": ["ddr"],
    
    # MAPK pathway
    "KRAS": ["mapk"],
    "NRAS": ["mapk"],
    "BRAF": ["mapk"],
    "MAP2K1": ["mapk"],
    "MAP2K2": ["mapk"],
    "MAPK1": ["mapk"],
    "MAPK3": ["mapk"],
    "RAF1": ["mapk"],
    
    # PI3K pathway
    "PIK3CA": ["pi3k"],
    "PIK3CB": ["pi3k"],
    "AKT1": ["pi3k"],
    "AKT2": ["pi3k"],
    "MTOR": ["pi3k"],
    "PTEN": ["pi3k", "tp53"],
    "TSC1": ["pi3k"],
    "TSC2": ["pi3k"],
    
    # VEGF pathway
    "VEGFA": ["vegf"],
    "VEGFR1": ["vegf"],
    "VEGFR2": ["vegf"],
    "HIF1A": ["vegf"],
    "FLT1": ["vegf"],
    "KDR": ["vegf"],
    
    # HER2 pathway
    "ERBB2": ["her2"],
    "HER2": ["her2"],
    
    # TP53 pathway
    "TP53": ["tp53"],
    "MDM2": ["tp53"],
    "MDM4": ["tp53"],
    "CDKN1A": ["tp53"],
    "CDKN2A": ["tp53"],
    "RB1": ["tp53"],
    "E2F1": ["tp53"],
}

# Pathway metadata (from template)
PATHWAY_METADATA = {
    "ddr": {
        "name": "DNA Damage Response & Homologous Recombination",
        "description": "BRCA1/2-mediated homologous recombination repair, DNA damage signaling, platinum sensitivity mechanism",
        "aliases": ["DNA repair", "HR", "HRD", "BRCA"],
        "genes_in_pathway": ["BRCA1", "BRCA2", "RAD51", "ATM", "ATR", "CHEK1", "CHEK2", "PALB2", "FANCA", "FANCD2"],
        "drugs_targeting": ["Platinum (cisplatin, carboplatin)", "PARP inhibitors (olaparib, niraparib, rucaparib)"]
    },
    "mapk": {
        "name": "MAPK/RAS/RAF Signaling",
        "description": "Oncogenic RAS→RAF→MEK→ERK cascade, resistance mechanism, growth signaling",
        "aliases": ["RAS pathway", "MAPK cascade", "ERK signaling"],
        "genes_in_pathway": ["KRAS", "NRAS", "BRAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3"],
        "drugs_targeting": ["MEK inhibitors (trametinib, cobimetinib)", "RAF inhibitors (dabrafenib, vemurafenib)"]
    },
    "pi3k": {
        "name": "PI3K/AKT/mTOR Pathway",
        "description": "Growth factor signaling, apoptosis resistance, metabolic reprogramming",
        "aliases": ["PI3K", "AKT", "mTOR", "PTEN"],
        "genes_in_pathway": ["PIK3CA", "PIK3CB", "AKT1", "AKT2", "MTOR", "PTEN", "TSC1", "TSC2"],
        "drugs_targeting": ["PI3K inhibitors (alpelisib)", "AKT inhibitors (capivasertib)", "mTOR inhibitors (everolimus)"]
    },
    "vegf": {
        "name": "VEGF/Angiogenesis",
        "description": "Tumor vascularization, hypoxia response, metastatic potential",
        "aliases": ["VEGF", "angiogenesis", "hypoxia"],
        "genes_in_pathway": ["VEGFA", "VEGFR1", "VEGFR2", "HIF1A", "FLT1", "KDR"],
        "drugs_targeting": ["Bevacizumab (anti-VEGF)", "VEGFR inhibitors (pazopanib, sorafenib)"]
    },
    "her2": {
        "name": "HER2/ERBB2 Signaling",
        "description": "HER2 amplification and signaling, growth factor receptor pathway",
        "aliases": ["HER2", "ERBB2", "HER2/neu"],
        "genes_in_pathway": ["ERBB2", "HER2"],
        "drugs_targeting": ["Trastuzumab", "Pertuzumab", "T-DM1"]
    },
    "tp53": {
        "name": "TP53/Cell Cycle Checkpoint",
        "description": "p53-mediated apoptosis, G1/S checkpoint, DNA damage response",
        "aliases": ["p53", "MDM2", "cell cycle"],
        "genes_in_pathway": ["TP53", "MDM2", "MDM4", "CDKN1A", "CDKN2A", "RB1", "E2F1"],
        "drugs_targeting": ["MDM2 inhibitors (experimental)"]
    }
}


def load_biomarker_results(biomarker_file: Path) -> List[Dict[str, Any]]:
    """Load biomarker analysis results."""
    with open(biomarker_file, 'r') as f:
        data = json.load(f)
    
    # Extract top features
    features = data.get("top_features", [])
    logger.info(f"Loaded {len(features)} features from biomarker analysis")
    return features


def load_cohort_data(cohort_file: Path) -> Dict[str, Any]:
    """Load SAE cohort data."""
    with open(cohort_file, 'r') as f:
        data = json.load(f)
    
    patients = data.get("patients", [])
    logger.info(f"Loaded {len(patients)} patients from cohort")
    return {"patients": patients, "metadata": data.get("provenance", {})}


def find_patients_with_feature(patients: List[Dict], feature_index: int, threshold_percentile: float = 90.0) -> List[str]:
    """Find patients where this feature is highly activated (top threshold_percentile)."""
    feature_values = []
    patient_values = []
    
    for patient in patients:
        # Get feature value for this patient
        # Check aggregated features first
        agg_features = patient.get("aggregated_features", {})
        mean_features = agg_features.get("mean_features", [])
        max_features = agg_features.get("max_features", [])
        
        # Use max feature value
        if feature_index < len(max_features):
            value = max_features[feature_index]
            feature_values.append(value)
            patient_values.append(patient.get("patient_id", "unknown"))
        elif feature_index < len(mean_features):
            value = mean_features[feature_index]
            feature_values.append(value)
            patient_values.append(patient.get("patient_id", "unknown"))
    
    if not feature_values:
        return []
    
    # Find threshold (top threshold_percentile)
    threshold = sorted(feature_values, reverse=True)[int(len(feature_values) * (100 - threshold_percentile) / 100)]
    
    # Return patient IDs above threshold
    high_activation_patients = [
        pid for pid, val in zip(patient_values, feature_values) if val >= threshold
    ]
    
    return high_activation_patients


def extract_genes_from_patients(patients: List[Dict], patient_ids: List[str]) -> Set[str]:
    """Extract unique genes from patients."""
    genes = set()
    
    patient_dict = {p.get("patient_id"): p for p in patients}
    
    for pid in patient_ids:
        patient = patient_dict.get(pid)
        if not patient:
            continue
        
        # Extract genes from variants
        variants = patient.get("variants", [])
        for variant in variants:
            var_data = variant.get("variant", {})
            gene = var_data.get("gene")
            if gene:
                genes.add(gene.upper())
    
    return genes


def infer_pathway_from_genes(genes: Set[str]) -> Dict[str, float]:
    """Infer pathway weights from genes using gene→pathway mapping."""
    pathway_counts = Counter()
    
    for gene in genes:
        pathways = GENE_TO_PATHWAY.get(gene, [])
        for pathway in pathways:
            pathway_counts[pathway] += 1
    
    # Normalize to weights (0-1)
    total = sum(pathway_counts.values())
    if total == 0:
        return {}
    
    pathway_weights = {pathway: count / total for pathway, count in pathway_counts.items()}
    
    # Only return pathways with weight >= 0.1 (at least 10% of genes)
    return {p: w for p, w in pathway_weights.items() if w >= 0.1}


def create_mapping_file(
    biomarker_features: List[Dict],
    cohort_data: Dict,
    output_file: Path,
    top_n: int = 100,
    min_correlation: float = 0.3,
    min_p_value: float = 0.05
) -> None:
    """Create feature→pathway mapping file."""
    
    # Filter top features
    filtered_features = [
        f for f in biomarker_features
        if abs(f.get("pearson_r", 0)) >= min_correlation
        and f.get("p_value", 1.0) <= min_p_value
    ][:top_n]
    
    logger.info(f"Filtering to top {len(filtered_features)} features (correlation >= {min_correlation}, p <= {min_p_value})")
    
    patients = cohort_data["patients"]
    
    # Group features by pathway
    pathway_features = defaultdict(list)
    
    for feature in filtered_features:
        feature_index = feature.get("feature_index")
        if feature_index is None:
            continue
        
        # Find patients with high activation
        high_activation_patients = find_patients_with_feature(patients, feature_index)
        
        if not high_activation_patients:
            logger.warning(f"Feature {feature_index}: No patients with high activation")
            continue
        
        # Extract genes from those patients
        genes = extract_genes_from_patients(patients, high_activation_patients)
        
        if not genes:
            logger.warning(f"Feature {feature_index}: No genes found in high-activation patients")
            continue
        
        # Infer pathway
        pathway_weights = infer_pathway_from_genes(genes)
        
        if not pathway_weights:
            logger.warning(f"Feature {feature_index}: No pathway inferred from genes {genes}")
            continue
        
        # Assign to primary pathway (highest weight)
        primary_pathway = max(pathway_weights.items(), key=lambda x: x[1])[0]
        
        pathway_features[primary_pathway].append({
            "index": feature_index,
            "weight": pathway_weights[primary_pathway],
            "genes": list(genes),
            "patient_count": len(high_activation_patients),
            "correlation": feature.get("pearson_r", 0),
            "p_value": feature.get("p_value", 1.0),
            "cohen_d": feature.get("cohen_d", 0)
        })
        
        logger.info(f"Feature {feature_index} → {primary_pathway} (weight: {pathway_weights[primary_pathway]:.2f}, genes: {len(genes)})")
    
    # Build mapping structure
    mapping = {
        "metadata": {
            "version": "v1",
            "last_updated": datetime.now().isoformat(),
            "source": "Sprint 1 biomarker correlation analysis + gene→pathway inference",
            "model": "Goodfire/Evo-2-Layer-26-Mixed",
            "layer": "blocks-26",
            "total_features": 32768,
            "mapped_features": sum(len(features) for features in pathway_features.values()),
            "status": "preliminary",
            "warnings": [
                "⚠️ Preliminary mapping based on gene→pathway inference",
                "⚠️ Requires validation on known cases (BRCA1→DDR, KRAS→MAPK)",
                "⚠️ Manager approval required before use in production"
            ]
        },
        "pathways": {}
    }
    
    # Populate pathways
    for pathway_name, features in pathway_features.items():
        pathway_meta = PATHWAY_METADATA.get(pathway_name, {})
        
        # Compute aggregate statistics
        correlations = [f["correlation"] for f in features]
        p_values = [f["p_value"] for f in features]
        cohen_ds = [f["cohen_d"] for f in features if f.get("cohen_d")]
        
        mapping["pathways"][pathway_name] = {
            "name": pathway_meta.get("name", pathway_name),
            "description": pathway_meta.get("description", ""),
            "aliases": pathway_meta.get("aliases", []),
            "feature_indices": [f["index"] for f in features],
            "correlation_with_platinum": "positive" if sum(correlations) > 0 else "negative",
            "mean_pearson_r": sum(correlations) / len(correlations) if correlations else None,
            "mean_cohen_d": sum(cohen_ds) / len(cohen_ds) if cohen_ds else None,
            "mean_p_value": sum(p_values) / len(p_values) if p_values else None,
            "confidence": "preliminary",
            "evidence": f"Gene→pathway inference from {len(features)} features",
            "biological_rationale": pathway_meta.get("description", ""),
            "genes_in_pathway": pathway_meta.get("genes_in_pathway", []),
            "drugs_targeting": pathway_meta.get("drugs_targeting", [])
        }
    
    mapping["provenance"] = {
        "creation_method": "gene_to_pathway_inference",
        "curator": "Zo (Lead Commander)",
        "timestamp": datetime.now().isoformat(),
        "sprint": "Sprint 2 Task 1",
        "validation_status": "preliminary",
        "manager_approval": "pending",
        "next_steps": [
            "Validate on known cases (BRCA1→DDR, KRAS→MAPK, HER2→HER2)",
            "Literature review for top features",
            "Manager review and approval"
        ]
    }
    
    # Save mapping file
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(mapping, f, indent=2)
    
    logger.info(f"✅ Mapping file created: {output_file}")
    logger.info(f"   Mapped {mapping['metadata']['mapped_features']} features to {len(mapping['pathways'])} pathways")


def main():
    parser = argparse.ArgumentParser(description="Create SAE feature→pathway mapping")
    parser.add_argument("--biomarkers", type=Path, required=True,
                        help="Biomarker analysis results JSON file")
    parser.add_argument("--cohort", type=Path, required=True,
                        help="SAE cohort file with patient features")
    parser.add_argument("--output", type=Path, required=True,
                        help="Output mapping JSON file")
    parser.add_argument("--top-n", type=int, default=100,
                        help="Top N features to map (default: 100)")
    parser.add_argument("--min-correlation", type=float, default=0.3,
                        help="Minimum correlation threshold (default: 0.3)")
    parser.add_argument("--min-p-value", type=float, default=0.05,
                        help="Maximum p-value threshold (default: 0.05)")
    
    args = parser.parse_args()
    
    logger.info("=" * 80)
    logger.info("🗺️  Creating SAE Feature→Pathway Mapping")
    logger.info("=" * 80)
    
    # Load data
    biomarker_features = load_biomarker_results(args.biomarkers)
    cohort_data = load_cohort_data(args.cohort)
    
    # Create mapping
    create_mapping_file(
        biomarker_features,
        cohort_data,
        args.output,
        top_n=args.top_n,
        min_correlation=args.min_correlation,
        min_p_value=args.min_p_value
    )
    
    logger.info("=" * 80)
    logger.info("✅ Mapping creation complete!")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()

