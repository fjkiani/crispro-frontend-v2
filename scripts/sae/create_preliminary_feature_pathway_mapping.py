#!/usr/bin/env python3
"""
Create Preliminary Feature→Pathway Mapping

Since biomarker analysis found 0 significant features (small sample size),
we create a preliminary mapping based on:
1. Top frequent features across cohort
2. Gene→pathway inference from patient mutations
3. Known pathway relationships

This mapping will be validated and refined when more data is available.
"""

import json
from collections import Counter, defaultdict
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Set, Any

# Gene→Pathway Mappings
GENE_TO_PATHWAY = {
    # DDR pathway
    "BRCA1": ["ddr"], "BRCA2": ["ddr"], "RAD51": ["ddr"], "RAD51C": ["ddr"],
    "RAD51D": ["ddr"], "ATM": ["ddr", "tp53"], "ATR": ["ddr", "tp53"],
    "CHEK2": ["ddr", "tp53"], "PALB2": ["ddr"], "FANCA": ["ddr"],
    "FANCD2": ["ddr"], "PARP1": ["ddr"], "XRCC1": ["ddr"], "MBD4": ["ddr"],
    
    # MAPK pathway
    "KRAS": ["mapk"], "NRAS": ["mapk"], "BRAF": ["mapk"], "MAP2K1": ["mapk"],
    "MAP2K2": ["mapk"], "MAPK1": ["mapk"], "MAPK3": ["mapk"], "RAF1": ["mapk"],
    
    # PI3K pathway
    "PIK3CA": ["pi3k"], "PIK3CB": ["pi3k"], "AKT1": ["pi3k"], "AKT2": ["pi3k"],
    "MTOR": ["pi3k"], "PTEN": ["pi3k", "tp53"], "TSC1": ["pi3k"], "TSC2": ["pi3k"],
    
    # VEGF pathway
    "VEGFA": ["vegf"], "VEGFR1": ["vegf"], "VEGFR2": ["vegf"], "HIF1A": ["vegf"],
    "FLT1": ["vegf"], "KDR": ["vegf"],
    
    # HER2 pathway
    "ERBB2": ["her2"], "HER2": ["her2"], "ERBB3": ["her2"],
    
    # TP53 pathway
    "TP53": ["tp53"], "MDM2": ["tp53"], "MDM4": ["tp53"], "CDKN1A": ["tp53"],
    "CDKN2A": ["tp53"], "RB1": ["tp53"], "E2F1": ["tp53"],
}

# Pathway metadata
PATHWAY_METADATA = {
    "ddr": {
        "name": "DNA Damage Response & Homologous Recombination",
        "description": "BRCA1/2-mediated homologous recombination repair, DNA damage signaling",
        "aliases": ["DNA repair", "HR", "HRD", "BRCA"],
        "genes_in_pathway": ["BRCA1", "BRCA2", "RAD51", "ATM", "ATR", "CHEK2", "PALB2", "MBD4"],
        "drugs_targeting": ["PARP inhibitors", "Platinum chemotherapy"]
    },
    "mapk": {
        "name": "MAPK/RAS/RAF Signaling",
        "description": "Oncogenic RAS→RAF→MEK→ERK cascade",
        "aliases": ["RAS pathway", "MAPK cascade"],
        "genes_in_pathway": ["KRAS", "NRAS", "BRAF", "MAP2K1", "MAP2K2"],
        "drugs_targeting": ["MEK inhibitors", "RAF inhibitors"]
    },
    "pi3k": {
        "name": "PI3K/AKT/mTOR Pathway",
        "description": "Growth factor signaling, apoptosis resistance",
        "aliases": ["PI3K", "AKT", "mTOR"],
        "genes_in_pathway": ["PIK3CA", "PTEN", "AKT1", "AKT2", "MTOR"],
        "drugs_targeting": ["PI3K inhibitors", "AKT inhibitors"]
    },
    "vegf": {
        "name": "VEGF/Angiogenesis",
        "description": "Tumor vascularization, hypoxia response",
        "aliases": ["VEGF", "angiogenesis"],
        "genes_in_pathway": ["VEGFA", "VEGFR1", "VEGFR2"],
        "drugs_targeting": ["Bevacizumab", "VEGFR inhibitors"]
    },
    "her2": {
        "name": "HER2/ERBB2 Signaling",
        "description": "HER2 amplification, growth signaling",
        "aliases": ["HER2", "ERBB2"],
        "genes_in_pathway": ["ERBB2", "ERBB3"],
        "drugs_targeting": ["Trastuzumab", "HER2 inhibitors"]
    },
    "tp53": {
        "name": "TP53/Cell Cycle Checkpoint",
        "description": "p53-mediated apoptosis, DNA damage response",
        "aliases": ["p53", "cell cycle"],
        "genes_in_pathway": ["TP53", "MDM2", "CDKN1A"],
        "drugs_targeting": ["MDM2 inhibitors (experimental)"]
    },
    "io": {
        "name": "Immune/Inflammatory Response",
        "description": "T-cell infiltration, immune checkpoint signaling",
        "aliases": ["Immune", "T-cell", "PD-L1"],
        "genes_in_pathway": ["CD8A", "CD274", "PDCD1"],
        "drugs_targeting": ["Checkpoint inhibitors"]
    }
}


def create_preliminary_mapping(cohort_file: Path, output_file: Path, top_n: int = 100):
    """Create preliminary feature→pathway mapping from cohort data."""
    
    # Load cohort
    with open(cohort_file) as f:
        data = json.load(f)
    
    patients = data.get("patients", [])
    print(f"Loaded {len(patients)} patients")
    
    # Collect feature frequencies and gene associations
    feature_counts = Counter()
    feature_to_genes = defaultdict(set)
    gene_to_features = defaultdict(set)
    
    for patient in patients:
        variants = patient.get("variants", [])
        for variant in variants:
            var_data = variant.get("variant", {})
            gene = var_data.get("gene", "").upper()
            
            top_features = variant.get("top_features", [])
            for feat in top_features:
                idx = feat.get("index")
                if idx is not None:
                    feature_counts[idx] += 1
                    if gene:
                        feature_to_genes[idx].add(gene)
                        gene_to_features[gene].add(idx)
    
    # Get top N most frequent features
    top_features = [idx for idx, count in feature_counts.most_common(top_n)]
    print(f"Top {len(top_features)} features identified")
    
    # Map features to pathways via genes
    pathway_features = defaultdict(list)
    
    for feat_idx in top_features:
        genes = feature_to_genes.get(feat_idx, set())
        if not genes:
            continue
        
        # Infer pathways from genes
        pathway_counts = Counter()
        for gene in genes:
            pathways = GENE_TO_PATHWAY.get(gene, [])
            for pathway in pathways:
                pathway_counts[pathway] += 1
        
        if not pathway_counts:
            continue
        
        # Normalize to weights
        total = sum(pathway_counts.values())
        pathway_weights = {p: count / total for p, count in pathway_counts.items()}
        
        # Assign to primary pathway
        primary_pathway = max(pathway_weights.items(), key=lambda x: x[1])[0]
        
        pathway_features[primary_pathway].append({
            "index": feat_idx,
            "weight": pathway_weights[primary_pathway],
            "genes": list(genes),
            "frequency": feature_counts[feat_idx],
            "pathway_weights": pathway_weights
        })
    
    # Build mapping structure
    mapping = {
        "metadata": {
            "version": "v1_preliminary",
            "last_updated": datetime.now().isoformat(),
            "source": "Preliminary mapping from cohort feature frequencies + gene→pathway inference",
            "model": "Goodfire/Evo-2-Layer-26-Mixed",
            "layer": "blocks-26",
            "total_features": 32768,
            "mapped_features": sum(len(features) for features in pathway_features.values()),
            "status": "preliminary",
            "warnings": [
                "⚠️ PRELIMINARY MAPPING - Not validated against biomarker analysis",
                "⚠️ Based on feature frequency and gene→pathway inference",
                "⚠️ Requires validation on known cases (BRCA1→DDR, KRAS→MAPK)",
                "⚠️ Manager approval required before use in production"
            ]
        },
        "pathways": {}
    }
    
    # Populate pathways
    for pathway_name, features in pathway_features.items():
        pathway_meta = PATHWAY_METADATA.get(pathway_name, {})
        
        mapping["pathways"][pathway_name] = {
            "name": pathway_meta.get("name", pathway_name),
            "description": pathway_meta.get("description", ""),
            "aliases": pathway_meta.get("aliases", []),
            "feature_indices": [f["index"] for f in features],
            "feature_count": len(features),
            "mean_frequency": sum(f["frequency"] for f in features) / len(features) if features else 0,
            "confidence": "preliminary",
            "evidence": f"Gene→pathway inference from {len(features)} features",
            "biological_rationale": pathway_meta.get("description", ""),
            "genes_in_pathway": pathway_meta.get("genes_in_pathway", []),
            "drugs_targeting": pathway_meta.get("drugs_targeting", []),
            "features": features[:20]  # Store top 20 for reference
        }
    
    mapping["provenance"] = {
        "creation_method": "preliminary_gene_to_pathway_inference",
        "curator": "Automated from cohort data",
        "timestamp": datetime.now().isoformat(),
        "cohort_size": len(patients),
        "top_features_analyzed": top_n,
        "validation_status": "preliminary",
        "manager_approval": "pending",
        "next_steps": [
            "Validate on known cases (BRCA1→DDR, KRAS→MAPK, HER2→HER2)",
            "Re-run biomarker analysis with larger cohort",
            "Literature review for top features",
            "Manager review and approval"
        ]
    }
    
    # Save mapping
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(mapping, f, indent=2)
    
    print(f"✅ Preliminary mapping created: {output_file}")
    print(f"   Mapped {mapping['metadata']['mapped_features']} features to {len(mapping['pathways'])} pathways")
    
    # Print summary
    print("\nPathway Summary:")
    for pathway, data in mapping["pathways"].items():
        print(f"  {pathway}: {data['feature_count']} features")


if __name__ == "__main__":
    import sys
    
    cohort_file = Path("data/validation/sae_cohort/sae_features_tcga_ov_platinum.json")
    output_file = Path("oncology-coPilot/oncology-backend-minimal/api/resources/sae_feature_mapping.json")
    
    create_preliminary_mapping(cohort_file, output_file, top_n=100)

