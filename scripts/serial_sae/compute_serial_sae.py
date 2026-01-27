#!/usr/bin/env python3
"""
Compute Serial SAE Pathway Scores
Compute pathway scores at baseline and progression timepoints
Calculate pathway kinetics (ΔSAE)
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional
import sys

# Add backend to path for SAE service
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "oncology-coPilot" / "oncology-backend-minimal"))

from api.services.sae_feature_service import SAEFeatureService
from api.services.pathway_to_mechanism_vector import convert_pathway_scores_to_mechanism_vector

DATA_DIR = Path("data/serial_sae")
RESULTS_DIR = DATA_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

def compute_pathway_scores_from_mutations(mutations: List[Dict], pathway_genes: Dict[str, List[str]]) -> Dict[str, float]:
    """
    Compute pathway scores from mutation list.
    
    Args:
        mutations: List of mutation dicts with 'gene' field
        pathway_genes: Dict mapping pathway names to gene lists
    
    Returns:
        Dict of pathway_name -> score (0.0-1.0)
    """
    pathway_scores = {}
    
    # Get unique mutated genes
    mutated_genes = set(m.get("gene", "").upper() for m in mutations if m.get("gene"))
    
    for pathway_name, genes in pathway_genes.items():
        # Count how many pathway genes are mutated
        pathway_mutated = len([g for g in genes if g.upper() in mutated_genes])
        # Score = proportion of pathway genes mutated (capped at 1.0)
        pathway_scores[pathway_name] = min(pathway_mutated / len(genes) if len(genes) > 0 else 0.0, 1.0)
    
    return pathway_scores

def compute_pathway_scores_from_expression(expression_data: Dict[str, float], pathway_genes: Dict[str, List[str]]) -> Dict[str, float]:
    """
    Compute pathway scores from expression data.
    
    Args:
        expression_data: Dict of gene -> expression value
        pathway_genes: Dict mapping pathway names to gene lists
    
    Returns:
        Dict of pathway_name -> score (0.0-1.0)
    """
    pathway_scores = {}
    
    for pathway_name, genes in pathway_genes.items():
        # Get expression values for pathway genes
        pathway_expr = [expression_data.get(g.upper(), 0.0) for g in genes if g.upper() in expression_data]
        
        if len(pathway_expr) > 0:
            # Score = mean expression normalized (assuming Z-scores or log-normalized)
            # For Z-scores: positive = high expression, negative = low
            # Convert to 0-1 scale: sigmoid or percentile
            mean_expr = np.mean(pathway_expr)
            # Simple normalization: sigmoid to 0-1
            pathway_scores[pathway_name] = 1.0 / (1.0 + np.exp(-mean_expr))
        else:
            pathway_scores[pathway_name] = 0.0
    
    return pathway_scores

def compute_mechanism_vector(pathway_scores: Dict[str, float]) -> List[float]:
    """Convert pathway scores to 7D mechanism vector."""
    # Map pathway scores to mechanism vector
    # [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
    
    mechanism_vector = [
        pathway_scores.get("ddr", 0.0) + pathway_scores.get("tp53", 0.0) * 0.5,  # DDR + TP53 contribution
        pathway_scores.get("mapk", 0.0),  # MAPK
        pathway_scores.get("pi3k", 0.0),  # PI3K
        pathway_scores.get("vegf", 0.0),  # VEGF
        pathway_scores.get("her2", 0.0),  # HER2
        0.0,  # IO (needs TMB/MSI - set separately)
        pathway_scores.get("efflux", 0.0)  # Efflux
    ]
    
    return mechanism_vector

def compute_pathway_kinetics(baseline_scores: Dict[str, float], progression_scores: Dict[str, float], time_interval_months: float) -> Dict[str, float]:
    """
    Compute pathway kinetics (rate of change).
    
    Args:
        baseline_scores: Pathway scores at baseline
        progression_scores: Pathway scores at progression
        time_interval_months: Time between samples in months
    
    Returns:
        Dict of pathway_name -> delta (change per month)
    """
    kinetics = {}
    
    for pathway in baseline_scores.keys():
        if pathway in progression_scores:
            delta = (progression_scores[pathway] - baseline_scores[pathway]) / time_interval_months if time_interval_months > 0 else 0.0
            kinetics[f"delta_{pathway}"] = delta
            kinetics[f"absolute_delta_{pathway}"] = progression_scores[pathway] - baseline_scores[pathway]
    
    return kinetics

def main():
    print("="*80)
    print("Serial SAE Pathway Score Computation")
    print("="*80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Pathway gene definitions
    PATHWAY_GENES = {
        "ddr": ["BRCA1", "BRCA2", "ATM", "CHEK2", "RAD51", "PALB2", "RAD51C", "RAD51D", "BARD1", "NBN", "TP53"],
        "mapk": ["KRAS", "NRAS", "BRAF", "MEK1", "MEK2", "ERK1", "ERK2", "RAF1"],
        "pi3k": ["PIK3CA", "PIK3CB", "PIK3CD", "PTEN", "AKT1", "AKT2", "AKT3", "MTOR"],
        "vegf": ["VEGFA", "VEGFR1", "VEGFR2", "KDR", "FLT1"],
        "her2": ["ERBB2", "HER2"],
        "efflux": ["ABCB1", "MDR1", "ABCC1", "ABCG2"]
    }
    
    # Load paired patients
    paired_file = DATA_DIR / "cbioportal_paired" / "paired_patients.json"
    if not paired_file.exists():
        print(f"❌ Error: {paired_file} not found. Run download_cbioportal_paired.py first.")
        return
    
    with open(paired_file) as f:
        paired_data = json.load(f)
    
    paired_patients = paired_data.get("paired_patients", {})
    print(f"📊 Processing {len(paired_patients)} paired patients")
    print()
    
    results = []
    
    for patient_id, patient_data in list(paired_patients.items())[:10]:  # Process first 10 for now
        print(f"🔬 Processing patient: {patient_id}")
        
        primary_samples = patient_data.get("primary_samples", [])
        recurrent_samples = patient_data.get("recurrent_samples", [])
        
        if not primary_samples or not recurrent_samples:
            print(f"  ⚠️ Skipping - missing primary or recurrent samples")
            continue
        
        # For now, use first primary and first recurrent
        primary_sample = primary_samples[0]["sample_id"]
        recurrent_sample = recurrent_samples[0]["sample_id"]
        
        # TODO: Download actual mutation/expression data for these samples
        # For now, create placeholder structure
        
        patient_result = {
            "patient_id": patient_id,
            "primary_sample_id": primary_sample,
            "recurrent_sample_id": recurrent_sample,
            "baseline_sae": {
                "pathway_scores": {},  # To be filled with actual data
                "mechanism_vector": [0.0] * 7
            },
            "progression_sae": {
                "pathway_scores": {},  # To be filled with actual data
                "mechanism_vector": [0.0] * 7
            },
            "pathway_kinetics": {},
            "time_interval_months": None,  # To be filled from clinical data
            "status": "placeholder - needs actual data download"
        }
        
        results.append(patient_result)
        print(f"  ✅ Added to results (placeholder)")
    
    # Save results
    output_file = RESULTS_DIR / "serial_sae_computation.json"
    with open(output_file, "w") as f:
        json.dump({
            "computation_date": datetime.now().isoformat(),
            "n_patients_processed": len(results),
            "results": results,
            "note": "Placeholder - needs actual mutation/expression data download"
        }, f, indent=2)
    
    print(f"\n💾 Results saved to: {output_file}")
    print(f"\n⚠️ NOTE: This is a placeholder. Next steps:")
    print(f"  1. Download mutations for paired samples")
    print(f"  2. Download expression data for paired samples")
    print(f"  3. Re-run this script with actual data")

if __name__ == "__main__":
    main()
