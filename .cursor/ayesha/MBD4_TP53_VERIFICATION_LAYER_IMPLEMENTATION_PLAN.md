# MBD4+TP53 Verification Layer Implementation Plan

**Date**: January 21, 2025  
**Purpose**: Implement automated verification layer for MBD4+TP53 analysis answers  
**Scope**: Focus on what we can control now (deterministic checks, formula validation, consistency checks)  
**Reference**: `.cursor/ayesha/MBD4_TP53_ANSWERS_CRITICAL_ANALYSIS.md`

---

## 🎯 OBJECTIVE

Build an automated verification layer that validates analysis answers against:
1. **Deterministic Sources** (ClinVar, COSMIC, KEGG, Reactome, FDA, NCCN) - 90-100% confidence
2. **Formula Correctness** (DNA repair capacity, mechanism vectors) - 75-90% confidence
3. **Consistency Checks** (pathway mapping, variant classification) - 85-90% confidence
4. **Biological Plausibility** (expected ranges, known biology) - 80-90% confidence

**What We CANNOT Control (Out of Scope)**:
- Clinical validation (requires patient outcomes)
- Prospective validation (requires time)
- LLM extraction quality (variable, manual review)

---

## 📋 PHASE 1: DETERMINISTIC VERIFICATION (HIGH PRIORITY)

### **Task 1.1: Variant Classification Verification**

**File**: `scripts/sae/verify_variant_classification.py` (NEW)

**Purpose**: Verify variant impact predictions against ClinVar and COSMIC

**Verification Methods** (from `MBD4_TP53_ANSWERS_CRITICAL_ANALYSIS.md` lines 31-35):

1. **ClinVar Check**:
   ```python
   def verify_clinvar_classification(variant: Dict) -> Dict:
       """
       Verify variant classification against ClinVar.
       
       Expected for MBD4 frameshift:
       - Classification: "Pathogenic"
       - Review status: "reviewed by expert panel" or "criteria provided"
       - Clinical significance: "Pathogenic" or "Likely pathogenic"
       """
       # Call /api/evidence/clinvar endpoint
       # Compare to expected classification
       # Return: {verified: bool, expected: str, actual: str, confidence: float}
   ```

2. **COSMIC Check**:
   ```python
   def verify_cosmic_hotspot(variant: Dict) -> Dict:
       """
       Verify hotspot mutation against COSMIC database.
       
       Expected for TP53 R175H:
       - In hotspot database (high frequency)
       - Mutation frequency > 0.01 (1% of TP53 mutations)
       """
       # Check against COSMIC hotspot database (local or API)
       # Compare to expected hotspot status
       # Return: {verified: bool, is_hotspot: bool, frequency: float}
   ```

3. **Evo2 Validation**:
   ```python
   def verify_evo2_scores(variant: Dict, analysis_result: Dict) -> Dict:
       """
       Verify Evo2 delta scores are highly negative (disruptive).
       
       Expected:
       - MBD4 frameshift: delta < -5.0 (highly disruptive)
       - TP53 R175H: delta < -3.0 (moderately disruptive)
       """
       # Extract delta scores from analysis_result
       # Compare to expected ranges
       # Return: {verified: bool, delta_score: float, expected_range: tuple}
   ```

**Integration**: Call after `run_mbd4_tp53_analysis.py` completes

**Output**: Verification report with pass/fail for each variant

---

### **Task 1.2: Pathway Mapping Verification**

**File**: `scripts/sae/verify_pathway_mapping.py` (NEW)

**Verification Methods** (from `MBD4_TP53_ANSWERS_CRITICAL_ANALYSIS.md` lines 83-89):

1. **KEGG Check**:
   ```python
   def verify_kegg_pathway(gene: str, expected_pathway: str) -> Dict:
       """
       Verify gene→pathway mapping against KEGG database.
       
       Expected for MBD4:
       - Pathway: "Base excision repair" (BER) → DDR pathway
       
       Expected for TP53:
       - Pathway: "p53 signaling pathway" → DDR pathway
       """
       # Query KEGG API or local database
       # Compare to expected pathway
       # Return: {verified: bool, expected: str, actual: List[str], confidence: float}
   ```

2. **Reactome Check**:
   ```python
   def verify_reactome_pathway(gene: str, expected_pathway: str) -> Dict:
       """
       Verify gene→pathway mapping against Reactome database.
       
       Expected for MBD4:
       - Pathway: "BER pathway" → DDR pathway
       
       Expected for TP53:
       - Pathway: "TP53 checkpoint pathway" → DDR pathway
       """
       # Query Reactome API or local database
       # Compare to expected pathway
       # Return: {verified: bool, expected: str, actual: List[str], confidence: float}
   ```

3. **Formula Check**:
   ```python
   def verify_dna_repair_formula(pathway_scores: Dict, computed_dna_repair: float) -> Dict:
       """
       Verify DNA repair capacity formula correctness.
       
       Formula: (0.6 × DDR) + (0.2 × HRR) + (0.2 × exon)
       
       Expected for MBD4+TP53:
       - DDR score: 0.70-0.90
       - DNA repair capacity: 0.75-0.90
       """
       # Recompute formula manually
       # Compare to computed value
       # Return: {verified: bool, expected: float, actual: float, difference: float}
   ```

4. **TCGA Validation**:
   ```python
   def verify_tcga_pathway_weights(pathway_scores: Dict, disease: str) -> Dict:
       """
       Verify pathway weights come from real TCGA mutation frequencies.
       
       Expected:
       - Pathway weights from TCGA-OV (ovarian cancer)
       - Weights sum to reasonable range (not all 0.0 or 1.0)
       """
       # Check pathway weights against TCGA data
       # Validate weights are from real mutation frequencies
       # Return: {verified: bool, source: str, weights: Dict}
   ```

**Integration**: Call after pathway scores computed in analysis

**Output**: Verification report with pass/fail for each pathway mapping

---

### **Task 1.3: Functional Annotation Verification**

**File**: `scripts/sae/verify_functional_annotation.py` (NEW)

**Verification Methods** (from `MBD4_TP53_ANSWERS_CRITICAL_ANALYSIS.md` lines 57-60):

1. **UniProt Check**:
   ```python
   def verify_uniprot_function(gene: str, expected_function: str) -> Dict:
       """
       Verify protein function against UniProt database.
       
       Expected for MBD4:
       - Function: "DNA glycosylase" (BER pathway)
       
       Expected for TP53:
       - Function: "Tumor suppressor" (checkpoint, apoptosis)
       """
       # Query UniProt API or local database
       # Compare to expected function
       # Return: {verified: bool, expected: str, actual: str, confidence: float}
   ```

2. **Insights Bundle Validation**:
   ```python
   def verify_insights_bundle(insights: Dict, expected_ranges: Dict) -> Dict:
       """
       Verify insights bundle scores are within expected ranges.
       
       Expected for MBD4 frameshift:
       - Functionality: 0.0-0.3 (loss-of-function)
       - Essentiality: 0.7+ (high)
       
       Expected for TP53 R175H:
       - Functionality: 0.0-0.3 OR ≥0.80 (hotspot boost)
       - Essentiality: 0.7+ (high)
       """
       # Compare insights scores to expected ranges
       # Return: {verified: bool, insights: Dict, expected_ranges: Dict}
   ```

**Integration**: Call after insights bundle computed

**Output**: Verification report with pass/fail for each functional annotation

---

### **Task 1.4: Eligibility & IO Verification**

**File**: `scripts/sae/verify_eligibility_io.py` (NEW)

**Verification Methods** (from `MBD4_TP53_ANSWERS_CRITICAL_ANALYSIS.md` lines 198-200):

1. **FDA Labels Check**:
   ```python
   def verify_fda_io_eligibility(tmb: float, msi_status: str) -> Dict:
       """
       Verify IO eligibility against FDA labels.
       
       Expected:
       - TMB ≥20 → IO eligible (FDA label)
       - MSI-H → IO eligible (FDA label)
       
       For MBD4+TP53:
       - TMB = 25.0 → IO eligible (TRUE)
       """
       # Check FDA labels (local database or API)
       # Compare to computed eligibility
       # Return: {verified: bool, expected: bool, actual: bool, criteria: str}
   ```

2. **NCCN Guidelines Check**:
   ```python
   def verify_nccn_drug_recommendation(drug: str, disease: str) -> Dict:
       """
       Verify drug recommendations against NCCN guidelines.
       
       Expected for ovarian cancer:
       - Carboplatin + Paclitaxel = First-line (95-100% confidence)
       - PARP inhibitors = HRD+ maintenance (FDA label)
       """
       # Check NCCN guidelines (local database)
       # Compare to computed recommendation
       # Return: {verified: bool, expected: str, actual: str, confidence: float}
   ```

**Integration**: Call after IO eligibility and drug recommendations computed

**Output**: Verification report with pass/fail for eligibility checks

---

## 📋 PHASE 2: FORMULA & CONSISTENCY VERIFICATION (HIGH PRIORITY)

### **Task 2.1: DNA Repair Capacity Formula Verification**

**File**: `scripts/sae/verify_dna_repair_formula.py` (NEW)

**Purpose**: Verify DNA repair capacity formula correctness and expected ranges

**Verification Methods**:

1. **Formula Correctness**:
   ```python
   def verify_formula_correctness(pathway_scores: Dict, computed_dna_repair: float) -> Dict:
       """
       Verify formula: (0.6 × DDR) + (0.2 × HRR) + (0.2 × exon)
       
       Expected for MBD4+TP53:
       - DDR: 0.70-0.90
       - HRR essentiality: 0.7+ (if HRR genes present)
       - Exon disruption: 0.7+ (if essentiality > 0.65)
       - DNA repair capacity: 0.75-0.90
       """
       # Recompute formula manually
       ddr = pathway_scores.get("ddr", 0.0)
       hrr_essentiality = insights.get("essentiality_hrr", 0.0)
       exon_disruption = insights.get("exon_disruption_score", 0.0)
       
       expected = (0.6 * ddr) + (0.2 * hrr_essentiality) + (0.2 * exon_disruption)
       
       # Compare to computed value
       difference = abs(computed_dna_repair - expected)
       
       return {
           "verified": difference < 0.01,  # Within 1% tolerance
           "expected": expected,
           "actual": computed_dna_repair,
           "difference": difference,
           "formula": "0.6×DDR + 0.2×HRR + 0.2×exon"
       }
   ```

2. **Expected Range Check**:
   ```python
   def verify_expected_range(computed_dna_repair: float, expected_range: tuple) -> Dict:
       """
       Verify DNA repair capacity is within expected range.
       
       Expected for MBD4+TP53 (DDR-high):
       - Range: 0.75-0.90 (very high)
       """
       min_val, max_val = expected_range
       in_range = min_val <= computed_dna_repair <= max_val
       
       return {
           "verified": in_range,
           "expected_range": expected_range,
           "actual": computed_dna_repair,
           "within_range": in_range
       }
   ```

**Integration**: Call after DNA repair capacity computed

**Output**: Verification report with formula correctness and range checks

---

### **Task 2.2: Mechanism Vector Verification**

**File**: `scripts/sae/verify_mechanism_vector.py` (NEW)

**Purpose**: Verify mechanism vector structure and pathway mapping

**Verification Methods**:

1. **Vector Structure Check**:
   ```python
   def verify_vector_structure(mechanism_vector: List[float], expected_dim: int) -> Dict:
       """
       Verify mechanism vector has correct structure.
       
       Expected:
       - 7D vector: [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
       - All values in range [0.0, 1.0]
       - Sum of values ≤ 7.0 (not normalized to 1.0)
       """
       # Check dimension
       dim_correct = len(mechanism_vector) == expected_dim
       
       # Check value ranges
       all_in_range = all(0.0 <= v <= 1.0 for v in mechanism_vector)
       
       # Check sum (should be ≤ dim, not normalized)
       sum_valid = sum(mechanism_vector) <= expected_dim
       
       return {
           "verified": dim_correct and all_in_range and sum_valid,
           "dimension": len(mechanism_vector),
           "expected_dim": expected_dim,
           "all_in_range": all_in_range,
           "sum_valid": sum_valid
       }
   ```

2. **Pathway Mapping Check**:
   ```python
   def verify_pathway_mapping(pathway_scores: Dict, mechanism_vector: List[float]) -> Dict:
       """
       Verify pathway scores correctly map to mechanism vector indices.
       
       Expected for MBD4+TP53:
       - DDR pathway → mechanism_vector[0] = 0.70-0.90
       - MAPK pathway → mechanism_vector[1] = 0.10-0.30 (low)
       - IO → mechanism_vector[5] = 1.0 (TMB ≥20)
       """
       # Check DDR mapping
       ddr_score = pathway_scores.get("ddr", 0.0)
       ddr_idx = 0  # DDR is index 0 in 7D vector
       ddr_mapped = abs(mechanism_vector[ddr_idx] - ddr_score) < 0.10
       
       # Check TP53 → DDR contribution (50% rule)
       tp53_score = pathway_scores.get("tp53", 0.0)
       expected_ddr_with_tp53 = ddr_score + (tp53_score * 0.5)
       tp53_contribution_correct = abs(mechanism_vector[ddr_idx] - expected_ddr_with_tp53) < 0.10
       
       return {
           "verified": ddr_mapped and tp53_contribution_correct,
           "ddr_mapping": ddr_mapped,
           "tp53_contribution": tp53_contribution_correct,
           "expected_ddr": expected_ddr_with_tp53,
           "actual_ddr": mechanism_vector[ddr_idx]
       }
   ```

**Integration**: Call after mechanism vector computed

**Output**: Verification report with structure and mapping checks

---

### **Task 2.3: Consistency Checks**

**File**: `scripts/sae/verify_consistency.py` (NEW)

**Purpose**: Verify consistency across different analysis components

**Verification Methods**:

1. **Pathway Score Consistency**:
   ```python
   def verify_pathway_consistency(efficacy_response: Dict, sae_features: Dict) -> Dict:
       """
       Verify pathway scores are consistent between efficacy and SAE features.
       
       Expected:
       - Efficacy pathway_scores match SAE pathway_burden_* fields
       - Mechanism vector derived from same pathway scores
       """
       efficacy_pathways = efficacy_response.get("provenance", {}).get("confidence_breakdown", {}).get("pathway_disruption", {})
       sae_pathways = {
           "ddr": sae_features.get("pathway_burden_ddr", 0.0),
           "mapk": sae_features.get("pathway_burden_mapk", 0.0),
           "pi3k": sae_features.get("pathway_burden_pi3k", 0.0),
           "vegf": sae_features.get("pathway_burden_vegf", 0.0),
           "her2": sae_features.get("pathway_burden_her2", 0.0)
       }
       
       # Compare pathway scores
       differences = {}
       for pathway, sae_score in sae_pathways.items():
           eff_score = efficacy_pathways.get(pathway, 0.0)
           differences[pathway] = abs(eff_score - sae_score)
       
       all_consistent = all(diff < 0.10 for diff in differences.values())
       
       return {
           "verified": all_consistent,
           "differences": differences,
           "efficacy_pathways": efficacy_pathways,
           "sae_pathways": sae_pathways
       }
   ```

2. **Variant Annotation Consistency**:
   ```python
   def verify_variant_consistency(mutations: List[Dict], analysis_result: Dict) -> Dict:
       """
       Verify variant annotations are consistent across analysis.
       
       Expected:
       - Input mutations match output variant annotations
       - HGVS notation consistent
       - Gene names consistent
       """
       input_genes = {m.get("gene") for m in mutations}
       output_genes = {v.get("gene") for v in analysis_result.get("variants", [])}
       
       genes_match = input_genes == output_genes
       
       return {
           "verified": genes_match,
           "input_genes": input_genes,
           "output_genes": output_genes,
           "match": genes_match
       }
   ```

**Integration**: Call after all analysis components complete

**Output**: Verification report with consistency checks

---

## 📋 PHASE 3: BIOLOGICAL PLAUSIBILITY VERIFICATION (MEDIUM PRIORITY)

### **Task 3.1: Expected Range Verification**

**File**: `scripts/sae/verify_biological_plausibility.py` (NEW)

**Purpose**: Verify analysis results are biologically plausible

**Verification Methods**:

1. **Pathway Score Ranges**:
   ```python
   def verify_pathway_ranges(pathway_scores: Dict, expected_ranges: Dict) -> Dict:
       """
       Verify pathway scores are within biologically plausible ranges.
       
       Expected for MBD4+TP53:
       - DDR: 0.70-0.90 (high, as expected)
       - MAPK: 0.10-0.30 (low, as expected - no MAPK activation)
       - PI3K: 0.10-0.30 (low, as expected)
       """
       results = {}
       for pathway, expected_range in expected_ranges.items():
           actual_score = pathway_scores.get(pathway, 0.0)
           min_val, max_val = expected_range
           in_range = min_val <= actual_score <= max_val
           results[pathway] = {
               "verified": in_range,
               "expected_range": expected_range,
               "actual": actual_score
           }
       
       all_plausible = all(r["verified"] for r in results.values())
       
       return {
           "verified": all_plausible,
           "pathway_checks": results
       }
   ```

2. **Drug Efficacy Ranges**:
   ```python
   def verify_drug_efficacy_ranges(drug_rankings: List[Dict], expected_ranges: Dict) -> Dict:
       """
       Verify drug efficacy scores are within expected ranges.
       
       Expected for MBD4+TP53:
       - PARP inhibitors: 0.70-0.90 (high, if HRD ≥42)
       - Platinum: 0.60-0.80 (high, for DDR-high tumors)
       - MEK inhibitors: 0.20-0.40 (low, no MAPK activation)
       """
       results = {}
       for drug in drug_rankings:
           drug_name = drug.get("name", "")
           drug_class = drug.get("drug_class", "")
           efficacy_score = drug.get("efficacy_score", 0.0)
           
           # Check against expected ranges
           expected_range = expected_ranges.get(drug_class, (0.0, 1.0))
           min_val, max_val = expected_range
           in_range = min_val <= efficacy_score <= max_val
           
           results[drug_name] = {
               "verified": in_range,
               "expected_range": expected_range,
               "actual": efficacy_score,
               "drug_class": drug_class
           }
       
       all_plausible = all(r["verified"] for r in results.values())
       
       return {
           "verified": all_plausible,
           "drug_checks": results
       }
   ```

**Integration**: Call after analysis complete

**Output**: Verification report with biological plausibility checks

---

## 📋 PHASE 4: INTEGRATION & AUTOMATION (HIGH PRIORITY)

### **Task 4.1: Create Unified Verification Script**

**File**: `scripts/sae/verify_mbd4_analysis.py` (NEW)

**Purpose**: Unified script that runs all verification checks

**Structure**:

```python
#!/usr/bin/env python3
"""
Unified verification script for MBD4+TP53 analysis.
Runs all verification checks and generates comprehensive report.
"""

import json
import sys
from pathlib import Path

# Import all verification modules
from verify_variant_classification import (
    verify_clinvar_classification,
    verify_cosmic_hotspot,
    verify_evo2_scores
)
from verify_pathway_mapping import (
    verify_kegg_pathway,
    verify_reactome_pathway,
    verify_dna_repair_formula,
    verify_tcga_pathway_weights
)
from verify_functional_annotation import (
    verify_uniprot_function,
    verify_insights_bundle
)
from verify_eligibility_io import (
    verify_fda_io_eligibility,
    verify_nccn_drug_recommendation
)
from verify_dna_repair_formula import (
    verify_formula_correctness,
    verify_expected_range
)
from verify_mechanism_vector import (
    verify_vector_structure,
    verify_pathway_mapping
)
from verify_consistency import (
    verify_pathway_consistency,
    verify_variant_consistency
)
from verify_biological_plausibility import (
    verify_pathway_ranges,
    verify_drug_efficacy_ranges
)

def run_all_verifications(analysis_result_path: str) -> Dict:
    """
    Run all verification checks on analysis results.
    
    Args:
        analysis_result_path: Path to analysis JSON output
        
    Returns:
        Comprehensive verification report
    """
    # Load analysis results
    with open(analysis_result_path, 'r') as f:
        analysis_result = json.load(f)
    
    verification_report = {
        "timestamp": datetime.now().isoformat(),
        "analysis_file": analysis_result_path,
        "checks": {}
    }
    
    # Phase 1: Deterministic Verification
    verification_report["checks"]["variant_classification"] = {
        "clinvar": verify_clinvar_classification(analysis_result["variants"][0]),
        "cosmic": verify_cosmic_hotspot(analysis_result["variants"][1]),
        "evo2": verify_evo2_scores(analysis_result["variants"], analysis_result)
    }
    
    verification_report["checks"]["pathway_mapping"] = {
        "kegg": verify_kegg_pathway("MBD4", "DDR"),
        "reactome": verify_reactome_pathway("TP53", "DDR"),
        "formula": verify_dna_repair_formula(
            analysis_result["pathway_scores"],
            analysis_result["sae_features"]["dna_repair_capacity"]
        ),
        "tcga": verify_tcga_pathway_weights(analysis_result["pathway_scores"], "ovarian_cancer")
    }
    
    verification_report["checks"]["functional_annotation"] = {
        "uniprot": verify_uniprot_function("MBD4", "DNA glycosylase"),
        "insights": verify_insights_bundle(
            analysis_result["insights"],
            {"functionality": (0.0, 0.3), "essentiality": (0.7, 1.0)}
        )
    }
    
    verification_report["checks"]["eligibility_io"] = {
        "fda": verify_fda_io_eligibility(
            analysis_result["tumor_context"]["tmb_score"],
            analysis_result["tumor_context"]["msi_status"]
        ),
        "nccn": verify_nccn_drug_recommendation("Carboplatin", "ovarian_cancer")
    }
    
    # Phase 2: Formula & Consistency Verification
    verification_report["checks"]["dna_repair_formula"] = {
        "correctness": verify_formula_correctness(
            analysis_result["pathway_scores"],
            analysis_result["sae_features"]["dna_repair_capacity"]
        ),
        "range": verify_expected_range(
            analysis_result["sae_features"]["dna_repair_capacity"],
            (0.75, 0.90)
        )
    }
    
    verification_report["checks"]["mechanism_vector"] = {
        "structure": verify_vector_structure(
            analysis_result["sae_features"]["mechanism_vector"],
            7
        ),
        "mapping": verify_pathway_mapping(
            analysis_result["pathway_scores"],
            analysis_result["sae_features"]["mechanism_vector"]
        )
    }
    
    verification_report["checks"]["consistency"] = {
        "pathway": verify_pathway_consistency(
            analysis_result["efficacy_response"],
            analysis_result["sae_features"]
        ),
        "variant": verify_variant_consistency(
            analysis_result["mutations"],
            analysis_result
        )
    }
    
    # Phase 3: Biological Plausibility
    verification_report["checks"]["biological_plausibility"] = {
        "pathway_ranges": verify_pathway_ranges(
            analysis_result["pathway_scores"],
            {"ddr": (0.70, 0.90), "mapk": (0.10, 0.30)}
        ),
        "drug_efficacy_ranges": verify_drug_efficacy_ranges(
            analysis_result["drug_rankings"],
            {"PARP": (0.70, 0.90), "Platinum": (0.60, 0.80), "MEK": (0.20, 0.40)}
        )
    }
    
    # Compute overall verification score
    all_checks = []
    for category, checks in verification_report["checks"].items():
        if isinstance(checks, dict):
            for check_name, check_result in checks.items():
                if isinstance(check_result, dict) and "verified" in check_result:
                    all_checks.append(check_result["verified"])
    
    verification_report["overall_score"] = {
        "total_checks": len(all_checks),
        "passed": sum(all_checks),
        "failed": len(all_checks) - sum(all_checks),
        "pass_rate": sum(all_checks) / len(all_checks) if all_checks else 0.0
    }
    
    return verification_report

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 verify_mbd4_analysis.py <analysis_result.json>")
        sys.exit(1)
    
    analysis_path = sys.argv[1]
    report = run_all_verifications(analysis_path)
    
    # Save verification report
    output_path = analysis_path.replace(".json", "_verification.json")
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"Verification complete. Report saved to: {output_path}")
    print(f"Overall pass rate: {report['overall_score']['pass_rate']:.1%}")
```

**Integration**: Call after `run_mbd4_tp53_analysis.py` completes

**Output**: Comprehensive verification report JSON

---

### **Task 4.2: Integrate Verification into Analysis Pipeline**

**File**: `scripts/sae/run_mbd4_tp53_analysis.py` (MODIFY)

**Changes**:

1. Add verification step after analysis completes:
   ```python
   # After analysis completes
   if args.verify:
       from verify_mbd4_analysis import run_all_verifications
       verification_report = run_all_verifications(output_path)
       
       # Add verification to analysis result
       analysis_result["verification"] = verification_report
       
       # Save updated result
       with open(output_path, 'w') as f:
           json.dump(analysis_result, f, indent=2)
   ```

2. Add `--verify` flag to command-line arguments

**Usage**:
```bash
python3 scripts/sae/run_mbd4_tp53_analysis.py --verify
```

---

### **Task 4.3: Create Verification Report Generator**

**File**: `scripts/sae/generate_verification_report.py` (NEW)

**Purpose**: Generate human-readable verification report from JSON

**Output Format**:

```markdown
# MBD4+TP53 Analysis Verification Report

**Date**: 2025-01-21  
**Analysis File**: mbd4_tp53_analysis.json  
**Overall Pass Rate**: 92.5% (37/40 checks passed)

---

## Phase 1: Deterministic Verification

### Variant Classification
- ✅ **ClinVar Check**: MBD4 frameshift = Pathogenic (VERIFIED)
- ✅ **COSMIC Check**: TP53 R175H = Hotspot (VERIFIED)
- ✅ **Evo2 Validation**: Delta scores highly negative (VERIFIED)

### Pathway Mapping
- ✅ **KEGG Check**: MBD4 → DDR pathway (VERIFIED)
- ✅ **Reactome Check**: TP53 → DDR pathway (VERIFIED)
- ✅ **Formula Check**: DNA repair capacity formula correct (VERIFIED)
- ✅ **TCGA Validation**: Pathway weights from real mutation frequencies (VERIFIED)

### Functional Annotation
- ✅ **UniProt Check**: MBD4 = DNA glycosylase (VERIFIED)
- ✅ **Insights Bundle**: Scores within expected ranges (VERIFIED)

### Eligibility & IO
- ✅ **FDA Labels**: TMB ≥20 → IO eligible (VERIFIED)
- ✅ **NCCN Guidelines**: Carboplatin first-line (VERIFIED)

---

## Phase 2: Formula & Consistency Verification

### DNA Repair Capacity Formula
- ✅ **Formula Correctness**: (0.6 × DDR) + (0.2 × HRR) + (0.2 × exon) = 0.82 (VERIFIED)
- ✅ **Expected Range**: 0.75-0.90 (VERIFIED)

### Mechanism Vector
- ✅ **Vector Structure**: 7D vector, all values in [0.0, 1.0] (VERIFIED)
- ✅ **Pathway Mapping**: DDR index = 0.88, TP53 contribution correct (VERIFIED)

### Consistency Checks
- ✅ **Pathway Consistency**: Efficacy and SAE pathway scores match (VERIFIED)
- ✅ **Variant Consistency**: Input mutations match output annotations (VERIFIED)

---

## Phase 3: Biological Plausibility

### Pathway Score Ranges
- ✅ **DDR Range**: 0.88 (expected: 0.70-0.90) (VERIFIED)
- ✅ **MAPK Range**: 0.12 (expected: 0.10-0.30) (VERIFIED)

### Drug Efficacy Ranges
- ✅ **PARP Inhibitors**: 0.85 (expected: 0.70-0.90) (VERIFIED)
- ✅ **Platinum**: 0.72 (expected: 0.60-0.80) (VERIFIED)
- ✅ **MEK Inhibitors**: 0.28 (expected: 0.20-0.40) (VERIFIED)

---

## Summary

**Total Checks**: 40  
**Passed**: 37  
**Failed**: 3  
**Pass Rate**: 92.5%

**Failed Checks**:
1. ⚠️ Food recommendation evidence quality (variable, LLM extraction)
2. ⚠️ Neoantigen prediction (heuristic, not sequence-based)
3. ⚠️ Dosage extraction (regex + LLM fallback, may be inaccurate)

**Note**: Failed checks are in Level 4 (Speculative/Heuristic) category, which is expected for research-use-only features.
```

**Integration**: Call after verification completes

**Output**: Human-readable Markdown report

---

## 📋 PHASE 5: DATA SOURCES & APIS (MEDIUM PRIORITY)

### **Task 5.1: Set Up External API Clients**

**File**: `scripts/sae/verification_clients.py` (NEW)

**Purpose**: Create clients for external verification sources

**APIs to Integrate**:

1. **ClinVar API**:
   ```python
   def query_clinvar(gene: str, variant: str) -> Dict:
       """
       Query ClinVar API for variant classification.
       
       API: https://www.ncbi.nlm.nih.gov/clinvar/docs/api/
       """
       # Implementation using requests or httpx
   ```

2. **COSMIC API** (if available) or Local Database:
   ```python
   def query_cosmic_hotspot(gene: str, hgvs_p: str) -> Dict:
       """
       Query COSMIC for hotspot mutation status.
       
       Note: COSMIC API may require authentication
       Fallback: Use local COSMIC hotspot database
       """
       # Implementation
   ```

3. **KEGG API**:
   ```python
   def query_kegg_pathway(gene: str) -> List[str]:
       """
       Query KEGG API for gene pathways.
       
       API: https://www.genome.jp/kegg/rest/keggapi.html
       """
       # Implementation
   ```

4. **Reactome API**:
   ```python
   def query_reactome_pathway(gene: str) -> List[str]:
       """
       Query Reactome API for gene pathways.
       
       API: https://reactome.org/ContentService/
       """
       # Implementation
   ```

5. **UniProt API**:
   ```python
   def query_uniprot_function(gene: str) -> Dict:
       """
       Query UniProt API for protein function.
       
       API: https://www.uniprot.org/help/api
       """
       # Implementation
   ```

**Integration**: Use in verification scripts

**Caching**: Cache API responses to avoid rate limits

---

### **Task 5.2: Create Local Verification Databases**

**File**: `data/verification/` (NEW DIRECTORY)

**Purpose**: Store local copies of verification data for offline use

**Files to Create**:

1. `cosmic_hotspots.json`: COSMIC hotspot mutations (KRAS, BRAF, TP53, etc.)
2. `kegg_pathways.json`: KEGG gene→pathway mappings (common genes)
3. `reactome_pathways.json`: Reactome gene→pathway mappings (common genes)
4. `uniprot_functions.json`: UniProt protein functions (common genes)
5. `fda_labels.json`: FDA drug labels (IO eligibility, PARP approvals)
6. `nccn_guidelines.json`: NCCN guideline recommendations (ovarian cancer)

**Update Frequency**: Monthly or as needed

**Usage**: Fallback when APIs unavailable

---

## 📋 PHASE 6: TESTING & VALIDATION (HIGH PRIORITY)

### **Task 6.1: Create Verification Test Suite**

**File**: `tests/test_verification_layer.py` (NEW)

**Test Cases**:

1. **Test Variant Classification Verification**:
   ```python
   def test_clinvar_verification():
       """Test ClinVar classification verification"""
       result = verify_clinvar_classification({
           "gene": "MBD4",
           "hgvs_p": "p.Ile413Serfs*2"
       })
       assert result["verified"] == True
       assert result["expected"] == "Pathogenic"
   ```

2. **Test Pathway Mapping Verification**:
   ```python
   def test_kegg_pathway_verification():
       """Test KEGG pathway mapping verification"""
       result = verify_kegg_pathway("MBD4", "DDR")
       assert result["verified"] == True
       assert "DDR" in result["actual"]
   ```

3. **Test Formula Verification**:
   ```python
   def test_dna_repair_formula_verification():
       """Test DNA repair capacity formula verification"""
       pathway_scores = {"ddr": 0.80, "mapk": 0.10}
       insights = {"essentiality_hrr": 0.75, "exon_disruption_score": 0.70}
       computed = 0.82
       
       result = verify_formula_correctness(pathway_scores, computed, insights)
       assert result["verified"] == True
       assert abs(result["expected"] - computed) < 0.01
   ```

**Coverage**: Test all verification functions

---

### **Task 6.2: Integration Test with MBD4+TP53 Analysis**

**File**: `tests/test_mbd4_verification_integration.py` (NEW)

**Purpose**: End-to-end test of verification layer with real analysis

**Test Flow**:

1. Run MBD4+TP53 analysis
2. Run verification layer
3. Assert all deterministic checks pass
4. Assert formula checks pass
5. Assert consistency checks pass

**Expected Results**:

- Deterministic checks: 100% pass rate
- Formula checks: 100% pass rate
- Consistency checks: 100% pass rate
- Biological plausibility: 90%+ pass rate

---

## 📋 DELIVERABLES

1. **Verification Scripts** (8 new files):
   - `verify_variant_classification.py`
   - `verify_pathway_mapping.py`
   - `verify_functional_annotation.py`
   - `verify_eligibility_io.py`
   - `verify_dna_repair_formula.py`
   - `verify_mechanism_vector.py`
   - `verify_consistency.py`
   - `verify_biological_plausibility.py`

2. **Unified Verification Script**:
   - `verify_mbd4_analysis.py` (runs all checks)

3. **Report Generator**:
   - `generate_verification_report.py` (human-readable output)

4. **Verification Clients**:
   - `verification_clients.py` (external API clients)

5. **Local Databases**:
   - `data/verification/` directory with local copies

6. **Test Suite**:
   - `tests/test_verification_layer.py`
   - `tests/test_mbd4_verification_integration.py`

7. **Documentation**:
   - Verification framework documentation
   - API integration guides
   - Expected results reference

---

## 🎯 SUCCESS CRITERIA

1. ✅ All deterministic checks automated (ClinVar, COSMIC, KEGG, Reactome, FDA, NCCN)
2. ✅ All formula checks automated (DNA repair capacity, mechanism vectors)
3. ✅ All consistency checks automated (pathway scores, variant annotations)
4. ✅ Biological plausibility checks automated (expected ranges)
5. ✅ Unified verification script runs all checks
6. ✅ Human-readable report generated
7. ✅ Test suite covers all verification functions
8. ✅ Integration test passes with MBD4+TP53 analysis
9. ✅ Verification integrated into analysis pipeline (optional `--verify` flag)

---

## ⚠️ LIMITATIONS & OUT OF SCOPE

**What We CANNOT Verify (Out of Scope)**:

1. ❌ **Clinical Validation**: Drug efficacy scores require patient outcomes
2. ❌ **Prospective Validation**: Resistance prediction requires 3-6 month follow-up
3. ❌ **LLM Extraction Quality**: Variable, requires manual review
4. ❌ **Neoantigen Prediction**: Heuristic (TMB proxy), not sequence-based
5. ❌ **Food/Supplement Recommendations**: Research use only, variable evidence quality

**What We CAN Verify (In Scope)**:

1. ✅ **Deterministic Sources**: ClinVar, COSMIC, KEGG, Reactome, FDA, NCCN
2. ✅ **Formula Correctness**: DNA repair capacity, mechanism vectors
3. ✅ **Consistency**: Pathway scores, variant annotations
4. ✅ **Biological Plausibility**: Expected ranges, known biology

---

## 📊 IMPLEMENTATION PRIORITY

**P0 (Immediate - This Week)**:
- Task 1.1: Variant Classification Verification
- Task 1.2: Pathway Mapping Verification
- Task 2.1: DNA Repair Capacity Formula Verification
- Task 4.1: Unified Verification Script

**P1 (High Priority - Next Week)**:
- Task 1.3: Functional Annotation Verification
- Task 1.4: Eligibility & IO Verification
- Task 2.2: Mechanism Vector Verification
- Task 2.3: Consistency Checks
- Task 4.2: Integration into Analysis Pipeline

**P2 (Medium Priority - Week 3)**:
- Task 3.1: Biological Plausibility Verification
- Task 4.3: Report Generator
- Task 5.1: External API Clients
- Task 5.2: Local Verification Databases

**P3 (Lower Priority - Week 4)**:
- Task 6.1: Verification Test Suite
- Task 6.2: Integration Test

---

## 🚀 EXECUTION TIMELINE

**Week 1**: P0 tasks (deterministic + formula verification)  
**Week 2**: P1 tasks (functional + eligibility + consistency + integration)  
**Week 3**: P2 tasks (plausibility + reports + APIs + databases)  
**Week 4**: P3 tasks (testing + integration tests)

**Total**: 4 weeks to complete verification layer

---

**DOCTRINE STATUS: ACTIVE** ⚔️  
**LAST UPDATED**: January 21, 2025  
**NEXT STEP**: Start with P0 tasks (Task 1.1, 1.2, 2.1, 4.1)

