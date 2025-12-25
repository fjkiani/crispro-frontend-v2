# 🧪 Synthetic Lethality Benchmark Plan

**Date:** January 28, 2025  
**Goal:** Validate accuracy of synthetic lethality predictions using 100 test cases  
**Status:** 📋 READY FOR IMPLEMENTATION  
**Version:** 2.0 (Manager Review - Enhanced)

**Note:** This is a detailed implementation plan. For overview, see `SYNTHETIC_LETHALITY_COMPLETE.md`

---

## 🔍 Manager Review Summary

**Issues Fixed:**
1. ✅ Aligned with actual `/api/guidance/synthetic_lethality` response format
2. ✅ Referenced existing benchmark patterns (`benchmark_sota_*.py`)
3. ✅ Added practical data collection strategy
4. ✅ Clarified DepMap integration approach
5. ✅ Added ablation studies for component validation
6. ✅ Included comparison to existing SOTA benchmarks
7. ✅ Added phased rollout (10 → 50 → 100 cases)
8. ✅ Specified realistic ground truth sources
9. ✅ Added missing implementation details

---

## 🎯 Objectives (Enhanced)

1. **Measure Prediction Accuracy** - How well do we predict synthetic lethality?
2. **Validate Drug Recommendations** - Are recommended drugs clinically appropriate?
3. **Assess Essentiality Scoring** - Do our scores match known gene essentiality?
4. **Test Pathway Detection** - Can we correctly identify broken/essential pathways?
5. **Compare to Ground Truth** - Benchmark against published data and clinical outcomes
6. **Ablation Studies** - Which components (S/P/E) contribute most to accuracy?
7. **Comparison to Existing Benchmarks** - How do we compare to ovarian/MM benchmarks?

---

## 📊 Validation Metrics (Same, Validated)

### 1. **Drug Recommendation Accuracy**
- **Metric:** Top-1, Top-3, Top-5 accuracy
- **Definition:** % of cases where correct drug appears in top N recommendations
- **Ground Truth:** Known effective drugs for specific mutations (e.g., PARP for BRCA1/BRCA2)

### 2. **Synthetic Lethality Detection**
- **Metric:** True Positive Rate (TPR), False Positive Rate (FPR), Precision, Recall
- **Definition:** Can we detect known SL pairs? (e.g., BRCA1 loss + PARP inhibition)
- **Ground Truth:** Published synthetic lethality relationships

### 3. **Essentiality Score Calibration**
- **Metric:** Pearson correlation with DepMap CRISPR scores
- **Definition:** Do our scores correlate with experimental gene essentiality?
- **Ground Truth:** DepMap 24Q4 public data (Achilles_gene_effect.csv)

### 4. **Pathway Prediction Accuracy**
- **Metric:** Jaccard similarity for pathway sets
- **Definition:** % of cases where we correctly identify broken/essential pathways
- **Ground Truth:** Known pathway disruptions from literature

### 5. **Clinical Relevance Score**
- **Metric:** FDA approval match rate
- **Definition:** % of top recommendations that are FDA-approved for the indication
- **Ground Truth:** FDA drug labels and NCCN guidelines

---

## 📁 Test Dataset Structure (REVISED - Phased Approach)

### Phase 1: 10 Cases (Pilot - Days 1-2)
- 5 known BRCA1/BRCA2 + PARP cases (high confidence)
- 3 negative controls (no known SL)
- 2 edge cases
- **Goal:** Validate benchmark infrastructure

### Phase 2: 50 Cases (Validation - Week 2)
- 20 known SL pairs (BRCA, MBD4, TP53, etc.)
- 10 negative controls
- 15 diverse cancer types
- 5 edge cases
- **Goal:** Get reliable metric estimates

### Phase 3: 100 Cases (Full Benchmark - Week 3)
- 40 known SL pairs
- 20 negative controls
- 30 diverse cancer types
- 10 edge cases
- **Goal:** Publication-ready results

---

## 🔧 Actual API Response Format (CRITICAL FIX)

**Current `/api/guidance/synthetic_lethality` Returns:**

```python
{
  "suggested_therapy": "platinum",  # Single therapy suggestion
  "damage_report": [
    {
      "variant": {...},
      "vep": {...},
      "functionality": {...}
    }
  ],
  "essentiality_report": [
    {
      "gene": "BRCA1",
      "result": {
        "essentiality_score": 0.85,
        "flags": {"truncation": true, "frameshift": false},
        "rationale": "...",
        "confidence": 0.70,
        "pathway_impact": "HR pathway NON-FUNCTIONAL"
      }
    }
  ],
  "guidance": {
    # Chemo guidance payload
  }
}
```

**⚠️ Issue:** No `recommended_therapies` array or `pathway_analysis` object!

**Solution:** Need to either:
1. Enhance `/api/guidance/synthetic_lethality` to return structured output
2. OR parse `suggested_therapy` + `essentiality_report` to infer rankings
3. OR call `/api/efficacy/predict` separately for ranked drugs

**Recommendation:** Use hybrid approach:
- Call `/api/guidance/synthetic_lethality` for SL detection
- Parse `essentiality_report` for essentiality scores
- Infer pathways from `damage_report` + `essentiality_report`
- Compare `suggested_therapy` to ground truth (binary match)

---

## 🗂️ Ground Truth Data Sources (PRACTICAL)

### 1. **DepMap Data (Free, Public)**
- **Source:** https://depmap.org/portal/download/all/
- **File:** `Achilles_gene_effect.csv` (24Q4 release)
- **Contains:** CRISPR knockout scores for ~18K genes × 1000+ cell lines
- **Usage:** Extract essentiality scores for ovarian/breast cell lines
- **Script:** `scripts/benchmark_sl/download_depmap.py`

### 2. **Known SL Pairs from Literature**
- **Source:** Manual curation from PubMed
- **Key Papers:**
  - Lord & Ashworth (2017) - PARP inhibitors in BRCA-deficient tumors
  - Nijman (2011) - Synthetic lethality and cancer
  - DepMap publications (2019-2024)
- **Format:** `known_sl_pairs.json`
- **Example:**
```json
[
  {
    "gene1": "BRCA1",
    "gene2": "PARP1",
    "drugs": ["Olaparib", "Niraparib", "Rucaparib"],
    "evidence": "FDA approved",
    "pmid": "28355133",
    "cancer_types": ["ovarian", "breast"]
  }
]
```

### 3. **FDA Drug Labels (Automated)**
- **Source:** DailyMed API (https://dailymed.nlm.nih.gov/dailymed/)
- **Script:** `scripts/benchmark_sl/scrape_fda_labels.py`
- **Extract:** Indications for PARP inhibitors, ATR inhibitors, WEE1 inhibitors
- **Cache:** `fda_drug_labels.json`

### 4. **TCGA Clinical Data (If Available)**
- **Source:** Existing `tools/benchmarks/hrd_tcga_ov_labeled_1k_results.json`
- **Contains:** 1000 ovarian cancer cases with platinum response
- **Reuse:** Can extract BRCA mutations and test SL predictions
- **Advantage:** Already processed and labeled!

### 5. **Negative Controls (Important!)**
- **Source:** ClinVar benign variants + DepMap non-essential genes
- **Strategy:**
  - Pick genes with LOW DepMap dependency (e.g., non-essential genes)
  - Pick benign/likely benign variants from ClinVar
  - Expected: Should NOT predict high SL or recommend aggressive therapies
- **Example genes:** Benign variants in OR genes, HLA genes, etc.

---

## 🔧 Implementation Plan (REVISED)

### Phase 1: Infrastructure & Pilot (Days 1-3)

**Task 1.1: Download Ground Truth Data**
- **File:** `scripts/benchmark_sl/download_depmap.py`
```python
#!/usr/bin/env python3
"""Download and process DepMap data for benchmarking."""
import pandas as pd
import requests

DEPMAP_URL = "https://depmap.org/portal/download/api/download/external?file_name=public_24Q4_Achilles_gene_effect.csv"

def download_depmap():
    """Download DepMap gene effect scores."""
    print("Downloading DepMap data...")
    df = pd.read_csv(DEPMAP_URL)
    
    # Filter to ovarian and breast cell lines
    ovarian_lines = df[df['cell_line_name'].str.contains('OVARY', case=False, na=False)]
    breast_lines = df[df['cell_line_name'].str.contains('BREAST', case=False, na=False)]
    
    # Extract key genes: BRCA1, BRCA2, MBD4, TP53, PARP1, ATR, WEE1
    key_genes = ['BRCA1', 'BRCA2', 'MBD4', 'TP53', 'PARP1', 'ATR', 'WEE1']
    
    # Save to JSON
    essentiality_scores = {}
    for gene in key_genes:
        if gene in df.columns:
            essentiality_scores[gene] = {
                "mean_score": df[gene].mean(),
                "ovarian_mean": ovarian_lines[gene].mean(),
                "breast_mean": breast_lines[gene].mean()
            }
    
    with open('depmap_essentiality.json', 'w') as f:
        json.dump(essentiality_scores, f, indent=2)
    
    print(f"Saved essentiality scores for {len(key_genes)} genes")
```

**Task 1.2: Create 10-Case Pilot Dataset**
- **File:** `scripts/benchmark_sl/create_pilot_dataset.py`
- **Strategy:** Start with Ayesha-like cases (known BRCA1/BRCA2)
- **Example Cases:**
  1. BRCA1 C61G (known SL with PARP)
  2. BRCA2 truncating (known SL with PARP)
  3. MBD4 frameshift + TP53 (Ayesha's case)
  4. TP53 hotspot alone (moderate SL)
  5. BRCA1 + TP53 (double-hit)
  6-8. Negative controls (benign variants)
  9-10. Edge cases

**Task 1.3: Adapt Existing Benchmark Pattern**
- **Reference:** `scripts/benchmark_sota_ovarian.py` (already working)
- **Adapt:** Change API call to `/api/guidance/synthetic_lethality`
- **File:** `scripts/benchmark_sl/benchmark_synthetic_lethality.py`
- **Enhancement:** Add parallel execution for faster throughput (10 concurrent requests)
- **Enhancement:** Add result caching to avoid re-running expensive API calls

```python
#!/usr/bin/env python3
"""
Synthetic Lethality Prediction Benchmark
Adapted from benchmark_sota_ovarian.py pattern
"""
import asyncio
import json
from pathlib import Path
from typing import List, Dict, Any
from datetime import datetime

import httpx
from sklearn.metrics import accuracy_score, precision_recall_fscore_support
import numpy as np

API_ROOT = "http://127.0.0.1:8000"

async def predict_sl(client: httpx.AsyncClient, case: Dict[str, Any]) -> Dict[str, Any]:
    """Call synthetic lethality API."""
    try:
        payload = {
            "disease": case["disease"],
            "mutations": case["mutations"]
        }
        
        resp = await client.post(
            f"{API_ROOT}/api/guidance/synthetic_lethality",
            json=payload,
            timeout=120.0
        )
        resp.raise_for_status()
        return resp.json()
    except Exception as e:
        print(f"Error for case {case.get('case_id')}: {e}")
        return None

async def run_benchmark(test_file: str, max_concurrent: int = 10, cache_file: str = None):
    """Run benchmark on all test cases with parallel execution."""
    with open(test_file, 'r') as f:
        test_cases = json.load(f)
    
    # Load cache if exists
    cache = {}
    if cache_file and Path(cache_file).exists():
        with open(cache_file, 'r') as f:
            cache = json.load(f)
        print(f"Loaded {len(cache)} cached results")
    
    results = []
    semaphore = asyncio.Semaphore(max_concurrent)
    
    async def process_case(case: Dict[str, Any]) -> Dict:
        """Process a single case with semaphore for concurrency control."""
        async with semaphore:
            case_id = case['case_id']
            
            # Check cache first
            if case_id in cache:
                print(f"Using cached result for {case_id}")
                prediction = cache[case_id]
            else:
                print(f"Processing case {case_id}")
                async with httpx.AsyncClient() as client:
                    prediction = await predict_sl(client, case)
                    if prediction and cache_file:
                        cache[case_id] = prediction
            
            if prediction is None:
                return None
            
            # Compare to ground truth
            gt = case["ground_truth"]
            comparison = {
                "case_id": case_id,
                "ground_truth": gt,
                "prediction": {
                    "suggested_therapy": prediction.get("suggested_therapy"),
                    "essentiality_scores": [
                        {
                            "gene": e["gene"],
                            "score": e["result"]["essentiality_score"]
                        }
                        for e in prediction.get("essentiality_report", [])
                    ]
                },
                "metrics": calculate_metrics(gt, prediction)
            }
            return comparison
    
    # Process all cases in parallel
    tasks = [process_case(case) for case in test_cases]
    comparisons = await asyncio.gather(*tasks)
    results = [c for c in comparisons if c is not None]
    
    # Save cache
    if cache_file:
        with open(cache_file, 'w') as f:
            json.dump(cache, f, indent=2)
        print(f"Saved cache with {len(cache)} entries")
    
    # Aggregate metrics
    aggregate = aggregate_metrics(results)
    
    # Save results
    output = {
        "date": datetime.now().isoformat(),
        "num_cases": len(results),
        "aggregate_metrics": aggregate,
        "results": results
    }
    
    output_file = Path("results") / f"benchmark_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    output_file.parent.mkdir(exist_ok=True)
    
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\n✅ Benchmark complete! Results saved to {output_file}")
    print(f"\nAggregate Metrics:")
    for metric, value in aggregate.items():
        print(f"  {metric}: {value:.3f}")

def calculate_metrics(gt: Dict, pred: Dict) -> Dict:
    """Calculate metrics for a single case."""
    metrics = {}
    
    # Drug match (binary)
    suggested = pred.get("suggested_therapy", "").lower()
    effective_drugs = [d.lower() for d in gt.get("effective_drugs", [])]
    metrics["drug_match"] = any(drug in suggested for drug in effective_drugs)
    
    # Essentiality correlation (if available)
    if "depmap_essentiality" in gt and "essentiality_report" in pred:
        gt_scores = gt["depmap_essentiality"]
        pred_scores = {
            e["gene"]: e["result"]["essentiality_score"]
            for e in pred["essentiality_report"]
        }
        
        # Match genes
        common_genes = set(gt_scores.keys()) & set(pred_scores.keys())
        if common_genes:
            gt_vals = [gt_scores[g] for g in common_genes]
            pred_vals = [pred_scores[g] for g in common_genes]
            metrics["essentiality_correlation"] = np.corrcoef(gt_vals, pred_vals)[0, 1]
    
    return metrics

def aggregate_metrics(results: List[Dict]) -> Dict:
    """Aggregate metrics across all cases with confidence intervals."""
    drug_matches = [r["metrics"].get("drug_match", False) for r in results]
    correlations = [r["metrics"].get("essentiality_correlation") for r in results if "essentiality_correlation" in r["metrics"]]
    
    # Calculate with confidence intervals (95% CI using bootstrap or normal approximation)
    drug_accuracy = sum(drug_matches) / len(drug_matches) if drug_matches else 0.0
    n = len(drug_matches)
    # Wilson score interval for proportion
    from scipy.stats import norm
    z = norm.ppf(0.975)  # 95% CI
    p = drug_accuracy
    ci_lower = (p + z**2/(2*n) - z*np.sqrt((p*(1-p) + z**2/(4*n))/n)) / (1 + z**2/n)
    ci_upper = (p + z**2/(2*n) + z*np.sqrt((p*(1-p) + z**2/(4*n))/n)) / (1 + z**2/n)
    
    mean_corr = np.mean(correlations) if correlations else 0.0
    std_corr = np.std(correlations) if correlations else 0.0
    corr_ci = (mean_corr - 1.96*std_corr/np.sqrt(len(correlations)), 
               mean_corr + 1.96*std_corr/np.sqrt(len(correlations))) if correlations else (0.0, 0.0)
    
    return {
        "drug_accuracy": drug_accuracy,
        "drug_accuracy_ci_95": (ci_lower, ci_upper),
        "mean_essentiality_correlation": mean_corr,
        "essentiality_correlation_ci_95": corr_ci,
        "num_cases_with_correlation": len(correlations),
        "total_cases": len(results)
    }

if __name__ == "__main__":
    import sys
    test_file = sys.argv[1] if len(sys.argv) > 1 else "test_cases_pilot.json"
    asyncio.run(run_benchmark(test_file))
```

### Phase 2: Validation with 50 Cases (Days 4-10)

**Task 2.1: Expand Dataset to 50 Cases**
- Reuse TCGA-OV data where possible
- Add manual cases for diversity
- Include more negative controls

**Task 2.2: Add Ablation Studies**
- Run with different ablation modes (S-only, P-only, E-only, SP, SE, PE, SPE)
- Measure contribution of each component
- **Script:** `scripts/benchmark_sl/run_ablations.py`
- **Note:** Ablations require API changes or separate endpoints - may need to call `/api/efficacy/predict` with different options
- **Alternative:** Run full pipeline and analyze which components (S/P/E) contributed most to correct predictions

**Task 2.3: Compare to Existing SOTA Benchmarks**
- **Ovarian Benchmark:** `benchmark_sota_ovarian.py` (AUROC ~0.75 target)
- **MM Benchmark:** `benchmark_sota_mm.py` (100% pathway accuracy)
- **Compare:** Is SL detection better/worse than general efficacy prediction?
- **Method:** Run same test cases through `/api/efficacy/predict` and compare drug recommendations
- **Metric:** Agreement rate between SL-specific and general efficacy predictions

### Phase 3: Full 100-Case Benchmark (Days 11-15)

**Task 3.1: Complete Dataset Curation**
- Reach 100 cases with full coverage
- Balance across categories
- Ensure diverse representation

**Task 3.2: Run Full Benchmark**
- Execute on all 100 cases
- Collect comprehensive metrics
- Generate visualizations

**Task 3.3: Generate Report**
- **File:** `scripts/benchmark_sl/generate_report.py`
- Compare to targets
- Identify failure modes
- Recommendations for improvement

---

## 📈 Expected Metrics & Targets (REVISED)

### Success Criteria (Realistic):

| Metric | Minimum | Target | Excellent |
|--------|---------|--------|-----------|
| **Drug Match (Binary)** | >50% | >70% | >85% |
| **Essentiality Correlation** | >0.5 | >0.7 | >0.85 |
| **SL Detection TPR** | >60% | >75% | >90% |
| **SL Detection FPR** | <30% | <20% | <10% |

**Note:** Lowered targets for initial benchmark since:
1. API response format is limited
2. Ground truth is incomplete
3. Negative controls are hard to define
4. First iteration - expect to iterate

### Baseline Comparison:

- **Random Baseline:** ~20% (1 in 5 drug match)
- **Rule-Based:** ~50% (if BRCA → PARP)
- **Our Target:** >70% (better than rules)

---

## 🔄 Continuous Validation (ADDED)

### Version Control for Benchmarks:
```
results/
├── v1_baseline/
│   └── benchmark_20250128.json
├── v2_enhanced_pathways/
│   └── benchmark_20250205.json
└── comparison_v1_v2.md
```

### Track Over Time:
- **Metric:** Drug match accuracy
- **Expected:** Should improve with each release
- **Alert:** If drops >5%, investigate regression

---

## ✅ Success Criteria (REVISED)

**Phase 1 Success (Pilot):**
- ✅ 10 cases run without errors
- ✅ Metrics calculated correctly
- ✅ At least 50% drug match accuracy (random > 20%)

**Phase 2 Success (Validation):**
- ✅ 50 cases run successfully
- ✅ Drug match accuracy >60%
- ✅ Essentiality correlation >0.5

**Phase 3 Success (Full Benchmark):**
- ✅ 100 cases complete
- ✅ Drug match accuracy >70%
- ✅ Essentiality correlation >0.7
- ✅ Report generated with recommendations

**Ready for Clinical Use:**
- ✅ Drug match accuracy >75%
- ✅ SL detection TPR >75%, FPR <20%
- ✅ Validated on diverse cancer types
- ✅ Error analysis complete

---

## 🚨 Known Limitations (TRANSPARENCY)

1. **API Response Format** - Limited to `suggested_therapy` (single drug, not ranked list)
2. **Pathway Detection** - Must infer from `essentiality_report` (no explicit pathway_analysis object)
3. **Ground Truth** - Incomplete for many SL pairs (literature gaps, especially for rare variants)
4. **Negative Controls** - Hard to define "no SL" definitively (absence of evidence ≠ evidence of absence)
5. **Clinical Outcomes** - Real patient data limited (mostly preclinical/cell line data)
6. **Sample Size** - 100 cases small for publication (aim for 500+ in future)
7. **Binary Drug Match** - Can't measure Top-3/Top-5 accuracy (only binary match/no-match)
8. **Fast-Path Mode** - API has `GUIDANCE_FAST` mode that short-circuits for DDR genes (may miss nuanced predictions)
9. **No Pathway Ranking** - Can't validate which pathways are most critical
10. **Temporal Drift** - DepMap data updates quarterly, ground truth may become stale

**Mitigation:**
- Be transparent about limitations in report
- Focus on known-positive cases first (high-confidence ground truth)
- Use DepMap as proxy for essentiality (validated experimental data)
- Plan for larger benchmark in future (500-1000 cases)
- Consider enhancing API to return ranked drug list (parallel track)
- Document fast-path behavior and test with `GUIDANCE_FAST=0` for comparison
- Use multiple ground truth sources (DepMap + literature + clinical trials) for triangulation
- Add statistical significance testing to distinguish signal from noise

---

## 📝 Next Steps (ACTIONABLE)

**Week 1 (Manager Approval Required):**
1. ☐ Review and approve this plan
2. ☐ Download DepMap data (`download_depmap.py`)
3. ☐ Create 10-case pilot dataset
4. ☐ Run pilot benchmark
5. ☐ Validate metrics make sense

**Week 2 (Execute):**
1. ☐ Expand to 50 cases
2. ☐ Run ablation studies
3. ☐ Compare to existing benchmarks

**Week 3 (Complete):**
1. ☐ Reach 100 cases
2. ☐ Generate final report
3. ☐ Present findings
4. ☐ Plan improvements

---

## 📊 Manager Decision Points

**Decision 1:** Start with 10-case pilot or jump to 50?
- **Recommendation:** Start with 10 (de-risk, validate infrastructure)

**Decision 2:** Enhance API to return structured output?
- **Recommendation:** Yes, but parallel track (don't block benchmark)

**Decision 3:** How much manual curation vs. automated?
- **Recommendation:** 60% automated (DepMap, TCGA), 40% manual (literature)

**Decision 4:** Target for "good enough" to proceed?
- **Recommendation:** >70% drug match accuracy on pilot

**Decision 5:** Enable fast-path mode or full analysis?
- **Recommendation:** Test both - fast-path for speed, full analysis for accuracy comparison

**Decision 6:** Parallel execution or sequential?
- **Recommendation:** Parallel (10 concurrent) for speed, but respect API rate limits

---

**Status:** ✅ **READY FOR MANAGER APPROVAL & IMPLEMENTATION**

**Key Improvements Over V1:**
1. Phased approach (10 → 50 → 100)
2. Aligned with actual API response format
3. Practical ground truth sources (DepMap, TCGA)
4. References existing benchmark patterns
5. Ablation studies included
6. Realistic targets based on limitations
7. Clear decision points for manager

**Estimated Effort:** 2-3 weeks (one person full-time)

