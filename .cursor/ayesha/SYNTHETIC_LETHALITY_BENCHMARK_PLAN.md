# 🧪 Synthetic Lethality Prediction Validation & Benchmark Plan

**Date:** January 28, 2025  
**Goal:** Validate accuracy of synthetic lethality predictions using 100 test cases  
**Status:** 📋 PLANNING PHASE

---

## 🎯 Objectives

1. **Measure Prediction Accuracy** - How well do we predict synthetic lethality?
2. **Validate Drug Recommendations** - Are recommended drugs clinically appropriate?
3. **Assess Essentiality Scoring** - Do our scores match known gene essentiality?
4. **Test Pathway Detection** - Can we correctly identify broken/essential pathways?
5. **Compare to Ground Truth** - Benchmark against published data and clinical outcomes

---

## 📊 Validation Metrics

### 1. **Drug Recommendation Accuracy**
- **Metric:** Top-1, Top-3, Top-5 accuracy
- **Definition:** % of cases where correct drug appears in top N recommendations
- **Ground Truth:** Known effective drugs for specific mutations (e.g., PARP for BRCA1/BRCA2)

### 2. **Synthetic Lethality Detection**
- **Metric:** True Positive Rate (TPR), False Positive Rate (FPR), Precision, Recall
- **Definition:** Can we detect known SL pairs? (e.g., BRCA1 loss + PARP inhibition)
- **Ground Truth:** Published synthetic lethality relationships

### 3. **Essentiality Score Calibration**
- **Metric:** Correlation with DepMap essentiality scores
- **Definition:** Do our scores correlate with experimental gene essentiality?
- **Ground Truth:** DepMap CRISPR knockout data

### 4. **Pathway Prediction Accuracy**
- **Metric:** Pathway identification accuracy
- **Definition:** % of cases where we correctly identify broken/essential pathways
- **Ground Truth:** Known pathway disruptions (e.g., BER broken → HR essential)

### 5. **Clinical Relevance Score**
- **Metric:** FDA approval match rate
- **Definition:** % of top recommendations that are FDA-approved for the indication
- **Ground Truth:** FDA drug labels and NCCN guidelines

---

## 📁 Test Dataset Structure

### Dataset: 100 Curated Cases

**Categories:**

1. **Known Synthetic Lethality (40 cases)**
   - BRCA1/BRCA2 + PARP inhibitors (10 cases)
   - MBD4/BER + PARP inhibitors (5 cases)
   - TP53 + ATR/WEE1 inhibitors (5 cases)
   - Other known SL pairs (20 cases)

2. **Negative Controls (20 cases)**
   - Mutations without known SL relationships
   - Should NOT predict high synthetic lethality

3. **Diverse Cancer Types (30 cases)**
   - Ovarian (10), Breast (10), Other (10)
   - Mix of known and unknown SL relationships

4. **Edge Cases (10 cases)**
   - Multiple mutations (double-hit)
   - Germline + somatic combinations
   - Rare variants

### Ground Truth Format:

```json
{
  "case_id": "SL_001",
  "mutations": [
    {"gene": "BRCA1", "hgvs_p": "p.C61G", "consequence": "missense"},
    {"gene": "TP53", "hgvs_p": "p.R175H", "consequence": "missense"}
  ],
  "disease": "ovarian_cancer",
  "subtype": "high_grade_serous",
  "stage": "IVB",
  
  "ground_truth": {
    "synthetic_lethality_detected": true,
    "known_sl_pairs": ["BRCA1+PARP", "TP53+ATR"],
    "effective_drugs": ["Olaparib", "Niraparib", "Ceralasertib"],
    "broken_pathways": ["HR", "G1/S Checkpoint"],
    "essential_pathways": ["NER", "ATR/CHK1"],
    "depmap_essentiality": {
      "BRCA1": 0.92,
      "TP53": 0.88
    },
    "clinical_evidence": {
      "olaparib_fda_approved": true,
      "olaparib_indication": "BRCA-mutated ovarian cancer",
      "response_rate": 0.65
    }
  }
}
```

---

## 🔧 Implementation Plan

### Phase 1: Dataset Creation (Week 1)

**Task 1.1: Curate Test Cases**
- **File:** `oncology-coPilot/oncology-backend-minimal/scripts/benchmark_sl/create_test_dataset.py`
- **Sources:**
  - Published synthetic lethality papers (PubMed)
  - DepMap essentiality data
  - FDA drug labels
  - Clinical trial results (ClinicalTrials.gov)
  - Known SL databases (SynLethDB, if available)

**Task 1.2: Ground Truth Annotation**
- **File:** `oncology-coPilot/oncology-backend-minimal/scripts/benchmark_sl/annotate_ground_truth.py`
- **Process:**
  1. For each case, identify known SL relationships
  2. Document effective drugs from literature/clinical trials
  3. Extract DepMap essentiality scores
  4. Note FDA approval status
  5. Record pathway disruptions

**Task 1.3: Dataset Validation**
- **File:** `oncology-coPilot/oncology-backend-minimal/scripts/benchmark_sl/validate_dataset.py`
- **Checks:**
  - All cases have complete ground truth
  - No duplicates
  - Balanced distribution across categories
  - Representative of real-world scenarios

### Phase 2: Benchmark Script (Week 2)

**Task 2.1: Create Benchmark Runner**
- **File:** `oncology-coPilot/oncology-backend-minimal/scripts/benchmark_sl/benchmark_synthetic_lethality.py`

**Structure:**
```python
"""
Synthetic Lethality Prediction Benchmark

Runs 100 test cases through our synthetic lethality analyzer
and compares predictions to ground truth.
"""

import json
import asyncio
from pathlib import Path
from typing import Dict, List, Any
import httpx

# Load test dataset
TEST_DATASET = Path(__file__).parent / "test_cases_100.json"

# API endpoint
API_BASE = "http://localhost:8000"

async def run_benchmark():
    """Run all test cases and collect predictions."""
    test_cases = load_test_cases()
    results = []
    
    for case in test_cases:
        prediction = await predict_synthetic_lethality(case)
        comparison = compare_to_ground_truth(case, prediction)
        results.append(comparison)
    
    return results

async def predict_synthetic_lethality(case: Dict) -> Dict:
    """Call our synthetic lethality API."""
    async with httpx.AsyncClient() as client:
        response = await client.post(
            f"{API_BASE}/api/guidance/synthetic_lethality",
            json={
                "disease": case["disease"],
                "mutations": case["mutations"]
            }
        )
        return response.json()

def compare_to_ground_truth(case: Dict, prediction: Dict) -> Dict:
    """Compare prediction to ground truth and calculate metrics."""
    gt = case["ground_truth"]
    
    return {
        "case_id": case["case_id"],
        "metrics": {
            "drug_top1_accuracy": check_drug_in_top_n(
                gt["effective_drugs"], 
                prediction["recommended_therapies"], 
                n=1
            ),
            "drug_top3_accuracy": check_drug_in_top_n(
                gt["effective_drugs"], 
                prediction["recommended_therapies"], 
                n=3
            ),
            "sl_detection_tpr": check_sl_detection(
                gt["synthetic_lethality_detected"],
                prediction["pathway_analysis"]["double_hit_detected"]
            ),
            "pathway_accuracy": check_pathway_match(
                gt["broken_pathways"],
                prediction["pathway_analysis"]["broken_pathways"]
            ),
            "essentiality_correlation": calculate_essentiality_correlation(
                gt["depmap_essentiality"],
                prediction["essentiality"]
            )
        }
    }
```

**Task 2.2: Metric Calculation Functions**
- **File:** `oncology-coPilot/oncology-backend-minimal/scripts/benchmark_sl/metrics.py`

**Functions:**
```python
def check_drug_in_top_n(ground_truth_drugs: List[str], 
                        predicted_drugs: List[Dict], 
                        n: int) -> bool:
    """Check if any ground truth drug appears in top N predictions."""
    top_n_drugs = [d["drug"] for d in predicted_drugs[:n]]
    return any(gt_drug in top_n_drugs for gt_drug in ground_truth_drugs)

def check_sl_detection(gt_detected: bool, pred_detected: bool) -> Dict:
    """Calculate TPR, FPR for SL detection."""
    return {
        "true_positive": gt_detected and pred_detected,
        "false_positive": not gt_detected and pred_detected,
        "true_negative": not gt_detected and not pred_detected,
        "false_negative": gt_detected and not pred_detected
    }

def calculate_essentiality_correlation(gt_scores: Dict[str, float],
                                       pred_scores: List[Dict]) -> float:
    """Calculate Pearson correlation between predicted and DepMap scores."""
    # Match genes and calculate correlation
    pass

def check_pathway_match(gt_pathways: List[str],
                       pred_pathways: List[str]) -> float:
    """Calculate Jaccard similarity for pathway sets."""
    gt_set = set(gt_pathways)
    pred_set = set(pred_pathways)
    intersection = len(gt_set & pred_set)
    union = len(gt_set | pred_set)
    return intersection / union if union > 0 else 0.0
```

**Task 2.3: Report Generation**
- **File:** `oncology-coPilot/oncology-backend-minimal/scripts/benchmark_sl/generate_report.py`

**Output:**
- Summary statistics
- Per-metric breakdown
- Confusion matrices
- ROC curves (if applicable)
- Case-by-case analysis
- Visualizations

### Phase 3: Analysis & Reporting (Week 3)

**Task 3.1: Statistical Analysis**
- Calculate aggregate metrics across all 100 cases
- Category-wise breakdown (known SL vs. negative controls)
- Confidence intervals
- Statistical significance testing

**Task 3.2: Error Analysis**
- Identify failure modes
- Analyze false positives/negatives
- Common prediction errors
- Edge case handling

**Task 3.3: Visualization**
- **File:** `oncology-coPilot/oncology-backend-minimal/scripts/benchmark_sl/visualize_results.py`
- **Charts:**
  - Drug recommendation accuracy by rank
  - SL detection ROC curve
  - Essentiality score correlation scatter plot
  - Pathway accuracy heatmap
  - Per-category performance comparison

---

## 📈 Expected Metrics & Targets

### Success Criteria:

| Metric | Target | Excellent |
|--------|--------|-----------|
| **Drug Top-1 Accuracy** | >60% | >75% |
| **Drug Top-3 Accuracy** | >80% | >90% |
| **SL Detection TPR** | >70% | >85% |
| **SL Detection FPR** | <20% | <10% |
| **Pathway Accuracy (Jaccard)** | >0.6 | >0.8 |
| **Essentiality Correlation** | >0.7 | >0.85 |
| **FDA Approval Match Rate** | >70% | >85% |

### Baseline Comparison:

- **Random Baseline:** ~10-20% (depending on number of drugs)
- **Rule-Based Baseline:** ~40-50% (simple if-then rules)
- **Our Target:** >70% (ML/AI-enhanced)

---

## 📁 File Structure

```
oncology-coPilot/oncology-backend-minimal/scripts/benchmark_sl/
├── README.md                          # Benchmark documentation
├── create_test_dataset.py             # Generate 100 test cases
├── annotate_ground_truth.py           # Add ground truth labels
├── validate_dataset.py                # Validate dataset quality
├── benchmark_synthetic_lethality.py   # Main benchmark runner
├── metrics.py                         # Metric calculation functions
├── generate_report.py                 # Generate analysis report
├── visualize_results.py               # Create visualizations
├── test_cases_100.json                # Test dataset (100 cases)
├── results/                           # Benchmark results
│   ├── benchmark_YYYYMMDD.json       # Raw results
│   ├── report_YYYYMMDD.md             # Human-readable report
│   └── figures/                       # Visualizations
└── ground_truth_sources/              # Source data (optional, can reference external)
    ├── depmap_essentiality.csv
    ├── known_sl_pairs.json
    ├── fda_drug_labels.json
    └── clinical_trial_results.json
```

---

## 🔍 Ground Truth Sources

### 1. **Synthetic Lethality Databases**
- **SynLethDB** - Database of synthetic lethal gene pairs
- **DepMap** - Cancer dependency map (gene essentiality)
- **GDSC** - Genomics of Drug Sensitivity in Cancer

### 2. **Clinical Evidence**
- **FDA Drug Labels** - Approved indications
- **NCCN Guidelines** - Treatment recommendations
- **ClinicalTrials.gov** - Trial results
- **PubMed** - Published SL relationships

### 3. **Experimental Data**
- **DepMap CRISPR Screens** - Gene essentiality scores
- **Cell Line Drug Response** - GDSC, CCLE
- **Functional Genomics** - Pathway disruption data

### 4. **Known SL Pairs (Examples)**
```python
KNOWN_SL_PAIRS = [
    {"gene1": "BRCA1", "gene2": "PARP1", "drug": "Olaparib", "evidence": "FDA approved"},
    {"gene1": "BRCA2", "gene2": "PARP1", "drug": "Niraparib", "evidence": "FDA approved"},
    {"gene1": "MBD4", "gene2": "PARP1", "drug": "Olaparib", "evidence": "Preclinical"},
    {"gene1": "TP53", "gene2": "ATR", "drug": "Ceralasertib", "evidence": "Clinical trials"},
    {"gene1": "TP53", "gene2": "WEE1", "drug": "Adavosertib", "evidence": "Clinical trials"},
    # ... more pairs
]
```

---

## 🚀 Execution Plan

### Step 1: Dataset Creation (Days 1-3)
1. Research and compile 100 test cases
2. Annotate with ground truth
3. Validate dataset quality
4. Save to `test_cases_100.json`

### Step 2: Benchmark Implementation (Days 4-7)
1. Create benchmark runner script
2. Implement metric calculation functions
3. Add error handling and logging
4. Test on small subset (10 cases)

### Step 3: Full Benchmark Run (Day 8)
1. Run all 100 cases
2. Collect predictions
3. Calculate metrics
4. Generate raw results JSON

### Step 4: Analysis & Reporting (Days 9-10)
1. Statistical analysis
2. Error analysis
3. Generate visualizations
4. Write comprehensive report

### Step 5: Iteration (Days 11+)
1. Identify improvement opportunities
2. Fix issues found
3. Re-run benchmark
4. Compare before/after

---

## 📊 Report Template

### Benchmark Report Structure:

```markdown
# Synthetic Lethality Prediction Benchmark Report

**Date:** YYYY-MM-DD
**Dataset:** 100 test cases
**Model Version:** v2.0

## Executive Summary
- Overall accuracy: X%
- Key findings
- Recommendations

## Detailed Metrics

### 1. Drug Recommendation Accuracy
- Top-1: X%
- Top-3: X%
- Top-5: X%
- Analysis: [What worked well, what didn't]

### 2. Synthetic Lethality Detection
- TPR: X%
- FPR: X%
- Precision: X%
- Recall: X%
- ROC-AUC: X.X

### 3. Essentiality Score Calibration
- Correlation with DepMap: X.X
- Scatter plot: [figure]
- Analysis: [How well calibrated]

### 4. Pathway Prediction
- Accuracy: X%
- Jaccard similarity: X.X
- Common errors: [list]

### 5. Clinical Relevance
- FDA approval match: X%
- NCCN guideline alignment: X%

## Category Breakdown
- Known SL pairs: X% accuracy
- Negative controls: X% accuracy
- Diverse cancers: X% accuracy
- Edge cases: X% accuracy

## Error Analysis
- False positives: [examples]
- False negatives: [examples]
- Common failure modes: [list]

## Recommendations
1. [Improvement 1]
2. [Improvement 2]
3. [Improvement 3]

## Appendix
- Full results table
- Case-by-case analysis
- Methodology details
```

---

## 🔄 Continuous Validation

### Automated Benchmarking:
- Run benchmark on every major release
- Track metrics over time
- Alert on performance degradation
- Compare to previous versions

### Integration:
- Add to CI/CD pipeline (optional)
- Run weekly/monthly
- Store results in database
- Generate trend reports

---

## ✅ Success Criteria

**Benchmark is successful if:**
1. ✅ All 100 test cases run without errors
2. ✅ Metrics calculated correctly
3. ✅ Report generated with visualizations
4. ✅ Performance meets or exceeds targets
5. ✅ Error analysis identifies improvement areas

**Ready for clinical use if:**
- Drug Top-3 accuracy >80%
- SL detection TPR >70% and FPR <20%
- Essentiality correlation >0.7
- FDA approval match >70%

---

## 📝 Next Steps

1. **Approve Plan** - Review and approve this benchmark plan
2. **Create Dataset** - Start curating 100 test cases
3. **Implement Scripts** - Build benchmark infrastructure
4. **Run Benchmark** - Execute validation
5. **Analyze Results** - Generate report and insights
6. **Iterate** - Improve based on findings

---

**Status:** Ready for implementation! 🚀

