# Mechanism-Based Trial Matching: Validation Plan

**Date:** January 28, 2025  
**Status:** 📋 **VALIDATION PLAN** - Ready to execute  
**Purpose:** Validate all claimed metrics from `mechanism_trial_matching_contribution.mdc`

---

## 🎯 Validation Objectives

**Goal:** Verify all claimed metrics are accurate and reproducible.

**Claims to Validate:**
1. ✅ **0.92 avg mechanism fit** for DDR-high patients
2. ✅ **96.6% trial match accuracy**
3. ✅ **Shortlist compression**: 50+ → 5-12 trials (60-65% reduction)
4. ✅ **Time-to-first-trial**: 60-65% reduction
5. ✅ **Combined score formula**: 0.7×eligibility + 0.3×mechanism_fit

---

## 📋 Validation Methods

### **Method 1: Run Existing Validation Scripts** ✅ **READY**

**Script 1: `validate_mechanism_trial_matching.py`**

**What it validates:**
- ✅ Trial data quality (47 MoA-tagged trials)
- ✅ Mechanism vector structure (7D vectors)
- ✅ Mechanism fit computation (cosine similarity)
- ✅ Combined score formula (α=0.7, β=0.3)
- ✅ Ranking accuracy (Top-3 accuracy, MRR)
- ✅ Pathway alignment (DDR-focused trials)
- ✅ Edge cases (thresholds, all-zero vectors)
- ✅ Consistency (deterministic results)

**How to run:**
```bash
cd oncology-coPilot/oncology-backend-minimal
python scripts/validation/validate_mechanism_trial_matching.py
```

**Expected output:**
- Validation report JSON file: `trial_matching_report_YYYYMMDD_HHMMSS.json`
- Console output with task-by-task results
- Metrics: Top-3 accuracy, MRR, pathway coverage

**Previous results (from `.cursor/ayesha/SYSTEM_ANALYSIS_AND_VALIDATION_MASTER.md`):**
- ✅ Top-3 Accuracy: **1.00** (MVP target: ≥0.70) ✅
- ✅ MRR: **0.75** (MVP target: ≥0.65) ✅
- ✅ 31 DDR-focused trials found

---

**Script 2: `validate_mbd4_tp53_mechanism_capabilities.py`**

**What it validates:**
- ✅ End-to-end integration (trial matching + resistance prediction)
- ✅ MBD4+TP53 patient profile (DDR burden: 0.88)
- ✅ Mechanism fit scores for DDR-high patients
- ✅ Trial ranking accuracy

**How to run:**
```bash
cd oncology-coPilot/oncology-backend-minimal
python scripts/validation/validate_mbd4_tp53_mechanism_capabilities.py
```

**Expected output:**
- Integration report with trial matching results
- Average mechanism fit score for MBD4+TP53 patient
- Top-ranked trials with mechanism fit scores

**Previous results (from `.cursor/ayesha/SYSTEM_ANALYSIS_AND_VALIDATION_MASTER.md`):**
- ✅ **Average mechanism fit: 0.99** (excellent)
- ✅ 20 trials ranked
- ✅ Top trial: NCT04284969 (score: 0.99)

---

### **Method 2: Manual Validation with MBD4+TP53 Patient** 🔴 **HIGH PRIORITY**

**Purpose:** Verify 0.92 mechanism fit claim with real patient profile.

**Test Patient:**
```python
patient_profile = {
    "mutations": [
        {"gene": "MBD4", "hgvs_p": "p.R361*", "type": "germline"},
        {"gene": "TP53", "hgvs_p": "p.R175H", "type": "somatic"}
    ],
    "disease": "ovarian_cancer_hgsoc",
    "stage": "IVB",
    "tmb": 8.5,  # Not TMB-high
    "msi_status": "MSS"
}

# Expected mechanism vector (DDR-high):
mechanism_vector = [0.88, 0.12, 0.05, 0.02, 0.0, 0.0, 0.0]
#                    DDR   MAPK  PI3K  VEGF  HER2 IO   Efflux
```

**Test Script:**
```python
#!/usr/bin/env python3
"""
Manual Validation: Verify 0.92 Mechanism Fit for DDR-High Patients
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

from api.services.mechanism_fit_ranker import MechanismFitRanker
from api.services.trials.trial_matching_agent import TrialMatchingAgent
import json

# Load trial MoA vectors
moa_path = "../../api/resources/trial_moa_vectors.json"
with open(moa_path, "r") as f:
    trial_moa_vectors = json.load(f)

# MBD4+TP53 patient mechanism vector (DDR-high)
patient_mechanism_vector = [0.88, 0.12, 0.05, 0.02, 0.0, 0.0, 0.0]

# Initialize ranker
ranker = MechanismFitRanker(alpha=0.7, beta=0.3)

# Prepare test trials (use 47 tagged trials)
test_trials = []
for nct_id, data in trial_moa_vectors.items():
    moa_dict = data.get("moa_vector", {})
    from api.services.pathway_to_mechanism_vector import convert_moa_dict_to_vector
    moa_vector = convert_moa_dict_to_vector(moa_dict, use_7d=True)
    
    # Estimate eligibility (simplified - use 0.85 for all)
    test_trials.append({
        "nct_id": nct_id,
        "title": data.get("title", "Unknown"),
        "eligibility_score": 0.85,  # Assume patient meets criteria
        "moa_vector": moa_vector
    })

# Rank by mechanism fit
ranked_scores = ranker.rank_trials(
    trials=test_trials,
    sae_mechanism_vector=patient_mechanism_vector,
    min_eligibility=0.60,
    min_mechanism_fit=0.50
)

# Analyze results
print("=" * 60)
print("MECHANISM FIT VALIDATION: MBD4+TP53 Patient")
print("=" * 60)
print(f"Patient Mechanism Vector: {patient_mechanism_vector}")
print(f"DDR Burden: {patient_mechanism_vector[0]:.2f}")
print()

print(f"Total Trials Ranked: {len(ranked_scores)}")
print()

# Filter DDR-focused trials (DDR > 0.5 in MoA vector)
ddr_trials = [
    s for s in ranked_scores 
    if trial_moa_vectors.get(s.nct_id, {}).get("moa_vector", {}).get("ddr", 0) > 0.5
]

print(f"DDR-Focused Trials: {len(ddr_trials)}")
print()

# Calculate average mechanism fit for DDR trials
if ddr_trials:
    avg_mechanism_fit = sum(s.mechanism_fit_score for s in ddr_trials) / len(ddr_trials)
    print(f"✅ Average Mechanism Fit (DDR trials): {avg_mechanism_fit:.2f}")
    print()
    
    # Show top 5 DDR trials
    print("Top 5 DDR-Focused Trials:")
    for i, score in enumerate(ddr_trials[:5], 1):
        print(f"{i}. {score.nct_id}: {score.title[:60]}")
        print(f"   Mechanism Fit: {score.mechanism_fit_score:.2f}")
        print(f"   Combined Score: {score.combined_score:.2f}")
        print(f"   Eligibility: {score.eligibility_score:.2f}")
        print()
else:
    print("⚠️ No DDR-focused trials found")

# Verify claim: 0.92 avg mechanism fit
if avg_mechanism_fit >= 0.90:
    print(f"✅ CLAIM VERIFIED: Average mechanism fit ({avg_mechanism_fit:.2f}) ≥ 0.90")
else:
    print(f"⚠️ CLAIM NOT VERIFIED: Average mechanism fit ({avg_mechanism_fit:.2f}) < 0.90")
```

**Save as:** `scripts/validation/validate_092_mechanism_fit_claim.py`

**How to run:**
```bash
cd oncology-coPilot/oncology-backend-minimal
python scripts/validation/validate_092_mechanism_fit_claim.py
```

**Expected output:**
- Average mechanism fit for DDR-focused trials
- Top 5 DDR trials with scores
- Verification of 0.92 claim

---

### **Method 3: Shortlist Compression Validation** 🟡 **MEDIUM PRIORITY**

**Purpose:** Verify 50+ → 5-12 trials compression claim.

**Test Script:**
```python
#!/usr/bin/env python3
"""
Shortlist Compression Validation: Verify 50+ → 5-12 trials compression
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

from api.services.trials.trial_matching_agent import TrialMatchingAgent
from api.services.autonomous_trial_agent import AutonomousTrialAgent
import asyncio

# MBD4+TP53 patient profile
patient_profile = {
    "mutations": [
        {"gene": "MBD4", "hgvs_p": "p.R361*", "type": "germline"},
        {"gene": "TP53", "hgvs_p": "p.R175H", "type": "somatic"}
    ],
    "disease": "ovarian_cancer_hgsoc",
    "stage": "IVB"
}

mechanism_vector = [0.88, 0.12, 0.05, 0.02, 0.0, 0.0, 0.0]

async def test_compression():
    # Step 1: Get all trials (generic search - no mechanism fit)
    agent = AutonomousTrialAgent()
    queries = await agent.generate_search_queries(patient_profile)
    
    # Use hybrid search to get ~50 trials
    from api.services.hybrid_trial_search import HybridTrialSearchService
    search_service = HybridTrialSearchService()
    
    all_trials = []
    for query in queries:
        results = await search_service.search_optimized(
            query=query,
            patient_context=patient_profile,
            top_k=20
        )
        all_trials.extend(results.get("found_trials", []))
    
    # Deduplicate
    seen = set()
    unique_trials = []
    for trial in all_trials:
        nct_id = trial.get("nct_id")
        if nct_id and nct_id not in seen:
            seen.add(nct_id)
            unique_trials.append(trial)
    
    print(f"Total trials (generic search): {len(unique_trials)}")
    
    # Step 2: Apply mechanism fit ranking
    matching_agent = TrialMatchingAgent()
    result = await matching_agent.match(
        patient_profile=patient_profile,
        biomarker_profile={},
        mechanism_vector=mechanism_vector,
        max_results=12
    )
    
    mechanism_aligned_trials = len(result.matches)
    print(f"Mechanism-aligned trials: {mechanism_aligned_trials}")
    
    # Calculate compression
    compression_ratio = mechanism_aligned_trials / len(unique_trials) if unique_trials else 0
    reduction_percent = (1 - compression_ratio) * 100
    
    print()
    print("=" * 60)
    print("SHORTLIST COMPRESSION VALIDATION")
    print("=" * 60)
    print(f"Generic Search: {len(unique_trials)} trials")
    print(f"Mechanism-Aligned: {mechanism_aligned_trials} trials")
    print(f"Compression Ratio: {compression_ratio:.2%}")
    print(f"Reduction: {reduction_percent:.1f}%")
    print()
    
    # Verify claim: 50+ → 5-12 trials (60-65% reduction)
    if len(unique_trials) >= 50 and 5 <= mechanism_aligned_trials <= 12:
        if reduction_percent >= 60:
            print(f"✅ CLAIM VERIFIED: {len(unique_trials)} → {mechanism_aligned_trials} trials ({reduction_percent:.1f}% reduction)")
        else:
            print(f"⚠️ PARTIAL: Compression works but reduction ({reduction_percent:.1f}%) < 60%")
    else:
        print(f"⚠️ CLAIM NOT VERIFIED: {len(unique_trials)} → {mechanism_aligned_trials} trials")

if __name__ == "__main__":
    asyncio.run(test_compression())
```

**Save as:** `scripts/validation/validate_shortlist_compression.py`

**How to run:**
```bash
cd oncology-coPilot/oncology-backend-minimal
python scripts/validation/validate_shortlist_compression.py
```

**Expected output:**
- Total trials from generic search
- Mechanism-aligned trials count
- Compression ratio and reduction percentage
- Verification of 60-65% reduction claim

---

### **Method 4: Accuracy Validation** 🟡 **MEDIUM PRIORITY**

**Purpose:** Verify 96.6% trial match accuracy claim.

**Challenge:** Need ground truth (expert-validated trial matches) to calculate accuracy.

**Approach:**
1. **Use existing validation results** (Top-3 accuracy: 1.00, MRR: 0.75)
2. **Calculate accuracy from validation script output**
3. **Document methodology** (how accuracy is measured)

**Test Script:**
```python
#!/usr/bin/env python3
"""
Accuracy Validation: Verify 96.6% trial match accuracy
"""
import sys
import os
import json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

# Load previous validation report
report_path = "scripts/validation/trial_matching_report_*.json"  # Most recent

# Parse validation report
with open(report_path, "r") as f:
    report = json.load(f)

# Extract accuracy metrics
top3_accuracy = report.get("metrics", {}).get("top3_accuracy", 0)
mrr = report.get("metrics", {}).get("mrr", 0)

# Calculate overall accuracy (weighted average)
# Top-3 accuracy: 70% weight
# MRR: 30% weight
overall_accuracy = (top3_accuracy * 0.7) + (mrr * 0.3)

print("=" * 60)
print("ACCURACY VALIDATION")
print("=" * 60)
print(f"Top-3 Accuracy: {top3_accuracy:.2%}")
print(f"MRR: {mrr:.2%}")
print(f"Overall Accuracy (weighted): {overall_accuracy:.2%}")
print()

# Verify claim: 96.6% accuracy
if overall_accuracy >= 0.966:
    print(f"✅ CLAIM VERIFIED: Overall accuracy ({overall_accuracy:.2%}) ≥ 96.6%")
else:
    print(f"⚠️ CLAIM NOT VERIFIED: Overall accuracy ({overall_accuracy:.2%}) < 96.6%")
    print(f"   Note: Previous validation showed Top-3: 1.00, MRR: 0.75")
    print(f"   Weighted: {(1.00 * 0.7) + (0.75 * 0.3):.2%} = 92.5%")
```

**Note:** 96.6% claim may need clarification - is it Top-3 accuracy, MRR, or weighted average?

---

## 📊 Validation Checklist

### **Phase 1: Core Functionality** ✅ **READY**

- [ ] Run `validate_mechanism_trial_matching.py`
  - [ ] Verify all 8 tasks pass
  - [ ] Document Top-3 accuracy and MRR
  - [ ] Verify MoA coverage (47 trials)

- [ ] Run `validate_mbd4_tp53_mechanism_capabilities.py`
  - [ ] Verify average mechanism fit ≥ 0.90
  - [ ] Document top-ranked trials
  - [ ] Verify DDR-focused trials rank first

### **Phase 2: Claim Verification** 🔴 **HIGH PRIORITY**

- [ ] Verify 0.92 mechanism fit claim
  - [ ] Run `validate_092_mechanism_fit_claim.py`
  - [ ] Calculate average mechanism fit for DDR trials
  - [ ] Document results

- [ ] Verify shortlist compression
  - [ ] Run `validate_shortlist_compression.py`
  - [ ] Measure compression ratio
  - [ ] Verify 60-65% reduction

- [ ] Verify accuracy claim
  - [ ] Calculate overall accuracy from validation results
  - [ ] Document methodology
  - [ ] Clarify if 96.6% is Top-3, MRR, or weighted

### **Phase 3: Documentation** 📝

- [ ] Create validation report
  - [ ] Document all metrics
  - [ ] Include test results
  - [ ] Compare claims vs reality

- [ ] Update contribution document
  - [ ] Add validation results
  - [ ] Update metrics if needed
  - [ ] Add methodology section

---

## 🎯 Expected Results

### **From Previous Validations:**

**Mechanism Fit Ranking:**
- ✅ Top-3 Accuracy: **1.00** (MVP target: ≥0.70)
- ✅ MRR: **0.75** (MVP target: ≥0.65)
- ✅ 31 DDR-focused trials found

**MBD4+TP53 Integration:**
- ✅ Average mechanism fit: **0.99** (excellent)
- ✅ 20 trials ranked
- ✅ Top trial: NCT04284969 (score: 0.99)

### **What We Need to Verify:**

1. **0.92 mechanism fit claim:**
   - Expected: Average mechanism fit for DDR trials ≥ 0.90
   - Previous result: 0.99 (exceeds claim)

2. **Shortlist compression:**
   - Expected: 50+ → 5-12 trials (60-65% reduction)
   - Need to test with real search

3. **Accuracy:**
   - Previous: Top-3 = 1.00, MRR = 0.75
   - Weighted: (1.00 × 0.7) + (0.75 × 0.3) = **0.925 (92.5%)**
   - Claim: 96.6% (need to verify calculation method)

---

## 🚀 Quick Start

**Run all validations:**
```bash
cd oncology-coPilot/oncology-backend-minimal

# 1. Core functionality
python scripts/validation/validate_mechanism_trial_matching.py

# 2. MBD4+TP53 integration
python scripts/validation/validate_mbd4_tp53_mechanism_capabilities.py

# 3. 0.92 mechanism fit claim (create script first)
python scripts/validation/validate_092_mechanism_fit_claim.py

# 4. Shortlist compression (create script first)
python scripts/validation/validate_shortlist_compression.py
```

**Expected time:** 30-60 minutes total

---

## 📝 Validation Report Template

After running validations, create a report:

```markdown
# Mechanism-Based Trial Matching: Validation Report

**Date:** [Date]
**Validated By:** [Name]
**Status:** ✅ VERIFIED / ⚠️ PARTIAL / ❌ NOT VERIFIED

## Metrics Verified

| Metric | Claimed | Actual | Status |
|--------|---------|--------|--------|
| Mechanism Fit (DDR-high) | 0.92 avg | [Result] | ✅/⚠️/❌ |
| Shortlist Compression | 50+ → 5-12 | [Result] | ✅/⚠️/❌ |
| Time Reduction | 60-65% | [Result] | ✅/⚠️/❌ |
| Accuracy | 96.6% | [Result] | ✅/⚠️/❌ |

## Test Results

[Include test outputs, screenshots, JSON reports]

## Conclusion

[Summary of validation results]
```

---

*Validation Plan Created: January 28, 2025*  
*Status: 📋 READY TO EXECUTE*


