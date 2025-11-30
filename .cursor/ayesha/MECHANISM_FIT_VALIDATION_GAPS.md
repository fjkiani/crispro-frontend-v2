# ⚔️ MECHANISM FIT VALIDATION - GAPS & QUESTIONS FOR MANAGER

**Date:** January 14, 2025  
**Status:** 🔍 **GAPS IDENTIFIED - NEED MANAGER CLARIFICATION**  
**Purpose:** Document verified gaps and questions before building validation script

---


Questions for manager:
Q1: Data type conversion — add conversion code, fix ranker, or fix schema?
Q2: Test patient construction — call API, call service directly, or manually construct vectors?
Q3: Expected rankings — which specific NCT IDs should be top-3 for BRCA1/KRAS/HER2?
Q4: Eligibility score — is match_score correct, or use confidence gate?
Q5: Success criteria — what metrics define "success"?
Created a gaps document: .cursor/ayesha/MECHANISM_FIT_VALIDATION_GAPS.md with details.

## ✅ Manager‑Ready Answers (Verified)

- Q1: Data type conversion
  - Convert Dict→List in the router before calling the ranker; keep the ranker typed to List[float].
  - Pathway order (both SAE and MoA): ["ddr","mapk","pi3k","vegf","her2","io","efflux"].
  - Convert both:
    - SAE: `[request.sae_mechanism_vector.get(k, 0.0) for k in order]`
    - MoA: `[trial["moa_vector"].get(k, 0.0) for k in order]`
  - Keep Pydantic schema as Dict[str,float] (self‑describing payloads); do conversion at call site.

- Q2: Test patient construction
  - Use Option B for the harness: call `SAEFeatureService.compute_sae_features()` directly with controlled inputs to produce the 7D vector; then run ranking offline.
  - After harness passes, add two E2E checks via `/api/ayesha/complete_care_v2` (BRCA1 and KRAS scenarios) to verify wiring.

- Q3: Expected rankings
  - Validate by mechanism class, not fixed NCT IDs:
    - BRCA1/HRD+ (DDR high): DDR/PARP trials in top‑3; mechanism_fit ≥ 0.75; MAPK/HER2 monotherapy not in top‑3 unless mixed MoA.
    - KRAS G12D with MAPK ≥ 0.40: MEK/RAF trials in top‑3; PARP not in top‑3; mechanism_fit ≥ 0.70.
    - HER2 amplified with HER2 burden ≥ 0.70: HER2 trials in top‑3; mechanism_fit ≥ 0.80.
    - ABCB1 high (efflux): taxane trials demoted (combined score drops vs eligibility‑only).
  - Use `api/resources/trial_moa_vectors.json` to determine whether a trial matches the expected mechanism (component ≥ 0.70 for that pathway).

- Q4: Eligibility score source
  - Use existing `match_score` (post hard filters + soft boosts) as `eligibility_score` for the ranker.
  - Keep thresholds: eligibility ≥ 0.60, mechanism_fit ≥ 0.50.
  - Do not substitute the eligibility checklist “confidence gate” for the ranker input.

- Q5: Success criteria for validation
  - Top‑3 accuracy by expected mechanism ≥ 80% across 20 scenarios.
  - Mean Reciprocal Rank (MRR) ≥ 0.75.
  - Mechanism_fit for expected mechanism ≥ 0.75; non‑expected pathways ≤ 0.50 (per scenario).
  - Combined score lift vs eligibility‑only ≥ 0.10 for matching‑mechanism trials.
  - Sanity checks:
    - BRCA1/DDR: DDR trials present in top‑3; PARP appears; MAPK‑only trials not top‑3.
    - KRAS/MAPK: MEK/RAF present in top‑3; PARP not in top‑3.
    - HER2: HER2 trials present in top‑3.

## ✅ **CODE VERIFICATION (FROM ACTUAL CODEBASE)**

### **Q1: Data Type Mismatch - VERIFIED BUG** ✅

**Evidence from code:**
- **`ayesha_trials.py` line 600:** Passes `request.sae_mechanism_vector` (Dict) directly to ranker
- **`mechanism_fit_ranker.py` line 257:** Function signature expects `sae_mechanism_vector: List[float]`
- **`ayesha_trials.py` line 594:** Stores `moa_vector` as Dict: `{"ddr": 0.95, "mapk": 0.0, ...}`
- **`mechanism_fit_ranker.py` line 111:** Expects `trial_moa_vector = trial.get("moa_vector", [])` (List)

**Conclusion:** ✅ **CONFIRMED BUG** - Dict→List conversion needed before ranker call

**Manager's Answer Verified:** ✅ Correct - convert in router before calling ranker

---

### **Q2: Test Patient Construction - VERIFIED APPROACH** ✅

**Evidence from code:**
- **`test_sae_phase2_services.py` line 109-112:** Shows calling `compute_sae_features()` directly with mock data:
  ```python
  result = compute_sae_features(
      insights_bundle=insights_bundle,
      pathway_scores=pathway_scores,
      tumor_context=tumor_context
  )
  ```
- **`test_sae_phase2_services.py` line 253:** Shows constructing SAE vectors as Lists for testing

**Conclusion:** ✅ **VERIFIED** - Previous agents/testers call `compute_sae_features()` directly

**Manager's Answer Verified:** ✅ Correct - Option B (call service directly) is the established pattern

---

### **Q3: Expected Rankings - PARTIALLY VERIFIED** ⚠️

**Evidence from code:**
- **`test_sae_phase2_services.py` line 250-267:** Shows HER2 test case with NCT06819007
  - Validates: `result[0]["mechanism_alignment"]["HER2"] > 0.5`
  - Does NOT validate specific top-3 rankings
- **`trial_moa_vectors.json`:** Contains 47 trials with MoA vectors, but no "expected ranking" metadata

**Conclusion:** ⚠️ **PARTIALLY VERIFIED** - Tests validate mechanism alignment, not specific top-3 rankings

**Manager's Answer Verified:** ✅ Reasonable - Validate by mechanism class (not fixed NCT IDs) aligns with test patterns

---

### **Q4: Eligibility Score Source - VERIFIED** ✅

**Evidence from code:**
- **`ayesha_trials.py` line 593:** `"eligibility_score": trial.get("match_score", 0.7)`
- **`mechanism_fit_ranker.py` line 103:** `eligibility_score = trial.get("eligibility_score", 0.0)`
- **`ayesha_trials.py` line 225:** `match_score` computed from soft boosts (capped at 1.0)

**Conclusion:** ✅ **VERIFIED** - `match_score` is used as `eligibility_score`

**Manager's Answer Verified:** ✅ Correct - Use `match_score` as eligibility_score

---

### **Q5: Success Criteria - NOT VERIFIED IN CODE** ⚠️

**Evidence from code:**
- **`test_sae_phase2_services.py` line 234-244:** Shows min threshold tests (eligibility ≥0.60, mechanism_fit ≥0.50)
- **`test_sae_phase2_services.py` line 250-267:** Shows mechanism alignment validation (>0.5)
- **No test files found** with top-3 accuracy, MRR, or combined score lift metrics

**Conclusion:** ⚠️ **NOT VERIFIED** - No existing validation tests with these metrics

**Manager's Answer Verified:** ✅ Reasonable - Success criteria align with standard ranking validation metrics

---

## ✅ **WHAT I VERIFIED (75-80% Confidence)**

### **1. Mechanism Fit Ranker Code:**
- ✅ **Location:** `api/services/mechanism_fit_ranker.py`
- ✅ **Formula:** `(0.7 × eligibility) + (0.3 × mechanism_fit)` (line 129)
- ✅ **Process:** L2-normalize vectors → cosine similarity → combine scores (lines 118-129)
- ✅ **Thresholds:** eligibility ≥0.60, mechanism_fit ≥0.50 (lines 106, 125)
- ✅ **Input:** `sae_mechanism_vector: List[float]` (7D: DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux) (line 257)

### **2. Trial Data:**
- ✅ **47 trials** with MoA vectors in `trial_moa_vectors.json` (verified via terminal command)
- ✅ **Structure:** `{nct_id: {moa_vector: {ddr: 0.95, mapk: 0.0, ...}, confidence, ...}}`
- ✅ **1,000 trials** in SQLite database
- ✅ **MoA vectors attached** to trials in `ayesha_trials.py` (line 568-594)

### **3. SAE Mechanism Vector:**
- ✅ **Source:** `sae_feature_service.py` → `mechanism_vector: List[float]` (line 203-212)
- ✅ **Order:** `[DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]` (line 203-212)
- ✅ **Computed from:** Pathway burdens + IO eligibility + cross-resistance risk

### **4. Eligibility Scoring:**
- ✅ **Hard filters:** Stage IV, first-line, recruiting, NYC metro (line 96-157)
- ✅ **Soft boosts:** Frontline (+0.25), stage IV specific (+0.15), carbo/pacli (+0.20), bevacizumab (+0.15), phase III (+0.10), multi-center (+0.05)
- ✅ **Final score:** `match_score` (capped at 1.0) used as `eligibility_score` (line 593)

---

## 🚨 **CRITICAL GAPS IDENTIFIED**

### **Gap #1: Data Type Mismatch (CRITICAL BUG)**

**Issue:**
- **Schema:** `request.sae_mechanism_vector: Optional[Dict[str, float]]` (line 78)
- **Ranker expects:** `sae_mechanism_vector: List[float]` (line 257)
- **Code passes:** `request.sae_mechanism_vector` (dict) directly to ranker (line 600)
- **Result:** ❌ **TYPE ERROR** - ranker will fail when called

**Also:**
- **MoA vector in file:** `{ddr: 0.95, mapk: 0.0, ...}` (dict)
- **Ranker expects:** `[0.95, 0.0, ...]` (list)
- **Code passes:** `trial.get("moa_vector")` (dict) directly (line 594)
- **Result:** ❌ **TYPE ERROR** - ranker will fail when called

**Question for Manager:**
> **Q1:** Is there conversion code I'm missing, or is this a bug that needs fixing?
> 
> **Expected conversion:**
> ```python
> # SAE vector: Dict → List
> pathway_order = ["ddr", "mapk", "pi3k", "vegf", "her2", "io", "efflux"]
> sae_list = [request.sae_mechanism_vector.get(pathway, 0.0) for pathway in pathway_order]
> 
> # MoA vector: Dict → List
> moa_list = [trial["moa_vector"].get(pathway, 0.0) for pathway in pathway_order]
> ```
> 
> **Should I:**
> - A) Add conversion code in `ayesha_trials.py` before calling ranker?
> - B) Fix the ranker to accept dicts?
> - C) Fix the schema to use List instead of Dict?

---

### **Gap #2: How to Get SAE Mechanism Vector for Test Patients**

**Issue:**
- I know SAE vector comes from `sae_feature_service.compute_sae_features()`
- But I need to know: **How do I construct test patient profiles?**

**Question for Manager:**
> **Q2:** For validation, should I:
> 
> **Option A:** Call `/api/ayesha/complete_care_v2` with test patient data?
> - Pros: Uses real API, tests end-to-end
> - Cons: Requires full patient context (mutations, tumor context, etc.)
> 
> **Option B:** Call `SAEFeatureService.compute_sae_features()` directly?
> - Pros: Faster, more control
> - Cons: Bypasses API layer
> 
> **Option C:** Manually construct SAE vectors for test scenarios?
> - Pros: Simplest, fastest
> - Cons: Doesn't test SAE computation
> 
> **Which approach do you prefer?**

---

### **Gap #3: Expected Rankings for Test Scenarios**

**Issue:**
- I know mechanism fit ranking should prioritize PARP for BRCA1, MEK/RAF for KRAS
- But I need: **Which specific trials should be in top-3 for each profile?**

**Question for Manager:**
> **Q3:** For validation test cases, what are the expected top-3 trials?
> 
> **Test Case 1: BRCA1 Patient (HRD+, DDR high)**
> - Expected: PARP trials (DDR mechanism) in top-3
> - Which specific NCT IDs should I expect? (e.g., NCT04284969, NCT04001023?)
> 
> **Test Case 2: KRAS G12D Patient (MAPK high)**
> - Expected: MEK/RAF trials (MAPK mechanism) in top-3
> - Which specific NCT IDs should I expect?
> 
> **Test Case 3: HER2 Amplified Patient**
> - Expected: HER2 trials in top-3
> - Which specific NCT IDs should I expect? (e.g., NCT06331130?)
> 
> **Or should I:**
> - Just validate that mechanism fit scores are higher for matching mechanisms?
> - Not worry about specific NCT IDs, just validate ranking logic?

---

### **Gap #4: Eligibility Score Source**

**Issue:**
- I see `eligibility_score = trial.get("match_score", 0.7)` (line 593)
- But `match_score` comes from soft boosts (line 225)
- **Question:** Is `match_score` the right source for eligibility, or should it be computed differently?

**Question for Manager:**
> **Q4:** For mechanism fit ranking, what should `eligibility_score` be?
> 
> **Current code uses:** `trial.get("match_score", 0.7)` (soft boost score)
> 
> **Is this correct, or should it be:**
> - A) Hard filter pass/fail (1.0 if passes, 0.0 if fails)?
> - B) Confidence gate from eligibility checklist (0.75-0.90)?
> - C) Something else?

---

### **Gap #5: Validation Success Criteria**

**Issue:**
- I know we want to validate mechanism fit ranking works
- But I need: **What metrics define "success"?**

**Question for Manager:**
> **Q5:** What are the success criteria for mechanism fit validation?
> 
> **Option A: Ranking Accuracy**
> - Top-ranked trial has mechanism fit >0.85?
> - Top-3 trials all have mechanism fit >0.70?
> - Mechanism fit correlates with expected mechanism (r ≥0.60)?
> 
> **Option B: Clinical Alignment**
> - BRCA1 patients → PARP trials ranked #1-3?
> - KRAS patients → MEK/RAF trials ranked #1-3?
> - HER2 patients → HER2 trials ranked #1-3?
> 
> **Option C: Score Improvement**
> - Mechanism fit ranking improves combined scores vs eligibility-only?
> - Top-ranked trial has combined score >0.80?
> 
> **What metrics should I report?**

---

## 📋 **WHAT I CAN BUILD NOW (WITH ASSUMPTIONS)**

### **Validation Script Structure (Pending Manager Answers):**

```python
# Test scenarios
test_cases = [
    {
        "name": "BRCA1 Patient (HRD+, DDR high)",
        "mutations": [{"gene": "BRCA1", "hgvs_p": "Q1756*"}],
        "tumor_context": {"hrd_score": 65, "tmb": 5, "msi_status": "MSS"},
        "expected_mechanism": "DDR",
        "expected_top_trials": ["NCT04284969", "NCT04001023", ...]  # NEED MANAGER INPUT
    },
    {
        "name": "KRAS G12D Patient (MAPK high)",
        "mutations": [{"gene": "KRAS", "hgvs_p": "p.G12D"}],
        "tumor_context": {"hrd_score": 25, "tmb": 8, "msi_status": "MSS"},
        "expected_mechanism": "MAPK",
        "expected_top_trials": [...]  # NEED MANAGER INPUT
    },
    # ... more test cases
]

# For each test case:
# 1. Compute SAE features → get mechanism_vector
# 2. Call /api/ayesha/trials/search with mechanism_vector
# 3. Validate rankings match expectations
# 4. Report metrics
```

**Blockers:**
- ❌ Data type conversion (Gap #1) - **CRITICAL**
- ❌ Test patient construction (Gap #2)
- ❌ Expected rankings (Gap #3)
- ❌ Success criteria (Gap #5)

---

## 🎯 **RECOMMENDED NEXT STEPS**

### **Before Building Validation Script:**

1. **Fix Gap #1 (CRITICAL):** Resolve data type mismatch
   - Either add conversion code or fix ranker/schema
   - **This will break if not fixed**

2. **Get Manager Answers:**
   - Q1: Data type conversion approach
   - Q2: Test patient construction method
   - Q3: Expected rankings for test cases
   - Q4: Eligibility score source
   - Q5: Success criteria

3. **Then Build Validation Script:**
   - Create test scenarios
   - Call APIs with proper data types
   - Validate rankings
   - Report metrics

---

## ✅ **FINAL VERIFICATION & REMAINING QUESTIONS**

### **✅ Fully Clear (100% Confidence):**

1. **Q1 (Data Type Conversion)**: ✅ **VERIFIED BUG**
   - **Issue**: `ayesha_trials.py` passes `Dict[str, float]` but ranker expects `List[float]`
   - **Solution**: Convert Dict→List in router before calling ranker
   - **Pathway Order**: `["ddr","mapk","pi3k","vegf","her2","io","efflux"]` (verified in `sae_feature_service.py` lines 204-212)
   - **Edge Cases**: Missing pathways → `.get(k, 0.0)` handles it; missing `moa_vector` → defaults to zero vector (line 594)

2. **Q2 (Test Patient Construction)**: ✅ **VERIFIED APPROACH**
   - **Solution**: Call `SAEFeatureService.compute_sae_features()` directly (Option B)
   - **Pattern Verified**: `test_sae_phase2_services.py` shows this pattern (lines 109-112, 253)
   - **Inputs Needed**: `insights_bundle`, `pathway_scores`, `tumor_context`, `treatment_history`, `ca125_intelligence`

3. **Q4 (Eligibility Score)**: ✅ **VERIFIED**
   - **Source**: `match_score` (post hard filters + soft boosts) used as `eligibility_score`
   - **Verified**: `ayesha_trials.py` line 593 uses `trial.get("match_score", 0.7)`

4. **Missing MoA Vector Handling**: ✅ **VERIFIED**
   - **Ranker Logic**: Lines 114-116 default to zero vector if missing/invalid
   - **Trial Default**: Line 594 defaults to all zeros if `moa_vector` missing
   - **Result**: No crashes, but mechanism fit will be 0.0 for untagged trials

---

### **⚠️ Partially Clear (80% Confidence - Need Minor Clarification):**

1. **Q3 (Expected Rankings)**: ⚠️ **STRATEGIC GOAL, NOT IMPLEMENTED**
   - **Manager's Answer**: Validate by mechanism class (not fixed NCT IDs)
   - **Verified**: Tests validate mechanism alignment but not top-3 rankings
   - **Question**: Should I build this validation metric from scratch, or is there existing code I'm missing?

2. **Q5 (Success Criteria)**: ⚠️ **STRATEGIC GOAL, NOT IMPLEMENTED**
   - **Manager's Answer**: Top-3 accuracy ≥80%, MRR ≥0.75, mechanism_fit ≥0.75, combined score lift ≥0.10
   - **Verified**: No existing tests with these metrics
   - **Question**: Should I implement these metrics from scratch, or is there existing validation infrastructure I should use?

---

### **❓ Remaining Questions (Need Manager Input):**

#### **Q6: Pathway Scores Source for Validation** ✅ **RESOLVED FROM CODEBASE PATTERNS**

**Answer (Based on `test_sae_phase2_services.py` Pattern):**

**✅ Use Option C: Manually Construct Realistic Test Values**

**Evidence from Codebase:**
- **`test_sae_phase2_services.py` lines 106, 123, 138, 165, 442**: All tests use manually constructed pathway scores
- **Example pattern** (line 442-448):
  ```python
  pathway_scores = {
      "ddr": 0.90,  # High DDR burden (BRCA1 biallelic)
      "mapk": 0.20,
      "pi3k": 0.15,
      "vegf": 0.40,
      "her2": 0.0  # Unknown HER2 status
  }
  ```

**Rationale:**
- ✅ **Fast**: No Evo2 service dependency
- ✅ **Controlled**: Precise test values for each scenario
- ✅ **Established Pattern**: Matches existing test suite
- ✅ **Focus**: Tests mechanism fit ranking logic, not pathway computation

**Implementation for Validation:**
```python
# Test Case 1: BRCA1 Patient (HRD+, DDR high)
pathway_scores = {
    "ddr": 0.90,  # High DDR burden
    "mapk": 0.10,
    "pi3k": 0.15,
    "vegf": 0.30,
    "her2": 0.0
}

# Test Case 2: KRAS G12D Patient (MAPK high)
pathway_scores = {
    "ddr": 0.20,
    "mapk": 0.85,  # High MAPK burden
    "pi3k": 0.25,
    "vegf": 0.40,
    "her2": 0.0
}

# Test Case 3: HER2 Amplified Patient
pathway_scores = {
    "ddr": 0.30,
    "mapk": 0.20,
    "pi3k": 0.15,
    "vegf": 0.50,
    "her2": 0.80  # High HER2 burden
}
```

**Note on Gene-to-Pathway Mapping:**
- For validation, use realistic pathway scores directly (don't need to map genes)
- The `drug_mapping.py` mapping (BRCA1 → tp53) is for efficacy orchestrator, not SAE
- SAE uses pathway scores directly from `pathway_scores` dict (line 181-185 in `sae_feature_service.py`)

---

#### **Q7: Edge Cases for Validation** ✅ **RESOLVED FROM CODEBASE PATTERNS**

**Answer (Based on `test_sae_phase2_services.py` Pattern):**

**✅ Test These Edge Cases (Following Existing Test Suite):**

**1. Minimum Thresholds** (Already tested in `test_minimum_thresholds`, line 233-248):
   - Eligibility <0.60 → filtered out
   - Mechanism fit <0.50 → filtered out
   - Both thresholds must pass

**2. L2 Normalization** (Already tested in `test_l2_normalization`, line 183-196):
   - Zero vector handling
   - Normalization correctness

**3. Cosine Similarity** (Already tested in `test_cosine_similarity`, line 198-214):
   - Identical vectors → cosine = 1.0
   - Orthogonal vectors → cosine = 0.0

**4. Pathway Alignment Breakdown** (Already tested in `test_pathway_alignment_breakdown`, line 269-285):
   - Per-pathway alignment computation
   - Product of SAE burden × trial MoA

**5. HER2 Pathway Integration** (Already tested in `test_her2_pathway_integration`, line 250-267):
   - 7D vector with HER2 component
   - HER2-targeted trial matching

**Additional Edge Cases for Validation Script:**
1. **Empty SAE vector (all zeros)**: Verify ranking falls back to eligibility-only
2. **Missing MoA vectors**: Verify untagged trials get mechanism_fit = 0.0, but still ranked if eligibility ≥0.60
3. **All trials below thresholds**: Verify empty result set (no crashes)
4. **Single trial above thresholds**: Verify it ranks #1
5. **Pathway mismatch**: Verify non-matching mechanisms get low scores (e.g., BRCA1 patient vs MAPK-only trial)

**Note**: The ranker already handles most edge cases gracefully (see code lines 106-127, 114-116, 189). Validation should verify these behaviors work as expected.

---

#### **Q8: Validation Output Format** ✅ **RESOLVED FROM CODEBASE PATTERNS**

**Answer (Based on `validate_sae_tcga.py` Pattern):**

**✅ Use JSON-Only Output (Following `validate_sae_tcga.py` Pattern)**

**Evidence from Codebase:**
- **`validate_sae_tcga.py` lines 431-467**: Uses JSON-only output (`validation_results.json`)
- **Structure**:
  ```python
  {
      "validation_date": "2025-01-14",
      "data_source": "...",
      "backend_url": "...",
      "metrics": {
          "auroc": 0.85,
          "auprc": 0.78,
          "sensitivity": 0.90,
          "specificity": 0.75,
          ...
      },
      "results": [
          {
              "patient_id": "...",
              "scenario": "BRCA1",
              "top_3_trials": [...],
              "metrics": {...}
          }
      ]
  }
  ```

**Console Output:**
- Print metrics to console with interpretation (lines 504-520)
- Human-readable summary with pass/fail status
- Interpretation guidance (e.g., "EXCELLENT: Model shows strong discrimination ability!")

**Implementation:**
```python
# Save JSON
output_path = OUTPUT_DIR / "mechanism_fit_validation_results.json"
with open(output_path, 'w') as f:
    json.dump(output, f, indent=2)

# Console output
print("\n⚔️ VALIDATION COMPLETE ⚔️\n")
print(f"📊 Top-3 Accuracy: {metrics['top3_accuracy']:.3f}")
print(f"📊 MRR: {metrics['mrr']:.3f}")
print(f"📊 Mechanism Fit (expected): {metrics['mechanism_fit_expected']:.3f}")

# Interpretation
if metrics['top3_accuracy'] >= 0.80:
    print("\n✅ EXCELLENT: Mechanism fit ranking shows strong alignment!")
elif metrics['top3_accuracy'] >= 0.70:
    print("\n✅ GOOD: Mechanism fit ranking shows acceptable alignment!")
else:
    print("\n⚠️ MODERATE: Mechanism fit ranking needs improvement.")
```

**Note**: No Markdown report needed - JSON + console output is the established pattern.

---

## 📊 **CURRENT CONFIDENCE LEVEL**

**Overall:** 85-90% (up from 75-80% after final verification)

**What I'm Confident About:**
- ✅ Mechanism fit ranker algorithm (cosine similarity, L2-normalization) - **100%**
- ✅ SAE mechanism vector computation (pathway burdens → 7D vector) - **100%**
- ✅ Trial data structure (47 MoA-tagged trials) - **100%**
- ✅ Eligibility scoring (hard filters + soft boosts) - **100%**
- ✅ Data type conversion approach (Dict→List in router) - **100%**
- ✅ Test patient construction (call service directly) - **100%**
- ✅ Missing MoA vector handling (defaults to zero vector) - **100%**

**What I Need Minor Clarification On:**
- ✅ Pathway scores source for validation (Q6) - **100% RESOLVED** (use manual construction, following `test_sae_phase2_services.py` pattern)
- ✅ Edge cases to test (Q7) - **100% RESOLVED** (follow existing test suite patterns)
- ✅ Validation output format (Q8) - **100% RESOLVED** (JSON-only, following `validate_sae_tcga.py` pattern)

**What Needs Implementation (Not Questions):**
- 🔨 Build top-3 accuracy metric (Q3) - **Clear what to build, just need to implement**
- 🔨 Build success criteria metrics (Q5) - **Clear what to build, just need to implement**

---

**Status:** ⚔️ **READY TO BUILD - ALL QUESTIONS RESOLVED** ⚔️

**Next Action:** 
1. ✅ **Q6-Q8 RESOLVED** - Answers found in codebase patterns (`test_sae_phase2_services.py`, `validate_sae_tcga.py`)
2. Fix Q1 bug (Dict→List conversion) - **CRITICAL** - **READY TO IMPLEMENT**
3. Build validation script with all answers (Q2-Q8) - **READY TO IMPLEMENT**
4. Implement Q3/Q5 metrics (top-3 accuracy, MRR, etc.) - **READY TO IMPLEMENT**

**Key Insights Added:**
- ✅ **Pathway Score Computation**: Understood full flow from Evo2 → aggregation → SAE features
- ✅ **Q6 Answer**: Use manual pathway score construction (following `test_sae_phase2_services.py` pattern)
- ✅ **Q7 Answer**: Follow existing test suite edge case patterns (thresholds, normalization, alignment)
- ✅ **Q8 Answer**: Use JSON-only output (following `validate_sae_tcga.py` pattern)
- ✅ **All Questions Resolved**: Found answers in codebase patterns, no manager input needed

---

**Last Updated:** January 14, 2025 (All Questions Resolved from Codebase Patterns)  
**Owner:** Zo  
**Reference:** This documents all gaps before building mechanism fit validation script

**Recent Updates:**
- ✅ **Q6 RESOLVED**: Use manual pathway score construction (following `test_sae_phase2_services.py` pattern)
- ✅ **Q7 RESOLVED**: Follow existing test suite edge case patterns (from `test_sae_phase2_services.py`)
- ✅ **Q8 RESOLVED**: Use JSON-only output (following `validate_sae_tcga.py` pattern)
- ✅ **All Questions Resolved**: Found definitive answers in codebase patterns - ready to build validation script

