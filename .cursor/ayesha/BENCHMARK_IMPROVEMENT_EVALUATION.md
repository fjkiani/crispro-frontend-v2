# Benchmark Improvement Evaluation: Gaps & Recommendations

**Date**: January 27, 2025  
**Purpose**: Evaluate benchmark findings against MBD4 plan capabilities and identify improvement opportunities  
**Status**: Ready for implementation

---

## Executive Summary

**Benchmark Findings**: Drug ranking is excellent (100% Top-5), but **outcome correlations are weak** (r=0.037-0.278) and **critical features (TMB/HRD/MSI) are not being utilized** in predictions.

**MBD4 Plan Status**: System **supports** TMB/HRD/MSI via tumor context and sporadic gates, but **benchmark scripts are not passing these values**.

**Key Gap**: **Infrastructure exists, but benchmark integration is missing**. The system can use TMB/HRD/MSI, but benchmarks aren't providing them.

---

## Critical Findings: What Needs Improvement

### 1. ⚠️ **CRITICAL**: TMB/HRD/MSI Not Being Used in Benchmarks

**Benchmark Report Finding**:
> "Missing features: TMB, HRD, MSI not included in predictions"

**MBD4 Plan Reality**:
- ✅ Tumor context schema supports TMB/HRD/MSI (`tumor_context.py`)
- ✅ Sporadic gates apply IO boosts based on TMB/MSI (`sporadic_gates.py`)
- ✅ PARP penalty/rescue based on HRD score (`sporadic_gates.py`)
- ✅ Disease priors can estimate TMB/HRD from disease type (`tumor_quick_intake.py`)

**Root Cause**: **Benchmark scripts are not extracting or passing TMB/HRD/MSI values** from TCGA dataset to `/api/efficacy/predict`.

**Impact**:
- IO boost not applied (TMB ≥20 or MSI-High should boost checkpoint inhibitors)
- PARP penalty/rescue not applied (HRD score affects PARP inhibitor efficacy)
- Efficacy scores missing critical biomarker signals
- **This explains weak correlations** - predictions are missing key prognostic features

**Priority**: **P0 - CRITICAL** (directly addresses weak correlation issue)

---

### 2. ⚠️ **HIGH**: Validation Set Bias (Lowest Mutation Counts)

**Benchmark Report Finding**:
> "Validation Mode: Selected patients with **lowest mutation counts** to avoid timeouts"

**Impact**:
- Low-mutation patients may have different outcomes (better prognosis)
- Not representative of typical ovarian cancer patients
- Could explain weak correlations (low-mutation patients may have different biology)

**MBD4 Plan Reality**:
- System can handle high mutation counts (MBD4+TP53 analysis worked)
- Timeout issues resolved with validation mode (2-18 mutations vs 60-80+)

**Recommendation**: 
- **Test with stratified sample** (mix of low/medium/high mutations)
- **Test with random sample** (not just lowest mutations)
- **Compare validation set vs full dataset** characteristics

**Priority**: **P1 - HIGH** (affects correlation validity)

---

### 3. ⚠️ **MEDIUM**: Top-1/Top-3 Accuracy is 0% (Top-5 is 100%)

**Benchmark Report Finding**:
> "Top-1 Accuracy: 0.000 (0%), Top-3 Accuracy: 0.000 (0%), Top-5 Accuracy: 1.000 (100%)"

**Interpretation**:
- System recommends clinically appropriate drugs (carboplatin, olaparib in top-5)
- But exact received drug is not ranked #1 (ranking precision issue)

**Possible Causes**:
1. **Drug name normalization**: "Carboplatin" vs "carboplatin" vs "CARBOPLATIN"
2. **Combination therapy**: Patient received "Carboplatin + Paclitaxel" but we rank single drugs
3. **Treatment line**: Drug may be appropriate for different line than received
4. **Efficacy score calibration**: Scores may need adjustment for ranking precision

**MBD4 Plan Reality**:
- Drug ranking works (100% Top-5)
- But ranking precision could be improved

**Recommendation**:
- Investigate drug name normalization
- Consider combination therapy matching
- Review efficacy score calibration for ranking

**Priority**: **P2 - MEDIUM** (Top-5 is excellent, but Top-1 would be better)

---

### 4. ⚠️ **MEDIUM**: Classification Cannot Be Assessed (PFS_STATUS Parsing)

**Benchmark Report Finding**:
> "Classification: AUC=0.000 (0 events, 0 censored - parsing issue?)"

**Possible Causes**:
1. PFS_STATUS field format not being parsed correctly
2. Validation set bias (low-mutation patients may have better outcomes)
3. Data quality issue (PFS_STATUS may be missing or incorrectly formatted)

**Recommendation**:
- Investigate PFS_STATUS field format in TCGA dataset
- Check actual PFS_STATUS values for test patients
- Fix parsing logic if needed
- Test with patients known to have progression events

**Priority**: **P2 - MEDIUM** (classification metrics are important but not blocking)

---

### 5. ⚠️ **LOW**: Small Sample Size (10-20 patients)

**Benchmark Report Finding**:
> "10 patients: Minimum for correlation (barely sufficient), 20 patients: Minimum for classification (barely sufficient)"

**Impact**:
- Low statistical power (weak correlations may be due to small sample)
- Need 50-100+ patients for robust statistics

**Recommendation**:
- Complete 50-patient test (survival analysis)
- Run 100-patient test for robust correlation analysis
- Use stratified sampling (not just lowest mutations)

**Priority**: **P3 - LOW** (infrastructure ready, just need to scale)

---

## Detailed Gap Analysis

### Gap 1: TMB/HRD/MSI Integration (CRITICAL)

**Current State**:
- ✅ System supports TMB/HRD/MSI via `tumor_context` parameter
- ✅ Sporadic gates apply boosts/penalties based on these values
- ❌ **Benchmark scripts do not extract or pass TMB/HRD/MSI from TCGA dataset**

**What Needs to Happen**:

1. **Extract TMB from TCGA Dataset**:
   ```python
   # In benchmark script
   tmb = patient_data.get("TMB", None)  # mutations per Mb
   # Or compute from mutation count: tmb = len(mutations) / genome_size_mb
   ```

2. **Extract HRD Score from TCGA Dataset**:
   ```python
   # In benchmark script
   hrd_score = patient_data.get("HRD_SCORE", None)
   # Or estimate from BRCA1/BRCA2 mutations + genomic scar signatures
   ```

3. **Extract MSI Status from TCGA Dataset**:
   ```python
   # In benchmark script
   msi_status = patient_data.get("MSI_STATUS", None)  # "MSI-H", "MSS", "MSI-L"
   ```

4. **Pass to API**:
   ```python
   # In benchmark API call
   tumor_context = {
       "disease": "ovarian_cancer",
       "tmb": tmb,
       "hrd_score": hrd_score,
       "msi_status": msi_status
   }
   
   response = await client.post("/api/efficacy/predict", {
       "mutations": mutations,
       "disease": "ovarian_cancer",
       "tumor_context": tumor_context
   })
   ```

**Expected Impact**:
- **IO Boost**: TMB ≥20 or MSI-High → checkpoint inhibitors get 1.3-1.35x boost
- **PARP Rescue**: HRD ≥42 → PARP inhibitors get full effect (no penalty)
- **Better Correlations**: Efficacy scores now include biomarker signals → should improve correlation with outcomes

**Files to Modify**:
- `scripts/benchmark/benchmark_small_test.py` - Extract TMB/HRD/MSI from TCGA data
- `scripts/benchmark/benchmark_common/api_client.py` - Pass tumor_context to API

**Estimated Effort**: 2-4 hours (data extraction + API integration)

---

### Gap 2: Validation Set Bias (HIGH)

**Current State**:
- ✅ Validation mode selects lowest mutation counts (2-18 mutations)
- ⚠️ This may not be representative of typical ovarian cancer patients

**What Needs to Happen**:

1. **Add Stratified Sampling**:
   ```python
   # In patient_selection.py
   def select_stratified_patients(patients, n, mutation_bins=[0, 10, 30, 100]):
       """Select patients stratified by mutation count."""
       low = [p for p in patients if mutation_bins[0] <= len(p.mutations) < mutation_bins[1]]
       medium = [p for p in patients if mutation_bins[1] <= len(p.mutations) < mutation_bins[2]]
       high = [p for p in patients if len(p.mutations) >= mutation_bins[2]]
       
       # Select equal numbers from each bin
       n_per_bin = n // 3
       return random.sample(low, min(n_per_bin, len(low))) + \
              random.sample(medium, min(n_per_bin, len(medium))) + \
              random.sample(high, min(n_per_bin, len(high)))
   ```

2. **Add Random Sampling Option**:
   ```python
   # In patient_selection.py
   def select_random_patients(patients, n):
       """Select random patients (not just lowest mutations)."""
       return random.sample(patients, min(n, len(patients)))
   ```

3. **Compare Validation Set vs Full Dataset**:
   ```python
   # In benchmark script
   validation_set_stats = {
       "mean_mutations": np.mean([len(p.mutations) for p in validation_patients]),
       "mean_tmb": np.mean([p.tmb for p in validation_patients if p.tmb]),
       "hrd_positive_rate": sum(1 for p in validation_patients if p.hrd_score >= 42) / len(validation_patients)
   }
   
   full_dataset_stats = {
       "mean_mutations": np.mean([len(p.mutations) for p in all_patients]),
       "mean_tmb": np.mean([p.tmb for p in all_patients if p.tmb]),
       "hrd_positive_rate": sum(1 for p in all_patients if p.hrd_score >= 42) / len(all_patients)
   }
   
   # Compare and report bias
   ```

**Expected Impact**:
- More representative sample → better correlation estimates
- Identifies if low-mutation bias explains weak correlations

**Files to Modify**:
- `scripts/benchmark/benchmark_common/patient_selection.py` - Add stratified/random sampling
- `scripts/benchmark/benchmark_small_test.py` - Add sampling mode option

**Estimated Effort**: 3-5 hours (sampling logic + comparison analysis)

---

### Gap 3: Drug Ranking Precision (MEDIUM)

**Current State**:
- ✅ Top-5 accuracy: 100% (excellent)
- ⚠️ Top-1 accuracy: 0% (needs improvement)

**What Needs to Happen**:

1. **Investigate Drug Name Normalization**:
   ```python
   # In drug_ranking.py
   def normalize_drug_name(name):
       """Normalize drug names for matching."""
       # Remove case sensitivity
       name = name.lower().strip()
       # Handle common variations
       variations = {
           "carboplatin": ["carboplatin", "carboplatinum", "paraplatin"],
           "olaparib": ["olaparib", "lynparza"],
           "paclitaxel": ["paclitaxel", "taxol", "abraxane"]
       }
       for canonical, variants in variations.items():
           if name in variants:
               return canonical
       return name
   ```

2. **Handle Combination Therapy**:
   ```python
   # In drug_ranking.py
   def match_combination_therapy(received_drugs, recommended_drugs):
       """Match combination therapy to single drug recommendations."""
       # If patient received "Carboplatin + Paclitaxel"
       # Check if both drugs are in top-5 recommendations
       received_normalized = [normalize_drug_name(d) for d in received_drugs]
       recommended_normalized = [normalize_drug_name(d["name"]) for d in recommended_drugs[:5]]
       
       # Match if any received drug is in top-5
       return any(d in recommended_normalized for d in received_normalized)
   ```

3. **Review Efficacy Score Calibration**:
   - Check if scores need adjustment for ranking precision
   - Consider treatment line appropriateness
   - Review confidence vs efficacy score relationship

**Expected Impact**:
- Improved Top-1 accuracy (exact drug match)
- Better ranking precision

**Files to Modify**:
- `scripts/benchmark/benchmark_common/metrics/drug_ranking.py` - Add normalization + combination matching
- Review `api/services/efficacy_orchestrator/drug_scorer.py` - Score calibration

**Estimated Effort**: 2-3 hours (normalization + combination matching)

---

### Gap 4: PFS_STATUS Parsing (MEDIUM)

**Current State**:
- ⚠️ Classification shows 0 events AND 0 censored (parsing issue)

**What Needs to Happen**:

1. **Investigate PFS_STATUS Format**:
   ```python
   # In benchmark script
   # Check actual PFS_STATUS values
   for patient in test_patients:
       print(f"Patient {patient.id}: PFS_STATUS = {patient.pfs_status}, PFS_MONTHS = {patient.pfs_months}")
   ```

2. **Fix Parsing Logic**:
   ```python
   # In classification.py
   def parse_pfs_status(pfs_status, pfs_months):
       """Parse PFS_STATUS field."""
       # Common formats:
       # - "0" = censored, "1" = event
       # - "Censored" = censored, "Progressed" = event
       # - "0: Censored" = censored, "1: Progressed" = event
       
       if pfs_status is None:
           # Try to infer from PFS_MONTHS
           return None  # Unknown
       
       pfs_status_str = str(pfs_status).upper()
       
       if "0" in pfs_status_str or "CENSORED" in pfs_status_str:
           return False  # Censored
       elif "1" in pfs_status_str or "PROGRESSED" in pfs_status_str or "EVENT" in pfs_status_str:
           return True  # Event
       else:
           return None  # Unknown
   ```

3. **Test with Known Progression Events**:
   - Select patients with known progression events
   - Validate parsing logic works correctly

**Expected Impact**:
- Classification metrics can be assessed
- Better understanding of progression prediction

**Files to Modify**:
- `scripts/benchmark/benchmark_common/metrics/classification.py` - Fix parsing logic
- `scripts/benchmark/benchmark_small_test.py` - Add PFS_STATUS inspection

**Estimated Effort**: 2-3 hours (parsing fix + validation)

---

## Implementation Priority & Timeline

### Phase 1: Critical Fixes (Week 1) - **P0**

**Goal**: Address weak correlation issue by integrating TMB/HRD/MSI

**Tasks**:
1. ✅ Extract TMB/HRD/MSI from TCGA dataset (2-4 hours)
2. ✅ Pass tumor_context to API in benchmark scripts (1-2 hours)
3. ✅ Validate sporadic gates are applied correctly (1 hour)
4. ✅ Re-run 20-patient test with TMB/HRD/MSI (1 hour)

**Expected Outcome**:
- Efficacy scores include biomarker signals
- IO boost applied for TMB-high patients
- PARP rescue applied for HRD-high patients
- **Correlation should improve** (r=0.037-0.278 → r=0.2-0.4+)

**Success Criteria**:
- TMB/HRD/MSI values extracted from TCGA dataset
- Tumor context passed to API in all benchmark calls
- Sporadic gates applied (IO boost, PARP rescue visible in scores)
- Correlation improves (r > 0.2, p < 0.1)

---

### Phase 2: Validation Set Bias (Week 2) - **P1**

**Goal**: Test with representative sample to validate correlation improvements

**Tasks**:
1. ✅ Add stratified sampling (mix of low/medium/high mutations) (2-3 hours)
2. ✅ Add random sampling option (1 hour)
3. ✅ Compare validation set vs full dataset characteristics (1-2 hours)
4. ✅ Re-run 50-patient test with stratified sampling (1 hour)

**Expected Outcome**:
- More representative sample
- Better correlation estimates
- Identify if low-mutation bias explains weak correlations

**Success Criteria**:
- Stratified sampling implemented
- Validation set characteristics match full dataset
- Correlation results are representative

---

### Phase 3: Ranking Precision & Classification (Week 3) - **P2**

**Goal**: Improve Top-1 accuracy and enable classification assessment

**Tasks**:
1. ✅ Fix drug name normalization (1-2 hours)
2. ✅ Handle combination therapy matching (1 hour)
3. ✅ Fix PFS_STATUS parsing (2-3 hours)
4. ✅ Re-run benchmarks with fixes (1 hour)

**Expected Outcome**:
- Top-1 accuracy improves (0% → 20-40%)
- Classification metrics can be assessed
- Better ranking precision

**Success Criteria**:
- Top-1 accuracy > 20%
- Classification metrics compute correctly
- PFS_STATUS parsing works

---

## Expected Improvements After Fixes

### Correlation Improvements

**Current**: r=0.037-0.278 (negligible to weak, not significant)

**After TMB/HRD/MSI Integration**:
- **Expected**: r=0.2-0.4 (weak to moderate, potentially significant)
- **Rationale**: Efficacy scores now include biomarker signals (TMB, HRD, MSI) that are known to correlate with outcomes

**After Stratified Sampling**:
- **Expected**: More stable correlation estimates
- **Rationale**: Representative sample reduces bias

**Target**: r > 0.3, p < 0.05 (moderate correlation, statistically significant)

---

### Drug Ranking Improvements

**Current**: Top-5: 100%, Top-1: 0%

**After Drug Name Normalization**:
- **Expected**: Top-1: 20-40%
- **Rationale**: Better matching between received and recommended drugs

**After Combination Therapy Matching**:
- **Expected**: Top-1: 30-50%
- **Rationale**: Handles "Carboplatin + Paclitaxel" combinations

**Target**: Top-1 > 40%, Top-3 > 60%

---

### Classification Improvements

**Current**: Cannot assess (0 events, 0 censored)

**After PFS_STATUS Parsing Fix**:
- **Expected**: Classification metrics compute correctly
- **Rationale**: Proper parsing enables event/censored detection

**Target**: AUC > 0.6 (better than random)

---

## Code Locations for Implementation

### TMB/HRD/MSI Integration

**Files to Modify**:
1. `scripts/benchmark/benchmark_small_test.py` (lines 200-300) - Extract TMB/HRD/MSI
2. `scripts/benchmark/benchmark_common/api_client.py` (lines 50-100) - Pass tumor_context
3. `scripts/benchmark/benchmark_common/data_loader.py` (lines 100-200) - Load TCGA biomarker data

**Key Functions**:
- `extract_tmb_from_patient(patient_data)` - Extract TMB
- `extract_hrd_from_patient(patient_data)` - Extract HRD score
- `extract_msi_from_patient(patient_data)` - Extract MSI status
- `build_tumor_context(tmb, hrd, msi)` - Build tumor context dict

---

### Stratified Sampling

**Files to Modify**:
1. `scripts/benchmark/benchmark_common/patient_selection.py` (lines 50-150) - Add sampling methods

**Key Functions**:
- `select_stratified_patients(patients, n, bins)` - Stratified sampling
- `select_random_patients(patients, n)` - Random sampling
- `compare_validation_set_stats(validation, full_dataset)` - Bias analysis

---

### Drug Ranking Precision

**Files to Modify**:
1. `scripts/benchmark/benchmark_common/metrics/drug_ranking.py` (lines 20-80) - Add normalization

**Key Functions**:
- `normalize_drug_name(name)` - Normalize drug names
- `match_combination_therapy(received, recommended)` - Handle combinations

---

### PFS_STATUS Parsing

**Files to Modify**:
1. `scripts/benchmark/benchmark_common/metrics/classification.py` (lines 30-70) - Fix parsing

**Key Functions**:
- `parse_pfs_status(pfs_status, pfs_months)` - Parse PFS_STATUS field
- `validate_pfs_status_parsing(patients)` - Validation

---

## Success Metrics

### Phase 1 Success Criteria

- ✅ TMB/HRD/MSI extracted from TCGA dataset (100% of patients)
- ✅ Tumor context passed to API (100% of API calls)
- ✅ Sporadic gates applied (IO boost, PARP rescue visible in scores)
- ✅ Correlation improves: r > 0.2 (from r=0.037-0.278)
- ✅ P-value improves: p < 0.1 (from p=0.4-0.9)

### Phase 2 Success Criteria

- ✅ Stratified sampling implemented
- ✅ Validation set characteristics match full dataset (within 20%)
- ✅ Correlation results are representative

### Phase 3 Success Criteria

- ✅ Top-1 accuracy > 20% (from 0%)
- ✅ Classification metrics compute correctly
- ✅ PFS_STATUS parsing works (events/censored detected)

---

## Risk Assessment

### Risk 1: TCGA Dataset May Not Have TMB/HRD/MSI

**Mitigation**:
- Check TCGA dataset schema first
- If missing, estimate from mutations (TMB = mutation count / genome size)
- Use disease priors as fallback (`tumor_quick_intake.py`)

**Probability**: Medium  
**Impact**: High  
**Mitigation Effort**: 2-4 hours

---

### Risk 2: Correlation May Not Improve Even With TMB/HRD/MSI

**Mitigation**:
- This is expected - correlation improvement is not guaranteed
- Focus on drug ranking (already excellent)
- Consider treatment-adjusted outcomes

**Probability**: Medium  
**Impact**: Medium  
**Mitigation Effort**: N/A (expected outcome)

---

### Risk 3: Stratified Sampling May Not Be Representative

**Mitigation**:
- Compare validation set vs full dataset characteristics
- Use random sampling as alternative
- Document limitations

**Probability**: Low  
**Impact**: Medium  
**Mitigation Effort**: 1-2 hours (comparison analysis)

---

## Conclusion

**Key Finding**: **Infrastructure exists, but benchmark integration is missing**. The system can use TMB/HRD/MSI, but benchmarks aren't providing them.

**Critical Action**: **Extract and pass TMB/HRD/MSI from TCGA dataset to API** - This is the #1 priority and should directly address weak correlation issue.

**Expected Impact**: 
- Correlation: r=0.037-0.278 → r=0.2-0.4+ (weak to moderate)
- Drug ranking: Top-1: 0% → 20-40% (with normalization)
- Classification: Can be assessed (with parsing fix)

**Timeline**: 3 weeks (Phase 1: 1 week, Phase 2: 1 week, Phase 3: 1 week)

---

**Report Prepared By**: Auto (Agent)  
**Date**: January 27, 2025  
**For Review By**: Manager  
**Status**: Ready for Implementation

