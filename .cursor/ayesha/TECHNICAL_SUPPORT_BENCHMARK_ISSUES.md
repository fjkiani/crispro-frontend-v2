# Benchmark Accuracy Issues - My Ownership & Deliverables

**Date**: January 27, 2025  
**Owner**: Auto (Agent with SPE/WIWFM expertise)  
**Report Reviewed**: `BENCHMARK_ACCURACY_REPORT.md`  
**Status**: 🔧 **READY TO IMPLEMENT** - Full ownership defined

---

## My Ownership: What I Will Deliver

### ✅ **PRIMARY OWNERSHIP**: Benchmark Integration Fixes (5 Issues)

**My Primary Scope**: Fix benchmark script bugs and integration issues. These are definite bugs I can fix.

**What I Will Deliver**:

1. **TMB/HRD/MSI Integration** (P0 - CRITICAL)
   - Extract biomarkers from TCGA dataset
   - Pass `tumor_context` to API
   - **Deliverable**: Biomarker extraction utility + API integration
   - **Acceptance**: 100% of API calls include tumor_context, sporadic gates applied

2. **PFS_STATUS Parser** (P0)
   - Fix classification parsing bug
   - **Deliverable**: Parser utility + fixed classification metrics
   - **Acceptance**: Classification metrics compute correctly (events/censored detected)

3. **Drug Name Normalizer** (P1)
   - Improve Top-1/Top-3 accuracy
   - **Deliverable**: Normalizer utility + fixed drug ranking metrics
   - **Acceptance**: Top-1 accuracy > 20% (from 0%)

4. **Patient Sampling** (P1)
   - Random/stratified sampling options
   - **Deliverable**: Sampling functions + bias comparison
   - **Acceptance**: Stratified sampling works, bias comparison reports correctly

5. **Data Quality Checker** (P1)
   - Validate dataset before running
   - **Deliverable**: Data quality validation script
   - **Acceptance**: Validates PFS_STATUS parsing, reports data issues

**Total Deliverables**: 5 utilities + 3 metric fixes = 8 components  
**Total Time**: 14-15 hours (2-3 days for P0, 1 week for all)

---

### ✅ **SECONDARY OWNERSHIP**: SPE/WIWFM Improvements (If Needed)

**My Secondary Scope**: After benchmark integration fixes, if correlations are still weak, I can investigate and improve SPE/WIWFM model.

**Why I Can Own This**: I have deep knowledge of:
- S/P/E framework architecture (`drug_scorer.py`, `orchestrator.py`)
- Efficacy score formula: `0.3 * seq_pct + 0.4 * path_pct + 0.3 * evidence + clinvar_prior`
- Confidence computation (`confidence_computation.py`)
- Pathway aggregation (`pathway/aggregation.py`)
- Sporadic gates (`sporadic_gates.py`)

**Potential SPE/WIWFM Improvements** (If Needed After Benchmark Fixes):

1. **Score Calibration** (If correlations still weak)
   - Efficacy score calibration for outcome prediction
   - **My Approach**: Review S/P/E weights, add outcome-based calibration
   - **Files I Know**: `drug_scorer.py`, `confidence_computation.py`

2. **Model Improvements** (If correlations still weak)
   - Review S/P/E formula weights (currently 0.3/0.4/0.3)
   - Consider treatment-adjusted outcomes
   - **My Approach**: Analyze correlation patterns, adjust weights if needed
   - **Files I Know**: `drug_scorer.py`, `orchestrator.py`

**Decision Point**: After benchmark integration fixes complete:
- **If correlation improves** (r > 0.2): ✅ Success, no SPE/WIWFM changes needed
- **If correlation still weak** (r < 0.2): Investigate SPE/WIWFM improvements (I can handle this)

**Priority**: Fix benchmark integration FIRST, then assess if SPE/WIWFM improvements needed.

---

## Issue 0: TMB/HRD/MSI Integration (P0 - CRITICAL - I OWN)

### Problem Statement

**Symptom**: Benchmark report says "Missing features: TMB, HRD, MSI not included in predictions"  
**Root Cause**: Benchmark scripts are not extracting or passing TMB/HRD/MSI from TCGA dataset to `/api/efficacy/predict`  
**Impact**: 
- IO boost not applied (TMB ≥20 or MSI-High should boost checkpoint inhibitors)
- PARP penalty/rescue not applied (HRD score affects PARP inhibitor efficacy)
- Efficacy scores missing critical biomarker signals
- **This explains weak correlations** - predictions are missing key prognostic features

### System Status

**✅ System SUPPORTS TMB/HRD/MSI** (MBD4 Plan Analysis):
- Tumor context schema supports TMB/HRD/MSI (`tumor_context.py`)
- Sporadic gates apply IO boosts based on TMB/MSI (`sporadic_gates.py`)
- PARP penalty/rescue based on HRD score (`sporadic_gates.py`)
- Disease priors can estimate TMB/HRD from disease type (`tumor_quick_intake.py`)

**❌ Benchmark Integration Missing**:
- Benchmark scripts do not extract TMB/HRD/MSI from TCGA dataset
- Benchmark scripts do not pass `tumor_context` to API

**Conclusion**: This is a **benchmark integration issue** (I own), not a SPE/WIWFM model issue.

---

### My Deliverables

#### Deliverable 1: Biomarker Extractor Utility

**File**: `scripts/benchmark/benchmark_common/utils/biomarker_extractor.py` (NEW)

**Functions I Will Create**:

```python
def extract_tmb_from_patient(patient_data: Dict) -> Optional[float]:
    """
    Extract TMB from TCGA patient data.
    
    Tries multiple field names:
    - TMB (mutations per Mb)
    - TMB_SCORE
    - Computed from mutation count if available
    
    Returns: TMB value (mutations per Mb) or None
    """
    # Try direct field
    tmb = patient_data.get("TMB") or patient_data.get("TMB_SCORE")
    if tmb is not None:
        return float(tmb)
    
    # Compute from mutation count if available
    mutations = patient_data.get("mutations", [])
    if mutations:
        # Genome size ~3.2 Gb, so ~3200 Mb
        tmb = len(mutations) / 3200.0
        return tmb
    
    return None

def extract_hrd_from_patient(patient_data: Dict) -> Optional[float]:
    """
    Extract HRD score from TCGA patient data.
    
    Tries multiple field names:
    - HRD_SCORE
    - HRD_SCORE_MYCHOICE
    - Estimate from BRCA1/BRCA2 mutations if available
    
    Returns: HRD score or None
    """
    hrd = patient_data.get("HRD_SCORE") or patient_data.get("HRD_SCORE_MYCHOICE")
    if hrd is not None:
        return float(hrd)
    
    # Estimate from BRCA mutations (if available)
    mutations = patient_data.get("mutations", [])
    brca_mutations = [m for m in mutations if m.get("gene", "").upper() in ["BRCA1", "BRCA2"]]
    if brca_mutations:
        # Rough estimate: BRCA mutation suggests HRD+ (score ≥42)
        return 45.0  # Conservative estimate
    
    return None

def extract_msi_from_patient(patient_data: Dict) -> Optional[str]:
    """
    Extract MSI status from TCGA patient data.
    
    Returns: "MSI-H", "MSS", "MSI-L", or None
    """
    msi = patient_data.get("MSI_STATUS") or patient_data.get("MSI")
    if msi is None:
        return None
    
    msi_str = str(msi).upper().strip()
    
    # Normalize to standard format
    if "MSI-H" in msi_str or "MSI-HIGH" in msi_str or msi_str == "1":
        return "MSI-H"
    elif "MSI-L" in msi_str or "MSI-LOW" in msi_str:
        return "MSI-L"
    elif "MSS" in msi_str or msi_str == "0":
        return "MSS"
    
    return None

def build_tumor_context(
    disease: str,
    tmb: Optional[float] = None,
    hrd_score: Optional[float] = None,
    msi_status: Optional[str] = None
) -> Dict[str, Any]:
    """
    Build tumor_context dict for API call.
    
    Returns: tumor_context dict ready for API
    """
    context = {
        "disease": disease
    }
    
    if tmb is not None:
        context["tmb"] = tmb
    
    if hrd_score is not None:
        context["hrd_score"] = hrd_score
    
    if msi_status is not None:
        context["msi_status"] = msi_status
    
    return context
```

**Acceptance Criteria**:
- ✅ Handles multiple field name formats (TMB, TMB_SCORE, etc.)
- ✅ Computes TMB from mutation count if field missing
- ✅ Estimates HRD from BRCA mutations if field missing
- ✅ Normalizes MSI to standard format
- ✅ Unit tests pass (test with sample TCGA data)

**Estimated Time**: 2 hours

---

#### Deliverable 2: API Client Update

**File**: `scripts/benchmark/benchmark_common/api_client.py` (UPDATE)

**Changes I Will Make**:

```python
async def call_efficacy_predict(
    client: httpx.AsyncClient,
    mutations: List[Dict],
    disease: str,
    tumor_context: Optional[Dict] = None,  # ← ADD THIS PARAMETER
    api_base: str = "http://localhost:8000"
) -> Dict:
    """
    Call /api/efficacy/predict with tumor_context support.
    """
    payload = {
        "mutations": mutations,
        "disease": disease
    }
    
    # Add tumor_context if provided ← ADD THIS LOGIC
    if tumor_context:
        payload["tumor_context"] = tumor_context
    
    response = await client.post(
        f"{api_base}/api/efficacy/predict",
        json=payload,
        timeout=60.0
    )
    
    return response.json()
```

**Acceptance Criteria**:
- ✅ `tumor_context` parameter added to function signature
- ✅ `tumor_context` included in API payload when provided
- ✅ Backward compatible (works without tumor_context)
- ✅ Logging added for biomarker values (debug mode)

**Estimated Time**: 1 hour

---

#### Deliverable 3: Benchmark Script Integration

**File**: `scripts/benchmark/benchmark_small_test.py` (UPDATE)

**Changes I Will Make**:

```python
from benchmark_common.utils.biomarker_extractor import (
    extract_tmb_from_patient,
    extract_hrd_from_patient,
    extract_msi_from_patient,
    build_tumor_context
)

async def process_patient(patient: Dict, disease: str, api_client):
    """
    Process single patient with biomarker extraction.
    """
    # Extract biomarkers ← ADD THIS
    tmb = extract_tmb_from_patient(patient)
    hrd_score = extract_hrd_from_patient(patient)
    msi_status = extract_msi_from_patient(patient)
    
    # Build tumor context ← ADD THIS
    tumor_context = build_tumor_context(
        disease=disease,
        tmb=tmb,
        hrd_score=hrd_score,
        msi_status=msi_status
    )
    
    # Call API with tumor_context ← UPDATE THIS
    response = await api_client.call_efficacy_predict(
        mutations=patient["mutations"],
        disease=disease,
        tumor_context=tumor_context  # ← ADD THIS
    )
    
    return response
```

**Acceptance Criteria**:
- ✅ Biomarkers extracted for each patient
- ✅ `tumor_context` built and passed to API
- ✅ 100% of API calls include tumor_context (when available)
- ✅ Validation: Check scores for IO boost/PARP rescue (TMB-high/HRD-high patients)

**Estimated Time**: 1-2 hours

---

#### Deliverable 4: Integration Validation

**File**: `scripts/benchmark/validate_biomarker_integration.py` (NEW - Optional)

**What I Will Create**:
- Script to validate sporadic gates are applied
- Check IO boost for TMB-high patients
- Check PARP rescue for HRD-high patients
- Report biomarker extraction success rate

**Acceptance Criteria**:
- ✅ IO boost visible in scores (TMB ≥20 patients)
- ✅ PARP rescue visible in scores (HRD ≥42 patients)
- ✅ Biomarker extraction success rate > 80%

**Estimated Time**: 1 hour

---

### Expected Impact

**After My Fixes**:
- **IO Boost**: TMB ≥20 or MSI-High → checkpoint inhibitors get 1.3-1.35x boost
- **PARP Rescue**: HRD ≥42 → PARP inhibitors get full effect (no penalty)
- **Better Correlations**: Efficacy scores now include biomarker signals → should improve correlation with outcomes
- **Expected Correlation Improvement**: r=0.037-0.278 → r=0.2-0.4+ (weak to moderate)

**Success Metrics**:
- TMB/HRD/MSI extracted from TCGA dataset (100% of patients where available)
- Tumor context passed to API (100% of API calls)
- Sporadic gates applied (IO boost, PARP rescue visible in scores)
- Correlation improves: r > 0.2 (from r=0.037-0.278)

---

### Risk Mitigation

**Risk**: TCGA dataset may not have TMB/HRD/MSI fields  
**My Mitigation**: 
- Check dataset schema first
- If missing, estimate from mutations (TMB = mutation count / genome size)
- Use disease priors as fallback (reference `tumor_quick_intake.py`)
- Document fallback logic in code comments

**Probability**: Medium  
**Impact**: High  
**Mitigation Effort**: 2-4 hours (fallback logic)

---

## Issue 1: Classification Parsing Bug (P0 - I OWN)

### Problem Statement

**Symptom**: Classification shows `0 events AND 0 censored` (impossible)  
**Expected**: Should sum to n patients  
**Impact**: Cannot assess classification accuracy

### My Deliverables

#### Deliverable 1: PFS_STATUS Parser Utility

**File**: `scripts/benchmark/benchmark_common/utils/pfs_status_parser.py` (NEW)

**Function I Will Create**:

```python
def parse_pfs_status(value: Any) -> Tuple[Optional[int], bool]:
    """
    Parse PFS_STATUS field from cBioPortal data.
    
    Handles multiple formats:
    - "0:DiseaseFree" or "1:Recurred/Progressed"
    - "DISEASEFREE", "RECURRED", "PROGRESSED"
    - Numeric (0 or 1)
    
    Returns:
        (event, is_valid): (1=progressed, 0=censored, None=invalid), bool
    """
    if value is None or pd.isna(value):
        return None, False
    
    value_str = str(value).strip().upper()
    
    # Pattern 1: "0:DiseaseFree" or "1:Recurred/Progressed"
    if ":" in value_str:
        code = value_str.split(":")[0]
        if code == "0":
            return 0, True  # Censored
        elif code == "1":
            return 1, True  # Event
        else:
            return None, False
    
    # Pattern 2: "DISEASEFREE", "RECURRED", "PROGRESSED"
    if any(word in value_str for word in ["DISEASEFREE", "CENSORED", "NO RECURRENCE"]):
        return 0, True  # Censored
    elif any(word in value_str for word in ["RECURRED", "PROGRESSED", "EVENT"]):
        return 1, True  # Event
    
    # Pattern 3: Numeric (0 or 1)
    try:
        code = int(value_str)
        if code == 0:
            return 0, True
        elif code == 1:
            return 1, True
    except ValueError:
        pass
    
    return None, False
```

**Acceptance Criteria**:
- ✅ Handles multiple field name formats (PFS_STATUS, DFS_STATUS, PFS_STATUS_STR)
- ✅ Handles multiple value formats (numeric, string, colon-separated)
- ✅ Returns (event, is_valid) tuple
- ✅ Unit tests pass (test with sample cBioPortal data)

**Estimated Time**: 1 hour

---

#### Deliverable 2: Validation Script

**File**: `scripts/benchmark/validate_pfs_status.py` (NEW)

**What I Will Create**:
- Script to inspect actual PFS_STATUS field values in dataset
- Report field name distribution
- Report value format distribution
- Test parser on real data

**Acceptance Criteria**:
- ✅ Inspects dataset field names
- ✅ Reports value distributions
- ✅ Tests parser on sample values
- ✅ Identifies parsing issues

**Estimated Time**: 1 hour

---

#### Deliverable 3: Fixed Classification Metrics

**File**: `scripts/benchmark/benchmark_common/metrics/classification.py` (FIX)

**Changes I Will Make**:

```python
from benchmark_common.utils.pfs_status_parser import parse_pfs_status

def compute_classification_metrics(patients: List[Dict]) -> Dict:
    """
    Fixed version with proper PFS_STATUS parsing.
    """
    events = []
    labels = []
    
    for patient in patients:
        # Try multiple field names ← FIX THIS
        pfs_status = (
            patient.get("PFS_STATUS") or 
            patient.get("DFS_STATUS") or 
            patient.get("PFS_STATUS_STR") or
            patient.get("pfs_status")
        )
        
        event, is_valid = parse_pfs_status(pfs_status)  # ← USE PARSER
        
        if is_valid:
            events.append(event)
            # Use efficacy score as prediction
            labels.append(patient.get("efficacy_score", 0.0))
    
    if len(events) == 0:
        return {
            "roc_auc": None,
            "n_events": 0,
            "n_censored": 0,
            "error": "No valid PFS_STATUS values found"
        }
    
    # Compute ROC-AUC
    # ... existing logic
```

**Acceptance Criteria**:
- ✅ Uses PFS_STATUS parser utility
- ✅ Handles multiple field names
- ✅ Error handling for missing/invalid data
- ✅ Logging for debugging
- ✅ Classification metrics compute correctly (events/censored detected)

**Estimated Time**: 1 hour

---

### Expected Impact

**After My Fixes**:
- Classification metrics can be assessed
- Events and censored patients correctly identified
- ROC-AUC computed correctly

**Success Metrics**:
- PFS_STATUS parsing works (events/censored detected)
- Classification metrics compute correctly
- No more "0 events AND 0 censored" error

---

## Issue 2: Drug Name Matching (P1 - I OWN)

### Problem Statement

**Symptom**: Top-1 accuracy = 0%, Top-3 accuracy = 0%, Top-5 accuracy = 100%  
**Interpretation**: System recommends correct drugs but not in exact order  
**Impact**: Lower precision in top recommendations

### My Deliverables

#### Deliverable 1: Drug Name Normalizer Utility

**File**: `scripts/benchmark/benchmark_common/utils/drug_normalizer.py` (NEW)

**Functions I Will Create**:

```python
DRUG_SYNONYMS = {
    "carboplatin": ["carboplatin", "paraplatin", "cbdca"],
    "paclitaxel": ["paclitaxel", "taxol", "abraxane"],
    "olaparib": ["olaparib", "lynparza"],
    "niraparib": ["niraparib", "zejula"],
    "rucaparib": ["rucaparib", "rubraca"],
    # ... more mappings for ovarian cancer drugs
}

def normalize_drug_name(name: str) -> str:
    """
    Normalize drug name for matching.
    
    Handles:
    - Case insensitivity
    - Synonyms (brand names, abbreviations)
    - Whitespace
    
    Returns: Canonical drug name
    """
    name_lower = name.lower().strip()
    
    # Check synonyms
    for canonical, synonyms in DRUG_SYNONYMS.items():
        if name_lower in synonyms or name_lower == canonical:
            return canonical
    
    return name_lower

def match_drug(received: str, recommended: List[str]) -> bool:
    """
    Check if received drug matches any recommended drug.
    
    Returns: True if match found
    """
    received_norm = normalize_drug_name(received)
    recommended_norm = [normalize_drug_name(d) for d in recommended]
    
    return received_norm in recommended_norm
```

**Acceptance Criteria**:
- ✅ Synonym dictionary for common ovarian cancer drugs
- ✅ Case-insensitive matching
- ✅ Handles brand names and abbreviations
- ✅ Unit tests pass

**Estimated Time**: 2 hours

---

#### Deliverable 2: Fixed Drug Ranking Metrics

**File**: `scripts/benchmark/benchmark_common/metrics/drug_ranking.py` (FIX)

**Changes I Will Make**:

```python
from benchmark_common.utils.drug_normalizer import normalize_drug_name

def compute_drug_ranking_accuracy(
    patients: List[Dict],
    predictions: List[Dict]
) -> Dict:
    """
    Fixed version with proper drug name matching.
    """
    top_1_matches = 0
    top_3_matches = 0
    top_5_matches = 0
    
    for patient, prediction in zip(patients, predictions):
        received_drugs = patient.get("treatments", [])
        recommended_drugs = prediction.get("top_drugs", [])
        
        # Normalize all drug names ← ADD THIS
        received_norm = [normalize_drug_name(d) for d in received_drugs]
        recommended_norm = [normalize_drug_name(d) for d in recommended_drugs]
        
        # Check Top-1 ← FIX THIS
        if received_norm and recommended_norm:
            if received_norm[0] == recommended_norm[0]:
                top_1_matches += 1
        
        # Check Top-3 ← FIX THIS
        if any(drug in recommended_norm[:3] for drug in received_norm):
            top_3_matches += 1
        
        # Check Top-5
        if any(drug in recommended_norm[:5] for drug in received_norm):
            top_5_matches += 1
    
    return {
        "top_1_accuracy": top_1_matches / len(patients),
        "top_3_accuracy": top_3_matches / len(patients),
        "top_5_accuracy": top_5_matches / len(patients),
    }
```

**Acceptance Criteria**:
- ✅ Uses drug normalizer utility
- ✅ Top-1/Top-3 matching logic fixed
- ✅ Logging for debugging
- ✅ Top-1 accuracy > 20% (from 0%)

**Estimated Time**: 1 hour

---

### Expected Impact

**After My Fixes**:
- Top-1 accuracy: 0% → 20-40%
- Top-3 accuracy: 0% → 40-60%
- Better ranking precision

**Success Metrics**:
- Top-1 accuracy > 20% (from 0%)
- Top-3 accuracy > 40% (from 0%)
- Drug name matching works correctly

---

## Issue 3: Patient Selection Bias (P1 - I OWN)

### Problem Statement

**Symptom**: Testing with lowest mutation counts (2-18 mutations)  
**Impact**: May not represent typical ovarian cancer patients  
**Concern**: Could explain weak correlations

### My Deliverables

#### Deliverable 1: Sampling Functions

**File**: `scripts/benchmark/benchmark_common/patient_selection.py` (UPDATE)

**Functions I Will Add**:

```python
def select_random_sample(
    patients: List[Dict],
    n: int,
    seed: int = 42
) -> List[Dict]:
    """
    Select random sample of patients.
    
    Returns: Random sample of n patients
    """
    import random
    random.seed(seed)
    return random.sample(patients, min(n, len(patients)))

def select_stratified_sample(
    patients: List[Dict],
    n: int,
    stratify_by: str = "mutation_count",
    seed: int = 42
) -> List[Dict]:
    """
    Select stratified sample (mix of low/medium/high mutations).
    
    Returns: Stratified sample of n patients
    """
    # Group by mutation count
    low = [p for p in patients if p.get("mutation_count", 0) < 20]
    medium = [p for p in patients if 20 <= p.get("mutation_count", 0) < 50]
    high = [p for p in patients if p.get("mutation_count", 0) >= 50]
    
    # Sample proportionally
    n_per_group = n // 3
    random.seed(seed)
    
    selected = []
    selected.extend(random.sample(low, min(n_per_group, len(low))))
    selected.extend(random.sample(medium, min(n_per_group, len(medium))))
    selected.extend(random.sample(high, min(n_per_group, len(high))))
    
    # Fill remaining with random
    remaining = n - len(selected)
    if remaining > 0:
        all_patients = low + medium + high
        selected.extend(random.sample(all_patients, min(remaining, len(all_patients))))
    
    return selected

def compare_validation_set_stats(
    validation_patients: List[Dict],
    full_dataset_patients: List[Dict]
) -> Dict:
    """
    Compare validation set vs full dataset characteristics.
    
    Returns: Comparison stats with bias detection
    """
    def compute_stats(patients):
        mutations = [len(p.get("mutations", [])) for p in patients]
        return {
            "mean_mutations": np.mean(mutations) if mutations else 0,
            "median_mutations": np.median(mutations) if mutations else 0,
            "min_mutations": min(mutations) if mutations else 0,
            "max_mutations": max(mutations) if mutations else 0,
            "n_patients": len(patients)
        }
    
    validation_stats = compute_stats(validation_patients)
    full_stats = compute_stats(full_dataset_patients)
    
    return {
        "validation": validation_stats,
        "full_dataset": full_stats,
        "bias_detected": abs(validation_stats["mean_mutations"] - full_stats["mean_mutations"]) > 5
    }
```

**Acceptance Criteria**:
- ✅ Random sampling works correctly
- ✅ Stratified sampling works correctly (mix of low/medium/high)
- ✅ Bias comparison function works
- ✅ Unit tests pass

**Estimated Time**: 1.5 hours

---

#### Deliverable 2: Benchmark Script Update

**File**: `scripts/benchmark/benchmark_small_test.py` (UPDATE)

**Changes I Will Make**:

```python
def select_patients(
    patients: List[Dict],
    n: int,
    mode: str = "random"  # ← ADD MODE PARAMETER
) -> List[Dict]:
    """
    Select patients based on mode.
    
    Modes:
    - "validation": Lowest mutations (for timeout avoidance)
    - "random": Random sample
    - "stratified": Mix of low/medium/high mutations
    """
    if mode == "validation":
        # Original: lowest mutations (for timeout avoidance)
        sorted_patients = sorted(patients, key=lambda p: p.get("mutation_count", 0))
        return sorted_patients[:n]
    elif mode == "random":
        return select_random_sample(patients, n)
    elif mode == "stratified":
        return select_stratified_sample(patients, n)
    else:
        raise ValueError(f"Unknown mode: {mode}")

# Add --mode argument to CLI
# Add bias comparison reporting
```

**Acceptance Criteria**:
- ✅ `--mode` argument added (random/stratified/validation)
- ✅ Patient selection logic updated
- ✅ Bias comparison reporting added
- ✅ Works with different modes

**Estimated Time**: 30 minutes

---

### Expected Impact

**After My Fixes**:
- More representative sample
- Better correlation estimates
- Identifies if low-mutation bias explains weak correlations

**Success Metrics**:
- Stratified sampling implemented
- Validation set characteristics match full dataset (within 20%)
- Bias comparison reports correctly

---

## Issue 4: Data Quality Validation (P1 - I OWN)

### Problem Statement

**Symptom**: Need to validate PFS_STATUS field exists and is parseable  
**Impact**: Prevents classification bug from recurring

### My Deliverables

#### Deliverable: Data Quality Checker

**File**: `scripts/benchmark/benchmark_common/data_quality.py` (NEW)

**Function I Will Create**:

```python
def validate_dataset_quality(dataset_path: str) -> Dict:
    """
    Validates benchmark dataset quality.
    
    Checks:
    - Required fields exist
    - PFS_STATUS field format
    - Outcome distributions
    - TMB/HRD/MSI field availability
    
    Returns: Validation report with issues
    """
    # Load dataset
    # Check required fields exist
    # Check PFS_STATUS field format
    # Check outcome distributions
    # Check TMB/HRD/MSI field availability
    # Report issues
```

**Acceptance Criteria**:
- ✅ Checks required fields
- ✅ Validates PFS_STATUS parsing
- ✅ Reports TMB/HRD/MSI field availability
- ✅ Reports outcome distributions
- ✅ Identifies data quality issues

**Estimated Time**: 1 hour

---

## Implementation Plan: My Ownership

### Week 1: Critical Fixes (P0) - **I OWN**

**Goal**: Address weak correlation issue by integrating TMB/HRD/MSI

**My Tasks**:
1. ✅ Create biomarker extractor utility (2 hours)
2. ✅ Update API client to pass tumor_context (1 hour)
3. ✅ Update benchmark script to extract and pass biomarkers (1-2 hours)
4. ✅ Create PFS_STATUS parser utility (1 hour)
5. ✅ Fix classification metrics (1 hour)
6. ✅ Validate integration (1 hour)

**My Deliverables**:
- `biomarker_extractor.py` utility
- Updated `api_client.py` with tumor_context support
- Updated `benchmark_small_test.py` with biomarker extraction
- `pfs_status_parser.py` utility
- Fixed `classification.py` metrics

**My Success Criteria**:
- TMB/HRD/MSI values extracted from TCGA dataset (100% of patients where available)
- Tumor context passed to API (100% of API calls)
- Sporadic gates applied (IO boost, PARP rescue visible in scores)
- Classification metrics compute correctly (events/censored detected)
- Correlation improves: r > 0.2 (from r=0.037-0.278)

**My Time**: 7-8 hours

---

### Week 2: Ranking Precision (P1) - **I OWN**

**Goal**: Fix parsing bugs and improve ranking precision

**My Tasks**:
1. ✅ Create drug normalizer utility (2 hours)
2. ✅ Fix drug ranking metrics (1 hour)
3. ✅ Create data quality checker (1 hour)
4. ✅ Re-run benchmarks with fixes (1 hour)

**My Deliverables**:
- `drug_normalizer.py` utility
- Fixed `drug_ranking.py` metrics
- `data_quality.py` checker

**My Success Criteria**:
- Top-1 accuracy > 20% (from 0%)
- Drug name matching works correctly
- Data quality checker validates datasets

**My Time**: 5 hours

---

### Week 3: Validation Set Bias (P1) - **I OWN**

**Goal**: Test with representative sample to validate correlation improvements

**My Tasks**:
1. ✅ Add stratified sampling functions (1.5 hours)
2. ✅ Update benchmark script with sampling modes (30 min)
3. ✅ Add bias comparison function (30 min)
4. ✅ Compare validation set vs full dataset (1-2 hours)
5. ✅ Re-run 50-patient test with stratified sampling (1 hour)

**My Deliverables**:
- Updated `patient_selection.py` with sampling functions
- Updated `benchmark_small_test.py` with `--mode` argument
- Bias comparison reporting

**My Success Criteria**:
- Stratified sampling implemented
- Validation set characteristics match full dataset (within 20%)
- Bias comparison reports correctly

**My Time**: 4.5-5.5 hours

---

## My Complete Deliverables List

### Utilities I Will Create (5 files)

1. **`biomarker_extractor.py`** (NEW)
   - `extract_tmb_from_patient()`
   - `extract_hrd_from_patient()`
   - `extract_msi_from_patient()`
   - `build_tumor_context()`

2. **`pfs_status_parser.py`** (NEW)
   - `parse_pfs_status()`

3. **`drug_normalizer.py`** (NEW)
   - `normalize_drug_name()`
   - `match_drug()`

4. **`data_quality.py`** (NEW)
   - `validate_dataset_quality()`

5. **`validate_pfs_status.py`** (NEW - Optional)
   - PFS_STATUS field validation script

### Files I Will Update (4 files)

1. **`api_client.py`** (UPDATE)
   - Add `tumor_context` parameter to `call_efficacy_predict()`

2. **`benchmark_small_test.py`** (UPDATE)
   - Add biomarker extraction
   - Add `--mode` argument for sampling
   - Add bias comparison reporting

3. **`classification.py`** (FIX)
   - Use PFS_STATUS parser
   - Fix event/censored detection

4. **`drug_ranking.py`** (FIX)
   - Use drug normalizer
   - Fix Top-1/Top-3 matching

5. **`patient_selection.py`** (UPDATE)
   - Add random sampling
   - Add stratified sampling
   - Add bias comparison function

---

## What I Will NOT Do (Handoff Points)

### ❌ **I DO NOT OWN**: SPE/WIWFM Model Issues

**After My Fixes, If Correlations Still Weak** (r < 0.2):

**Handoff to SPE/WIWFM Team**:
1. **Score Calibration**
   - Efficacy score calibration for outcome prediction
   - **Why**: This is SPE/WIWFM model logic, not benchmark integration

2. **Model Improvements**
   - SPE/WIWFM prediction model enhancements
   - **Why**: This is model architecture, not benchmark integration

**Handoff Criteria**: 
- After TMB/HRD/MSI integration complete
- After all benchmark integration fixes complete
- If correlation still weak (r < 0.2), handoff to SPE/WIWFM team

---

## My Success Metrics

### Phase 1 Success (Week 1)

**My Deliverables**:
- ✅ TMB/HRD/MSI extracted from TCGA dataset (100% of patients where available)
- ✅ Tumor context passed to API (100% of API calls)
- ✅ Sporadic gates applied (IO boost, PARP rescue visible in scores)
- ✅ Classification metrics compute correctly (events/censored detected)

**My Success Criteria**:
- Correlation improves: r > 0.2 (from r=0.037-0.278)
- P-value improves: p < 0.1 (from p=0.4-0.9)
- Classification metrics work (no more "0 events AND 0 censored")

---

### Phase 2 Success (Week 2)

**My Deliverables**:
- ✅ Drug name normalizer working
- ✅ Data quality checker working

**My Success Criteria**:
- Top-1 accuracy > 20% (from 0%)
- Top-3 accuracy > 40% (from 0%)
- Data quality checker validates datasets

---

### Phase 3 Success (Week 3)

**My Deliverables**:
- ✅ Stratified sampling implemented
- ✅ Bias comparison reporting

**My Success Criteria**:
- Validation set characteristics match full dataset (within 20%)
- Bias comparison reports correctly

---

## Risk Assessment: My Ownership

### Risk 1: TCGA Dataset May Not Have TMB/HRD/MSI

**My Mitigation**:
- Check TCGA dataset schema first
- If missing, estimate from mutations (TMB = mutation count / genome size)
- Use disease priors as fallback (reference `tumor_quick_intake.py`)
- Document fallback logic in code comments

**Probability**: Medium  
**Impact**: High  
**My Mitigation Effort**: 2-4 hours (fallback logic)

---

### Risk 2: Correlation May Not Improve Even With TMB/HRD/MSI

**My Response**:
- This is expected - correlation improvement is not guaranteed
- My job is to integrate biomarkers (benchmark integration)
- If no improvement, indicates SPE/WIWFM model issue (handoff to SPE/WIWFM team)
- I will document results and handoff criteria

**Probability**: Medium  
**Impact**: Medium  
**My Action**: Document results, handoff to SPE/WIWFM team if needed

---

### Risk 3: Stratified Sampling May Not Be Representative

**My Mitigation**:
- Compare validation set vs full dataset characteristics
- Use random sampling as alternative
- Document limitations in code comments

**Probability**: Low  
**Impact**: Medium  
**My Mitigation Effort**: 1-2 hours (comparison analysis)

---

## My Timeline

### Week 1: Critical Fixes (P0)
- **Monday-Tuesday**: TMB/HRD/MSI integration (5-6 hours)
- **Wednesday**: PFS_STATUS parser (3 hours)
- **Thursday**: Integration validation (1 hour)
- **Friday**: Re-run benchmarks, document results

**Total**: 9-10 hours

---

### Week 2: Ranking Precision (P1)
- **Monday**: Drug normalizer (3 hours)
- **Tuesday**: Data quality checker (1 hour)
- **Wednesday**: Re-run benchmarks, document results

**Total**: 4 hours

---

### Week 3: Validation Set Bias (P1)
- **Monday**: Patient sampling functions (2 hours)
- **Tuesday**: Bias comparison (1-2 hours)
- **Wednesday**: Re-run benchmarks, document results

**Total**: 3-4 hours

---

## My Complete File Manifest

```
scripts/benchmark/
├── benchmark_common/
│   ├── utils/
│   │   ├── biomarker_extractor.py      # NEW - I CREATE
│   │   ├── pfs_status_parser.py        # NEW - I CREATE
│   │   └── drug_normalizer.py          # NEW - I CREATE
│   ├── patient_selection.py            # UPDATE - I MODIFY
│   ├── data_quality.py                 # NEW - I CREATE
│   ├── api_client.py                   # UPDATE - I MODIFY
│   └── metrics/
│       ├── classification.py           # FIX - I MODIFY
│       └── drug_ranking.py             # FIX - I MODIFY
└── validate_pfs_status.py              # NEW - I CREATE (Optional)
```

**Total Files**: 5 new utilities + 4 updates = 9 files I will create/modify

---

## Conclusion: My Ownership

**What I Own**:
- ✅ Benchmark integration fixes (5 issues)
- ✅ 5 utility modules (biomarker extractor, PFS parser, drug normalizer, data quality, validation)
- ✅ 4 metric/script fixes (API client, classification, drug ranking, patient selection)
- ✅ Integration validation and testing

**What I Do NOT Own**:
- ❌ SPE/WIWFM model improvements
- ❌ Score calibration (after my fixes, if still needed)
- ❌ Model architecture changes

**My Deliverables**: 9 files (5 new, 4 updates)  
**My Timeline**: 3 weeks (16-18 hours total)  
**My Success**: Correlation improves (r > 0.2), classification works, Top-1 accuracy > 20%

---

**Status**: ✅ **READY TO IMPLEMENT** - Clear ownership, deliverables, and success criteria defined

**Estimated Completion**: 16-18 hours (3 weeks, phased approach)
