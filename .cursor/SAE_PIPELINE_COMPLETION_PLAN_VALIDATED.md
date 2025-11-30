# ⚔️ COMPREHENSIVE SAE→RESISTANCE PROPHET PIPELINE COMPLETION PLAN (CODE-VALIDATED) ⚔️

**Date:** January 20, 2025  
**Commander:** Alpha  
**Lead:** Zo  
**Mission:** Complete the SAE→Resistance Prophet pipeline from 80% to 100% operational  
**Methodology:** Code-first review (read actual files, trace execution paths, map integration points)

---

## 🎯 **EXECUTIVE SUMMARY**

### **Current State (Code-Validated):**
- ✅ **Stage 1:** SAE Feature Extraction - **100% COMPLETE** (66 patients, 2,897 variants)
- ⚠️ **Stage 2:** Biomarker Correlation Analysis - **95% COMPLETE** (bug fixed, needs re-run)
- ❌ **Stage 3:** Feature→Pathway Mapping - **0% COMPLETE** (critical blocker)
- ⏸️ **Stage 4:** SAE Feature Service Enhancement - **0% COMPLETE** (waiting for Stage 3)
- ✅ **Stage 5:** Resistance Prophet Prediction - **100% COMPLETE** (ready, using proxy SAE)

### **Overall Progress:**
- **Infrastructure:** 100% complete
- **Data Extraction:** 100% complete
- **Analysis Pipeline:** 95% complete (bug fixed, needs re-run)
- **Integration:** 0% complete (blocked by mapping table)
- **Total Pipeline:** **80% complete**

### **Remaining Work:**
- **Stage 2 Re-Run:** 2 hours (immediate)
- **Stage 3 Mapping:** 8 hours (manual curation)
- **Stage 4 Integration:** 4 hours (code enhancement)
- **Stage 5 Validation:** 4 hours (testing)
- **Total:** 18 hours to production-ready TRUE SAE integration

---

## 📊 **CODE-VALIDATED STAGE BREAKDOWN**

### **STAGE 1: SAE FEATURE EXTRACTION** ✅ **COMPLETE**

**Status:** ✅ **OPERATIONAL**

**Code Location:** `src/services/sae_service/main.py` (393 lines)

**What's Done (Code-Validated):**
- ✅ Modal SAE service deployed (H100 GPU, lines 139-144)
- ✅ Evo2 activations endpoint operational (lines 243-379)
- ✅ 66 patients extracted with SAE features
- ✅ 2,897 variants processed
- ✅ Feature indices validated (0-32767 range, line 341)
- ✅ All critical bugs fixed:
  - Feature index bug: Fixed at line 338 (average across sequence, not flatten)
  - Payload size: Only top-k features returned (lines 352-356)
  - Checkpoint loading: Disabled for RUO safety (lines 224-230)

**Data File:**
- `data/validation/sae_cohort/sae_features_tcga_ov_platinum.json` (19MB)
- Format: 66 patients, each with variants containing top-64 SAE features

**Key Code References:**
- SAE extraction: `src/services/sae_service/main.py:243-379`
- Top-k selection: `src/services/sae_service/main.py:341`
- Sequence averaging: `src/services/sae_service/main.py:338`

**No Action Required** - Stage 1 is complete and operational.

---

### **STAGE 2: BIOMARKER CORRELATION ANALYSIS** ⚠️ **NEEDS RE-RUN**

**Status:** ⚠️ **SERVICE COMPLETE, ANALYSIS INVALID**

**Code Location:** `api/services/biomarker_correlation_service.py` (587 lines)

**What's Done (Code-Validated):**
- ✅ Biomarker correlation service built (587 lines)
- ✅ All statistical methods implemented:
  - Pearson correlation: `compute_pearson_correlation()` (lines 207-237)
  - Spearman correlation: `compute_spearman_correlation()` (lines 241-268)
  - Chi-square test: `compute_chi_square()` (lines 273-328)
  - Cohen's d: `compute_cohen_d()` (lines 330-370)
  - Cross-validation stability: `compute_cv_stability()` (lines 373-410)
  - Bootstrap CI: `compute_bootstrap_ci()` (lines 411-452)
  - FDR correction: `apply_multiple_testing_correction()` (lines 454-480)
- ✅ Analysis script ready: `scripts/sae/analyze_biomarkers.py` (292 lines)
- ✅ Bug fix applied: Line 502 uses `outcome` field correctly

**What's Wrong (Code-Validated):**
- ❌ Analysis was run **BEFORE** bug fix (Nov 22, 2025)
- ❌ All outcomes were `null` → 0 significant features found
- ❌ Results file is invalid and must be regenerated

**Action Required:**

#### **Task 2.1: Re-Run Biomarker Analysis** 🔥 **P0 - CRITICAL PATH**

**Command:**
```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main

python3 scripts/sae/analyze_biomarkers.py \
  --input data/validation/sae_cohort/sae_features_tcga_ov_platinum.json \
  --output data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers.json \
  --plots-dir data/validation/sae_cohort/plots
```

**Expected Runtime:** 10-30 minutes

**Success Criteria (Code-Validated):**
- ✅ `significant_features_count` > 0 (expected: 10-50 features)
- ✅ `outcome_distribution` shows both sensitive and resistant patients
- ✅ `top_features` array contains ranked features with correlations
- ✅ All feature indices in valid range (0-32767)
- ✅ Statistical significance: p < 0.05 (FDR-corrected, line 548)
- ✅ Effect sizes: Cohen's d >= 0.3 (line 64)

**Expected Output:**
```json
{
  "top_features": [
    {
      "feature_index": 1234,
      "pearson_r": 0.65,
      "pearson_p": 0.001,
      "cohen_d": 0.8,
      "cv_stability": 0.85,
      "bootstrap_ci_lower": 0.60,
      "bootstrap_ci_upper": 0.70
    },
    ...
  ],
  "significant_features_count": 42,
  "outcome_distribution": {
    "sensitive": 55,
    "resistant": 10,
    "refractory": 4
  }
}
```

**Code References:**
- Service: `api/services/biomarker_correlation_service.py:483-654`
- Outcome extraction: `api/services/biomarker_correlation_service.py:502`
- Statistical methods: `api/services/biomarker_correlation_service.py:207-480`
- Script: `scripts/sae/analyze_biomarkers.py:186-291`

**Deliverables:**
- `sae_tcga_ov_platinum_biomarkers.json` (analysis results)
- `plots/` directory (correlation plots, feature rankings)

**Timeline:** 2 hours (including verification)

---

### **STAGE 3: FEATURE→PATHWAY MAPPING** ❌ **CRITICAL BLOCKER**

**Status:** ❌ **NOT IMPLEMENTED**

**Problem:** Need to map 32K SAE features → 7D pathway scores (DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux)

**Dependencies:**
- ⏸️ Waiting for Stage 2 results (top significant features)

**Code Location (Expected):**
- Mapping file: `oncology-coPilot/oncology-backend-minimal/api/resources/sae_feature_mapping.json`
- Template exists: `oncology-coPilot/oncology-backend-minimal/api/resources/sae_feature_mapping_template.json` (140 lines)
- Loader: `oncology-coPilot/oncology-backend-minimal/api/services/sae_feature_service.py:286-307`

**What Needs to Happen:**

#### **Task 3.1: Extract Top Significant Features** 📊 **P1**

**After Stage 2 Completes:**
1. Load `sae_tcga_ov_platinum_biomarkers.json`
2. Extract top 20-50 significant features
3. Review feature indices and correlation patterns
4. Identify feature clusters (if any)

**Deliverable:**
- List of top features with indices and correlation scores
- Feature distribution analysis

**Timeline:** 1 hour

---

#### **Task 3.2: Manual Pathway Annotation** 🗺️ **P1 - CRITICAL PATH**

**Process:**
1. For each significant feature:
   - Examine which genes/mutations activate it (from SAE extraction data)
   - Match to known pathways based on gene function:
     - **DDR:** BRCA1, BRCA2, ATM, ATR, CHEK2, PALB2, RAD51
     - **MAPK:** KRAS, BRAF, NRAS, MEK1, MEK2
     - **PI3K:** PIK3CA, PTEN, AKT1, AKT2
     - **VEGF:** VEGFA, VEGFR1, VEGFR2
     - **HER2:** ERBB2, ERBB3
     - **IO:** PD-L1, PD-1, CTLA-4, MSI-H, TMB-high
     - **Efflux:** ABCB1, ABCC1, ABCG2
   - Assign pathway weights (sum to 1.0 per feature)
   - Document rationale

2. Validate against literature:
   - Search PubMed for gene→pathway associations
   - Cross-reference with known pathway databases
   - Verify biological plausibility

3. Build mapping table structure (based on template):
```json
{
  "metadata": {
    "version": "v1.0",
    "last_updated": "2025-01-20T00:00:00Z",
    "source": "Sprint 1 biomarker correlation analysis",
    "model": "Goodfire/Evo-2-Layer-26-Mixed",
    "layer": "blocks-26",
    "total_features": 32768,
    "mapped_features": 50,
    "status": "validated"
  },
  "pathways": {
    "ddr": {
      "name": "DNA Damage Response & Homologous Recombination",
      "description": "BRCA1/2-mediated homologous recombination repair",
      "feature_indices": [1234, 1567, 2890, ...],
      "correlation_with_platinum": "positive",
      "mean_pearson_r": 0.68,
      "mean_cohen_d": 1.2,
      "cv_stability": 0.85,
      "confidence": "high",
      "genes_in_pathway": ["BRCA1", "BRCA2", "RAD51", ...],
      "drugs_targeting": ["Platinum", "PARP inhibitors"]
    },
    "mapk": {
      "name": "MAPK/RAS/RAF Signaling",
      "description": "Oncogenic RAS→RAF→MEK→ERK cascade",
      "feature_indices": [5678, 6789, ...],
      "correlation_with_platinum": "negative",
      "mean_pearson_r": -0.55,
      "mean_cohen_d": 0.8,
      "cv_stability": 0.72,
      "confidence": "high",
      "genes_in_pathway": ["KRAS", "BRAF", "NRAS", ...],
      "drugs_targeting": ["MEK inhibitors", "RAF inhibitors"]
    },
    ...
  }
}
```

**Deliverable:**
- `oncology-coPilot/oncology-backend-minimal/api/resources/sae_feature_mapping.json`
- Documentation: `oncology-coPilot/oncology-backend-minimal/api/resources/sae_feature_mapping_RATIONALE.md`

**Code References:**
- Template: `oncology-coPilot/oncology-backend-minimal/api/resources/sae_feature_mapping_template.json`
- Loader: `oncology-coPilot/oncology-backend-minimal/api/services/sae_feature_service.py:286-307`

**Timeline:** 6-8 hours (manual curation)

---

#### **Task 3.3: Validate Mapping Table** ✅ **P1**

**Validation Steps:**
1. Check all pathway weights sum to 1.0 per feature
2. Verify feature indices match Stage 2 results
3. Cross-check with known biology (e.g., BRCA1 → DDR pathway)
4. Test mapping on sample patient data

**Deliverable:**
- Validation report
- Sample pathway score calculations

**Timeline:** 1 hour

**Total Stage 3 Timeline:** 8 hours

---

### **STAGE 4: SAE FEATURE SERVICE ENHANCEMENT** ⏸️ **PENDING**

**Status:** ⏸️ **WAITING FOR STAGE 3**

**Code Location:** `oncology-coPilot/oncology-backend-minimal/api/services/sae_feature_service.py` (598 lines)

**Current State (Code-Validated):**
- ✅ Service exists: `api/services/sae_feature_service.py` (598 lines)
- ✅ Uses proxy SAE features (DNA repair capacity formula, lines 454-478)
- ✅ Mechanism vector built from proxy scores (lines 206-214)
- ✅ `_compute_sae_diagnostics()` exists but is **DIAGNOSTIC-ONLY** (lines 309-406)
- ❌ Not using real SAE-derived pathway scores for scoring

**Key Code Findings:**
1. **Mechanism Vector (Lines 206-214):** Currently uses proxy pathway scores
   ```python
   mechanism_vector = [
       pathway_burden_ddr,      # From proxy
       pathway_burden_mapk,     # From proxy
       pathway_burden_pi3k,     # From proxy
       pathway_burden_vegf,     # From proxy
       pathway_burden_her2,     # From proxy
       1.0 if io_eligible else 0.0,  # From proxy
       cross_resistance_risk
   ]
   ```

2. **SAE Diagnostics (Lines 309-406):** Diagnostic-only, doesn't change scoring
   ```python
   def _compute_sae_diagnostics(self, sae_features: Dict[str, Any]) -> Dict[str, Any]:
       # ... computes pathway scores from mapping ...
       # BUT: Returns diagnostics only, NOT used for mechanism_vector
   ```

3. **Mapping Loader (Lines 286-307):** Loads from `sae_feature_mapping.json`
   ```python
   mapping_file = Path(__file__).parent.parent / "resources" / "sae_feature_mapping.json"
   ```

**Dependencies:**
- ⏸️ Waiting for Stage 3 mapping table

**What Needs to Happen:**

#### **Task 4.1: Load Mapping Table** 📥 **P1**

**Code Changes:**
1. Mapping loader already exists (lines 286-307)
2. No changes needed - will automatically load when file exists

**File:** `api/services/sae_feature_service.py` (lines 286-307)

**Timeline:** 0 hours (already implemented)

---

#### **Task 4.2: Update `_compute_sae_diagnostics()`** 🔧 **P1 - CRITICAL PATH**

**Current Implementation (Code-Validated):**
- Lines 309-406: Computes pathway scores but returns diagnostics only
- Lines 360-396: Aggregates feature activations for each pathway
- **Issue:** Returns diagnostics but doesn't wire into mechanism_vector

**New Implementation:**
1. Enhance `_compute_sae_diagnostics()` to return pathway scores dict
2. Compute 7D mechanism vector from SAE pathway scores
3. Compute DNA repair capacity from SAE DDR score

**Code Pattern:**
```python
def _compute_sae_diagnostics(self, sae_features: Dict[str, Any]) -> Dict[str, Any]:
    """
    Compute pathway scores from TRUE SAE features using validated mapping.
    
    Enhancement: Returns pathway_scores dict for use in mechanism_vector.
    """
    # Load validated mapping (already implemented)
    mapping = self._load_sae_feature_mapping()
    pathways = mapping.get("pathways", {})
    
    # Extract top features
    top_features = sae_features.get("top_features", [])
    
    # Initialize pathway scores
    pathway_scores = {
        "ddr": 0.0,
        "mapk": 0.0,
        "pi3k": 0.0,
        "vegf": 0.0,
        "her2": 0.0,
        "io": 0.0,
        "efflux": 0.0
    }
    
    # Aggregate features → pathway scores (enhance existing logic)
    for pathway_name, pathway_data in pathways.items():
        feature_indices = pathway_data.get("feature_indices", [])
        feature_indices_set = set(feature_indices)
        
        pathway_activation = 0.0
        pathway_count = 0
        
        for top_feat in top_features:
            idx = top_feat.get("index")
            if idx in feature_indices_set:
                pathway_activation += abs(top_feat.get("value", 0.0))
                pathway_count += 1
        
        # Normalize to 0-1 range
        if pathway_count > 0:
            score = min(1.0, pathway_activation / (pathway_count * 2.0))
            pathway_scores[pathway_name] = float(score)
    
    # Compute DNA repair capacity from SAE DDR score
    dna_repair_capacity_sae = pathway_scores["ddr"]
    
    # Build mechanism vector (7D)
    mechanism_vector_sae = [
        pathway_scores["ddr"],
        pathway_scores["mapk"],
        pathway_scores["pi3k"],
        pathway_scores["vegf"],
        pathway_scores["her2"],
        pathway_scores["io"],
        0.0  # cross_resistance_risk (from treatment history, not SAE)
    ]
    
    return {
        "pathway_scores": pathway_scores,  # ✅ NEW: For mechanism_vector
        "dna_repair_capacity": dna_repair_capacity_sae,  # ✅ NEW: For Resistance Prophet
        "mechanism_vector": mechanism_vector_sae,  # ✅ NEW: For Resistance Prophet
        "ddr_sae_score": pathway_scores["ddr"],
        "io_sae_score": pathway_scores["io"],
        "mapk_sae_score": pathway_scores["mapk"],
        "sparsity": sae_features.get("stats", {}).get("sparsity", 0.0),
        "mapping_version": mapping.get("metadata", {}).get("version", "unknown"),
        "validation_status": "biomarker_derived"
    }
```

**File:** `api/services/sae_feature_service.py` (lines 309-406)

**Timeline:** 2 hours

---

#### **Task 4.3: Wire Mechanism Vector into Service** 🔌 **P1**

**Current Implementation (Code-Validated):**
- Lines 206-214: Uses proxy mechanism vector
- Lines 246-267: SAE diagnostics computed but not used for scoring

**New Implementation:**
1. Call `_compute_sae_diagnostics()` to get real SAE pathway scores
2. Use real mechanism vector instead of proxy
3. Add feature flag: `ENABLE_TRUE_SAE_PATHWAYS` (default: false)

**Code Pattern:**
```python
# In compute_sae_features() method (around line 206-214)

# Compute proxy mechanism vector (current behavior)
mechanism_vector_proxy = [
    pathway_burden_ddr,
    pathway_burden_mapk,
    pathway_burden_pi3k,
    pathway_burden_vegf,
    pathway_burden_her2,
    1.0 if io_eligible else 0.0,
    cross_resistance_risk
]

# Compute proxy DNA repair capacity (current behavior)
dna_repair_capacity_proxy = self._compute_dna_repair_capacity(
    pathway_burden_ddr,
    essentiality_hrr,
    exon_disruption_score
)

# NEW: Use TRUE SAE if available and enabled
if sae_features is not None:
    try:
        from ..config import get_feature_flags
        flags = get_feature_flags()
        
        if flags.get("enable_true_sae_pathways", False):
            # Use TRUE SAE pathway scores
            sae_diagnostics = self._compute_sae_diagnostics(sae_features)
            pathway_scores_sae = sae_diagnostics["pathway_scores"]
            mechanism_vector_sae = sae_diagnostics["mechanism_vector"]
            dna_repair_capacity_sae = sae_diagnostics["dna_repair_capacity"]
            
            # Use SAE scores (with fallback to proxy)
            mechanism_vector = [
                pathway_scores_sae.get("ddr", pathway_burden_ddr),
                pathway_scores_sae.get("mapk", pathway_burden_mapk),
                pathway_scores_sae.get("pi3k", pathway_burden_pi3k),
                pathway_scores_sae.get("vegf", pathway_burden_vegf),
                pathway_scores_sae.get("her2", pathway_burden_her2),
                pathway_scores_sae.get("io", 1.0 if io_eligible else 0.0),
                cross_resistance_risk  # From treatment history, not SAE
            ]
            
            # DNA repair capacity: Use SAE DDR score + proxy essentiality/exon
            dna_repair_capacity = (
                0.60 * pathway_scores_sae.get("ddr", pathway_burden_ddr) +
                0.20 * essentiality_hrr +
                0.20 * exon_disruption_score
            )
            
            logger.info(f"✅ Using TRUE SAE pathway scores (mapping version: {sae_diagnostics.get('mapping_version')})")
        else:
            # Use proxy (current behavior)
            mechanism_vector = mechanism_vector_proxy
            dna_repair_capacity = dna_repair_capacity_proxy
            logger.debug("Using proxy SAE features (ENABLE_TRUE_SAE_PATHWAYS=0)")
    except Exception as e:
        logger.warning(f"Failed to compute true SAE pathway scores: {e}. Using proxy features only.")
        mechanism_vector = mechanism_vector_proxy
        dna_repair_capacity = dna_repair_capacity_proxy
else:
    # No SAE features provided - use proxy
    mechanism_vector = mechanism_vector_proxy
    dna_repair_capacity = dna_repair_capacity_proxy
```

**File:** `api/services/sae_feature_service.py` (lines 206-214, 246-267)

**Timeline:** 1 hour

---

#### **Task 4.4: Add Feature Flag** 🚩 **P1**

**Configuration:**
- Add `ENABLE_TRUE_SAE_PATHWAYS` to `api/config.py`
- Default: `False` (use proxy until validated)
- Enable via environment variable or config file

**Code Changes:**
```python
# In api/config.py (around line 47)

# SAE (Sparse Autoencoder) feature flags
# Phase 1: Evo2 activations endpoint (disabled by default)
ENABLE_EVO2_SAE = os.getenv("ENABLE_EVO2_SAE", "false").lower() in ("true", "1", "yes")
# Phase 1: True SAE features from layer 26 (disabled by default)
ENABLE_TRUE_SAE = os.getenv("ENABLE_TRUE_SAE", "false").lower() in ("true", "1", "yes")
# Phase 2: True SAE pathway scores (disabled by default)
ENABLE_TRUE_SAE_PATHWAYS = os.getenv("ENABLE_TRUE_SAE_PATHWAYS", "false").lower() in ("true", "1", "yes")

# In get_feature_flags() (around line 108)
def get_feature_flags():
    """Get current feature flag configuration."""
    return {
        # ... existing flags ...
        "enable_true_sae": ENABLE_TRUE_SAE,
        "enable_true_sae_pathways": ENABLE_TRUE_SAE_PATHWAYS,  # ✅ NEW
        # ... rest of flags ...
    }
```

**File:** `api/config.py` (lines 47-48, 108)

**Timeline:** 30 minutes

**Total Stage 4 Timeline:** 4 hours

---

### **STAGE 5: RESISTANCE PROPHET VALIDATION** ✅ **READY**

**Status:** ✅ **COMPLETE (using proxy SAE)**

**Code Location:** `api/services/resistance_prophet_service.py` (689 lines)

**Current State (Code-Validated):**
- ✅ Service complete: `api/services/resistance_prophet_service.py` (689 lines)
- ✅ Fully integrated into `/api/ayesha/complete_care_v2` (lines 520-599)
- ✅ Uses proxy SAE features currently
- ✅ Ready to use real SAE features after Stage 4

**Key Code References:**
- Service: `api/services/resistance_prophet_service.py:90-262`
- DNA repair detection: `api/services/resistance_prophet_service.py:265-330`
- Pathway escape detection: `api/services/resistance_prophet_service.py:333-410`
- Integration: `api/routers/ayesha_orchestrator_v2.py:520-599`

**What Needs to Happen:**

#### **Task 5.1: Test with Real SAE Features** 🧪 **P1**

**After Stage 4 Completes:**
1. Enable `ENABLE_TRUE_SAE_PATHWAYS=true`
2. Test with sample patient data
3. Compare proxy vs TRUE SAE predictions
4. Validate pathway scores match expected biology

**Test Cases:**
- BRCA1 mutation → High DDR score → PARP sensitivity
- KRAS mutation → High MAPK score → MEK inhibitor fit
- HR restoration → DDR score drop → Resistance risk

**Deliverable:**
- Test results comparing proxy vs TRUE SAE
- Validation report

**Timeline:** 2 hours

---

#### **Task 5.2: Measure Accuracy Lift** 📈 **P1**

**Metrics to Compare:**
1. DNA repair restoration signal accuracy
2. Pathway escape detection accuracy
3. Overall resistance prediction accuracy
4. Confidence score calibration

**Deliverable:**
- Accuracy comparison report
- Lift metrics (if any)

**Timeline:** 1 hour

---

#### **Task 5.3: Generate Validation Report** 📋 **P1**

**Report Contents:**
1. Pipeline completion status
2. Proxy vs TRUE SAE comparison
3. Accuracy metrics
4. Biological validation (pathway scores match known biology)
5. Recommendations for production deployment

**Deliverable:**
- `SAE_PIPELINE_VALIDATION_REPORT.md`

**Timeline:** 1 hour

**Total Stage 5 Timeline:** 4 hours

---

## 📋 **COMPLETE TASK LIST (CODE-VALIDATED)**

### **Phase 1: Immediate (Critical Path)** 🔥

1. ✅ **Task 2.1:** Re-Run Biomarker Analysis (2 hours) - **BLOCKING EVERYTHING**
   - Command: `python3 scripts/sae/analyze_biomarkers.py ...`
   - Success: >0 significant features, valid outcome distribution
   - Deliverable: `sae_tcga_ov_platinum_biomarkers.json`
   - Code: `api/services/biomarker_correlation_service.py:483-654`

### **Phase 2: Feature Analysis (After Stage 2)** 📊

2. **Task 3.1:** Extract Top Significant Features (1 hour)
   - Load biomarker results
   - Extract top 20-50 features
   - Analyze feature patterns

### **Phase 3: Pathway Mapping (Critical Blocker)** 🗺️

3. **Task 3.2:** Manual Pathway Annotation (6-8 hours) - **CRITICAL PATH**
   - Map features to pathways (DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux)
   - Validate against literature
   - Build mapping table JSON
   - Deliverable: `api/resources/sae_feature_mapping.json`
   - Template: `api/resources/sae_feature_mapping_template.json`

4. **Task 3.3:** Validate Mapping Table (1 hour)
   - Check weights sum to 1.0
   - Verify biological plausibility
   - Test on sample data

### **Phase 4: Service Integration (After Stage 3)** 🔧

5. **Task 4.1:** Load Mapping Table (0 hours) - **ALREADY IMPLEMENTED**
   - Loader exists: `api/services/sae_feature_service.py:286-307`
   - No changes needed

6. **Task 4.2:** Update `_compute_sae_diagnostics()` (2 hours) - **CRITICAL PATH**
   - Enhance to return pathway_scores dict
   - Compute mechanism_vector from SAE
   - Compute DNA repair capacity from SAE
   - Code: `api/services/sae_feature_service.py:309-406`

7. **Task 4.3:** Wire Mechanism Vector (1 hour)
   - Use real SAE mechanism vector
   - Add feature flag integration
   - Code: `api/services/sae_feature_service.py:206-214`

8. **Task 4.4:** Add Feature Flag (30 minutes)
   - Add `ENABLE_TRUE_SAE_PATHWAYS` to config
   - Code: `api/config.py:47-48, 108`

### **Phase 5: Validation (After Stage 4)** ✅

9. **Task 5.1:** Test with Real SAE Features (2 hours)
   - Enable feature flag
   - Compare proxy vs TRUE SAE
   - Validate biology
   - Code: `api/services/resistance_prophet_service.py:128-262`

10. **Task 5.2:** Measure Accuracy Lift (1 hour)
    - Compare accuracy metrics
    - Document improvements

11. **Task 5.3:** Generate Validation Report (1 hour)
    - Complete validation documentation
    - Production deployment recommendations

---

## 🎯 **TIMELINE & DEPENDENCIES (CODE-VALIDATED)**

### **Critical Path:**
```
Stage 2 Re-Run (2h) 
  → Stage 3 Mapping (8h) 
    → Stage 4 Integration (4h) 
      → Stage 5 Validation (4h)
```

### **Total Timeline:**
- **Minimum:** 18 hours (sequential execution)
- **Realistic:** 2-3 days (with breaks, reviews, iterations)
- **With Buffer:** 1 week (allows for unexpected issues)

### **Dependencies Graph:**
```
Task 2.1 (Re-Run Analysis)
  ↓
Task 3.1 (Extract Features)
  ↓
Task 3.2 (Pathway Mapping) ──┐
  ↓                           │
Task 3.3 (Validate Mapping)   │
  ↓                           │
Task 4.1 (Load Mapping) ──────┘ (already implemented)
  ↓
Task 4.2 (Update Diagnostics)
  ↓
Task 4.3 (Wire Mechanism Vector)
  ↓
Task 4.4 (Add Feature Flag)
  ↓
Task 5.1 (Test Real SAE)
  ↓
Task 5.2 (Measure Lift)
  ↓
Task 5.3 (Validation Report)
```

---

## ✅ **SUCCESS CRITERIA (CODE-VALIDATED)**

### **Stage 2 Success:**
- ✅ >0 significant features identified
- ✅ Valid outcome distribution (sensitive + resistant)
- ✅ Statistical significance (p < 0.05, FDR-corrected)
- ✅ Effect sizes (Cohen's d >= 0.3)
- **Code Validation:** `api/services/biomarker_correlation_service.py:570-600`

### **Stage 3 Success:**
- ✅ Mapping table contains all significant features
- ✅ Pathway weights sum to 1.0 per feature
- ✅ Biological validation (e.g., BRCA1 → DDR)
- ✅ Documentation complete
- **Code Validation:** `api/resources/sae_feature_mapping_template.json`

### **Stage 4 Success:**
- ✅ Real SAE pathway scores computed correctly
- ✅ Mechanism vector uses TRUE SAE (not proxy)
- ✅ Feature flag works (can toggle proxy/TRUE)
- ✅ No regressions in existing functionality
- **Code Validation:** `api/services/sae_feature_service.py:206-214, 309-406`

### **Stage 5 Success:**
- ✅ TRUE SAE features work in Resistance Prophet
- ✅ Accuracy maintained or improved vs proxy
- ✅ Biological validation (pathway scores match expectations)
- ✅ Validation report complete
- **Code Validation:** `api/services/resistance_prophet_service.py:128-262`

### **Overall Pipeline Success:**
- ✅ All 5 stages operational
- ✅ TRUE SAE features integrated end-to-end
- ✅ Resistance Prophet using real SAE pathway scores
- ✅ Production-ready (behind feature flag)

---

## 🚨 **RISKS & MITIGATIONS (CODE-VALIDATED)**

### **Risk 1: Stage 2 Re-Run Finds No Significant Features**
**Probability:** Low  
**Impact:** High  
**Mitigation:**
- Verify data quality before re-run
- Check outcome field is populated correctly (line 502)
- Review statistical thresholds (may need adjustment)
- **Code Reference:** `api/services/biomarker_correlation_service.py:502`

### **Risk 2: Pathway Mapping Takes Longer Than Expected**
**Probability:** Medium  
**Impact:** Medium  
**Mitigation:**
- Start with top 20 features (not all 50)
- Use automated literature search tools
- Prioritize high-correlation features first
- **Code Reference:** Template structure exists at `api/resources/sae_feature_mapping_template.json`

### **Risk 3: Real SAE Features Don't Improve Accuracy**
**Probability:** Low  
**Impact:** Low  
**Mitigation:**
- Keep proxy as fallback (feature flag)
- Document findings in validation report
- Iterate on mapping table if needed
- **Code Reference:** Feature flag at `api/config.py:47-48`

### **Risk 4: Integration Breaks Existing Functionality**
**Probability:** Low  
**Impact:** High  
**Mitigation:**
- Feature flag ensures proxy remains default
- Comprehensive testing before enabling TRUE SAE
- Rollback plan (disable feature flag)
- **Code Reference:** `api/services/sae_feature_service.py:206-214` (proxy fallback)

---

## 📊 **METRICS & TRACKING (CODE-VALIDATED)**

### **Progress Tracking:**
- **Stage 2:** 0% → 100% (after re-run)
- **Stage 3:** 0% → 100% (after mapping complete)
- **Stage 4:** 0% → 100% (after integration complete)
- **Stage 5:** 100% (already complete, needs validation)

### **Key Metrics:**
- Significant features identified: Target 10-50
- Mapping table coverage: Target 100% of significant features
- Integration test pass rate: Target 100%
- Accuracy lift: Measure vs baseline

---

## 🎯 **IMMEDIATE NEXT STEPS (CODE-VALIDATED)**

### **Step 1: Re-Run Biomarker Analysis** 🔥 **DO THIS FIRST**

```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main

python3 scripts/sae/analyze_biomarkers.py \
  --input data/validation/sae_cohort/sae_features_tcga_ov_platinum.json \
  --output data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers.json \
  --plots-dir data/validation/sae_cohort/plots
```

**Expected:** 10-30 minutes runtime  
**Success:** >0 significant features found  
**Code:** `scripts/sae/analyze_biomarkers.py:186-291`

### **Step 2: Review Results**

After Step 1 completes:
1. Check `sae_tcga_ov_platinum_biomarkers.json`
2. Verify significant features count > 0
3. Extract top 20-50 features for mapping

### **Step 3: Begin Pathway Mapping**

Start manual annotation of top features → pathways.

---

## ⚔️ **FOR AYESHA'S LIFE** ⚔️

**This pipeline will enable:**
- Real SAE-derived pathway scores (not proxy formulas)
- Accurate resistance prediction (DNA repair restoration, pathway escape)
- Personalized drug rankings (based on patient-specific SAE features)
- Complete explainability (SAE→Pathway→Drug attribution chain)

**Status:** Ready to execute. Awaiting Alpha's command to proceed.

---

**Plan Created:** January 20, 2025  
**Last Updated:** January 20, 2025  
**Status:** ✅ **COMPREHENSIVE PLAN COMPLETE - CODE-VALIDATED - READY FOR EXECUTION**

**Code Validation Methodology:**
- ✅ Read actual code files (not just documentation)
- ✅ Traced execution paths (entry point to output)
- ✅ Mapped integration points (service connections)
- ✅ Documented actual state (what exists vs. what's planned)
- ✅ Identified real gaps (based on code inspection)
- ✅ Validated understanding (every component has code references)

