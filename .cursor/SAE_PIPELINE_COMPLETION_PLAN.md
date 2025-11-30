# ⚔️ COMPREHENSIVE SAE→RESISTANCE PROPHET PIPELINE COMPLETION PLAN ⚔️

**Date:** January 20, 2025  
**Commander:** Alpha  
**Lead:** Zo  
**Mission:** Complete the SAE→Resistance Prophet pipeline from 80% to 100% operational

---

## 🎯 **EXECUTIVE SUMMARY**

### **Current State:**
- ✅ **Stage 1:** SAE Feature Extraction - **100% COMPLETE** (66 patients, 2,897 variants)
- ⚠️ **Stage 2:** Biomarker Correlation Analysis - **95% COMPLETE** (needs re-run with fixed code)
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

## 📊 **DETAILED STAGE BREAKDOWN**

### **STAGE 1: SAE FEATURE EXTRACTION** ✅ **COMPLETE**

**Status:** ✅ **OPERATIONAL**

**What's Done:**
- ✅ Modal SAE service deployed (H100 GPU)
- ✅ Evo2 activations endpoint operational
- ✅ 66 patients extracted with SAE features
- ✅ 2,897 variants processed
- ✅ Feature indices validated (0-32767 range)
- ✅ All critical bugs fixed (feature index, payload size, checkpoint loading)

**Data File:**
- `data/validation/sae_cohort/sae_features_tcga_ov_platinum.json` (19MB)
- Format: 66 patients, each with variants containing top-64 SAE features

**Key Files:**
- `src/services/sae_service/main.py` (Modal service)
- `scripts/sae/extract_sae_features_cohort.py` (extraction script)

**No Action Required** - Stage 1 is complete and operational.

---

### **STAGE 2: BIOMARKER CORRELATION ANALYSIS** ⚠️ **NEEDS RE-RUN**

**Status:** ⚠️ **SERVICE COMPLETE, ANALYSIS INVALID**

**What's Done:**
- ✅ Biomarker correlation service built (587 lines)
- ✅ All statistical methods implemented (Pearson, Spearman, Chi-square, Cohen's d, FDR, Bootstrap, CV)
- ✅ Analysis script ready (`scripts/sae/analyze_biomarkers.py`)
- ✅ Bug fix applied (line 490: uses `outcome` field correctly)

**What's Wrong:**
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

**Success Criteria:**
- ✅ `significant_features_count` > 0 (expected: 10-50 features)
- ✅ `outcome_distribution` shows both sensitive and resistant patients
- ✅ `top_features` array contains ranked features with correlations
- ✅ All feature indices in valid range (0-32767)
- ✅ Statistical significance: p < 0.05 (FDR-corrected)
- ✅ Effect sizes: Cohen's d >= 0.3

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
      "bootstrap_ci": [0.60, 0.70]
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

3. Build mapping table structure:
```json
{
  "feature_1234": {
    "ddr": 0.6,
    "mapk": 0.1,
    "pi3k": 0.0,
    "vegf": 0.0,
    "her2": 0.0,
    "io": 0.2,
    "efflux": 0.1,
    "rationale": "Feature activated by BRCA1/BRCA2 mutations, associated with DNA repair deficiency",
    "validated": true,
    "pmid": ["12345678", "87654321"]
  },
  ...
}
```

**Deliverable:**
- `api/resources/sae_feature_to_pathway_mapping.json`
- Documentation: `api/resources/sae_feature_to_pathway_mapping_RATIONALE.md`

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

**Current State:**
- ✅ Service exists: `api/services/sae_feature_service.py` (448 lines)
- ✅ Uses proxy SAE features (DNA repair capacity formula)
- ❌ Not using real SAE-derived pathway scores

**Dependencies:**
- ⏸️ Waiting for Stage 3 mapping table

**What Needs to Happen:**

#### **Task 4.1: Load Mapping Table** 📥 **P1**

**Code Changes:**
1. Add mapping table loader to `sae_feature_service.py`
2. Load `sae_feature_to_pathway_mapping.json` on service initialization
3. Cache mapping table in memory

**File:** `api/services/sae_feature_service.py` (add to `__init__`)

**Timeline:** 30 minutes

---

#### **Task 4.2: Update `_compute_sae_diagnostics()`** 🔧 **P1 - CRITICAL PATH**

**Current Implementation (lines 246-267):**
- Uses proxy formula: `0.6×pathway_ddr + 0.2×essentiality_hrr + 0.2×exon_disruption`

**New Implementation:**
1. Extract patient SAE features (from variants)
2. For each feature, look up pathway weights in mapping table
3. Aggregate pathway scores: `pathway_score = Σ(feature_value × pathway_weight)`
4. Compute 7D mechanism vector: `[ddr, mapk, pi3k, vegf, her2, io, efflux]`
5. Compute DNA repair capacity: `ddr_score` (from mechanism vector)

**Code Pattern:**
```python
def _compute_sae_diagnostics(self, patient_sae_features: List[Dict]) -> Dict:
    # Load mapping table (cached)
    mapping = self._load_feature_pathway_mapping()
    
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
    
    # Aggregate features → pathway scores
    for variant in patient_sae_features:
        for feat in variant.get("top_features", []):
            feat_idx = feat["index"]
            feat_value = feat["value"]
            
            if str(feat_idx) in mapping:
                weights = mapping[str(feat_idx)]
                for pathway, weight in weights.items():
                    if pathway != "rationale" and pathway != "validated" and pathway != "pmid":
                        pathway_scores[pathway] += feat_value * weight
    
    # Normalize pathway scores
    total = sum(pathway_scores.values())
    if total > 0:
        pathway_scores = {k: v/total for k, v in pathway_scores.items()}
    
    # Compute DNA repair capacity (DDR score)
    dna_repair_capacity = pathway_scores["ddr"]
    
    # Build mechanism vector (7D)
    mechanism_vector = [
        pathway_scores["ddr"],
        pathway_scores["mapk"],
        pathway_scores["pi3k"],
        pathway_scores["vegf"],
        pathway_scores["her2"],
        pathway_scores["io"],
        pathway_scores["efflux"]
    ]
    
    return {
        "dna_repair_capacity": dna_repair_capacity,
        "mechanism_vector": mechanism_vector,
        "pathway_scores": pathway_scores,
        "sae_method": "true_sae"  # vs "proxy"
    }
```

**File:** `api/services/sae_feature_service.py` (lines 246-267)

**Timeline:** 2 hours

---

#### **Task 4.3: Wire Mechanism Vector into Service** 🔌 **P1**

**Current Implementation (lines 206-214):**
- Uses proxy mechanism vector from pathway aggregation

**New Implementation:**
1. Call `_compute_sae_diagnostics()` to get real SAE pathway scores
2. Use real mechanism vector instead of proxy
3. Add feature flag: `ENABLE_TRUE_SAE_PATHWAYS` (default: false)

**Code Pattern:**
```python
# In compute_sae_features() method
if self.config.get("ENABLE_TRUE_SAE_PATHWAYS", False):
    # Use real SAE pathway scores
    sae_diagnostics = self._compute_sae_diagnostics(patient_sae_features)
    mechanism_vector = sae_diagnostics["mechanism_vector"]
    dna_repair_capacity = sae_diagnostics["dna_repair_capacity"]
else:
    # Use proxy (current implementation)
    mechanism_vector = self._compute_proxy_mechanism_vector(...)
    dna_repair_capacity = self._compute_proxy_dna_repair_capacity(...)
```

**File:** `api/services/sae_feature_service.py` (lines 206-214)

**Timeline:** 1 hour

---

#### **Task 4.4: Add Feature Flag** 🚩 **P1**

**Configuration:**
- Add `ENABLE_TRUE_SAE_PATHWAYS` to `api/config.py`
- Default: `False` (use proxy until validated)
- Enable via environment variable or config file

**File:** `api/config.py`

**Timeline:** 30 minutes

**Total Stage 4 Timeline:** 4 hours

---

### **STAGE 5: RESISTANCE PROPHET VALIDATION** ✅ **READY**

**Status:** ✅ **COMPLETE (using proxy SAE)**

**Current State:**
- ✅ Service complete: `api/services/resistance_prophet_service.py` (689 lines)
- ✅ Fully integrated into `/api/ayesha/complete_care_v2`
- ✅ Uses proxy SAE features currently
- ✅ Ready to use real SAE features after Stage 4

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

## 📋 **COMPLETE TASK LIST**

### **Phase 1: Immediate (Critical Path)** 🔥

1. ✅ **Task 2.1:** Re-Run Biomarker Analysis (2 hours) - **BLOCKING EVERYTHING**
   - Command: `python3 scripts/sae/analyze_biomarkers.py ...`
   - Success: >0 significant features, valid outcome distribution
   - Deliverable: `sae_tcga_ov_platinum_biomarkers.json`

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
   - Deliverable: `sae_feature_to_pathway_mapping.json`

4. **Task 3.3:** Validate Mapping Table (1 hour)
   - Check weights sum to 1.0
   - Verify biological plausibility
   - Test on sample data

### **Phase 4: Service Integration (After Stage 3)** 🔧

5. **Task 4.1:** Load Mapping Table (30 minutes)
   - Add loader to `sae_feature_service.py`
   - Cache in memory

6. **Task 4.2:** Update `_compute_sae_diagnostics()` (2 hours) - **CRITICAL PATH**
   - Implement real SAE pathway score computation
   - Replace proxy formula with mapping table lookup

7. **Task 4.3:** Wire Mechanism Vector (1 hour)
   - Use real SAE mechanism vector
   - Add feature flag integration

8. **Task 4.4:** Add Feature Flag (30 minutes)
   - Add `ENABLE_TRUE_SAE_PATHWAYS` to config

### **Phase 5: Validation (After Stage 4)** ✅

9. **Task 5.1:** Test with Real SAE Features (2 hours)
   - Enable feature flag
   - Compare proxy vs TRUE SAE
   - Validate biology

10. **Task 5.2:** Measure Accuracy Lift (1 hour)
    - Compare accuracy metrics
    - Document improvements

11. **Task 5.3:** Generate Validation Report (1 hour)
    - Complete validation documentation
    - Production deployment recommendations

---

## 🎯 **TIMELINE & DEPENDENCIES**

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
Task 4.1 (Load Mapping) ──────┘
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

## ✅ **SUCCESS CRITERIA**

### **Stage 2 Success:**
- ✅ >0 significant features identified
- ✅ Valid outcome distribution (sensitive + resistant)
- ✅ Statistical significance (p < 0.05, FDR-corrected)
- ✅ Effect sizes (Cohen's d >= 0.3)

### **Stage 3 Success:**
- ✅ Mapping table contains all significant features
- ✅ Pathway weights sum to 1.0 per feature
- ✅ Biological validation (e.g., BRCA1 → DDR)
- ✅ Documentation complete

### **Stage 4 Success:**
- ✅ Real SAE pathway scores computed correctly
- ✅ Mechanism vector uses TRUE SAE (not proxy)
- ✅ Feature flag works (can toggle proxy/TRUE)
- ✅ No regressions in existing functionality

### **Stage 5 Success:**
- ✅ TRUE SAE features work in Resistance Prophet
- ✅ Accuracy maintained or improved vs proxy
- ✅ Biological validation (pathway scores match expectations)
- ✅ Validation report complete

### **Overall Pipeline Success:**
- ✅ All 5 stages operational
- ✅ TRUE SAE features integrated end-to-end
- ✅ Resistance Prophet using real SAE pathway scores
- ✅ Production-ready (behind feature flag)

---

## 🚨 **RISKS & MITIGATIONS**

### **Risk 1: Stage 2 Re-Run Finds No Significant Features**
**Probability:** Low  
**Impact:** High  
**Mitigation:**
- Verify data quality before re-run
- Check outcome field is populated correctly
- Review statistical thresholds (may need adjustment)

### **Risk 2: Pathway Mapping Takes Longer Than Expected**
**Probability:** Medium  
**Impact:** Medium  
**Mitigation:**
- Start with top 20 features (not all 50)
- Use automated literature search tools
- Prioritize high-correlation features first

### **Risk 3: Real SAE Features Don't Improve Accuracy**
**Probability:** Low  
**Impact:** Low  
**Mitigation:**
- Keep proxy as fallback (feature flag)
- Document findings in validation report
- Iterate on mapping table if needed

### **Risk 4: Integration Breaks Existing Functionality**
**Probability:** Low  
**Impact:** High  
**Mitigation:**
- Feature flag ensures proxy remains default
- Comprehensive testing before enabling TRUE SAE
- Rollback plan (disable feature flag)

---

## 📊 **METRICS & TRACKING**

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

## 🎯 **IMMEDIATE NEXT STEPS**

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
**Status:** ✅ **COMPREHENSIVE PLAN COMPLETE - READY FOR EXECUTION**

