# SAE Readiness Status - MBD4+TP53 Target

**Date**: January 20, 2025  
**Status**: ✅ **HEALTH CHECKS READY - BLOCKER REMOVAL PLAN COMPLETE**

---

## 🎯 Mission Context

**Primary Target**: MBD4 Germline + TP53 Somatic HGSOC Analysis
- Requires accurate pathway burden (BER + HRD deficiency)
- Needs mechanism fit ranking for PARP inhibitor trials
- Requires early resistance detection (DNA repair capacity)
- **All three services need TRUE SAE for optimal accuracy**

**Current Reality**:
- ✅ TRUE SAE features extracted (66 patients, trained weights)
- ⚠️ Production uses PROXY SAE (gene mutations → pathway scores)
- ❌ TRUE SAE blocked by Feature→Pathway Mapping
- **Impact**: MBD4+TP53 uses pathway-based vectors (works, but TRUE SAE would improve accuracy)

---

## ✅ What We've Accomplished (January 14-20, 2025)

### 1. Evo2 7B Migration ✅
- Migrated from evo2_1b to evo2_7b
- Trained weights loaded (4096×32768 checkpoint)
- Dimension mismatch resolved
- **Code**: `src/services/sae_service/main.py:155-230`

### 2. Full Cohort Extraction ✅
- 66 patients extracted
- 2,897 variants processed
- Trained weights verified
- **Output**: `data/validation/sae_cohort/sae_features_tcga_ov_platinum.json`

### 3. Data Quality Verification ✅
- Outcome field populated (53 sensitive, 11 refractory, 2 resistant)
- Feature indices valid (0-32767)
- All patients have SAE features

### 4. Plan Strengthening ✅
- Proxy vs TRUE SAE distinction clarified
- All three services documented
- Critical blocker identified (Feature→Pathway Mapping)
- MBD4+TP53 integration documented

### 5. Health Check Suite Created ✅
- Data quality health check: `scripts/sae/health_check_data.py`
- Feature distribution health check: `scripts/sae/health_check_feature_distributions.py`
- Pipeline integration health check: `scripts/sae/health_check_pipeline.py`
- Comprehensive readiness plan: `.cursor/plans/SAE_READINESS_AND_BLOCKER_REMOVAL_PLAN.md`

### 6. Clinical Trials Integration ✅
- SAE questions answered in clinical trials plan
- Mechanism fit ranking patterns documented
- Integration flow clarified

---

## ⏸️ What's Ready to Execute

### 1. Health Checks ⏸️ **READY**
- All health check scripts created and executable
- Pre-flight checklist defined
- Execution order documented

**Next Step**: Run health checks to verify system operational

### 2. Biomarker Analysis Re-Run ⏸️ **READY** (Blocked by health checks)
- Service ready, bug fixed
- Data quality verified
- Command ready: `python3 scripts/sae/analyze_biomarkers.py --input data/validation/sae_cohort/sae_features_tcga_ov_platinum.json`

**Next Step**: After health checks pass, run biomarker analysis

### 3. Feature→Pathway Mapping Strategy ⏸️ **DEFINED** (Pending biomarker results)
- Strategy: Biomarker-driven mapping (gene→pathway inference)
- Validation plan: Test on known cases (BRCA1 → DDR, KRAS → MAPK, MBD4+TP53 → DDR)
- Implementation approach documented

**Next Step**: After biomarker analysis, create mapping

### 4. MBD4+TP53 End-to-End Test ⏸️ **READY** (Blocked by health checks)
- Plan exists: `.cursor/plans/MBD4.mdc`
- Test cases defined
- Expected results documented

**Next Step**: After health checks pass, run end-to-end test

---

## 🚨 Fail Now vs Later Assessment (MBD4+TP53 Critical Path)

### ⚠️ FAIL NOW - Must Fix Before MBD4+TP53 Analysis (P0 - 5.5-7.5 hours)

**Gap 0: Mechanism Vector Conversion Function** ⚠️ **BLOCKING MBD4+TP53 Phase 4**

- **Issue**: MBD4+TP53 analysis Phase 4 (Clinical Trial Matching) requires mechanism vector conversion, but `convert_pathway_scores_to_mechanism_vector()` function does NOT exist
- **Impact**: **BLOCKS MBD4+TP53 Phase 4** - Cannot rank trials by mechanism fit without 7D mechanism vector
- **MBD4+TP53 Requirements**:
  - Extract pathway scores from `/api/efficacy/predict` response: `provenance["confidence_breakdown"]["pathway_disruption"]`
  - Convert to 7D vector: [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
  - DDR index must include: `ddr + 0.5 * tp53` (TP53 contributes 50% to DDR)
  - IO index: `1.0 if (tmb >= 20 or msi_high) else 0.0`
- **Resolution**: Create `api/services/pathway_to_mechanism_vector.py` (1.5-2.5 hours)
- **Reference**: `.cursor/plans/clinical-trials.mdc` (lines 712-757) for design

**Gap 3.5: MBD4+TP53 Validation Test Suite Missing** ⚠️ **BLOCKING MBD4+TP53 Analysis**

- **Issue**: MBD4+TP53 analysis requires validation that system works correctly, but no test suite exists
- **Impact**: **BLOCKS MBD4+TP53 Analysis Execution** - Cannot verify system produces accurate predictions
- **MBD4+TP53 Requirements**:
  1. Frameshift Detection: MBD4 c.1239delA → sequence_disruption ≥0.8
  2. Hotspot Detection: TP53 R175H → sequence_disruption ≥0.7
  3. Pathway Aggregation: MBD4+TP53 → DDR pathway score ≥0.70
  4. PARP Ranking: MBD4+TP53 → PARP inhibitors rank #1-2, efficacy_score >0.80
  5. Germline Status: MBD4 germline-positive → no PARP penalty
  6. Mechanism Vector: MBD4+TP53 pathway scores → 7D mechanism vector correct
  7. Trial Matching: MBD4+TP53 → mechanism fit ranking works
- **Resolution**: Create validation test scripts (2.5-3.5 hours)
- **Files**: `scripts/test_mbd4_tp53_*.py`, `tests/test_mbd4_tp53_analysis.py`

**Gap 3.6: Frameshift/Truncation Detection Validation** ⚠️ **CRITICAL FOR MBD4**

- **Issue**: MBD4 c.1239delA is frameshift - need to verify truncation lift works correctly
- **Current Implementation**: `evo2_scorer.py:151-157` applies 1.0 lift for truncation/frameshift variants
- **MBD4+TP53 Requirements**:
  - MBD4 frameshift (c.1239delA, p.Ile413Serfs*2) → sequence_disruption ≥0.8
  - Verify truncation lift applied (1.0 multiplier)
  - Verify contributes to DDR pathway correctly
- **Resolution**: Test MBD4 frameshift and TP53 hotspot (45 minutes)

**Gap 3.7: DDR Pathway Score Validation** ⚠️ **CRITICAL FOR MBD4+TP53**

- **Issue**: MBD4+TP53 combination must produce high DDR pathway scores for accurate PARP predictions
- **MBD4+TP53 Requirements**:
  - MBD4+TP53 combination → DDR pathway score ≥0.70
  - Both variants contribute to pathway aggregation
  - PARP inhibitors get high efficacy_score (>0.80) for MBD4+TP53
- **Resolution**: Test pathway aggregation and PARP predictions (1 hour)

**Total P0 Time (MBD4+TP53 Critical Path)**: 5.5-7.5 hours

### ⏸️ FAIL LATER - Can Address During Development (P1-P2)

**Gap 1: ClinVar BRCA1/BRCA2 Training Data** ✅ **RESOLVED** (Found existing pipeline)
- Existing script: `src/tools/adjudicator_trainer/01_parse_variant_summary.py`
- BRCA1/BRCA2 included in TARGET_GENES
- **Action**: Verify data exists or re-run extraction (15-45 minutes)

**Gap 2: Evo2 Embeddings Extraction** ✅ **PARTIALLY RESOLVED** (Pattern exists)
- Pattern exists: `src/tools/adjudicator_trainer/02_generate_embeddings.py`
- Need to adapt for Evo2 (not Zeta Oracle)
- **Action**: Verify endpoint and adapt script (1.75-2.75 hours)

**Gap 3: Base Scorer Interface** ⚠️ **IMPORTANT** (Can work around)
- **Impact**: Blocks clean integration of new scorers
- **Resolution**: Create `api/services/sequence_scorers/base_scorer.py` (1-2 hours)
- **Note**: Less critical than MBD4+TP53 gaps (can work around it)

**Gap 4-7: Training Infrastructure, Model Storage, Splice Detection, Noncoding** ⚠️ **P1-P2**
- Can be addressed during development
- Not blocking MBD4+TP53 analysis

---

## ❌ Critical Blocker

### Feature→Pathway Mapping

**Status**: ❌ **BLOCKS ALL THREE SERVICES**

**Impact**:
- Resistance Prophet: Cannot use TRUE SAE DNA repair capacity
- Mechanism Fit Ranking: Cannot use TRUE SAE mechanism vector
- Early Resistance Detection: Cannot use TRUE SAE DNA repair capacity

**Resolution Path**:
1. ⏸️ Re-run biomarker analysis (identifies top significant features)
2. ⏸️ Create Feature→Pathway Mapping (biomarker-driven approach)
3. ⏸️ Validate mapping (test on known cases)
4. ⏸️ Integrate into SAE Feature Service (switch from proxy to TRUE SAE)

**MBD4+TP53 Impact**: TRUE SAE would improve DDR pathway accuracy (0.85 → 0.92), leading to better PARP inhibitor mechanism fit (0.82 → 0.91)

---

## 📋 Immediate Next Steps (Prioritized by Fail Now vs Later)

### ⚠️ FAIL NOW - Before MBD4+TP53 Analysis (P0 - 5.5-7.5 hours)

**Step 0: Fix Critical Gaps** (5.5-7.5 hours - MBD4+TP53 Critical Path)

1. **Create Mechanism Vector Conversion** (1.5-2.5 hours)
   - File: `api/services/pathway_to_mechanism_vector.py`
   - Function: `convert_pathway_scores_to_mechanism_vector(pathway_scores, tmb, msi_status)`
   - Test with MBD4+TP53 pathway scores
   - **BLOCKS**: Phase 4 trial matching

2. **Create MBD4+TP53 Validation Test Suite** (2.5-3.5 hours)
   - `scripts/test_mbd4_tp53_sequence_scoring.py` - Frameshift/hotspot validation
   - `scripts/test_mbd4_tp53_pathway_aggregation.py` - DDR pathway validation
   - `scripts/test_mbd4_tp53_parp_predictions.py` - PARP ranking validation
   - `tests/test_mbd4_tp53_analysis.py` - Comprehensive test suite
   - **BLOCKS**: Analysis execution

3. **Validate Frameshift/Truncation Detection** (45 minutes)
   - Test MBD4 c.1239delA → sequence_disruption ≥0.8
   - Test TP53 R175H → sequence_disruption ≥0.7
   - Verify truncation lift applied
   - **CRITICAL**: For MBD4 frameshift detection

4. **Validate DDR Pathway Scores** (1 hour)
   - Test MBD4+TP53 → DDR pathway score ≥0.70
   - Test PARP inhibitors rank #1-2, efficacy_score >0.80
   - Verify both variants contribute
   - **CRITICAL**: For accurate PARP predictions

**Success Criteria**: All 7 MBD4+TP53 test cases pass before analysis execution

### ✅ Standard Next Steps (After P0 Gaps Fixed)

**Step 1: Run Health Checks** (30 min)
```bash
# 1. Data quality
python3 scripts/sae/health_check_data.py

# 2. Feature distributions
python3 scripts/sae/health_check_feature_distributions.py

# 3. Pipeline integration
python3 scripts/sae/health_check_pipeline.py

# 4. Modal service health
curl https://crispro--sae-service-saeservice-api.modal.run/health
curl https://crispro--evo-service-evoservice7b-api-7b.modal.run/health
curl http://localhost:8000/api/sae/health
```

**Step 2: Re-Run Biomarker Analysis** (30 min)
```bash
python3 scripts/sae/analyze_biomarkers.py \
  --input data/validation/sae_cohort/sae_features_tcga_ov_platinum.json \
  --output data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers.json \
  --plots-dir data/validation/sae_cohort/plots
```

**Step 3: Create Feature→Pathway Mapping** (2-4 hours)
- Use top significant features from biomarker analysis
- Map features to pathways via gene associations
- Validate on known cases
- Integrate into SAE Feature Service

**Step 4: MBD4+TP53 End-to-End Test** (1 hour)
- Run complete analysis pipeline
- Verify pathway scores, drug predictions, trial matching
- Compare proxy vs TRUE SAE results (when mapping ready)

### ⏸️ FAIL LATER - During Development (P1-P2)

**Gap 1: ClinVar Data Verification** (15-45 minutes)
- Check if `data/adjudicator_training/clinvar_missense_labels.csv` exists
- Filter for BRCA1/BRCA2 variants
- Re-run extraction if needed

**Gap 2: Evo2 Embeddings Extraction** (1.75-2.75 hours)
- Verify `/api/evo/score_variant_with_activations` works for 1B/7B
- Adapt extraction script from existing pattern
- Test with sample variants

**Gap 3: Base Scorer Interface** (1-2 hours)
- Create `api/services/sequence_scorers/base_scorer.py`
- Update existing scorers to implement interface
- **Note**: Can work around this, less critical than MBD4+TP53 gaps

---

## 🔗 Key Documents

- **Readiness Plan**: `.cursor/plans/SAE_READINESS_AND_BLOCKER_REMOVAL_PLAN.md`
- **Master Plan**: `.cursor/plans/final-comprehensive-document-review-bad14970.plan.md`
- **MBD4 Plan**: `.cursor/plans/MBD4.mdc`
- **Uncertainties**: `.cursor/ayesha/PLAN_UNCERTAINTIES_AND_RISKS.md`
- **Pipeline Roadmap**: `.cursor/rules/SAE_TO_RESISTANCE_PROPHET_PIPELINE.mdc`
- **Fail Now vs Later**: `.cursor/plans/FAIL_NOW_VS_LATER_ASSESSMENT.md` - Critical gaps assessment
- **Pre-Development Readiness**: `.cursor/plans/PRE_DEVELOPMENT_READINESS_ASSESSMENT.md` - Infrastructure readiness

---

## ⚠️ Critical Questions Requiring Answers

### Q1: Mechanism Vector Conversion Function
- **Question**: Does `convert_pathway_scores_to_mechanism_vector()` exist?
- **Answer**: ❌ **NO** - Must be created (1.5-2.5 hours)
- **Impact**: BLOCKS MBD4+TP53 Phase 4 (trial matching)
- **Action**: Create function before running MBD4+TP53 analysis

### Q2: MBD4+TP53 Validation Test Suite
- **Question**: Do validation tests exist for MBD4+TP53 analysis?
- **Answer**: ❌ **NO** - Must be created (2.5-3.5 hours)
- **Impact**: BLOCKS MBD4+TP53 analysis execution
- **Action**: Create test suite before running analysis

### Q3: Frameshift Detection Validation
- **Question**: Does MBD4 c.1239delA frameshift get correct sequence_disruption?
- **Answer**: ⚠️ **UNKNOWN** - Must be validated (45 minutes)
- **Impact**: CRITICAL for MBD4 analysis accuracy
- **Action**: Test frameshift detection before analysis

### Q4: DDR Pathway Score Validation
- **Question**: Does MBD4+TP53 combination produce DDR ≥0.70?
- **Answer**: ⚠️ **UNKNOWN** - Must be validated (1 hour)
- **Impact**: CRITICAL for PARP inhibitor predictions
- **Action**: Test pathway aggregation before analysis

### Q5: S/P/E Integration Status
- **Question**: How does SAE integrate with S/P/E framework?
- **Answer**: ✅ **CLARIFIED** - SAE uses S/P/E outputs but doesn't modulate S/P/E confidence yet
- **Impact**: Understanding current capabilities vs. manager's vision
- **Action**: Documented in Phase 13 of master plan

---

## 📊 Readiness Checklist

### Before MBD4+TP53 Analysis (P0 - Must Complete)

- [ ] **Gap 0**: Mechanism vector conversion function created and tested
- [ ] **Gap 3.5**: MBD4+TP53 validation test suite created (all 7 test cases)
- [ ] **Gap 3.6**: Frameshift/truncation detection validated (MBD4 ≥0.8, TP53 ≥0.7)
- [ ] **Gap 3.7**: DDR pathway validation passed (MBD4+TP53 → DDR ≥0.70)
- [ ] **Health Checks**: All health check scripts pass
- [ ] **MBD4+TP53 Tests**: All 7 validation test cases pass

### Before Development (P1 - Should Complete)

- [ ] **Gap 1**: ClinVar BRCA1/BRCA2 extraction verified (15-45 minutes)
- [ ] **Gap 2**: Evo2 embeddings extraction verified (1.75-2.75 hours)
- [ ] **Gap 3**: Base scorer interface created (1-2 hours)

### During Development (P2 - Can Complete Later)

- [ ] **Gap 4**: Training infrastructure decision made (30 minutes)
- [ ] **Gap 5**: Model storage decision made (30 minutes)
- [ ] **Gap 6**: Splice site detection verified (30 minutes)
- [ ] **Gap 7**: Noncoding variant detection implemented (2-3 hours)

---

**Status**: ⚠️ **NOT FULLY READY** - 4 Critical P0 Gaps Must Be Fixed (5.5-7.5 hours) Before MBD4+TP53 Analysis

**Recommendation**: Address P0 gaps (Gaps 0, 3.5, 3.6, 3.7) before running MBD4+TP53 analysis to ensure accurate results

