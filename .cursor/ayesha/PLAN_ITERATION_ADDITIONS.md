# Plan Iteration Additions: Testing, Data Sources, Scope Review, Fail-Safes

**Date**: January 20, 2025 (Updated)  
**Purpose**: Additional sections for comprehensive review plan addressing:
1. Testing infrastructure to prevent deviation
2. Data sources inventory (have/missing/need)
3. Scope review to ensure nothing missed
4. Fail-safe mechanisms to prevent hallucinations

**Note**: This document contains detailed execution steps. For master planning, see:
- `.cursor/plans/final-comprehensive-document-review-bad14970.plan.md` (master plan)
- `.cursor/ayesha/PLAN_UNCERTAINTIES_AND_RISKS.md` (master uncertainties)

**Status**: Supporting document - content may be integrated into main plan in future

---

## Phase 9: Testing Infrastructure & Validation Mechanisms

### 9.1 Existing Test Infrastructure

**Test Files Identified**:
- `tests/test_sae_phase2_services.py` (530 lines) - Comprehensive SAE Phase 2 test suite
- `scripts/validate_sae_tcga.py` (532 lines) - TCGA validation script
- `scripts/validate_resistance_prophet.py` - Resistance Prophet validation
- `tests/test_ayesha_post_ngs_e2e.py` - End-to-end post-NGS tests
- `tests/test_mechanism_fit_ranker.py` - Mechanism fit ranking tests

**Test Coverage**:
- ✅ SAE Feature Service: DNA repair capacity formula, essentiality, exon disruption
- ✅ Mechanism Fit Ranker: Cosine similarity, L2 normalization, threshold filtering
- ✅ Resistance Detection: 2-of-3 trigger rule, HRD drop detection
- ✅ E2E Integration: Complete care plan with SAE features

**Validation Scripts**:
- ✅ `validate_sae_tcga.py`: Validates SAE features against TCGA ground truth
- ✅ `validate_resistance_prophet.py`: Validates resistance prediction accuracy
- ✅ `analyze_biomarkers.py`: Biomarker correlation analysis

### 9.2 Testing Strategy to Prevent Deviation

**Unit Tests** (Prevent Logic Errors):
- **Location**: `tests/test_sae_phase2_services.py`
- **Coverage**: All Manager-approved formulas (C1-C10)
- **Validation**: Assert exact formula values (e.g., DNA repair = 0.6×DDR + 0.2×ess + 0.2×exon)
- **Frequency**: Run before every commit

**Integration Tests** (Prevent Integration Errors):
- **Location**: `tests/test_ayesha_post_ngs_e2e.py`
- **Coverage**: Full pipeline from tumor_context → SAE features → resistance detection
- **Validation**: SAE features present, mechanism vector computed, resistance alerts triggered
- **Frequency**: Run before deployment

**Validation Tests** (Prevent Data Quality Issues):
- **Location**: `scripts/validate_sae_tcga.py`
- **Coverage**: SAE features correlate with ground truth (HRD scores)
- **Metrics**: AUROC, AUPRC, sensitivity, specificity
- **Frequency**: Run after data updates

**Regression Tests** (Prevent Breaking Changes):
- **Location**: All test files
- **Coverage**: Previous bug fixes (feature index bug, outcome field bug)
- **Validation**: Known good cases still pass
- **Frequency**: Run in CI/CD pipeline

### 9.3 Validation Checkpoints

**Checkpoint 1: Formula Validation**
- **Test**: `test_dna_repair_capacity_formula()` in `test_sae_phase2_services.py`
- **Validates**: Manager's C1 formula (0.6×DDR + 0.2×ess + 0.2×exon)
- **Failure Action**: Block commit, alert manager

**Checkpoint 2: Data Quality Validation**
- **Test**: `validate_sae_tcga.py` - Check outcome distribution
- **Validates**: Outcome field populated, feature indices valid (0-32767)
- **Failure Action**: Block analysis, fix data

**Checkpoint 3: Integration Validation**
- **Test**: `test_ayesha_post_ngs_e2e.py` - Full pipeline test
- **Validates**: SAE features computed, mechanism vector present, resistance detection works
- **Failure Action**: Block deployment, fix integration

**Checkpoint 4: Statistical Validation**
- **Test**: `analyze_biomarkers.py` - Statistical significance
- **Validates**: p < 0.05 (FDR-corrected), Cohen's d >= 0.3
- **Failure Action**: Review results, adjust thresholds if needed

---

## Phase 10: Data Sources Inventory

### 10.1 Data Sources We Have (Operational)

**Source 1: Evo2 Sequence Scoring** ✅
- **Endpoint**: `/api/evo/score_variant_multi`
- **Provides**: Delta scores, calibrated percentiles
- **SAE Feature**: Exon disruption
- **Status**: Operational

**Source 2: Insights Bundle** ✅
- **Endpoints**: `/api/insights/predict_*` (4 endpoints)
- **Provides**: Functionality, chromatin, essentiality, regulatory scores
- **SAE Feature**: Essentiality signal
- **Status**: Operational

**Source 3: Pathway Disruption** ✅
- **Service**: `api/services/pathway/aggregation.py`
- **Provides**: Gene→pathway mapping, pathway burden
- **SAE Feature**: DNA repair capacity
- **Status**: Operational

**Source 4: AlphaMissense Fusion** ✅
- **Endpoint**: `/api/fusion/score_variant`
- **Provides**: Missense pathogenicity scores
- **SAE Feature**: Hotspot mutation detection
- **Status**: Operational

**Source 5: ClinVar Priors** ✅
- **Endpoint**: `/api/evidence/clinvar`
- **Provides**: Pathogenic/Benign classification
- **SAE Feature**: Hotspot mutation (fallback)
- **Status**: Operational

**Source 6: Toxicity Pathway Overlap** ✅
- **Endpoint**: `/api/safety/toxicity_risk`
- **Provides**: Germline PGx detection, pathway overlap
- **SAE Feature**: DNA repair capacity, PGx flags
- **Status**: Operational

**Source 7: Off-Target Heuristics** ✅
- **Endpoint**: `/api/safety/off_target_preview`
- **Provides**: GC content, homopolymer detection
- **SAE Feature**: Seed region quality
- **Status**: Operational

**Source 8: Evidence & Literature** ✅
- **Endpoint**: `/api/evidence/deep_analysis`
- **Provides**: Evidence tier, citation count
- **SAE Feature**: Literature evidence strength
- **Status**: Operational

**Source 9: TCGA-OV Cohort Data** ✅
- **File**: `data/validation/sae_cohort/sae_features_tcga_ov_platinum.json`
- **Provides**: 66 patients, 2,897 variants with SAE features
- **Status**: Extracted, ready for biomarker analysis

### 10.2 Data Sources We're Missing

**Missing 1: Cohort Signals (Real-World Validation)** ⚠️
- **Planned Endpoint**: `/api/datasets/extract_and_benchmark`
- **Provides**: Cohort coverage, response rates, outcome data
- **SAE Feature**: Cohort overlap
- **Status**: Stub (future, Week 2)
- **Impact**: Cannot validate against real-world outcomes yet

**Missing 2: Feature→Pathway Mapping Table** ❌
- **Expected File**: `api/resources/sae_feature_mapping.json`
- **Provides**: Maps 32K SAE features → 7D pathway scores
- **Status**: Not implemented (critical blocker)
- **Impact**: Cannot use TRUE SAE features for pathway scoring

**Missing 3: Validated Biomarker Results** ⚠️
- **Expected File**: `data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers.json`
- **Provides**: Top significant features with correlations
- **Status**: Needs re-run (bug fixed, analysis invalid)
- **Impact**: Cannot identify which SAE features are predictive

**Missing 4: HRD Scores for TCGA-OV** ⚠️
- **Expected**: HRD scores (0-100) for validation cohort
- **Status**: Jr2 assigned but not executed
- **Impact**: Cannot validate DNA repair capacity formula

### 10.3 Data Sources We Need

**Need 1: Top Significant SAE Features** (After Stage 2 Re-Run)
- **Source**: Biomarker correlation analysis
- **Format**: List of feature indices with correlation scores
- **Timeline**: 2 hours (re-run analysis)
- **Priority**: P0 (critical path)

**Need 2: Feature→Pathway Mapping Table** (After Stage 3)
- **Source**: Manual pathway annotation
- **Format**: JSON mapping file
- **Timeline**: 8 hours (manual curation)
- **Priority**: P0 (critical blocker)

**Need 3: HRD Scores for Validation** (Jr2 Mission)
- **Source**: cBioPortal extraction
- **Format**: TCGA patient IDs → HRD scores
- **Timeline**: 1-2 days (Jr2 execution)
- **Priority**: P1 (validation enhancement)

**Need 4: Real-World Outcome Data** (Future)
- **Source**: Yale retrospective data, published trials
- **Format**: Patient outcomes (response, OS, PFS)
- **Timeline**: TBD (external collaboration)
- **Priority**: P2 (future validation)

---

## Phase 11: Scope Review & Completeness Check

### 11.1 Scope Coverage Verification

**SAE Theory & Implementation** ✅
- ✅ BatchTopKTiedSAE architecture understood
- ✅ 32K features, layer 26 activations documented
- ✅ Current proxy implementation vs true SAE documented
- ✅ Dimension mismatch issue documented

**WIWFM S/P/E Framework** ✅
- ✅ Sequence scoring flow traced
- ✅ Pathway aggregation flow traced
- ✅ Evidence gathering flow traced
- ✅ Drug scoring flow traced
- ✅ SAE integration gap documented

**Integration Points** ✅
- ✅ SAE → WIWFM (display-only, not integrated)
- ✅ SAE → Resistance Detection (integrated)
- ✅ SAE → Mechanism Fit Ranker (integrated, bug documented)
- ✅ SAE → Biomarker Correlation (operational)

**Pipeline Stages** ✅
- ✅ Stage 1: SAE Feature Extraction (100% complete)
- ✅ Stage 2: Biomarker Correlation (95% complete, needs re-run)
- ✅ Stage 3: Feature→Pathway Mapping (0% complete, blocker)
- ✅ Stage 4: Service Enhancement (0% complete, waiting)
- ✅ Stage 5: Resistance Prophet (100% complete)

**Advanced Care Plan Enablement** ✅
- ✅ All 8 features analyzed
- ✅ SAE role for each feature documented
- ✅ Current state vs future state documented
- ✅ Path to full vision documented

### 11.2 Potential Scope Gaps

**Gap 1: Evo2 SAE Notebook Integration** ⚠️
- **Question**: Have we reviewed the official Evo2 SAE notebook?
- **File**: `scripts/evo2/evo2/notebooks/sparse_autoencoder/sparse_autoencoder.ipynb`
- **Status**: Not reviewed in detail
- **Action**: Review notebook to understand official SAE implementation pattern

**Gap 2: Modal SAE Service Deployment** ⚠️
- **Question**: Is the Modal service actually deployed and accessible?
- **File**: `src/services/sae_service/main.py`
- **Status**: Code exists, deployment status unclear
- **Action**: Verify Modal service is live and responding

**Gap 3: Frontend Integration** ⚠️
- **Question**: How does frontend consume SAE features?
- **Files**: Frontend components (not reviewed)
- **Status**: Backend reviewed, frontend not reviewed
- **Action**: Review frontend integration if needed

**Gap 4: Manager Policy Documents** ⚠️
- **Question**: Are all Manager policy documents reviewed?
- **Files**: `MANAGER_ANSWERS_TO_ZO_SAE_QUESTIONS.md` and related
- **Status**: Referenced but not fully reviewed
- **Action**: Review all Manager policy documents

### 11.3 Scope Completeness Checklist

- [x] SAE theory and architecture
- [x] Current implementation (proxy vs true)
- [x] WIWFM S/P/E framework
- [x] Integration points mapped
- [x] Pipeline stages documented
- [x] Data sources inventoried
- [x] Testing infrastructure identified
- [x] Error handling patterns documented
- [x] Advanced Care Plan enablement
- [x] Mechanism fit validation bug
- [ ] Evo2 SAE notebook reviewed (potential gap)
- [ ] Modal service deployment verified (potential gap)
- [ ] Frontend integration reviewed (if needed)
- [ ] All Manager policy documents reviewed (potential gap)

---

## Phase 12: Fail-Safe Mechanisms to Prevent Hallucinations

### 12.1 Code-First Validation

**Principle**: Never assume - always verify with code

**Implementation**:
1. **Code References Required**: Every claim must have code file + line number
2. **Execution Path Tracing**: Follow code from entry point to output
3. **Integration Point Verification**: Verify service connections with actual code
4. **Data Source Verification**: Check data files exist and are accessible

**Examples**:
- ✅ "SAE features are proxy" → Code reference: `sae_feature_service.py:243`
- ✅ "Mechanism vector uses proxy" → Code reference: `sae_feature_service.py:206-214`
- ✅ "Biomarker analysis needs re-run" → Code reference: `biomarker_correlation_service.py:502`

### 12.2 Data Validation Checkpoints

**Checkpoint 1: File Existence**
- **Before Analysis**: Verify data files exist
- **Code**: `Path(file).exists()` check
- **Failure**: Raise `FileNotFoundError`, don't proceed

**Checkpoint 2: Data Structure Validation**
- **Before Processing**: Validate JSON structure
- **Code**: Schema validation or structure checks
- **Failure**: Raise `ValueError`, log error, don't proceed

**Checkpoint 3: Outcome Field Validation**
- **Before Correlation**: Verify outcome field populated
- **Code**: Check `outcome` field not null
- **Failure**: Block analysis, alert user

**Checkpoint 4: Feature Index Validation**
- **Before Processing**: Verify indices in valid range (0-32767)
- **Code**: `assert 0 <= idx <= 32767`
- **Failure**: Raise `ValueError`, skip invalid features

### 12.3 Formula Validation

**Manager-Approved Formulas** (Must Match Exactly):
- **DNA Repair Capacity**: `0.6×DDR + 0.2×ess + 0.2×exon`
- **Mechanism Fit**: `0.7×eligibility + 0.3×mechanism_fit`
- **Resistance Detection**: 2-of-3 trigger rule

**Validation**:
- **Unit Tests**: Assert exact formula values
- **Code**: `test_dna_repair_capacity_formula()` in `test_sae_phase2_services.py`
- **Failure**: Block commit, alert manager

### 12.4 Error Handling Patterns

**Pattern 1: Graceful Degradation**
- **Code**: `try/except` blocks with fallback values
- **Example**: `sae_feature_service.py:265-266` - Falls back to proxy if true SAE fails
- **Principle**: Never crash, always provide partial results

**Pattern 2: Provenance Tracking**
- **Code**: Log all fallbacks and errors in provenance
- **Example**: `response.provenance["fallback"] = "service_timeout"`
- **Principle**: Transparency - user knows what failed

**Pattern 3: Validation Before Use**
- **Code**: Check data exists and is valid before processing
- **Example**: `biomarker_correlation_service.py:502` - Validates outcome field
- **Principle**: Fail fast, don't process invalid data

**Pattern 4: Circuit Breakers**
- **Code**: Disable service after N failures
- **Example**: SAE service circuit breaker (not yet implemented)
- **Principle**: Prevent cascading failures

### 12.5 Anti-Hallucination Rules

**Rule 1: Code References Required**
- **Enforcement**: Every claim must cite code file + line
- **Example**: "SAE uses proxy" → `sae_feature_service.py:243`
- **Violation**: Mark as unverified, require code reference

**Rule 2: No Assumptions**
- **Enforcement**: Verify with code, don't assume
- **Example**: Don't assume mapping table exists - check `Path(file).exists()`
- **Violation**: Mark as assumption, require verification

**Rule 3: Data Validation First**
- **Enforcement**: Validate data before processing
- **Example**: Check outcome field populated before correlation
- **Violation**: Block processing, require data fix

**Rule 4: Test-Driven Validation**
- **Enforcement**: Write tests for all formulas and logic
- **Example**: `test_dna_repair_capacity_formula()` validates Manager's formula
- **Violation**: Block commit until tests pass

**Rule 5: Fail Gracefully**
- **Enforcement**: Never crash, always provide partial results
- **Example**: Return error dict instead of raising exception
- **Violation**: Refactor to use graceful degradation

---

## Integration into Main Plan

These sections should be added as:
- **Phase 9**: Testing Infrastructure & Validation Mechanisms
- **Phase 10**: Data Sources Inventory
- **Phase 11**: Scope Review & Completeness Check
- **Phase 12**: Fail-Safe Mechanisms to Prevent Hallucinations

Each phase includes:
- Code references for all claims
- Validation checkpoints
- Failure actions
- Anti-hallucination rules

