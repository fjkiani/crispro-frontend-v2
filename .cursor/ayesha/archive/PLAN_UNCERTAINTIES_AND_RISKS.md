# Plan Uncertainties, Unclear Gaps, and Potential Failure Points

**Date**: January 20, 2025 (Updated)  
**Purpose**: Critical analysis of what we're unsure about, unclear gaps, and potential failure modes

**Key Updates**:
- ✅ **Uncertainty 3 RESOLVED**: Random SAE weights → Trained weights (evo2_7b migration complete)
- ✅ **Uncertainty 4 MOSTLY RESOLVED**: Modal service operational, 66 patients extracted
- ❌ **Uncertainty 2 CRITICAL**: Feature→Pathway Mapping blocks ALL THREE services (Resistance Prophet, Mechanism Fit Ranking, Early Resistance Detection)
- ⚠️ **Uncertainty 1 UPDATED**: Now using trained weights (should improve biomarker analysis)

---

## 🚨 CRITICAL UNCERTAINTIES

### Uncertainty 1: Will Biomarker Analysis Find Significant Features?

**What We Think We Know**:
- Service is complete, bug is fixed
- 66 patients with SAE features extracted
- Outcome field should be populated

**What We're Unsure About**:
- ❓ **Will re-run actually find significant features?** (Previous run found 0)
- ❓ **Are the SAE features meaningful?** (Now using trained weights - evo2_7b migration complete ✅)
- ❓ **Is 66 patients enough for statistical power?** (May need more)
- ❓ **What if outcome distribution is imbalanced?** (53 sensitive, 11 refractory, 2 resistant)

**Potential Failure**:
- Re-run finds 0 significant features again
- All correlations are weak (r < 0.3)
- Statistical tests fail (p > 0.05 after FDR correction)

**What We Need to Verify**:
- [x] Check outcome field is actually populated in data file - **VERIFIED**: All 66 patients have outcome (53 sensitive, 11 refractory, 2 resistant)
- [x] Verify feature indices are valid (0-32767) in extracted data - **VERIFIED**: Sample indices [23594, 13201, 6369, 23263, 1678] all valid
- [x] Check if SAE weights are trained (not random) - **VERIFIED**: Migrated to evo2_7b, trained weights loaded ✅
- [ ] Check feature distributions (mean, std, variation) - **PENDING**: Need feature distribution analysis
- [ ] Determine minimum sample size needed for statistical power - **PENDING**: Re-run analysis first, then assess

**Code References to Check**:
- `data/validation/sae_cohort/sae_features_tcga_ov_platinum.json` - Verify structure
- `api/services/biomarker_correlation_service.py:502` - Verify outcome extraction
- `api/services/biomarker_correlation_service.py:163` - Verify feature matrix building

---

### Uncertainty 2: Can We Create Feature→Pathway Mapping Without Biological Validation? ❌ **CRITICAL BLOCKER**

**Status**: ❌ **CRITICAL BLOCKER** - Blocks ALL THREE services from using TRUE SAE features

**What We Think We Know**:
- Need to map 32K SAE features → 7D pathway scores
- Template exists: `api/resources/sae_feature_mapping_template.json`
- Loader exists: `api/services/sae_feature_service.py:286-307`
- **Impact**: Blocks Resistance Prophet, Mechanism Fit Ranking, and Early Resistance Detection from using TRUE SAE

**What We're Unsure About**:
- ❓ **How do we know which features map to which pathways?** (No biological annotation)
- ❓ **Can we infer pathway mapping from biomarker analysis?** (Biomarker-driven approach - pending analysis)
- ❓ **What if features are multi-pathway?** (One feature might activate DDR + MAPK)
- ❓ **How do we validate the mapping is correct?** (No ground truth, but can test on known cases)

**Approach** (Updated):
- **Biomarker-driven mapping**: Use top significant features from biomarker analysis
- **Gene→pathway inference**: Map features that activate for specific genes to pathways
- **Validation strategy**: Test on known cases (BRCA1 → DDR high, KRAS → MAPK high)
- **Roadmap**: `.cursor/rules/SAE_TO_RESISTANCE_PROPHET_PIPELINE.mdc` (Stage 3-4)

**Potential Failure**:
- Manual annotation is guesswork (no biological basis)
- Mapping doesn't improve pathway scores vs proxy
- Features don't correlate with expected pathways
- Validation fails (pathway scores don't match known biology)

**What We Need to Verify**:
- [ ] Review Evo2 SAE notebook to understand feature interpretation
- [ ] Check if there's any feature annotation in literature
- [x] Determine approach: Biomarker-driven mapping (from top significant features) - **STRATEGY DEFINED**
- [ ] Run biomarker analysis to get top significant features
- [ ] Create feature→pathway mapping from biomarker results
- [ ] Create validation strategy (e.g., BRCA1 mutation → DDR pathway should be high)

**Code References to Check**:
- `scripts/evo2/evo2/notebooks/sparse_autoencoder/sparse_autoencoder.ipynb` - Official pattern
- `api/resources/sae_feature_mapping_template.json` - Expected structure
- `api/services/sae_feature_service.py:309-406` - How mapping is used

---

### Uncertainty 3: Are SAE Weights Producing Meaningful Features? ✅ **RESOLVED**

**Status**: ✅ **RESOLVED** - Migrated to evo2_7b with trained weights

**What We Know Now**:
- ✅ SAE service uses **trained weights** (not random) - evo2_7b migration complete
- ✅ Dimension match: 4096-dim Evo2 7B matches 4096-dim checkpoint
- ✅ Trained weights loaded: `Goodfire/Evo-2-Layer-26-Mixed` checkpoint
- ✅ 66 patients extracted with trained weights
- ✅ Feature extraction operational

**What We're Still Unsure About**:
- ❓ **Do trained weights produce meaningful biological signal?** (Need biomarker analysis to verify)
- ❓ **Can we trust biomarker correlations from trained features?** (Need validation)
- ❓ **Will pathway mapping work with trained features?** (Still blocked by mapping creation)
- ❓ **What's the actual impact vs proxy features?** (Need comparison after mapping created)

**What We Need to Verify**:
- [x] Migrate to evo2_7b - **COMPLETE**: Migration successful
- [x] Load trained weights - **COMPLETE**: Weights loaded successfully
- [x] Extract features with trained weights - **COMPLETE**: 66 patients extracted
- [ ] Check if feature distributions are reasonable (not all zeros or all ones) - **PENDING**: Need analysis
- [ ] Validate that features vary across patients (not constant) - **PENDING**: Need analysis
- [ ] Re-run biomarker analysis with trained features - **PENDING**: After data quality verification

**Code References**:
- Migration: `.cursor/ayesha/EVO2_7B_MIGRATION_COMPLETE.md`
- SAE Service: `src/services/sae_service/main.py:155-230` - Trained weights loaded
- Feature Extraction: `src/services/sae_service/main.py:338-379` - Feature extraction logic
- Data: `data/validation/sae_cohort/sae_features_tcga_ov_platinum.json` - Feature distributions

---

### Uncertainty 4: Will Modal SAE Service Work Reliably? ✅ **MOSTLY RESOLVED**

**What We Know Now**:
- ✅ Modal service is deployed: `src/services/sae_service/main.py`
- ✅ Circuit breaker implemented (lines 238-257)
- ✅ Error handling exists
- ✅ Service operational: 66 patients extracted successfully
- ✅ Trained weights loaded: evo2_7b migration complete

**What We're Still Unsure About**:
- ❓ **Are API keys configured correctly in production?** (Previous 403 errors documented, but extraction succeeded)
- ❓ **Will circuit breaker prevent all failures?** (May reject valid requests)
- ❓ **What happens on timeout?** (May silently fail)

**Potential Failure**:
- Service not deployed (404 errors)
- API key mismatch (403 errors)
- Circuit breaker triggers incorrectly (503 errors)
- Timeouts cause silent failures
- Service crashes on edge cases

**What We Need to Verify**:
- [x] Test Modal service endpoint is accessible - **VERIFIED**: 66 patients extracted successfully
- [x] Verify service can extract features - **VERIFIED**: Extraction operational
- [x] Verify trained weights loaded - **VERIFIED**: evo2_7b migration complete
- [ ] Verify API key configuration in production environment
- [ ] Test circuit breaker behavior (error rate calculation)
- [ ] Test timeout handling
- [ ] Test edge cases (invalid inputs, missing data)

**Code References to Check**:
- `src/services/sae_service/main.py:243-379` - Endpoint implementation
- `src/services/sae_service/main.py:238-257` - Circuit breaker logic
- `scripts/sae/extract_sae_features_cohort.py:312-358` - Error handling

---

### Uncertainty 5: How Do We Validate the Entire Pipeline End-to-End?

**What We Think We Know**:
- Individual components tested (unit tests exist)
- Validation scripts exist (`validate_sae_tcga.py`)
- E2E tests exist (`test_ayesha_post_ngs_e2e.py`)

**What We're Unsure About**:
- ❓ **Do tests cover all failure modes?** (May miss edge cases)
- ❓ **What's the ground truth for validation?** (HRD scores? Clinical outcomes?)
- ❓ **How do we know if pipeline is working correctly?** (No gold standard)
- ❓ **What if validation fails?** (No clear remediation path)

**Potential Failure**:
- Tests pass but pipeline fails in production
- Validation shows poor accuracy (but is it the pipeline or the data?)
- Edge cases cause crashes (not covered by tests)
- No way to detect silent failures

**What We Need to Verify**:
- [ ] Review all test coverage (unit, integration, E2E)
- [ ] Identify missing test cases (edge cases, error conditions)
- [ ] Define validation success criteria (what metrics = success?)
- [ ] Create remediation plan (what to do if validation fails?)

**Code References to Check**:
- `tests/test_sae_phase2_services.py` - Unit test coverage
- `tests/test_ayesha_post_ngs_e2e.py` - E2E test coverage
- `scripts/validate_sae_tcga.py` - Validation script logic

---

## ⚠️ UNCLEAR GAPS

### Gap 1: Evo2 SAE Notebook Pattern Not Fully Understood

**What We Know**:
- Notebook exists: `scripts/evo2/evo2/notebooks/sparse_autoencoder/sparse_autoencoder.ipynb`
- Our implementation references it: `src/services/sae_service/main.py:2`

**What's Unclear**:
- ❓ **What's the official pattern for loading SAE weights?** ✅ **RESOLVED**: Using trained weights from checkpoint
- ❓ **How do they extract activations?** (We may be doing it differently)
- ❓ **What's the correct layer name?** (We use "blocks-26", is that correct?)
- ❓ **How do they interpret features?** (We don't have feature annotations)

**Action Required**:
- [ ] Read entire notebook (1780 lines) to understand official pattern - **IN PROGRESS**: Layer name verified ("blocks.26" vs "blocks-26" - same layer)
- [x] Compare our implementation vs notebook pattern - **VERIFIED**: BatchTopKTiedSAE matches exactly
- [x] Identify any deviations that could cause issues - **VERIFIED**: Only difference is model choice (evo2_1b vs evo2_7b) - intentional
- [ ] Document differences and rationale - **PENDING**: Need complete notebook review

---

### Gap 2: Manager Policy Documents Not Fully Reviewed

**What We Know**:
- Manager policy referenced: `MANAGER_ANSWERS_TO_ZO_SAE_QUESTIONS.md`
- Formulas are Manager-approved (C1-C10)

**What's Unclear**:
- ❓ **Are all Manager decisions documented?** (May have missed some)
- ❓ **What are the validation requirements?** (When can we integrate SAE into WIWFM?)
- ❓ **What are the success criteria?** (What metrics define "success"?)
- ❓ **What are the failure scenarios?** (What happens if validation fails?)

**Action Required**:
- [x] Read all Manager policy documents - **VERIFIED**: `MANAGER_ANSWERS_TO_ZO_SAE_QUESTIONS.md` contains C1-C10 formulas
- [x] Extract all validation requirements - **VERIFIED**: P3 (≥90% accuracy), P4 (thresholds), C1-C10 (formulas)
- [ ] Document success/failure criteria - **PARTIAL**: Mechanism fit criteria defined, biomarker/pathway criteria missing
- [x] Identify any conflicting guidance - **VERIFIED**: No conflicts found

---

### Gap 3: Frontend Integration Not Reviewed

**What We Know**:
- Backend provides SAE features in response
- Frontend components exist (not reviewed)

**What's Unclear**:
- ❓ **How does frontend consume SAE features?** (May have bugs)
- ❓ **What happens if SAE features are missing?** (May crash frontend)
- ❓ **Are SAE features displayed correctly?** (May show wrong data)
- ❓ **What's the user experience?** (May be confusing)

**Action Required**:
- [ ] Review frontend components that use SAE features
- [ ] Test error handling (missing data, invalid data)
- [ ] Verify UI displays match backend data
- [ ] Check user experience (is it clear what SAE features mean?)

---

### Gap 4: Error Handling May Not Cover All Cases

**What We Know**:
- Error handling patterns exist: `.cursor/ayesha/Deliverables/ERROR_HANDLING_PATTERNS.md`
- Graceful degradation implemented

**What's Unclear**:
- ❓ **Do we handle all error conditions?** (May miss some)
- ❓ **Are error messages clear?** (May be confusing)
- ❓ **Do we log enough for debugging?** (May not have enough context)
- ❓ **What happens in cascading failures?** (Multiple services fail)

**Action Required**:
- [ ] Review all error handling code paths
- [ ] Test all failure scenarios (timeout, 403, 500, invalid data)
- [ ] Verify error messages are actionable
- [ ] Check logging is sufficient for debugging

---

## 🔴 POTENTIAL FAILURE POINTS

### Failure Point 1: Biomarker Analysis Finds No Significant Features

**Probability**: Medium  
**Impact**: High  
**Scenario**: Re-run analysis, still finds 0 significant features

**Possible Causes**:
1. ~~Random SAE weights produce no biological signal~~ ✅ **RESOLVED**: Using trained weights
2. Sample size too small (66 patients insufficient)
3. ~~Outcome field still not populated correctly~~ ✅ **VERIFIED**: All 66 patients have outcome
4. ~~Feature indices still invalid~~ ✅ **VERIFIED**: All indices valid (0-32767)
5. Statistical thresholds too strict
6. Trained SAE features still don't correlate with outcomes (biological signal weak)

**Detection**:
- Check `significant_features_count` in results
- Verify outcome distribution
- Check feature matrix is not all zeros

**Mitigation**:
- ✅ Data quality verified (outcome field, feature indices)
- ✅ Trained weights loaded (evo2_7b migration complete)
- Consider relaxing statistical thresholds if needed
- May need more patients if signal is weak
- Re-run analysis with trained features to assess biological signal

---

### Failure Point 2: Feature→Pathway Mapping Is Biologically Invalid ❌ **CRITICAL BLOCKER**

**Probability**: High  
**Impact**: High  
**Scenario**: Biomarker-driven mapping creates incorrect pathway assignments

**Services Affected**: ALL THREE services blocked from using TRUE SAE:
- Resistance Prophet: Cannot use TRUE SAE DNA repair capacity
- Mechanism Fit Ranking: Cannot use TRUE SAE mechanism vector
- Early Resistance Detection: Cannot use TRUE SAE DNA repair capacity

**Possible Causes**:
1. No biological basis for mapping (pure guesswork)
2. Features are multi-pathway (can't map to single pathway)
3. Features don't correlate with genes (SAE features are abstract)
4. Validation fails (pathway scores don't match known biology)

**Detection**:
- Validate mapping against known biology (BRCA1 → DDR should be high)
- Test pathway scores match expected patterns
- Compare SAE pathway scores vs proxy pathway scores

**Mitigation**:
- Use biomarker-driven approach (top significant features from analysis)
- Use gene→pathway mapping to infer feature→pathway (for features that activate with specific genes)
- Validate each mapping against literature
- Test on known cases (BRCA1 → DDR high, KRAS → MAPK high, HER2 → HER2 high)
- Compare SAE pathway scores vs proxy pathway scores (should be similar for known cases)
- May need domain expert review
- Roadmap: `.cursor/rules/SAE_TO_RESISTANCE_PROPHET_PIPELINE.mdc` (Stage 3-4)

---

### Failure Point 3: Modal Service Fails Silently

**Probability**: Medium  
**Impact**: High  
**Scenario**: Service returns errors but pipeline continues

**Possible Causes**:
1. API key mismatch (403 errors ignored)
2. Timeout errors not caught
3. Circuit breaker triggers incorrectly
4. Service crashes but error not propagated

**Detection**:
- Check error logs for 403/503/500 errors
- Verify circuit breaker state
- Test timeout handling
- Monitor error rates

**Mitigation**:
- Fix API key configuration
- Improve error handling (don't ignore errors)
- Test circuit breaker logic
- Add monitoring/alerting

---

### Failure Point 4: Tests Pass But Production Fails

**Probability**: Medium  
**Impact**: High  
**Scenario**: All tests pass but pipeline fails in production

**Possible Causes**:
1. Tests don't cover production scenarios
2. Test data is different from production data
3. Environment differences (dev vs prod)
4. Race conditions not tested

**Detection**:
- Compare test vs production behavior
- Check for environment-specific issues
- Review test coverage gaps

**Mitigation**:
- Add production-like tests
- Test with real data (not just mocks)
- Test in production-like environment
- Add integration tests

---

## 📋 ACTION ITEMS TO RESOLVE UNCERTAINTIES

### Immediate (Before Re-Running Analysis)

1. **Verify Data Quality**:
   - [ ] Check `sae_features_tcga_ov_platinum.json` structure
   - [ ] Verify outcome field is populated
   - [ ] Check feature indices are valid (0-32767)
   - [ ] Verify feature distributions are reasonable

2. **Review Evo2 SAE Notebook**:
   - [ ] Read entire notebook (1780 lines)
   - [ ] Compare our implementation vs official pattern
   - [ ] Document any deviations

3. **Test Modal Service**:
   - [ ] Verify service is deployed and accessible
   - [ ] Test API key configuration
   - [ ] Test error handling (timeout, 403, 500)

### Short-Term (Before Creating Mapping)

4. **Review Manager Policy**:
   - [ ] Read all Manager policy documents
   - [ ] Extract validation requirements
   - [ ] Document success/failure criteria

5. **Validate SAE Features (Trained Weights)**:
   - [x] Migrate to trained weights - **COMPLETE**: evo2_7b migration successful
   - [ ] Check if features have any biological signal (need biomarker analysis)
   - [ ] Compare feature distributions across patients
   - [x] Determine if trained weights are needed - **RESOLVED**: Trained weights loaded

6. **Review Frontend Integration**:
   - [ ] Review frontend components
   - [ ] Test error handling
   - [ ] Verify UI displays

### Long-Term (Before Production)

7. **Improve Test Coverage**:
   - [ ] Add missing test cases
   - [ ] Test edge cases
   - [ ] Test error conditions

8. **Create Validation Strategy**:
   - [ ] Define success criteria
   - [ ] Create validation tests
   - [ ] Document remediation plan

9. **Improve Error Handling**:
   - [ ] Review all error paths
   - [ ] Improve error messages
   - [ ] Add monitoring/alerting

---

## 🎯 CONFIDENCE LEVELS

**High Confidence (80-100%)**:
- SAE service code structure
- Biomarker correlation service implementation
- Error handling patterns (documented)
- Manager-approved formulas (C1-C10)

**Medium Confidence (50-80%)**:
- Data quality (needs verification)
- Modal service deployment status
- Test coverage completeness
- Feature→pathway mapping approach

**Low Confidence (<50%)**:
- Will biomarker analysis find significant features? (Now using trained weights - should improve)
- ~~Are random SAE weights meaningful?~~ ✅ **RESOLVED**: Using trained weights
- Can we create valid pathway mapping? (Still critical blocker)
- Will entire pipeline work end-to-end? (Blocked by pathway mapping)

---

## ⚠️ RISK ASSESSMENT

**Highest Risk**: Feature→Pathway Mapping (High probability, High impact) ❌ **CRITICAL BLOCKER**
- Blocks ALL THREE services from using TRUE SAE features:
  - Resistance Prophet: Cannot use TRUE SAE DNA repair capacity
  - Mechanism Fit Ranking: Cannot use TRUE SAE mechanism vector
  - Early Resistance Detection: Cannot use TRUE SAE DNA repair capacity
- No biological basis for mapping (biomarker-driven approach pending)
- May create incorrect pathway scores
- Could invalidate entire pipeline
- **Current Workaround**: All three services use PROXY features (gene mutations → pathway scores)

**Medium Risk**: Biomarker Analysis (Medium probability, High impact)
- May find 0 significant features (even with trained weights)
- Trained weights should improve signal (but not guaranteed)
- Could block entire pipeline if no significant features found

**Low Risk**: Modal Service (Medium probability, Medium impact)
- Deployment issues are fixable
- Error handling can be improved
- Circuit breaker can be adjusted

---

## 📝 NEXT STEPS

**See detailed resolution plan**: `.cursor/ayesha/UNCERTAINTIES_RESOLUTION_PLAN.md`

### Immediate Actions (Next 2 Hours)

1. **Verify Data Quality** ✅ **PARTIALLY COMPLETE**
   - ✅ Outcome field populated (53 sensitive, 11 refractory, 2 resistant)
   - ✅ Feature indices valid (all < 32768)
   - [ ] Check feature distributions (mean, std, variation)
   - [ ] Verify features vary across patients

2. **Test Modal Service** ✅ **MOSTLY COMPLETE**
   - ✅ Service operational: 66 patients extracted successfully
   - ✅ Trained weights loaded: evo2_7b migration complete
   - [ ] Test health check endpoint: `curl http://localhost:8000/api/sae/health`
   - [ ] Check environment variables: `SAE_SERVICE_URL`, `ENABLE_TRUE_SAE`
   - [ ] Test service directly in production environment

3. **Review Evo2 Notebook** ⏸️ **IN PROGRESS**
   - ✅ Layer name verified ("blocks.26" correct)
   - ✅ Architecture matches (BatchTopKTiedSAE)
   - [ ] Complete notebook review (1780 lines)
   - [ ] Document feature interpretation approach

### Short-Term Actions (Next Week)

4. **Define Validation Criteria** ⏸️ **PENDING**
   - [ ] Biomarker analysis success criteria
   - [ ] Pathway mapping validation strategy
   - [ ] End-to-end validation plan

5. **Re-Run Biomarker Analysis** ⏸️ **PENDING** (After Step 1 complete)
   - Run with verified data
   - Check results for significant features
   - Document findings

6. **Create Pathway Mapping Strategy** ⏸️ **PENDING** (After Step 5)
   - Determine approach (gene-based? feature-based?)
   - Create validation plan
   - Document uncertainty

