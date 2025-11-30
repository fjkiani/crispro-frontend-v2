# Uncertainties Resolution Plan - Findings & Action Items

**Date**: January 20, 2025 (Updated)  
**Purpose**: Address each uncertainty from PLAN_UNCERTAINTIES_AND_RISKS.md with concrete findings and action items

**Note**: This document provides action items for each uncertainty. For master risk analysis, see:
- `.cursor/ayesha/PLAN_UNCERTAINTIES_AND_RISKS.md` (master uncertainties document)

**Status**: Supporting document - may be redundant with master uncertainties document

---

## ✅ UNCERTAINTY 1: Will Biomarker Analysis Find Significant Features?

### Findings (Verified)

**Data Quality - VERIFIED**:
- ✅ **Outcome field populated**: All 66 patients have outcome field
  - Distribution: 53 sensitive, 11 refractory, 2 resistant
  - **Code verification**: `data/validation/sae_cohort/sae_features_tcga_ov_platinum.json`
  - **Command result**: `Counter({'sensitive': 53, 'refractory': 11, 'resistant': 2})`

- ✅ **Feature indices valid**: Sample indices [23594, 13201, 6369, 23263, 1678] all < 32768
  - **Range**: All indices in valid range (0-32767)
  - **Structure**: Each variant has `top_features` array with `index` and `value`

- ✅ **Bug fix verified**: Line 502 uses `outcome` field correctly
  - **Code**: `self.outcome_labels = [p.get("platinum_response") or p.get("outcome") for p in patients]`
  - **Status**: Fix is in place, ready for re-run

**Remaining Concerns**:
- ⚠️ **Sample size**: 66 patients (53 sensitive, 11 refractory, 2 resistant)
  - **Imbalanced**: 80% sensitive vs 20% resistant/refractory
  - **Statistical power**: May need more resistant patients for robust correlations
  - **Action**: Re-run analysis first, then assess if more data needed

- ✅ **SAE weights**: Now using trained weights (evo2_7b migration complete) - see Uncertainty 3

### Action Items

**Immediate (Before Re-Run)**:
1. ✅ **Data quality verified** - Outcome field populated, indices valid
2. [ ] **Check feature distributions**: Verify features vary across patients (not constant)
3. [ ] **Check feature values**: Verify values are reasonable (not all zeros/ones)

**After Re-Run**:
4. [ ] **If 0 significant features found**: Investigate other causes (sample size, outcome imbalance, statistical thresholds)
5. [ ] **If weak correlations**: Consider sample size increase or feature selection
6. [ ] **If imbalanced outcomes**: Consider stratified analysis or resampling

### Code Verification Commands

```bash
# Verify outcome distribution
python3 -c "import json; data = json.load(open('data/validation/sae_cohort/sae_features_tcga_ov_platinum.json')); from collections import Counter; print(Counter([p.get('outcome') for p in data['patients']]))"

# Verify feature indices range
python3 -c "import json; data = json.load(open('data/validation/sae_cohort/sae_features_tcga_ov_platinum.json')); indices = [f['index'] for p in data['patients'] for v in p['variants'] for f in v.get('top_features', [])]; print(f'Min: {min(indices)}, Max: {max(indices)}, All valid: {all(0 <= i < 32768 for i in indices)}')"

# Check feature value distributions
python3 -c "import json, numpy as np; data = json.load(open('data/validation/sae_cohort/sae_features_tcga_ov_platinum.json')); values = [f['value'] for p in data['patients'] for v in p['variants'] for f in v.get('top_features', [])]; print(f'Mean: {np.mean(values):.2e}, Std: {np.std(values):.2e}, Min: {np.min(values):.2e}, Max: {np.max(values):.2e}')"
```

---

## ⚠️ UNCERTAINTY 2: Can We Create Feature→Pathway Mapping Without Biological Validation?

### Findings (Partially Resolved)

**Evo2 SAE Notebook Pattern - VERIFIED**:
- ✅ **Layer name**: We use `"blocks.26"`, notebook uses `"blocks-26"` (different separator, same layer)
  - **Code**: `src/services/sae_service/main.py:164` uses `"blocks.26"`
  - **Notebook**: Uses `"blocks-26"` (line 1516)
  - **Status**: Both refer to layer 26, implementation is correct

- ✅ **SAE Architecture**: Our `BatchTopKTiedSAE` matches notebook exactly
  - **Code**: `src/services/sae_service/main.py:79-136` matches notebook pattern
  - **Key methods**: `encoder_pre`, `encode`, `_batch_topk`, `decode`, `forward` all match

**Remaining Concerns**:
- ❓ **Feature interpretation**: Notebook doesn't provide feature→pathway mapping
  - **Finding**: Notebook shows how to extract features, not how to interpret them
  - **Gap**: No biological annotation for 32K features
  - **Action**: Need alternative approach (gene→pathway inference)

- ❓ **Mapping strategy**: No clear biological basis
  - **Option 1**: Gene→pathway mapping → infer feature→pathway (assumes features correlate with genes)
  - **Option 2**: Literature mining (search PubMed for feature associations)
  - **Option 3**: Evo2 activation analysis (which pathways activate when feature is active)
  - **Risk**: All options are indirect and may be inaccurate

### Action Items

**Immediate**:
1. [ ] **Review notebook completely**: Read all 1780 lines to understand feature interpretation
2. [ ] **Check literature**: Search for SAE feature annotation papers
3. [ ] **Design mapping strategy**: Choose approach (gene-based vs feature-based)

**After Stage 2 (Biomarker Results)**:
4. [ ] **Analyze significant features**: Which genes/mutations activate them?
5. [ ] **Create gene→pathway mapping**: Use known gene→pathway associations
6. [ ] **Infer feature→pathway**: Map features based on activating genes
7. [ ] **Validate mapping**: Test on known cases (BRCA1 → DDR should be high)

### Validation Strategy

**Test Cases for Mapping Validation**:
1. **BRCA1 mutation** → Should map to DDR pathway (high weight)
2. **KRAS mutation** → Should map to MAPK pathway (high weight)
3. **HER2 amplification** → Should map to HER2 pathway (high weight)
4. **TP53 mutation** → May map to multiple pathways (DDR, others)

**Success Criteria**:
- Pathway scores match expected biology (BRCA1 → high DDR)
- Pathway scores improve vs proxy (or at least match)
- Validation on known cases passes

---

## ✅ UNCERTAINTY 3: Are SAE Weights Producing Meaningful Features? **RESOLVED**

### Findings (RESOLVED ✅)

**Trained Weights Now Used**:
- ✅ **Migration complete**: Migrated to evo2_7b model
- ✅ **Trained weights loaded**: `Goodfire/Evo-2-Layer-26-Mixed` checkpoint
- ✅ **Dimension match**: 4096-dim Evo2 7B matches 4096-dim checkpoint
- ✅ **66 patients extracted**: All with trained weights
- ✅ **Code verified**: `src/services/sae_service/main.py:155-230` - Trained weights loaded successfully

**Status**: ✅ **RESOLVED** - No longer using random weights

**Remaining Questions**:
- ❓ **Do trained weights produce biological signal?** (Need biomarker analysis to verify)
  - **Status**: Trained weights should produce meaningful signal (vs random)
  - **Action**: Re-run biomarker analysis with trained features
  - **Impact**: Should improve biomarker discovery vs random weights

**What We Can Verify**:
- [ ] **Feature distributions**: Are they reasonable? (not all zeros/ones)
- [ ] **Feature variation**: Do features vary across patients? (not constant)
- [ ] **Feature activation patterns**: Do certain features activate for certain mutations?

### Action Items

**Immediate (Before Re-Run)**:
1. [ ] **Analyze feature distributions**: Check if features are reasonable
   ```bash
   # Check feature value distributions across all patients
   python3 -c "import json, numpy as np; data = json.load(open('data/validation/sae_cohort/sae_features_tcga_ov_platinum.json')); values = [f['value'] for p in data['patients'] for v in p['variants'] for f in v.get('top_features', [])]; print(f'Stats: mean={np.mean(values):.2e}, std={np.std(values):.2e}, min={np.min(values):.2e}, max={np.max(values):.2e}')"
   ```

2. [ ] **Check feature variation**: Verify features differ across patients
   ```bash
   # Check if same features activate for same mutations
   python3 -c "import json; data = json.load(open('data/validation/sae_cohort/sae_features_tcga_ov_platinum.json')); # Analyze feature patterns"
   ```

3. [ ] **Compare with mock data**: If we have mock data with known signals, compare patterns

**After Re-Run**:
4. [ ] **If 0 significant features**: Investigate other causes (sample size, outcome imbalance, statistical thresholds)
5. [ ] **If weak correlations**: Consider sample size increase or feature selection
6. [ ] **If strong correlations**: Validate against known biology, proceed with pathway mapping

### Risk Assessment (Updated)

**Status**: ✅ **RISK MITIGATED** - Trained weights now in use

**Remaining Scenarios**:
- **High Risk Scenario**: Trained weights still produce no biological signal
  - **Mitigation**: Investigate Evo2 activations, SAE model, or data quality
- **Medium Risk Scenario**: Trained weights produce weak signal
  - **Mitigation**: May need more patients, different SAE model, or feature selection
- **Low Risk Scenario**: Trained weights produce strong signal
  - **Mitigation**: Validate against known biology, proceed with pathway mapping

---

## ⚠️ UNCERTAINTY 4: Will Modal SAE Service Work Reliably?

### Findings (Partially Resolved)

**Service Code - VERIFIED**:
- ✅ **Code exists**: `src/services/sae_service/main.py` (393 lines)
- ✅ **Circuit breaker**: Implemented (lines 238-257)
- ✅ **Error handling**: Exists (try/except blocks)
- ✅ **Health check endpoint**: Exists in router (`api/routers/sae.py:188-219`)

**Deployment Status - UNCLEAR**:
- ❓ **Is service deployed?**: No clear evidence of deployment
  - **Code exists**: ✅
  - **Deployment status**: ❓ Unknown
  - **Previous issues**: 403 errors documented (API key mismatch)

**Health Check Endpoint**:
- ✅ **Code exists**: `api/routers/sae.py:188-219`
- ✅ **Checks**: Service URL configured, reachable, health endpoint
- **Action**: Test this endpoint to verify deployment

### Action Items

**Immediate**:
1. [ ] **Test health check endpoint**: 
   ```bash
   curl http://localhost:8000/api/sae/health
   # Or if backend running:
   # Check response for service status
   ```

2. [ ] **Check environment variables**:
   ```bash
   # Check if SAE_SERVICE_URL is set
   echo $SAE_SERVICE_URL
   # Check if ENABLE_TRUE_SAE is set
   echo $ENABLE_TRUE_SAE
   ```

3. [ ] **Test Modal service directly** (if URL known):
   ```bash
   curl -X POST $SAE_SERVICE_URL/extract_features \
     -H "Content-Type: application/json" \
     -d '{"chrom": "17", "pos": 43044295, "ref": "T", "alt": "G", "assembly": "GRCh38", "window": 8192}'
   ```

**If Service Not Deployed**:
4. [ ] **Deploy Modal service**:
   ```bash
   cd src/services/sae_service
   modal deploy main.py
   # Copy URL and set SAE_SERVICE_URL
   ```

**If Service Deployed But Failing**:
5. [ ] **Check API keys**: Verify `SAE_API_KEY` matches between backend and Modal
6. [ ] **Check circuit breaker**: Verify error rate calculation
7. [ ] **Test timeout handling**: Verify timeouts are caught and handled

### Code References

- **Health check**: `api/routers/sae.py:188-219`
- **Circuit breaker**: `src/services/sae_service/main.py:238-257`
- **Error handling**: `scripts/sae/extract_sae_features_cohort.py:312-358`

---

## ⚠️ UNCERTAINTY 5: How Do We Validate the Entire Pipeline End-to-End?

### Findings (Partially Resolved)

**Manager Policy - VERIFIED**:
- ✅ **Validation requirements found**: `MANAGER_ANSWERS_TO_ZO_SAE_QUESTIONS.md`
  - **P3**: Gemini trial tagging validation (≥90% accuracy)
  - **P4**: Mechanism fit thresholds (eligibility ≥0.60, mechanism_fit ≥0.50)
  - **C1-C10**: All formulas and thresholds documented

**Success Criteria - PARTIALLY DEFINED**:
- ✅ **Mechanism fit**: Top-3 accuracy ≥80%, MRR ≥0.75 (from MECHANISM_FIT_VALIDATION_GAPS.md)
- ❓ **Biomarker analysis**: Not explicitly defined in Manager policy
- ❓ **Pathway mapping**: No validation criteria defined
- ❓ **End-to-end**: No overall pipeline validation criteria

**Test Coverage - VERIFIED**:
- ✅ **Unit tests**: `tests/test_sae_phase2_services.py` (530 lines)
- ✅ **E2E tests**: `tests/test_ayesha_post_ngs_e2e.py`
- ✅ **Validation scripts**: `scripts/validate_sae_tcga.py` (532 lines)

**Remaining Concerns**:
- ❓ **Ground truth**: What's the gold standard?
  - **HRD scores**: Available for some patients (Jr2 mission)
  - **Clinical outcomes**: Available (platinum response)
  - **Pathway scores**: No ground truth (proxy is best we have)

- ❓ **Validation success criteria**: Not fully defined
  - **Biomarker analysis**: What metrics = success?
  - **Pathway mapping**: How do we know it's correct?
  - **End-to-end**: What accuracy is acceptable?

### Action Items

**Immediate**:
1. [ ] **Define biomarker analysis success criteria**:
   - Minimum significant features: ≥10?
   - Minimum correlation strength: r ≥ 0.3?
   - Minimum effect size: Cohen's d ≥ 0.5?
   - Statistical significance: p < 0.05 (FDR-corrected)?

2. [ ] **Define pathway mapping validation**:
   - Known cases test: BRCA1 → DDR should be high
   - Pathway score improvement: SAE scores better than proxy?
   - Biological plausibility: Do scores match known biology?

3. [ ] **Define end-to-end validation**:
   - Resistance prediction accuracy: ≥70%?
   - Drug ranking improvement: SAE improves vs baseline?
   - Clinical alignment: Do recommendations match guidelines?

**After Each Stage**:
4. [ ] **Stage 2 validation**: Re-run analysis, check for significant features
5. [ ] **Stage 3 validation**: Test mapping on known cases
6. [ ] **Stage 4 validation**: Compare SAE vs proxy pathway scores
7. [ ] **Stage 5 validation**: Test resistance prediction accuracy

### Code References

- **Manager policy**: `.cursor/ayesha/MANAGER_ANSWERS_TO_ZO_SAE_QUESTIONS.md`
- **Mechanism fit criteria**: `.cursor/ayesha/MECHANISM_FIT_VALIDATION_GAPS.md:45-54`
- **Validation scripts**: `scripts/validate_sae_tcga.py`

---

## 📋 SUMMARY OF FINDINGS

### ✅ Resolved (High Confidence)

1. **Data Quality**: Outcome field populated, feature indices valid
2. **Evo2 SAE Pattern**: Layer name correct, architecture matches notebook
3. **Manager Policy**: Formulas and thresholds documented
4. **Test Infrastructure**: Unit tests, E2E tests, validation scripts exist

### ⚠️ Partially Resolved (Medium Confidence)

1. **Modal Service Deployment**: Code exists, deployment status unclear
2. **Validation Criteria**: Some defined (mechanism fit), others missing (biomarker, pathway)
3. **Error Handling**: Patterns exist, coverage unclear

### 🔴 Unresolved (Low Confidence)

1. ~~**Random SAE Weights**: Unknown if meaningful~~ ✅ **RESOLVED**: Using trained weights
2. **Feature→Pathway Mapping**: No biological basis (needs strategy) - **CRITICAL BLOCKER**
3. **Biomarker Analysis Success**: Will it find features? (needs re-run with trained weights)

---

## 🎯 IMMEDIATE ACTION PLAN

### Step 1: Verify Data Quality (30 minutes)

```bash
# Run verification commands
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main

# 1. Verify outcome distribution
python3 -c "import json; from collections import Counter; data = json.load(open('data/validation/sae_cohort/sae_features_tcga_ov_platinum.json')); print('Outcomes:', Counter([p.get('outcome') for p in data['patients']]))"

# 2. Verify feature indices
python3 -c "import json; data = json.load(open('data/validation/sae_cohort/sae_features_tcga_ov_platinum.json')); indices = [f['index'] for p in data['patients'] for v in p['variants'] for f in v.get('top_features', [])]; print(f'Indices valid: {all(0 <= i < 32768 for i in indices)}, Range: {min(indices)}-{max(indices)}')"

# 3. Check feature distributions
python3 -c "import json, numpy as np; data = json.load(open('data/validation/sae_cohort/sae_features_tcga_ov_platinum.json')); values = [f['value'] for p in data['patients'] for v in p['variants'] for f in v.get('top_features', [])]; print(f'Feature stats: mean={np.mean(values):.2e}, std={np.std(values):.2e}, unique={len(set(indices))}')"
```

### Step 2: Test Modal Service (30 minutes)

```bash
# 1. Check health endpoint
curl http://localhost:8000/api/sae/health

# 2. Check environment variables
env | grep SAE

# 3. If service URL exists, test directly
curl -X POST $SAE_SERVICE_URL/extract_features -H "Content-Type: application/json" -d '{"chrom": "17", "pos": 43044295, "ref": "T", "alt": "G", "assembly": "GRCh38", "window": 8192}'
```

### Step 3: Review Evo2 Notebook (2 hours)

- [ ] Read entire notebook (1780 lines)
- [ ] Document feature interpretation approach
- [ ] Compare with our implementation
- [ ] Identify any gaps

### Step 4: Define Validation Criteria (1 hour)

- [ ] Biomarker analysis success criteria
- [ ] Pathway mapping validation strategy
- [ ] End-to-end validation plan

### Step 5: Re-Run Biomarker Analysis (After Steps 1-4)

```bash
python3 scripts/sae/analyze_biomarkers.py \
  --input data/validation/sae_cohort/sae_features_tcga_ov_platinum.json \
  --output data/validation/sae_cohort/sae_tcga_ov_platinum_biomarkers.json \
  --plots-dir data/validation/sae_cohort/plots
```

---

## 🚨 CRITICAL RISKS REMAINING

1. ~~**Random SAE Weights May Be Noise**~~ ✅ **RESOLVED**: Migrated to evo2_7b with trained weights

2. **Feature→Pathway Mapping Has No Biological Basis** (High Impact)
   - **Mitigation**: Use gene→pathway mapping to infer feature→pathway
   - **Fallback**: Keep proxy pathway scores until biological validation

3. **Biomarker Analysis May Find 0 Features** (High Impact)
   - **Mitigation**: Verify data quality, check feature distributions (trained weights now in use)
   - **Fallback**: May need more patients or different analysis approach

---

## 📝 NEXT ITERATION

After completing immediate action items, update this document with:
- Actual data quality results
- Modal service deployment status
- Feature distribution analysis
- Validation criteria definitions
- Re-run biomarker analysis results

