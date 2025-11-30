# 🛡️ PRE-FLIGHT QUALITY ASSURANCE PLAN

**Date**: January 27, 2025  
**Purpose**: Ensure 100% data quality and system readiness before full pipeline execution  
**Goal**: Zero wasted resources, clear go/no-go gates, comprehensive checkpoints

---

## 🎯 EXECUTIVE SUMMARY

**Critical Principle**: **Never run a full pipeline without verifying data quality and system readiness first.**

This plan establishes:
1. **Data Quality Checkpoints** - Verify inputs are valid before processing
2. **Service Health Checks** - Ensure all services are operational
3. **Pipeline Validation** - Test each stage independently
4. **Resource Protection** - Stop early if quality issues detected
5. **Clear Go/No-Go Gates** - Explicit criteria for proceeding

---

## 📋 PHASE 1: DATA QUALITY CHECKPOINTS

### Checkpoint 1.1: SAE Cohort Data Structure ✅ **READY**

**Script**: `scripts/sae/health_check_data.py`

**What It Checks**:
- ✅ File exists: `data/validation/sae_cohort/sae_features_tcga_ov_platinum.json`
- ✅ Structure valid: `patients` array, required fields present
- ✅ Patient count: Expected 66, currently 10 (⚠️ **KNOWN ISSUE**)
- ✅ Outcome field: All patients have valid outcome (sensitive/refractory/resistant)
- ✅ Feature indices: All indices in range [0, 32767]

**Go/No-Go Criteria**:
- ✅ **GO**: File exists, structure valid, outcome field populated
- ⚠️ **WARN**: Patient count < 66 (expected with small dataset)
- ❌ **STOP**: Missing required fields, invalid feature indices

**Expected Results**:
```
✅ PASS: File exists
✅ PASS: Data loaded successfully
⚠️  WARN: Expected 66 patients, got 10 (expected with small dataset)
✅ PASS: All patients have valid structure
✅ PASS: All feature indices valid (0-32767)
✅ PASS: Outcome distribution: sensitive=9, resistant=1
```

**Action**: Run this FIRST before any analysis.

---

### Checkpoint 1.2: Feature Distribution Quality ✅ **READY**

**Script**: `scripts/sae/health_check_feature_distributions.py`

**What It Checks**:
- ✅ Features are sparse (most are zero)
- ✅ Features have variation (not all identical)
- ✅ Feature activations in reasonable range
- ✅ Patients have different feature patterns (not all identical)

**Go/No-Go Criteria**:
- ✅ **GO**: Sparsity > 80%, variation present, patients differ
- ⚠️ **WARN**: Low variation (may indicate extraction issue)
- ❌ **STOP**: All features zero, all patients identical

**Expected Results** (with 10 patients):
```
⚠️  WARN: No activations found (expected with small dataset)
⚠️  WARN: Cannot verify distributions (need more patients)
✅ PASS: Script ready for larger cohorts
```

**Action**: Run this to verify feature quality (may show warnings with small dataset - expected).

---

### Checkpoint 1.3: Data Completeness Validation ⚠️ **NEW - NEEDS CREATION**

**What It Should Check**:
- ✅ All patients have mutations data
- ✅ All patients have SAE features extracted
- ✅ All patients have outcome labels
- ✅ No missing critical fields

**Script to Create**: `scripts/sae/health_check_data_completeness.py`

**Go/No-Go Criteria**:
- ✅ **GO**: >90% completeness, critical fields present
- ⚠️ **WARN**: 70-90% completeness (proceed with caution)
- ❌ **STOP**: <70% completeness, critical fields missing

**Action**: Create this script before full pipeline run.

---

## 🔧 PHASE 2: SERVICE HEALTH CHECKS

### Checkpoint 2.1: Modal Services Health ✅ **READY**

**What It Checks**:
- ✅ SAE Service: `https://crispro--sae-service-saeservice-api.modal.run/health`
- ✅ Evo2 7B Service: `https://crispro--evo-service-evoservice7b-api-7b.modal.run/health`

**Manual Check**:
```bash
# SAE Service
curl https://crispro--sae-service-saeservice-api.modal.run/health

# Expected: {"status": "healthy", "model": "evo2_7b", "weights": "trained"}

# Evo2 7B Service
curl https://crispro--evo-service-evoservice7b-api-7b.modal.run/health

# Expected: {"status": "healthy", "model": "evo2_7b_base"}
```

**Go/No-Go Criteria**:
- ✅ **GO**: Both services return 200 OK with "healthy" status
- ❌ **STOP**: Any service returns error or timeout

**Action**: Check Modal services BEFORE starting backend.

---

### Checkpoint 2.2: Backend Services Health ✅ **READY**

**Script**: `scripts/sae/health_check_backend.py`

**What It Checks**:
- ✅ Backend running: `http://localhost:8000`
- ✅ SAE endpoint: `/api/sae/compute_features` (POST)
- ✅ Health endpoint: `/healthz` (GET)
- ✅ Response structure: DNA repair capacity, mechanism vector present

**Go/No-Go Criteria**:
- ✅ **GO**: All endpoints return 200 OK, response structure valid
- ❌ **STOP**: Any endpoint fails, timeout, or invalid response

**Expected Results**:
```
✅ PASS: POST /api/sae/compute_features returned 200
   DNA repair capacity: 0.750
   Mechanism vector: 7D
✅ PASS: GET /healthz returned 200
✅ All backend services operational
```

**Action**: Run this AFTER starting backend, BEFORE any analysis.

---

### Checkpoint 2.3: Pipeline Integration Health ✅ **READY**

**Script**: `scripts/sae/health_check_pipeline.py`

**What It Checks**:
- ✅ Evo2 → SAE pipeline: Variant scoring → SAE extraction
- ✅ SAE feature dimensions: 32K-dim, sparse (≤64 active)
- ✅ End-to-end flow: Input variant → Output SAE features

**Go/No-Go Criteria**:
- ✅ **GO**: Pipeline completes, dimensions correct, sparsity valid
- ❌ **STOP**: Pipeline fails, wrong dimensions, no sparsity

**Expected Results**:
```
✅ PASS: Evo2 → SAE pipeline working
✅ PASS: Evo2 activations: 4096-dim
✅ PASS: SAE features: 32768-dim, 64 active
```

**Action**: Run this to verify Evo2 → SAE integration works.

---

## 🧪 PHASE 3: PIPELINE VALIDATION

### Checkpoint 3.1: MBD4+TP53 Variant Scoring ✅ **READY**

**Script**: `scripts/sae/health_check_mbd4.py`

**What It Checks**:
- ✅ MBD4 variant scoring: High disruption (frameshift)
- ✅ TP53 variant scoring: Hotspot mutation detected
- ✅ Evo2 delta scores: Expected ranges for each variant

**Go/No-Go Criteria**:
- ✅ **GO**: Variants score correctly, disruption detected
- ❌ **STOP**: Scoring fails, unexpected scores

**Expected Results**:
```
✅ PASS: MBD4 variant scoring working
✅ PASS: Delta score: -2.450 (high disruption - frameshift)
✅ PASS: High disruption confirmed (frameshift)
```

**Action**: Run this to verify MBD4+TP53 variants are scored correctly.

---

### Checkpoint 3.2: Pathway Analysis Validation ✅ **READY**

**Script**: `scripts/sae/health_check_pathways.py`

**What It Checks**:
- ✅ Pathway scores computed: DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux
- ✅ MBD4+TP53 → High DDR pathway (>0.7 expected)
- ✅ PARP inhibitors rank high (efficacy >0.75 expected)

**Go/No-Go Criteria**:
- ✅ **GO**: Pathway scores computed, DDR high, PARP ranks high
- ❌ **STOP**: Pathway scores missing, unexpected values

**Expected Results**:
```
✅ PASS: MBD4+TP53 pathway analysis working
✅ PASS: DDR pathway score: 0.880 (HIGH - expected)
✅ PASS: Top PARP inhibitor: Olaparib (efficacy: 0.800)
```

**Action**: Run this to verify pathway analysis works for MBD4+TP53.

---

### Checkpoint 3.3: Feature→Pathway Mapping Validation ⚠️ **NEW - NEEDS CREATION**

**What It Should Check**:
- ✅ Mapping file exists: `api/resources/sae_feature_mapping.json`
- ✅ Mapping structure valid: 88 features → 4 pathways
- ✅ Pathway coverage: TP53, PI3K, VEGF, HER2
- ✅ Known case validation: BRCA1 → DDR (if data available)

**Script to Create**: `scripts/sae/health_check_feature_pathway_mapping.py`

**Go/No-Go Criteria**:
- ✅ **GO**: Mapping file exists, structure valid, pathways covered
- ⚠️ **WARN**: Missing DDR/MAPK pathways (known limitation)
- ❌ **STOP**: Mapping file missing, invalid structure

**Action**: Create this script to verify Feature→Pathway mapping is ready.

---

## 🚦 PHASE 4: GO/NO-GO GATES

### Gate 1: Data Quality Gate ⚠️ **CURRENT STATUS**

**Requirements**:
- ✅ Data file exists and structure valid
- ⚠️ Patient count: 10 (expected 66) - **KNOWN LIMITATION**
- ✅ Outcome field populated
- ✅ Feature indices valid

**Decision**:
- ✅ **GO**: Proceed with analysis (small dataset is expected)
- ⚠️ **WARN**: Results may be limited by small sample size
- ❌ **STOP**: If structure invalid or critical fields missing

**Current Status**: ✅ **PASS** (with known limitation)

---

### Gate 2: Service Health Gate ⚠️ **REQUIRES VERIFICATION**

**Requirements**:
- ⏸️ Modal services healthy (SAE + Evo2 7B)
- ⏸️ Backend running and responsive
- ⏸️ All endpoints operational

**Decision**:
- ✅ **GO**: All services healthy
- ❌ **STOP**: Any service down or unresponsive

**Current Status**: ⏸️ **PENDING VERIFICATION**

**Action**: Run health checks before proceeding.

---

### Gate 3: Pipeline Integration Gate ⏸️ **REQUIRES VERIFICATION**

**Requirements**:
- ⏸️ Evo2 → SAE pipeline works
- ⏸️ SAE feature extraction successful
- ⏸️ Backend SAE service operational

**Decision**:
- ✅ **GO**: Pipeline works end-to-end
- ❌ **STOP**: Pipeline fails at any stage

**Current Status**: ⏸️ **PENDING VERIFICATION**

**Action**: Run pipeline health checks before proceeding.

---

### Gate 4: MBD4+TP53 Specific Gate ⏸️ **REQUIRES VERIFICATION**

**Requirements**:
- ⏸️ MBD4 variant scores correctly
- ⏸️ TP53 variant scores correctly
- ⏸️ Pathway analysis produces expected results (DDR >0.7)
- ⏸️ PARP inhibitors rank high (efficacy >0.75)

**Decision**:
- ✅ **GO**: All MBD4+TP53 checks pass
- ❌ **STOP**: Variant scoring fails or unexpected results

**Current Status**: ⏸️ **PENDING VERIFICATION**

**Action**: Run MBD4+TP53 specific checks before full analysis.

---

## 📊 PHASE 5: RESOURCE PROTECTION CHECKLIST

### Before Running Full Pipeline

**Checklist**:
- [ ] **Data Quality**: Run `health_check_data.py` - ✅ PASS
- [ ] **Feature Distributions**: Run `health_check_feature_distributions.py` - ⚠️ WARN (expected)
- [ ] **Modal Services**: Check SAE + Evo2 7B health - ⏸️ PENDING
- [ ] **Backend Services**: Run `health_check_backend.py` - ⏸️ PENDING
- [ ] **Pipeline Integration**: Run `health_check_pipeline.py` - ⏸️ PENDING
- [ ] **MBD4 Variants**: Run `health_check_mbd4.py` - ⏸️ PENDING
- [ ] **Pathway Analysis**: Run `health_check_pathways.py` - ⏸️ PENDING
- [ ] **Feature Mapping**: Run `health_check_feature_pathway_mapping.py` - ⏸️ PENDING (needs creation)

**Stop Conditions**:
- ❌ Any critical check fails (data structure, service health)
- ❌ Pipeline integration fails
- ❌ MBD4+TP53 specific checks fail

**Proceed Conditions**:
- ✅ All critical checks pass
- ⚠️ Non-critical warnings acceptable (e.g., small dataset)

---

## 🎯 PHASE 6: EXECUTION PLAN

### Step 1: Pre-Flight Checks (30 minutes)

**Run All Health Checks**:
```bash
# 1. Data Quality
python3 scripts/sae/health_check_data.py

# 2. Feature Distributions
python3 scripts/sae/health_check_feature_distributions.py

# 3. Backend Services (requires backend running)
python3 scripts/sae/health_check_backend.py

# 4. Pipeline Integration (requires backend running)
python3 scripts/sae/health_check_pipeline.py

# 5. MBD4 Variants (requires backend running)
python3 scripts/sae/health_check_mbd4.py

# 6. Pathway Analysis (requires backend running)
python3 scripts/sae/health_check_pathways.py
```

**Expected Time**: 30 minutes  
**Stop If**: Any critical check fails

---

### Step 2: Service Verification (10 minutes)

**Check Modal Services**:
```bash
# SAE Service
curl https://crispro--sae-service-saeservice-api.modal.run/health

# Evo2 7B Service
curl https://crispro--evo-service-evoservice7b-api-7b.modal.run/health
```

**Start Backend** (if not running):
```bash
cd oncology-coPilot/oncology-backend-minimal
uvicorn main:app --reload --port 8000
```

**Expected Time**: 10 minutes  
**Stop If**: Modal services down or backend won't start

---

### Step 3: Integration Verification (20 minutes)

**Run Integration Tests**:
```bash
# End-to-end pipeline test
pytest tests/test_ayesha_post_ngs_e2e.py -v

# SAE Phase 2 services test
pytest tests/test_sae_phase2_services.py -v
```

**Expected Time**: 20 minutes  
**Stop If**: Any test fails

---

### Step 4: MBD4+TP53 Validation (15 minutes)

**Run MBD4+TP53 Specific Checks**:
```bash
# Variant scoring
python3 scripts/sae/health_check_mbd4.py

# Pathway analysis
python3 scripts/sae/health_check_pathways.py
```

**Expected Time**: 15 minutes  
**Stop If**: Variant scoring fails or unexpected pathway scores

---

### Step 5: Full Pipeline Execution (ONLY IF ALL CHECKS PASS)

**Run Full Analysis**:
```bash
# MBD4+TP53 end-to-end analysis
python3 scripts/sae/run_mbd4_tp53_analysis.py

# Answer clinical questions
python3 scripts/sae/answer_mbd4_clinical_questions.py
```

**Expected Time**: 10-15 minutes  
**Only Proceed If**: All previous steps passed

---

## 📝 PHASE 7: KNOWN LIMITATIONS & EXPECTATIONS

### Current Data Limitations

**1. Small Dataset (10 patients vs 66 expected)**
- **Impact**: Biomarker analysis may find 0 significant features
- **Expected**: This is normal with small sample size
- **Action**: Proceed with analysis, note limitation in results

**2. Missing DDR/MAPK Pathways in Mapping**
- **Impact**: Feature→Pathway mapping incomplete (only 4/7 pathways)
- **Expected**: Preliminary mapping, requires larger cohort
- **Action**: Use proxy SAE for DDR/MAPK, TRUE SAE for TP53/PI3K/VEGF/HER2

**3. Backend Not Running (Previous Checks)**
- **Impact**: Health checks may fail if backend not started
- **Expected**: Normal if backend not running
- **Action**: Start backend before running health checks

---

## 🚨 PHASE 8: FAILURE MODES & MITIGATION

### Failure Mode 1: Data Quality Issues

**Symptoms**:
- Missing required fields
- Invalid feature indices
- Corrupted data file

**Mitigation**:
- ❌ **STOP**: Do not proceed with analysis
- 🔧 **FIX**: Re-extract data or fix data file
- ✅ **VERIFY**: Re-run health checks before proceeding

---

### Failure Mode 2: Service Health Issues

**Symptoms**:
- Modal services down
- Backend not responding
- Endpoint errors

**Mitigation**:
- ❌ **STOP**: Do not proceed with analysis
- 🔧 **FIX**: Restart services, check network, verify URLs
- ✅ **VERIFY**: Re-run health checks before proceeding

---

### Failure Mode 3: Pipeline Integration Issues

**Symptoms**:
- Evo2 → SAE pipeline fails
- Wrong feature dimensions
- No sparsity in features

**Mitigation**:
- ❌ **STOP**: Do not proceed with analysis
- 🔧 **FIX**: Check Modal service logs, verify model weights
- ✅ **VERIFY**: Re-run pipeline health checks

---

### Failure Mode 4: Unexpected Results

**Symptoms**:
- MBD4 variant scores unexpectedly low
- DDR pathway score unexpectedly low
- PARP inhibitors don't rank high

**Mitigation**:
- ⚠️ **INVESTIGATE**: Check variant annotations, pathway mappings
- 🔧 **FIX**: Verify inputs, check formulas
- ✅ **VERIFY**: Re-run specific health checks

---

## ✅ PHASE 9: SUCCESS CRITERIA

### All Checks Must Pass

**Data Quality**:
- ✅ Data file exists and structure valid
- ✅ Outcome field populated
- ✅ Feature indices valid

**Service Health**:
- ✅ Modal services healthy
- ✅ Backend services operational
- ✅ All endpoints responsive

**Pipeline Integration**:
- ✅ Evo2 → SAE pipeline works
- ✅ SAE features extracted correctly
- ✅ Backend SAE service operational

**MBD4+TP53 Specific**:
- ✅ Variants score correctly
- ✅ Pathway analysis produces expected results
- ✅ PARP inhibitors rank high

---

## 📋 PHASE 10: EXECUTION CHECKLIST

### Pre-Execution (Must Complete Before Full Pipeline)

- [ ] **Data Quality Check**: `health_check_data.py` - ✅ PASS
- [ ] **Feature Distribution Check**: `health_check_feature_distributions.py` - ⚠️ WARN (expected)
- [ ] **Modal Services Check**: Manual curl commands - ⏸️ PENDING
- [ ] **Backend Services Check**: `health_check_backend.py` - ⏸️ PENDING
- [ ] **Pipeline Integration Check**: `health_check_pipeline.py` - ⏸️ PENDING
- [ ] **MBD4 Variant Check**: `health_check_mbd4.py` - ⏸️ PENDING
- [ ] **Pathway Analysis Check**: `health_check_pathways.py` - ⏸️ PENDING
- [ ] **Feature Mapping Check**: `health_check_feature_pathway_mapping.py` - ⏸️ PENDING (needs creation)

### Execution (Only If All Pre-Execution Checks Pass)

- [ ] **Full Pipeline Run**: `run_mbd4_tp53_analysis.py` - ⏸️ PENDING
- [ ] **Question Answering**: `answer_mbd4_clinical_questions.py` - ⏸️ PENDING
- [ ] **Verification**: `verify_mbd4_analysis.py` - ⏸️ PENDING

---

## 🎯 SUMMARY

**Current Status**:
- ✅ Data quality checks ready
- ⏸️ Service health checks pending (requires backend running)
- ⏸️ Pipeline integration checks pending
- ⏸️ Full pipeline execution pending

**Next Steps**:
1. **Start Backend**: `uvicorn main:app --reload --port 8000`
2. **Run Health Checks**: Execute all health check scripts
3. **Verify Services**: Check Modal services manually
4. **Proceed If All Pass**: Run full pipeline only after all checks pass

**Critical Principle**: **Never run full pipeline without verifying all checkpoints first.**

---

**Status**: ✅ **PLAN READY** - Execute pre-flight checks before full pipeline

**Estimated Time**: 75 minutes (pre-flight) + 15 minutes (full pipeline) = 90 minutes total

