# 🛡️ PRE-FLIGHT QUALITY ASSURANCE - MASTER DOCUMENT

**Date**: January 27, 2025  
**Status**: ✅ **PLAN COMPLETE** - Ready for execution  
**Last Updated**: January 28, 2025

**Consolidated From:**
- `.cursor/ayesha/PRE_FLIGHT_SUMMARY.md` (executive summary)
- `.cursor/ayesha/PRE_FLIGHT_QUALITY_ASSURANCE_PLAN.md` (comprehensive plan)
- `.cursor/ayesha/PRE_FLIGHT_QUICK_CHECKLIST.md` (quick reference)

---

## 🎯 EXECUTIVE SUMMARY

**Critical Principle**: **Never run a full pipeline without verifying data quality and system readiness first.**

This master document establishes:
1. **Data Quality Checkpoints** - Verify inputs are valid before processing
2. **Service Health Checks** - Ensure all services are operational
3. **Pipeline Validation** - Test each stage independently
4. **Resource Protection** - Stop early if quality issues detected
5. **Clear Go/No-Go Gates** - Explicit criteria for proceeding

**Total Time**: 75 minutes (pre-flight) + 15 minutes (pipeline) = 90 minutes

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

**Action**: Run this to verify feature quality (may show warnings with small dataset - expected).

---

### Checkpoint 1.3: Feature→Pathway Mapping Validation ✅ **READY**

**Script**: `scripts/sae/health_check_feature_pathway_mapping.py` (NEW)

**What It Checks**:
- ✅ Mapping file exists: `api/resources/sae_feature_mapping.json`
- ✅ Mapping structure valid: 88 features → 4 pathways
- ✅ Pathway coverage: TP53, PI3K, VEGF, HER2
- ✅ Known case validation: BRCA1 → DDR (if data available)

**Go/No-Go Criteria**:
- ✅ **GO**: Mapping file exists, structure valid, pathways covered
- ⚠️ **WARN**: Missing DDR/MAPK pathways (known limitation)
- ❌ **STOP**: Mapping file missing, invalid structure

**Action**: Run this to verify Feature→Pathway mapping is ready.

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

# Evo2 7B Service
curl https://crispro--evo-service-evoservice7b-api-7b.modal.run/health
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

**Action**: Run this to verify pathway analysis works for MBD4+TP53.

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

### Gate 2: Service Health Gate ⏸️ **REQUIRES VERIFICATION**

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

## 🚦 EXECUTION FLOW

```
START
  ↓
[1] Data Quality (5 min)
  ├─ health_check_data.py
  └─ health_check_feature_distributions.py
  ↓
[2] Service Health (10 min)
  ├─ Modal services (curl)
  ├─ Start backend
  └─ health_check_backend.py
  ↓
[3] Pipeline Integration (20 min)
  ├─ health_check_pipeline.py
  └─ pytest tests
  ↓
[4] MBD4+TP53 Validation (15 min)
  ├─ health_check_mbd4.py
  ├─ health_check_pathways.py
  └─ health_check_feature_pathway_mapping.py
  ↓
[5] Full Pipeline (15 min) - ONLY IF ALL PASS
  ├─ run_mbd4_tp53_analysis.py
  ├─ answer_mbd4_clinical_questions.py
  └─ verify_mbd4_analysis.py
```

**Total Time**: 75 minutes (pre-flight) + 15 minutes (pipeline) = 90 minutes

---

## 🚨 STOP CONDITIONS

**STOP IMMEDIATELY IF**:
- ❌ Data file missing or corrupted
- ❌ Modal services down
- ❌ Backend won't start
- ❌ Pipeline integration fails
- ❌ MBD4 variants don't score correctly
- ❌ Unexpected pathway scores
- ❌ Mapping file invalid

**DO NOT PROCEED** until all issues resolved.

---

## ✅ SUCCESS CRITERIA

**All Must Pass**:
- ✅ Data quality: File exists, structure valid, outcome populated
- ✅ Service health: Modal + Backend operational
- ✅ Pipeline: Evo2 → SAE works end-to-end
- ✅ MBD4+TP53: Variants score, pathways correct, mapping valid

**Only Then**: Run full pipeline

---

## 📊 CURRENT STATUS

| Component | Status | Script | Action |
|-----------|--------|--------|--------|
| Data Quality | ✅ READY | `health_check_data.py` | Run first |
| Feature Distributions | ✅ READY | `health_check_feature_distributions.py` | Run second |
| Modal Services | ⏸️ PENDING | Manual curl | Check before backend |
| Backend Services | ⏸️ PENDING | `health_check_backend.py` | Run after backend starts |
| Pipeline Integration | ⏸️ PENDING | `health_check_pipeline.py` | Run after backend starts |
| MBD4 Variants | ⏸️ PENDING | `health_check_mbd4.py` | Run after backend starts |
| Pathway Analysis | ⏸️ PENDING | `health_check_pathways.py` | Run after backend starts |
| Feature Mapping | ✅ READY | `health_check_feature_pathway_mapping.py` | Run anytime |
| Full Pipeline | ⏸️ BLOCKED | `run_mbd4_tp53_analysis.py` | Wait for all checks |

---

## 📋 EXECUTION CHECKLIST

### Pre-Execution (Must Complete Before Full Pipeline)

- [ ] **Data Quality Check**: `health_check_data.py` - ✅ PASS
- [ ] **Feature Distribution Check**: `health_check_feature_distributions.py` - ⚠️ WARN (expected)
- [ ] **Modal Services Check**: Manual curl commands - ⏸️ PENDING
- [ ] **Backend Services Check**: `health_check_backend.py` - ⏸️ PENDING
- [ ] **Pipeline Integration Check**: `health_check_pipeline.py` - ⏸️ PENDING
- [ ] **MBD4 Variant Check**: `health_check_mbd4.py` - ⏸️ PENDING
- [ ] **Pathway Analysis Check**: `health_check_pathways.py` - ⏸️ PENDING
- [ ] **Feature Mapping Check**: `health_check_feature_pathway_mapping.py` - ⏸️ PENDING

### Execution (Only If All Pre-Execution Checks Pass)

- [ ] **Full Pipeline Run**: `run_mbd4_tp53_analysis.py` - ⏸️ PENDING
- [ ] **Question Answering**: `answer_mbd4_clinical_questions.py` - ⏸️ PENDING
- [ ] **Verification**: `verify_mbd4_analysis.py` - ⏸️ PENDING

---

## 🚨 FAILURE MODES & MITIGATION

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

## 📝 KEY PRINCIPLES

1. **Never run full pipeline without pre-flight checks**
2. **Stop early if any critical check fails**
3. **Document all warnings and limitations**
4. **Verify each stage independently**
5. **Protect resources by catching issues early**

---

## 🎯 NEXT STEPS

1. **Review Plans**: Read this master document for details
2. **Start Execution**: Follow the execution order (Step 1 → Step 5)
3. **Stop If Any Fail**: Do not proceed until all checks pass

---

**Status**: ✅ **READY FOR EXECUTION**

**Estimated Time**: 90 minutes total (75 min pre-flight + 15 min pipeline)

**Critical**: Do not skip pre-flight checks. They prevent wasted resources and ensure quality results.

