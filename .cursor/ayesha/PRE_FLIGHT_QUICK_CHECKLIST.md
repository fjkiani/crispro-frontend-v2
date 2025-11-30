# ⚡ PRE-FLIGHT QUICK CHECKLIST

**Date**: January 27, 2025  
**Purpose**: Quick reference for pre-flight checks before full pipeline execution

---

## 🚦 GO/NO-GO DECISION TREE

```
START
  ↓
[1] Data Quality Check
  ├─ ✅ PASS → Continue
  └─ ❌ FAIL → STOP (fix data)
  ↓
[2] Service Health Check
  ├─ ✅ PASS → Continue
  └─ ❌ FAIL → STOP (fix services)
  ↓
[3] Pipeline Integration Check
  ├─ ✅ PASS → Continue
  └─ ❌ FAIL → STOP (fix pipeline)
  ↓
[4] MBD4+TP53 Specific Check
  ├─ ✅ PASS → Continue
  └─ ❌ FAIL → STOP (investigate)
  ↓
[5] Full Pipeline Execution
  └─ ✅ PROCEED
```

---

## 📋 EXECUTION ORDER (75 minutes)

### Step 1: Data Quality (5 min)
```bash
python3 scripts/sae/health_check_data.py
python3 scripts/sae/health_check_feature_distributions.py
```
**Expected**: ✅ PASS (warnings OK for small dataset)

### Step 2: Service Health (10 min)
```bash
# Check Modal services
curl https://crispro--sae-service-saeservice-api.modal.run/health
curl https://crispro--evo-service-evoservice7b-api-7b.modal.run/health

# Start backend (if not running)
cd oncology-coPilot/oncology-backend-minimal
uvicorn main:app --reload --port 8000

# Check backend
python3 scripts/sae/health_check_backend.py
```
**Expected**: ✅ All services healthy

### Step 3: Pipeline Integration (20 min)
```bash
python3 scripts/sae/health_check_pipeline.py
pytest tests/test_ayesha_post_ngs_e2e.py -v
pytest tests/test_sae_phase2_services.py -v
```
**Expected**: ✅ All tests pass

### Step 4: MBD4+TP53 Validation (15 min)
```bash
python3 scripts/sae/health_check_mbd4.py
python3 scripts/sae/health_check_pathways.py
python3 scripts/sae/health_check_feature_pathway_mapping.py
```
**Expected**: ✅ Variants score correctly, DDR pathway high, mapping file valid

### Step 5: Full Pipeline (15 min) - ONLY IF ALL CHECKS PASS
```bash
python3 scripts/sae/run_mbd4_tp53_analysis.py
python3 scripts/sae/answer_mbd4_clinical_questions.py
python3 scripts/sae/verify_mbd4_analysis.py
```
**Expected**: ✅ Complete analysis with all 8 questions answered

---

## 🚨 STOP CONDITIONS

**STOP IMMEDIATELY IF**:
- ❌ Data file missing or corrupted
- ❌ Modal services down
- ❌ Backend won't start
- ❌ Pipeline integration fails
- ❌ MBD4 variants don't score correctly
- ❌ Unexpected pathway scores

**DO NOT PROCEED** until all issues resolved.

---

## ✅ SUCCESS CRITERIA

**All Must Pass**:
- ✅ Data quality: File exists, structure valid
- ✅ Service health: All services operational
- ✅ Pipeline: Evo2 → SAE works
- ✅ MBD4+TP53: Variants score, pathways correct

**Only Then**: Run full pipeline

---

## 📊 CURRENT STATUS

| Checkpoint | Status | Action |
|------------|--------|--------|
| Data Quality | ✅ READY | Run `health_check_data.py` |
| Feature Distributions | ✅ READY | Run `health_check_feature_distributions.py` |
| Modal Services | ⏸️ PENDING | Check manually |
| Backend Services | ⏸️ PENDING | Start backend, run check |
| Pipeline Integration | ⏸️ PENDING | Run after backend starts |
| MBD4 Variants | ⏸️ PENDING | Run after backend starts |
| Pathway Analysis | ⏸️ PENDING | Run after backend starts |
| Feature Mapping | ✅ READY | Run `health_check_feature_pathway_mapping.py` |
| Full Pipeline | ⏸️ BLOCKED | Wait for all checks to pass |

---

**Next Action**: Start backend, then run all health checks in order.

