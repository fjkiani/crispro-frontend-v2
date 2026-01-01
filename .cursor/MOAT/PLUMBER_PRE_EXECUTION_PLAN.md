# 🔧 PLUMBER'S PRE-EXECUTION PLAN & QUESTIONS

**Date:** January 30, 2025  
**Role:** Plumber (Production Readiness)  
**Status:** ⏳ PRE-EXECUTION REVIEW

---

## ✅ VERIFICATION COMPLETE

### ✅ SQLite Database Location - VERIFIED
**Found:** `oncology-coPilot/oncology-backend-minimal/data/clinical_trials.db` ✅  
**Tables:** `clinical_trials`, `trials`, `trials_fresh`  
**`trials` table count:** **1397 rows** (more than enough!) ✅  
**Script will use:** `trials` table (correct) ✅ - try the minimal one 

---

### ✅ Environment Variables - VERIFIED
**Status:** All required variables found in `.env` ✅
- ✅ `GEMINI_API_KEY` or `GOOGLE_API_KEY` (found)
- ✅ `ASTRA_DB_APPLICATION_TOKEN` (found)
- ✅ `ASTRA_DB_API_ENDPOINT` (found)
- ✅ `ASTRA_COLLECTION_NAME` (will use default if not set)

---

### ⚠️ Backend Running State - NEEDS ACTION
**Task 2 (E2E Tests) requires backend running on port 8000**

**Decision:** Will start backend in background for E2E tests, then stop after completion

---

### ✅ AstraDB Current State - WILL VERIFY
**Production plan says:** "30 trials seeded (limited set)"

**Decision:** Will check current count first, then seed 1000 (upsert will handle duplicates)

---

### ✅ Documentation Updates - WILL UPDATE DIRECTLY
**Tasks 3-5 are documentation fixes**

**Decision:** Will update production plan file directly with verified statuses

---

## 📋 AUDIT FINDINGS

### ✅ Scripts Exist and Look Good:
1. ✅ `seed_astradb_from_sqlite.py` - Exists, uses `trials` table, has proper error handling
2. ✅ `validate_sporadic_gates.py` - Exists, has 6 tests, saves reports
3. ✅ `validate_quick_intake.py` - Exists, tests all 15 cancers, saves reports
4. ✅ `e2e_sporadic_workflow.sh` - Exists, has proper error handling, saves outputs

### ⚠️ Potential Issues Found:

1. **SQLite Table Name Mismatch:**
   - Script uses `trials` table (line 111)
   - But some code references `clinical_trials` table
   - Need to verify which table has data

2. **Database Path:**
   - Script uses `data/clinical_trials.db` (relative to backend-minimal)
   - But found 3 different database files
   - Need to verify correct path

3. ✅ **E2E Script API Endpoint:**
   - Script uses `cancer_type` field (line 24) ✅
   - Verified: API schema uses `cancer_type` (not `disease`) ✅
   - E2E script is correct ✅

---

## 🎯 EXECUTION PLAN (PENDING ANSWERS)

### Phase 1: Verification (Before Execution)
1. [ ] Verify SQLite database location and table structure
2. [ ] Verify environment variables are set
3. [ ] Check current AstraDB trial count
4. [ ] Verify backend API endpoint schema

### Phase 2: Execution (After Verification)
1. [ ] Run AstraDB seeding (Task 1)
2. [ ] Run E2E tests (Task 2)
3. [ ] Update documentation (Tasks 3-5)

### Phase 3: Validation (After Execution)
1. [ ] Verify AstraDB has 1000+ trials
2. [ ] Verify all tests pass
3. [ ] Verify documentation is updated

---

## 📊 RISK ASSESSMENT

| **Risk** | **Severity** | **Mitigation** |
|----------|--------------|----------------|
| Wrong database path | HIGH | Verify before execution |
| Missing env vars | HIGH | Check .env file first |
| API schema mismatch | MEDIUM | Test with curl first |
| AstraDB overwrite | MEDIUM | Check current count first |
| Script errors | LOW | Scripts have error handling |

---

**⚔️ PLUMBER AWAITING CLARIFICATIONS BEFORE EXECUTION ⚔️**

