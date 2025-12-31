# 🔧 PLUMBER'S PRE-EXECUTION PLAN & QUESTIONS

**Date:** January 30, 2025  
**Role:** Plumber (Production Readiness)  
**Status:** ⏳ PRE-EXECUTION REVIEW

---

## ❓ CLARIFYING QUESTIONS (NEVER ASSUME)

### Question 1: SQLite Database Location
**Issue:** Found 3 different SQLite databases:
- `oncology-coPilot/oncology-backend-minimal/data/clinical_trials.db`
- `oncology-coPilot/oncology-backend/data/clinical_trials.db`
- `oncology-coPilot/oncology-backend/backend/data/clinical_trials.db`

**Script expects:** `oncology-coPilot/oncology-backend-minimal/data/clinical_trials.db` (uses `trials` table)

**Question:** Which database has the 1000 trials? Should I verify before seeding? - try the minimal one 

---

### Question 2: Environment Variables
**Required:**
- `GEMINI_API_KEY` or `GOOGLE_API_KEY` (for embeddings)
- `ASTRA_DB_APPLICATION_TOKEN` (for AstraDB)
- `ASTRA_DB_API_ENDPOINT` (for AstraDB)
- `ASTRA_COLLECTION_NAME` (optional, defaults to "clinical_trials_eligibility")

**Question:** Are these set in `.env` file? Should I verify before running scripts?

---

### Question 3: Backend Running State
**Task 2 (E2E Tests) requires backend running on port 8000**

**Question:** Should I start the backend, or is it already running? What's the preferred approach?

---

### Question 4: AstraDB Current State
**Production plan says:** "30 trials seeded (limited set)"

**Question:** Should I:
- A) Clear existing 30 trials and seed fresh 1000?
- B) Keep existing 30 and add 970 more?
- C) Verify current count first?

---

### Question 5: Documentation Updates
**Tasks 3-5 are documentation fixes**

**Question:** Should I:
- A) Update the production plan file directly?
- B) Create a separate update document?
- C) Just verify the fixes are correct?

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

3. **E2E Script API Endpoint:**
   - Script uses `cancer_type` field (line 24)
   - But production plan shows `disease` field
   - Need to verify correct API schema

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

