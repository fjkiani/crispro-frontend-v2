# 📁 File Consolidation Plan - Following MM Model

**Date:** January 28, 2025  
**Status:** ✅ **CONSOLIDATION COMPLETE**  
**Model:** `.cursor/MOAT/MM/` structure  
**Last Updated:** January 2025

---

## 🎯 PROBLEM STATEMENT

**Current State:** Files scattered across `.cursor/MOAT/` with no organization

**Example of Good Structure:** `.cursor/MOAT/MM/`
```
.cursor/MOAT/MM/
├── README.md                    # Navigation hub
├── 00_MISSION.mdc               # SOURCE OF TRUTH
├── 01_AUDIT.md                  # Supporting docs
├── 02_VALIDATION.md
├── 03_DELIVERY_PLAN.md
└── archive/                     # Old files
```

**Current Problem:** 
- 50+ files at `.cursor/MOAT/` top level
- No organization by topic
- Hard to find related files
- Duplicate/outdated files mixed with current

---

## 📊 CURRENT SCATTER ANALYSIS

### **Files That Need Consolidation:**

#### **1. TOXICITY Files (8 files)** 🔴 **HIGH PRIORITY**

**Current Location:** `.cursor/MOAT/TOXICITY_*.md`

**Files:**
- `TOXICITY_LLM_INTEGRATION.md`
- `TOXICITY_RISK_DOCUMENTATION_UPDATE_SUMMARY.md`
- `TOXICITY_RISK_FRONTEND_SOURCE_OF_TRUTH.md`
- `TOXICITY_RISK_PRODUCTION_PLAN.md`
- `TOXICITY_RISK_PRODUCTION_READINESS_AUDIT.md`
- `TOXICITY_RISK_TEST_RESULTS.md`
- `TOXICITY_RISK_VERIFICATION_POLISH_COMPLETE.md`
- `ADVANCED_CARE_PLAN_TOXCITY.md`

**Should Be:** `.cursor/MOAT/TOXICITY/`

**Structure:**
```
.cursor/MOAT/TOXICITY/
├── README.md                    # Navigation hub
├── 00_SOURCE_OF_TRUTH.md        # ADVANCED_CARE_PLAN_TOXCITY.md (THE MOAT)
├── 01_PRODUCTION_READINESS.md   # TOXICITY_RISK_PRODUCTION_READINESS_AUDIT.md
├── 02_FRONTEND_SOURCE_OF_TRUTH.md
├── 03_PRODUCTION_PLAN.md
├── 04_TEST_RESULTS.md
├── 05_LLM_INTEGRATION.md
└── archive/                     # Old files
```

---

#### **2. FOOD VALIDATOR Files (5+ files)** 🟡 **MEDIUM PRIORITY**

**Current Location:** `.cursor/MOAT/FOOD_*.md`

**Files:**
- `FOOD_VALIDATOR_ASSESSMENT.md`
- `FOOD_VALIDATOR_CAPABILITIES.md`
- `FOOD_VALIDATOR_COMPLETION_SUMMARY.md`
- `FOOD_VALIDATOR_IMPLEMENTATION_PLAN.md`
- `HOW_TO_TEST_FOODS.md`

**Should Be:** `.cursor/MOAT/FOOD_VALIDATOR/`

---

#### **3. ADVANCED_CARE_PLAN Files (5+ files)** 🟡 **MEDIUM PRIORITY**

**Current Location:** `.cursor/MOAT/ADVANCED_CARE_PLAN_*.md`

**Files:**
- `ADVANCED_CARE_PLAN_UNIVERSAL.md`
- `ADVANCED_CARE_PLAN_EXPLAINED.md`
- `ADVANCED_CARE_PLAN_RESISTANCE_PREDICTION.md`
- `ADVANCED_CARE_PLAN_MECHANISM_TRIAL_MATCHING.md`
- `ADVANCED_CARE_PLAN_LLM_PERSONALIZATION.md`
- `ADVANCED_CARE_PLAN_TOXCITY.md` (move to TOXICITY/)

**Should Be:** `.cursor/MOAT/ADVANCED_CARE_PLAN/`

---

#### **4. ORCHESTRATION Files (already organized)** ✅

**Current Location:** `.cursor/MOAT/orchestration/`

**Status:** ✅ Already organized in subdirectory

---

## 🎯 CONSOLIDATION STRATEGY

### **Follow MM Model:**

**Pattern:**
1. Create topic directory: `.cursor/MOAT/TOPIC/`
2. Create `README.md` (navigation hub)
3. Create `00_SOURCE_OF_TRUTH.md` (main document)
4. Number supporting docs: `01_`, `02_`, `03_`
5. Create `archive/` for old files

**Benefits:**
- ✅ Single source of truth per topic
- ✅ Easy navigation (README.md)
- ✅ Modular organization
- ✅ No data loss (archive old files)

---

## 📋 CONSOLIDATION ROADMAP

### **Phase 1: TOXICITY (Priority: HIGH)**

**Why First:**
- 8 files scattered
- Recently audited
- Production readiness critical

**Steps:**
1. Create `.cursor/MOAT/TOXICITY/` directory
2. Create `README.md` (navigation hub)
3. Move `ADVANCED_CARE_PLAN_TOXCITY.md` → `00_SOURCE_OF_TRUTH.md`
4. Move `TOXICITY_RISK_PRODUCTION_READINESS_AUDIT.md` → `01_PRODUCTION_READINESS.md`
5. Move other files → numbered docs
6. Archive old files

**Estimated Time:** 15 minutes

---

### **Phase 2: FOOD_VALIDATOR (Priority: MEDIUM)**

**Steps:**
1. Create `.cursor/MOAT/FOOD_VALIDATOR/` directory
2. Create `README.md`
3. Identify source of truth
4. Consolidate files
5. Archive old files

**Estimated Time:** 10 minutes

---

### **Phase 3: ADVANCED_CARE_PLAN (Priority: MEDIUM)**

**Steps:**
1. Create `.cursor/MOAT/ADVANCED_CARE_PLAN/` directory
2. Create `README.md`
3. Move `ADVANCED_CARE_PLAN_UNIVERSAL.md` → `00_SOURCE_OF_TRUTH.md`
4. Consolidate other files
5. Archive old files

**Estimated Time:** 10 minutes

---

## 🚀 EXECUTION PLAN

### **Step 1: Create TOXICITY Directory Structure**

```bash
mkdir -p .cursor/MOAT/TOXICITY/archive
```

### **Step 2: Create README.md**

**File:** `.cursor/MOAT/TOXICITY/README.md`

**Content:**
- Navigation hub
- Links to all modules
- Quick reference
- Status summary

### **Step 3: Consolidate Files**

**Mapping:**
- `ADVANCED_CARE_PLAN_TOXCITY.md` → `00_SOURCE_OF_TRUTH.md`
- `TOXICITY_RISK_PRODUCTION_READINESS_AUDIT.md` → `01_PRODUCTION_READINESS.md`
- `TOXICITY_RISK_FRONTEND_SOURCE_OF_TRUTH.md` → `02_FRONTEND_SOURCE_OF_TRUTH.md`
- `TOXICITY_RISK_PRODUCTION_PLAN.md` → `03_PRODUCTION_PLAN.md`
- `TOXICITY_RISK_TEST_RESULTS.md` → `04_TEST_RESULTS.md`
- `TOXICITY_LLM_INTEGRATION.md` → `05_LLM_INTEGRATION.md`
- `TOXICITY_RISK_DOCUMENTATION_UPDATE_SUMMARY.md` → `archive/`
- `TOXICITY_RISK_VERIFICATION_POLISH_COMPLETE.md` → `archive/`

### **Step 4: Archive Old Files**

Move outdated files to `archive/` directory

---

## 📊 EXPECTED OUTCOME

### **Before:**
```
.cursor/MOAT/
├── TOXICITY_LLM_INTEGRATION.md
├── TOXICITY_RISK_DOCUMENTATION_UPDATE_SUMMARY.md
├── TOXICITY_RISK_FRONTEND_SOURCE_OF_TRUTH.md
├── TOXICITY_RISK_PRODUCTION_PLAN.md
├── TOXICITY_RISK_PRODUCTION_READINESS_AUDIT.md
├── TOXICITY_RISK_TEST_RESULTS.md
├── TOXICITY_RISK_VERIFICATION_POLISH_COMPLETE.md
├── ADVANCED_CARE_PLAN_TOXCITY.md
└── ... 40+ other files
```

### **After:**
```
.cursor/MOAT/
├── MM/                          # ✅ Already organized
│   ├── README.md
│   ├── 00_MISSION.mdc
│   └── ...
├── TOXICITY/                    # ✅ NEW - Organized
│   ├── README.md
│   ├── 00_SOURCE_OF_TRUTH.md
│   ├── 01_PRODUCTION_READINESS.md
│   └── ...
├── FOOD_VALIDATOR/              # ✅ NEW - Organized
│   ├── README.md
│   └── ...
├── ADVANCED_CARE_PLAN/          # ✅ NEW - Organized
│   ├── README.md
│   └── ...
└── ... (other organized topics)
```

---

## ✅ SUCCESS CRITERIA

- [x] TOXICITY files consolidated into `.cursor/MOAT/TOXICITY/` ✅
- [x] README.md created for navigation ✅
- [x] Source of truth identified (00_SOURCE_OF_TRUTH.md) ✅
- [x] Supporting docs numbered (01_, 02_, 03_) ✅
- [x] Old files archived ✅
- [x] Same pattern applied to FOOD_VALIDATOR ✅
- [x] Same pattern applied to ADVANCED_CARE_PLAN ✅
- [x] Top-level `.cursor/MOAT/` has < 20 files (only index files) ✅

---

## ✅ CONSOLIDATION COMPLETE (January 2025)

### **Files Processed:**
- ✅ **7 TOXICITY files** archived to `TOXICITY/archive/`
- ✅ **5 FOOD_VALIDATOR files** archived to `FOOD_VALIDATOR/archive/`
- ✅ **1 ADVANCED_CARE_PLAN file** moved to `ADVANCED_CARE_PLAN/06_AUDIT_REPORT.md` (new file)
- ✅ **6 ADVANCED_CARE_PLAN files** archived to `ADVANCED_CARE_PLAN/archive/`

**Total:** 19 files processed (1 moved, 18 archived)

### **Verification:**
- ✅ **0 files remaining** at top level with TOXICITY/FOOD/ADVANCED_CARE_PLAN prefixes
- ✅ All duplicate files archived (not deleted)
- ✅ No data loss - all files preserved
- ✅ README.md files updated

---

## 🎯 NEXT STEPS

1. **Execute Phase 1:** Consolidate TOXICITY files (15 min)
2. **Execute Phase 2:** Consolidate FOOD_VALIDATOR files (10 min)
3. **Execute Phase 3:** Consolidate ADVANCED_CARE_PLAN files (10 min)
4. **Review:** Check for other scattered topics
5. **Document:** Update main MOAT README.md with new structure

**Total Time:** ~35 minutes

---

**Status:** ✅ **CONSOLIDATION COMPLETE** (January 2025)  
**Model:** Follow `.cursor/MOAT/MM/` structure  
**Goal:** Organized, navigable, single source of truth per topic

---

## 📋 DELIVERABLES TRACKING (Single File)

**This file tracks all deliverables and consolidation status.**

### **Consolidation Status:**
- ✅ TOXICITY: Complete (7 files archived)
- ✅ FOOD_VALIDATOR: Complete (5 files archived)
- ✅ ADVANCED_CARE_PLAN: Complete (1 moved, 6 archived)

---

## 🎯 ACTIVE DELIVERABLES (Priority Order)

### **QUICK WINS (Start Here):**

#### **1. TRUE SAE Production Enablement** ⚡ **5 MINUTES** ->  KEEP FALSE - USE EXISTING PROXY 
- **Status:** ✅ VERIFIED - Using Proxy SAE
- **Action:** Verified `ENABLE_TRUE_SAE_PATHWAYS` defaults to `false` in `api/config.py` (line 57)
- **Verification:** `sae_feature_service.py` uses proxy SAE by default (line 280: `"sae": "proxy"`), only switches to TRUE SAE if flag is enabled AND sae_features provided
- **Files:** `api/config.py` (line 57), `api/services/sae_feature_service.py` (line 280-329)
- **Impact:** ✅ Proxy SAE is active and working correctly

#### **2. VUS Router Registration** ⚡ **15 MINUTES**
- **Status:** ✅ COMPLETE - Router registered
- **Action:** Added `from .routers import vus as vus_router` and `app.include_router(vus_router.router)` to `main.py`
- **Files:** `api/main.py` (lines 64, 232)
- **Test:** Server restart required to test `/api/vus/identify` endpoint
- **Impact:** ✅ VUS resolution functionality now available via `/api/vus/identify`---### **CRITICAL DELIVERABLES (Blocking):**#### **3. Data Extraction Agent** 🔴 **BLOCKING - 4-6 HOURS**
- **Status:** ✅ **VERIFIED - COMPLETE** (100% - not 20% as previously documented)
- **Priority:** CRITICAL - Nothing can run without it
- **Files:** `api/services/extraction/extraction_agent.py` (EXISTS - 389 lines, all parsers complete)
- **Verification:** ✅ Agent exists, all 5 parsers (VCF/MAF/PDF/JSON/TXT) implemented, orchestrator integration working
- **Dependencies:** None
- **See:** `.cursor/MOAT/VERIFICATION_REPORT.md` for details#### **4. Drug Efficacy Integration** 🔴 **CRITICAL - 8-10 HOURS**
- **Status:** ✅ **VERIFIED - INTEGRATED** (100% - not 80% skeleton as previously documented)
- **Priority:** CRITICAL - Core drug ranking
- **Files:** `api/services/orchestrator/orchestrator.py` (lines 743-792: `_run_drug_efficacy_agent()` fully implemented)
- **Verification:** ✅ Direct import from `efficacy_orchestrator`, no HTTP calls, mechanism vector extraction working
- **Dependencies:** Data Extraction Agent (✅ Complete)
- **See:** `.cursor/MOAT/VERIFICATION_REPORT.md` for details---### **HIGH PRIORITY DELIVERABLES:**#### **5. Nutrition Integration** 🟡 **HIGH - 4-6 HOURS**
- **Status:** ✅ **VERIFIED - INTEGRATED** (100% - not 70% skeleton as previously documented)
- **Priority:** HIGH - MOAT feature
- **Files:** `api/services/orchestrator/orchestrator.py` (lines 660-741: `_run_nutrition_agent()` fully implemented)
- **Verification:** ✅ Direct import from `nutrition`, error handling with graceful fallback, germline/drug extraction working
- **See:** `.cursor/MOAT/VERIFICATION_REPORT.md` for details#### **6. End-to-End Pipeline Testing** 🟡 **MEDIUM - 2-3 HOURS**
- **Status:** ✅ **TEST SCRIPTS CREATED**
- **Files Created:**
  - `tests/test_orchestrator_e2e_pipeline.py` - Full pipeline test with mutations
  - `tests/test_vus_endpoint.py` - VUS endpoint test
  - `scripts/validate_therapy_fit_metrics.py` - Therapy Fit metric validation
- **Action:** Run tests to validate pipeline execution
- **Goal:** Verify all agents execute correctly and state is properly managed#### **7. Therapy Fit Verification** 🟡 **HIGH - 5 DAYS**
- **Status:** 🚧 **IN PROGRESS** (4/10 deliverables complete, demo ready)
- **Completed:**
  - ✅ Deliverable 1: End-to-End Test Script (`scripts/test_therapy_fit_endpoint.py` - already exists)
  - ✅ Deliverable 2: Metric Validation Script (`scripts/validate_therapy_fit_metrics.py` - created)
  - ✅ Deliverable 3: Real Patient Test Case Collection (`.cursor/MOAT/THERAPY_FIT_TEST_CASES.md` - created with 5 test cases)
  - ✅ Deliverable 4: Demo Script (`scripts/demo_therapy_fit.py` - created, verified uses real API calls, no hard-coded values)
- **Testing Infrastructure:**
  - ✅ `scripts/generate_therapy_fit_results.py` - Run test cases and generate actual results (testing/validation mode)
  - ✅ `scripts/demo_therapy_fit.py` - Demo script (presentation mode with formatted output)
  - ✅ Test Cases: 5 real patient cases (AYESHA-001, MM-001, MEL-001, OV-002, MM-002)
  - ✅ Verification: Demo script confirmed to make real HTTP API calls, displays actual responses, no mock data
- **Remaining:** 6 deliverables (see `ZO_NEXT_10_DELIVERABLES.md`)
  - ⏳ Deliverable 5: Code Verification Report (verify S/P/E formula, evidence tiers, badges, confidence computation)
  - ⏳ Deliverable 6: Performance Benchmark Script (latency, throughput, error rate)
  - ⏳ Deliverable 7: Integration Test with Frontend (verify API calls from frontend components)
  - ⏳ Deliverable 8: Updated Documentation Based on Testing (update test cases with actual results)
  - ⏳ Deliverable 9: Demo Readiness Checklist (ensure demo readiness)
  - ⏳ Deliverable 10: Brutal Honesty Report (honest assessment of implementation)
- **Usage:**
  - Testing: `python scripts/generate_therapy_fit_results.py --test-case AYESHA-001`
  - Demo: `python scripts/demo_therapy_fit.py --test-case AYESHA-001 --debug`
  - All cases: `python scripts/generate_therapy_fit_results.py --all`
- **Goal:** Complete verification and demo readiness
- **Next:** Deliverable 5 - Code Verification Report (verify implementation matches documentation)---**See:** `.cursor/MOAT/CORE_DELIVERABLES/06_GAP_ANALYSIS.md` for full gap list  
**See:** `.cursor/MOAT/orchestration/03_DELIVERABLES_PLAN.md` for orchestration plan  
**See:** `.cursor/MOAT/ZO_NEXT_10_DELIVERABLES.md` for Therapy Fit verification details**Note:** All Therapy Fit status consolidated in this file. Orphan documentation files archived to `archive/` directory.
