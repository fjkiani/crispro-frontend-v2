# 📁 File Consolidation Plan - Following MM Model

**Date:** January 28, 2025  
**Status:** 🔄 **CONSOLIDATION PLAN**  
**Model:** `.cursor/MOAT/MM/` structure

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

- [ ] TOXICITY files consolidated into `.cursor/MOAT/TOXICITY/`
- [ ] README.md created for navigation
- [ ] Source of truth identified (00_SOURCE_OF_TRUTH.md)
- [ ] Supporting docs numbered (01_, 02_, 03_)
- [ ] Old files archived
- [ ] Same pattern applied to FOOD_VALIDATOR
- [ ] Same pattern applied to ADVANCED_CARE_PLAN
- [ ] Top-level `.cursor/MOAT/` has < 20 files (only index files)

---

## 🎯 NEXT STEPS

1. **Execute Phase 1:** Consolidate TOXICITY files (15 min)
2. **Execute Phase 2:** Consolidate FOOD_VALIDATOR files (10 min)
3. **Execute Phase 3:** Consolidate ADVANCED_CARE_PLAN files (10 min)
4. **Review:** Check for other scattered topics
5. **Document:** Update main MOAT README.md with new structure

**Total Time:** ~35 minutes

---

**Status:** 🔄 **READY TO EXECUTE**  
**Model:** Follow `.cursor/MOAT/MM/` structure  
**Goal:** Organized, navigable, single source of truth per topic

