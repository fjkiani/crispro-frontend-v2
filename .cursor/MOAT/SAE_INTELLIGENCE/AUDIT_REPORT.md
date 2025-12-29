# SAE_INTELLIGENCE Directory - Systematic Audit Report

**Date:** January 28, 2025  
**Auditor:** Zo  
**Scope:** Complete systematic review of `.cursor/MOAT/SAE_INTELLIGENCE/` directory  
**Purpose:** Identify inconsistencies, gaps, contradictions, and alignment issues

---

## 🔍 EXECUTIVE SUMMARY

### **Critical Findings:**

1. **✅ DELIVERABLES EXIST — BUT THEY SUPPORT DIFFERENT CLAIMS**
   - **External, publishable platinum resistance signal (expression)**: **MFAP4 AUROC = 0.763** on **GSE63885** (n=101; 34 resistant / 67 sensitive). See `Publication-1/SAE_RESISTANCE/VALIDATION_SUMMARY_FINAL.md` and `Publication-1/SAE_RESISTANCE/WHY_WE_ARE_OFF.md`.
   - **Internal Tier-3 “TRUE SAE” signal (exploratory, contract-specific)**: mean CV-AUROC ≈ **0.783 ± 0.100** on **149** patients (pos = resistant+refractory = 24) using **29 selected SAE features**. This is *not external validation*; treat as internal until replicated with a clean contract and confound controls.
   - **DDR_bin on TCGA-style platinum labels fails**: AUROC(resistant) ≈ **0.517** (≈ random) under TCGA-OV platinum-response contract.
   - **DDR_bin appears prognostic (OS)**: HR ≈ **0.62** (p≈0.013) and Spearman(DDR_bin, OS) ≈ **0.252** (p≈0.0013) in TCGA-OV v2 analyses — but see confounding notes below.

2. **🚨 VALIDATION CONTRACT DRIFT + CONFOUNDING IS THE CENTRAL RISK**
   - Multiple cohorts/label contracts exist and disagree. Any “resistance biomarker” claim must always cite: cohort file + label mapping + endpoint (predictive vs prognostic).
   - DDR_bin is strongly confounded by **extraction coverage / variant count** in the v2 linked cohort (see `Publication-1/SAE_RESISTANCE/WHY_WE_ARE_OFF.md`).

3. **✅ FEATURE SETS ARE REAL; “BLOCKED” VS “MAPPED” NEEDS PRECISION**
   - TRUE SAE extraction exists (32K sparse features), but **full 32K→pathway mapping is not complete**.
   - A **partial mapping** does exist for the **9 “diamond” features** into a DDR_bin aggregate; this enables the DDR_bin and Tier-3 baselines. Any doc claiming “TRUE SAE is blocked solely by mapping” should be interpreted as “full mapping incomplete,” not “no mapping exists.”

4. **⚠️ DOC/DATE INCONSISTENCIES ARE SECONDARY BUT REAL**
   - There are mixed timestamps and mixed narratives across `.cursor/MOAT/SAE_INTELLIGENCE/` vs `Publication-1/SAE_RESISTANCE/`. The canonical publication truth should be taken from the publication folder’s postmortems + validation summaries.

5. **⚠️ DEMO-NARRATIVE DRIFT (MECHANISM FIT)**
   - Some demo readiness docs say “mechanism fit not wired,” while `complete_care_universal` does import SAE + mechanism-fit machinery. This is likely “which endpoint path” ambiguity, not a pure contradiction.

---

## 📊 SYSTEMATIC AUDIT FINDINGS

### **1. DATE INCONSISTENCIES** ⚠️ **MEDIUM PRIORITY**

**Issue:** Mixed dates across files suggest different update cadences

| File | Date | Status | Issue |
|------|------|--------|-------|
| `README.md` | December 23, 2025 | ✅ Most recent | Latest TRUE SAE validation |
| `03_GENERALS_BATTLE_MAP.mdc` | December 23, 2025 | ✅ Consistent | TRUE SAE breakthrough |
| `06_STRATEGIC_VISION.md` | December 23, 2025 | ✅ Consistent | Strategic vision |
| `07_TRUE_SAE_DIAMONDS_EXCAVATION.md` | 2025-12-24 | ✅ Consistent | PLUMBER workpack |
| `00_MISSION.mdc` | January 28, 2025 | ⚠️ Older | Needs update |
| `01_SAE_SYSTEM_DEBRIEF.mdc` | January 13, 2025 | ⚠️ Older | Pre-validation |
| `02_SAE_ARCHITECTURE.md` | January 28, 2025 | ✅ Current | Architecture |
| `04_INTELLIGENCE_FLOW.md` | January 28, 2025 | ✅ Current | Intelligence flow |
| `05_SAE_CAPABILITIES.md` | January 28, 2025 | ✅ Current | Capabilities |
| `07_STRATEGIC_DELIVERABLES_PLAN.md` | January 28, 2025 | ✅ Current | Strategic plan |

**Analysis:**
- Mixed dates appear to reflect different update cadences and/or copy/restore artifacts.
- The **canonical “what we can claim”** for the resistance publication lives in `Publication-1/SAE_RESISTANCE/` (Dec 2024), not in generalized SAE_intelligence docs.
- Treat “AUROC 0.783” as **Tier-3 internal contract** performance, not “validated ovarian platinum resistance” in the general sense.

**Recommendation:**
- Add a short note in SAE_INTELLIGENCE core docs clarifying: **internal Tier-3 exploratory signal vs external validation** (MFAP4).
- Ensure any mention of “0.783” explicitly states: **Tier-3 internal contract**, 29 features, 5-fold CV, and caveats.

---

### **2. FEATURE COUNT CONFUSION** ⚠️ **HIGH PRIORITY**

**Issue:** Inconsistent reporting of "29 features" vs "9 diamonds"

**Current State:**
- `README.md`: "29 features" for TRUE SAE AUROC 0.783
- `03_GENERALS_BATTLE_MAP.mdc`: "9 diamond features" map to DDR_bin
- `07_TRUE_SAE_DIAMONDS_EXCAVATION.md`: "29 candidate features" for multi-feature baseline
- `true_sae_diamonds_baseline.v1.json`: Confirms 29 features (9 diamonds + 20 additional)

**Reality Check:**
- ✅ **29 features total** = 9 diamond features + 20 additional top features
- ✅ **9 diamonds** = Large-effect features (d > 0.5, higher in resistant)
- ✅ **Multi-feature baseline** = Logistic regression on all 29 features
- ✅ **Single-feature AUROC** = Individual diamond features (~0.62-0.65)

**Gap:**
- Documentation doesn't clearly explain: 29 total features = 9 diamonds + 20 others
- Some files say "29 features" without clarifying composition
- Some files say "9 diamonds" without mentioning the 20 additional features

**Recommendation:**
- Add clarification to `README.md`: "29 features (9 diamonds + 20 additional top features)"
- Update `03_GENERALS_BATTLE_MAP.mdc` to clarify: "9 diamond features (of 29 total)"
- Ensure all files consistently explain the 29 vs 9 relationship

---

### **3. REFRACTORY/RESISTANT CLASSIFICATION** ⚠️ **HIGH PRIORITY**

**Issue:** Critical label definition only documented in workpack, missing from other files

**Current State:**
- ✅ `07_TRUE_SAE_DIAMONDS_EXCAVATION.md`: Clearly documents `refractory + resistant = resistant` (pos class = 24)
- ✅ `true_sae_diamonds_baseline.v1.json`: Confirms positive_class = ["resistant", "refractory"]
- ✅ `validate_true_sae_diamonds.py`: Code correctly implements `y.append(1 if out in ("resistant", "refractory") else 0)`
- ❌ **Other files**: No mention of refractory classification

**Impact:**
- If someone reads only `README.md` or `03_GENERALS_BATTLE_MAP.mdc`, they won't know about refractory classification
- Could lead to confusion when interpreting validation results
- Missing from strategic vision and mission documents

**Recommendation:**
- Add to `README.md` Quick Reference: "Label definition: refractory + resistant = resistant (pos class = 24)"
- Add to `03_GENERALS_BATTLE_MAP.mdc`: Note about label definition in validation section
- Add to `00_MISSION.mdc`: Document label definition in validated metrics section

---

### **4. PLUMBER WORKPACK STATUS** ✅ **COMPLETE BUT NOT DOCUMENTED**

**Issue:** Deliverables appear to exist but status not reflected in documentation

**Current State:**
- ✅ `sae_feature_mapping.true_sae_diamonds.v1.json` EXISTS
  - Location: `oncology-coPilot/oncology-backend-minimal/api/resources/`
  - Contains: All 9 diamond features mapped to DDR_bin with high confidence
  - Evidence: Top variants, top patients, pathway distribution
  
- ✅ `true_sae_diamonds_baseline.v1.json` EXISTS
  - Location: `oncology-coPilot/oncology-backend-minimal/data/validation/sae_cohort/checkpoints/`
  - Contains: 29 features, 5-fold CV results, mean AUROC ≈ 0.78
  - Analysis date: 2025-01-28

- ✅ Validation script EXISTS
  - `validate_true_sae_diamonds.py` checks both deliverables
  - Verifies label contract, feature mapping, baseline reproducibility

**Gap:**
- `07_TRUE_SAE_DIAMONDS_EXCAVATION.md` says "READY" but doesn't indicate completion
- Strategic deliverables plan doesn't mention this workpack
- README doesn't indicate workpack status

**Recommendation:**
- Update `07_TRUE_SAE_DIAMONDS_EXCAVATION.md` status to "✅ COMPLETE" or "⚠️ PARTIAL" (need to verify)
- Add to `07_STRATEGIC_DELIVERABLES_PLAN.md`: "Deliverable 0: TRUE SAE Diamonds Mapping (COMPLETE)"
- Update `README.md` to indicate workpack deliverables exist

**Update (re-audit):**
- ✅ `07_STRATEGIC_DELIVERABLES_PLAN.md` now includes **Deliverable 0: TRUE SAE Diamonds Mapping** as complete.
- ✅ `README.md` now links to the workpack and calls out the label definition.

---

### **5. MECHANISM VECTOR DIMENSION CONSISTENCY** ✅ **CONSISTENT**

**Finding:** All files consistently use 7D mechanism vector

**Current State:**
- ✅ All files: 7D = [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
- ✅ Consistent across: `00_MISSION.mdc`, `01_SAE_SYSTEM_DEBRIEF.mdc`, `02_SAE_ARCHITECTURE.md`, `05_SAE_CAPABILITIES.md`
- ✅ Workpack question: "7D mechanism vector dims (DDR/MAPK/PI3K/VEGF/HER2/IO/Efflux) or start narrower (DDR/MAPK/PI3K/OTHER)?"

**Status:** ✅ **NO ISSUE** - Consistent documentation

---

### **6. STEERABILITY DOCUMENTATION** ⚠️ **PARTIAL**

**Issue:** Steerability V1 mentioned but workpack connection not clear

**Current State:**
- ✅ `03_GENERALS_BATTLE_MAP.mdc`: Steerability V1 roadmap (DDR_bin interventions)
- ✅ `06_STRATEGIC_VISION.md`: Steerability V1 roadmap (clamp_ddr_bin function)
- ✅ `07_TRUE_SAE_DIAMONDS_EXCAVATION.md`: Mentions "justifies Steerability V1"
- ⚠️ **Gap:** Connection between workpack deliverables and Steerability V1 not explicit

**Recommendation:**
- Add to Steerability sections: "Workpack deliverables (feature mapping + baseline) required for Steerability V1"
- Update workpack status to indicate it's a prerequisite for Steerability V1

---

### **7. CROSS-REFERENCE ACCURACY** ✅ **MOSTLY GOOD**

**Finding:** Most cross-references are accurate

**Verified:**
- ✅ `README.md` → All numbered files (00-07)
- ✅ `00_MISSION.mdc` → `01_SAE_SYSTEM_DEBRIEF.mdc`, `03_GENERALS_BATTLE_MAP.mdc`
- ✅ `03_GENERALS_BATTLE_MAP.mdc` → `07_TRUE_SAE_DIAMONDS_EXCAVATION.md`
- ✅ `06_STRATEGIC_VISION.md` → `07_TRUE_SAE_DIAMONDS_EXCAVATION.md`

**Missing:**
- ⚠️ `07_STRATEGIC_DELIVERABLES_PLAN.md` doesn't reference workpack
- ⚠️ `01_SAE_SYSTEM_DEBRIEF.mdc` doesn't reference TRUE SAE validation (dated before validation)

**Recommendation:**
- Add workpack reference to strategic deliverables plan
- Add TRUE SAE validation section to `01_SAE_SYSTEM_DEBRIEF.mdc`

---

### **8. TECHNICAL CLAIM VALIDATION** ✅ **VERIFIED**

**Finding:** Artifacts exist and numbers are reproducible, but **interpretation depends on contract and confounds**.

**Verified Claims:**
- ✅ **Tier-3 internal TRUE SAE baseline**: mean CV-AUROC ≈ 0.783 ± 0.100 (29 features; pos = resistant+refractory = 24/149) → See Tier-3 artifacts referenced in `Publication-1/SAE_RESISTANCE/ERRATA.md`.
- ✅ **DDR_bin under TCGA-style platinum labels**: AUROC(resistant) ≈ 0.517 (≈ random) → See `Publication-1/SAE_RESISTANCE/ERRATA.md` and linked report.json.
- ✅ **DDR_bin is prognostic (OS) in v2 analyses**: HR ≈ 0.62 (p≈0.013), Spearman rho ≈ 0.252 (p≈0.0013) → See `Publication-1/SAE_RESISTANCE/VALIDATION_RESULTS_FINAL.md`.
- ✅ **MFAP4 external platinum resistance validation**: AUROC = 0.763 on GSE63885 → See `Publication-1/SAE_RESISTANCE/VALIDATION_SUMMARY_FINAL.md` and `WHY_WE_ARE_OFF.md`.
- ✅ **Label definition in Tier-3**: refractory + resistant = positive class → documented in workpack and referenced in publication errata.

**Update (re-audit):**
- ✅ The *existence* of the key artifacts is verified (mapping JSON, baseline JSON, validator outputs).
- ⚠️ However, the workspace contains **multiple cohort/label contracts** that yield **different** answers for “platinum resistance.” This is not a bug in the files; it is a **contract drift risk** if claims are made without citing the contract.
 - 🚨 Additionally, the publication postmortem identifies **coverage/variant-count confounding** for DDR_bin in linked cohorts; this must be stated whenever OS/prognostic associations are discussed.

---

## 🧾 CLAIMS → RECEIPTS TRUTH TABLE (CANONICAL LINKS)

### **A) Tier-3 internal contract (149 patients; pos = refractory+resistant = 24)**
- **Tier-3 cohort file**: `oncology-coPilot/oncology-backend-minimal/data/validation/sae_cohort/checkpoints/Tier3_validation_cohort.json`
- **TRUE SAE multi-feature baseline**: `oncology-coPilot/oncology-backend-minimal/data/validation/sae_cohort/checkpoints/true_sae_diamonds_baseline.v1.json`
  - Claim supported: **mean CV-AUROC ≈ 0.783 ± 0.100** (29 features = 9 diamonds + 20 additional).
- **DDR_bin response benchmark (Tier-3)**: `oncology-coPilot/oncology-backend-minimal/scripts/validation/out/ddr_bin_tcga_ov/report.json`
  - Claim supported: **DDR_bin AUROC(resistant) = 0.698** on Tier-3.

### **B) Platinum-merged contract (different cohort; TCGA-style platinum labels)**
- **Cohort file**: `oncology-coPilot/oncology-backend-minimal/data/validation/sae_cohort/checkpoints/OV_PLATINUM_TRUE_SAE_cohort.v2.json`
- **DDR_bin report**: `oncology-coPilot/oncology-backend-minimal/scripts/validation/out/ddr_bin_ov_platinum_TRUE_SAE_v2/report.json`
  - Claim supported: **DDR_bin AUROC(resistant) ≈ 0.517 (≈ random)** under this contract.
  - Also reports: **Spearman(DDR_bin, OS) ≈ 0.252 (p≈0.0013)** (interpret carefully; see confounding notes in publication package).

### **C) “p = 0.0020” source (needs explicit citation when used)**
- The “DDR_bin distinguishes resistant vs sensitive (p = 0.0020)” claim is produced by the **publication plotting scripts** running Mann–Whitney U on a specific cohort/encoding:
  - Script: `oncology-coPilot/oncology-backend-minimal/scripts/publication/generate_ddr_bin_distribution.py`
  - Output: publication figure annotations (and the publication folder references).
  - **Rule:** if `p=0.0020` is stated, the cohort + script output must be cited to avoid ambiguity with other contracts.

### **D) External validation (non-SAE, expression)**
- **MFAP4 external validation**: `Publication-1/SAE_RESISTANCE/VALIDATION_SUMMARY_FINAL.md` + figures in `Publication-1/SAE_RESISTANCE/figures/`
  - Claim supported: **MFAP4 AUROC = 0.763** on **GSE63885** (expression biomarker).

---

## ⚠️ DEMO NARRATIVE ALIGNMENT NOTE (MECHANISM FIT “WIRED” VS “NOT WIRED”)

There are two different “trial” surfaces discussed across docs:
- **Universal orchestration path**: `oncology-coPilot/oncology-backend-minimal/api/routers/complete_care_universal.py`
  - Imports `compute_sae_features` and `rank_trials_by_mechanism` (i.e., the mechanism-fit machinery is present in this orchestrator).
- **Autonomous trial agent path** (separate): demo-readiness docs may refer to this router not being wired for mechanism fit yet.

**Actionable guidance:** demo/readme claims should specify which endpoint/path is being referenced to avoid “it’s wired / it’s not wired” contradictions.

---

### **9. MISSING INFORMATION** ⚠️ **MEDIUM PRIORITY**

**Gaps Identified:**

1. **Open Questions Not Answered:**
   - Q1: Should "refractory" be treated as resistant globally or only in Tier-3?
   - Q2: Mapping bins alignment (7D vs narrower)?
   - **Status:** Questions remain open in workpack

2. **Steerability V1 Implementation Status:**
   - Roadmap exists but no implementation status
   - `clamp_ddr_bin()` function not mentioned as implemented
   - **Status:** Unclear if Steerability V1 is ready or pending

3. **Publication Status:**
   - Manuscript draft mentioned but submission status unclear
   - bioRxiv submission status not documented
   - **Status:** Needs update

4. **MM Expansion Status:**
   - Multiple mentions of "MM TRUE SAE extraction pending"
   - No clear status on MMRF data acquisition
   - **Status:** Needs tracking

**Recommendation:**
- Create status tracking for open questions
- Add implementation status for Steerability V1
- Update publication status in README
- Add MM expansion tracking to strategic plan

---

### **10. FILE ORGANIZATION** ✅ **GOOD**

**Finding:** File organization is consistent and logical

**Structure:**
- ✅ Numbered files (00-07) show progression
- ✅ README.md as navigation hub
- ✅ Archive directory for old files
- ✅ Cross-references properly maintained

**Status:** ✅ **NO ISSUES** - Well organized

---

## 🎯 PRIORITY ACTION ITEMS

### **🔴 HIGH PRIORITY (Fix Immediately):**

1. **Clarify Feature Count (29 vs 9)**
   - Update `README.md` to explain: "29 features (9 diamonds + 20 additional)"
   - Update `03_GENERALS_BATTLE_MAP.mdc` to clarify relationship
   - **Impact:** Prevents confusion about TRUE SAE composition

2. **Document Refractory Classification**
   - Add to `README.md` Quick Reference
   - Add to `03_GENERALS_BATTLE_MAP.mdc` validation section
   - Add to `00_MISSION.mdc` validated metrics
   - **Impact:** Critical for understanding validation results

3. **Update Workpack Status**
   - Verify deliverables are complete (mapping + baseline exist)
   - Update `07_TRUE_SAE_DIAMONDS_EXCAVATION.md` status
   - Add to strategic deliverables plan
   - **Impact:** Prevents duplicate work, clarifies what's done

4. **Lock the publication claim set (predictive vs prognostic)**
   - Explicitly distinguish:
     - **Predictive platinum response**: DDR_bin fails under TCGA-style platinum labels; MFAP4 succeeds on GSE63885.
     - **Prognostic survival**: DDR_bin shows OS association but is vulnerable to extraction confounds.
   - **Impact:** Prevents manuscript “narrative oscillation” between incompatible claims.

5. **Add confound controls as a hard requirement for any SAE→resistance claim**
   - At minimum: stratify/adjust for `ddr_bin_num_variants` and `ddr_bin_coverage` where present; report results with and without adjustment.
   - **Impact:** Converts “interesting correlation” into something defensible.

### **🟡 MEDIUM PRIORITY (Fix Soon):**

4. **Resolve Date Inconsistencies**
   - Update `00_MISSION.mdc` to reflect TRUE SAE validation (December 23)
   - Update `01_SAE_SYSTEM_DEBRIEF.mdc` to include validation results
   - Add note explaining date discrepancy
   - **Impact:** Prevents confusion about timeline

5. **Answer Open Questions**
   - Q1: Refractory classification (global vs Tier-3 only)
   - Q2: Mapping bins alignment (7D vs narrower)
   - Document answers in workpack
   - **Impact:** Unblocks Steerability V1 implementation

6. **Integrate Workpack into Strategic Plan**
   - Add "Deliverable 0: TRUE SAE Diamonds Mapping" to strategic plan
   - Mark as complete or in-progress based on verification
   - **Impact:** Aligns strategic planning with actual work

### **🟢 LOW PRIORITY (Nice to Have):**

7. **Add Implementation Status Tracking**
   - Steerability V1 implementation status
   - Publication submission status
   - MM expansion tracking
   - **Impact:** Better visibility into progress

---

## ✅ VERIFIED CORRECT

**These are confirmed accurate and need no changes:**

1. ✅ **MFAP4 predicts platinum resistance externally**: AUROC 0.763 on GSE63885 (n=101; 34 resistant / 67 sensitive).
2. ✅ **DDR_bin does not discriminate baseline platinum response under TCGA-style labels**: baseline sensitive vs resistant DDR_bin difference is null (p≈0.80 in publication summary).
3. ✅ **DDR_bin is prognostic for OS in v2 analyses**: HR≈0.62 (p≈0.013), rho≈0.252 (p≈0.0013) — interpretation must note confounding risk.
4. ✅ **Tier-3 internal contract** exists and is consistently defined: n=149; pos = resistant+refractory = 24.
5. ✅ **29 vs 9 relationship** is real: 29 selected features include 9 “diamonds” (DDR_bin aggregate is built from the 9).
6. ✅ **7D mechanism vector** is consistently documented ([DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]).
7. ✅ File organization is strong; cross-references are mostly accurate (publication folder is the canonical claim source for resistance).

---

## 📋 AUDIT SUMMARY

| Category | Status | Issues Found | Priority |
|----------|--------|--------------|----------|
| **Date Consistency** | ⚠️ | 2 files need update | MEDIUM |
| **Feature Count Clarity** | ⚠️ | Needs clarification | HIGH |
| **Refractory Classification** | ⚠️ | Missing from key files | HIGH |
| **Workpack Status** | ⚠️ | Deliverables exist but not documented | HIGH |
| **Technical Claims** | ✅ | All verified | N/A |
| **Cross-References** | ✅ | Mostly accurate | LOW |
| **File Organization** | ✅ | Well-structured | N/A |
| **Open Questions** | ⚠️ | 2 unanswered | MEDIUM |

**Overall Assessment:** Documentation is **GOOD** but needs **3 HIGH priority fixes** and **3 MEDIUM priority updates** to be fully accurate and complete.

---

## 🔗 NEXT STEPS

1. **Immediate:** Fix HIGH priority items (feature count, refractory classification, workpack status)
2. **This Week:** Resolve date inconsistencies and answer open questions
3. **This Month:** Integrate workpack into strategic plan, add status tracking

---

*Audit Completed: January 28, 2025*  
*Next Audit: After HIGH priority fixes are complete*


