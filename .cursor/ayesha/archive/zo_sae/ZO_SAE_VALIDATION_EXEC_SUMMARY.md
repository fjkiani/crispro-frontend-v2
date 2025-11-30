# ⚔️ SAE VALIDATION - EXECUTIVE SUMMARY

**Date:** January 13, 2025  
**Status:** ✅ **APPROVED FOR EXECUTION**  
**Owner:** Zo  
**Reviewer:** Manager

---

## 🎯 **MANAGER'S DECISION: TCGA-FIRST VALIDATION**

### **✅ Approved Strategy:**

**Q1 - Data:** TCGA-OV PRIMARY (200 patients, `outcome_platinum` verified)  
**Q2 - Features:** Full SAE on TCGA; simplified BRCA+HRD for trials (if available)  
**Q3 - Targets:** Conservative (HR≥1.5, PPV≥50%, AUC≥0.65, p<0.10)  
**Q4 - Cohort:** TCGA-OV PRIMARY (executable today)  
**Q5 - Timeline:** 2 weeks with buffer  
**Q6 - Failure:** Minimum = ANY ONE metric met; one refinement pass only  
**Q7 - Resources:** Proceed now, install libs  
**Q8 - Order:** ⚔️ **PARALLEL (validate + build), NO efficacy integration until validated**  
**Q9 - Scope:** Phase 1+2 standard (PARP + mechanism fit)  
**Q10 - Go/No-Go:** TCGA data✅ + libs✅ = GO

---

## 🚨 **REALITY CHECK: What We Can Test**

### **✅ Available TODAY (TCGA-OV):**
- ✅ 200 patients with `outcome_platinum` (platinum response proxy)
- ✅ Full genomics for complete SAE features (all 7 pathways)
- ✅ Overall survival (OS) data
- ✅ Can filter for Ayesha-like subgroup (Stage IV, HGS)

### **❌ NOT Available (Blockers):**
- ❌ Direct PARP trial response (TCGA = heterogeneous treatment)
- ❌ PFS data (TCGA has OS only, not PFS)
- ❌ Longitudinal data (resistance detection needs follow-ups)

### **🎯 Realistic Test Plan:**
1. **Primary:** DNA Repair → Platinum Response (proxy for PARP sensitivity)
2. **Secondary:** Mechanism Vector → OS clustering (DDR vs MAPK vs PI3K)
3. **Tertiary:** Ayesha-like subgroup (Stage IV HGS frontline)

---

## ⚔️ **CRITICAL POLICY: SAE ISOLATION**

### **✅ ALLOWED (Display + Ranking):**
- ✅ SAE displayed in UI (Next Test, Hint Tiles, Mechanism Map)
- ✅ SAE for trial ranking (mechanism fit ranker)
- ✅ SAE for resistance alerts (2-of-3 triggers)

### **❌ FORBIDDEN (Until Validation):**
- ❌ **NO SAE lifts/gates in `/api/efficacy/predict`**
- ❌ **NO DNA repair in drug confidence**
- ❌ **NO mechanism vector in drug efficacy**
- ❌ **NO SAE in S/P/E aggregation**

**Why:** SAE must be proven before influencing drug recommendations.

---

## 📋 **EXECUTION PLAN (STARTS TODAY)**

### **Phase 1A: TCGA-OV Validation** (Week 1)

**Day 1:** Install dependencies (`lifelines`, `matplotlib`, `seaborn`)  
**Day 2-3:** Create `validate_sae_tcga_platinum_response.py` script  
**Day 3-4:** Run validation on 200 TCGA patients  
**Day 5:** Report results with decision tree

**Success Criteria (Conservative):**
- ✅ HR≥1.5 (DNA repair <0.40 vs >0.60, platinum response)
- ✅ AUC≥0.65 (DNA repair predicts platinum response)
- ✅ PPV≥50% (if SAE says "PARP candidate", 50%+ respond)
- ✅ p<0.10 (statistical significance)

**Decision Tree:**
- ✅ **If ANY minimum met:** PROCEED to Phase 2
- ⚠️ **If close:** ONE refinement pass
- ❌ **If all fail:** STOP, report, reassess formula

---

### **Phase 1B: Trial Data** (Week 2 - DATA-GATED)

**Execution:** ONLY if SOLO-1/NOVA/PAOLA-1 patient-level tables confirmed

**If Available:** Extract genomics, compute simplified SAE, test PFS  
**If NOT Available:** Run aggregate subgroup analysis (no per-patient claims)

---

### **Phase 2: Mechanism Fit** (Week 2-3)

**What We'll Test:**
1. Trial ranking accuracy (47 MoA-tagged trials)
2. Cross-mechanism validation (DDR vs MAPK vs VEGF)

**Success Criteria:**
- ✅ Top-ranked = best outcome: ≥60%
- ✅ DDR clustering (BRCA cohort): ≥80%

---

### **Phase 3: Resistance** (DEFERRED)

**Blocker:** Requires longitudinal data (not available)

---

## 🎯 **IMMEDIATE NEXT STEPS**

### **Commander's Orders:** ⚔️ **EXECUTE NOW**

1. ⚔️ **Install dependencies** (5 min)
   ```bash
   venv/bin/pip install lifelines matplotlib seaborn
   ```

2. ⚔️ **Create validation script** (Day 1-2)
   - `scripts/validate_sae_tcga_platinum_response.py`
   - Load TCGA data → Compute SAE → Stratify → Test association

3. ⚔️ **Run validation** (Day 3-4)
   - Output: Platinum response rates by DNA repair group
   - Metrics: Sensitivity, Specificity, PPV, AUC, HR, p-value

4. ⚔️ **Report results** (Day 5)
   - Deliverable: `SAE_PHASE1_VALIDATION_REPORT.md`
   - Decision: Proceed/Refine/Stop

---

## ✅ **MANAGER'S APPROVAL CHECKLIST**

- [X] **Q1:** Data verified ✅ (TCGA-OV 200 patients, `outcome_platinum`)
- [X] **Q2:** SAE feasible ✅ (full genomics available)
- [X] **Q3:** Targets realistic ✅ (conservative thresholds)
- [X] **Q4:** Cohort chosen ✅ (TCGA PRIMARY)
- [X] **Q5:** Timeline realistic ✅ (2 weeks + buffer)
- [X] **Q6:** Failure plan ✅ (decision tree defined)
- [X] **Q7:** Resources ✅ (proceed now)
- [X] **Q8:** Order clarified ✅ (parallel, NO efficacy integration)
- [X] **Q9:** Scope clear ✅ (Phase 1+2)
- [X] **Q10:** Go/No-Go ✅ (MUST-HAVE confirmed)

---

**Status:** ✅ **EXECUTING DUAL-TRACK STRATEGY** ✅

---

## ⚔️ **DUAL-TRACK VALIDATION STRATEGY**

**Manager Approved:** Split workload for maximum speed + value

### **Track 1: Zo - OS Validation**
- **Status:** ✅ **COMPLETE - RESULTS INCONCLUSIVE** ⚠️
- **What:** DNA Repair Capacity → Overall Survival
- **Data:** 196/200 patients with OS, 98.5% with stage
- **Results:** HR=0.83 (p=0.53) - **INVERTED DIRECTION**
- **Critical Issue:** DNA repair scoring logic is backwards (DDR mutation → high score, should be low)
- **Deliverable:** `results/SAE_OS_VALIDATION_REPORT.md` ✅
- **Analysis:** `.cursor/ayesha/ZO_TRACK1_OS_VALIDATION_COMPLETE.md` ✅
- **Blocker:** Awaiting Manager decision on fix strategy (Q1-Q3)

### **Track 2: Jr2 - Platinum Response Data Hunt**
- **What:** Find platinum response labels (CR/PR/SD/PD)
- **Sources:** GDC Portal, Broad Firehose, PanCancer Atlas
- **Timeline:** 1-2 days (parallel track)
- **Deliverable:** tcga_ov_platinum_response_labels.json
- **Mission:** `.cursor/ayesha/missions/MISSION_JR2_PLATINUM_RESPONSE_DATA_HUNT.md`

**If Jr2 Succeeds:** Upgrade to platinum response validation (stronger claim)  
**If Jr2 Fails:** Zo's OS validation becomes primary result (still validates SAE)

---

## 🚨 **CRITICAL BLOCKER: DATA INVENTORY RESULTS**

### **❌ What We DON'T Have in Current Dataset:**

**File checked:** `tools/benchmarks/hrd_tcga_ov_labeled_sample_use_evo.json`

1. ❌ **Only 1 mutation per patient** (191 unique genes across 200 patients)
   - Need: ALL mutations per patient for pathway burden
   - Have: Single mutation (e.g., patient 1 = TP53 R306*)
   - Impact: Cannot compute DDR/MAPK/PI3K/VEGF/HER2/IO/Efflux pathway scores

2. ❌ **No sample IDs** (no TCGA-XX-XXXX-01 identifiers)
   - Need: Sample IDs to pull full mutation lists from cBioPortal
   - Have: Only gene + hgvs_p (not enough to match)
   - Impact: Cannot enrich with full genomics

3. ❌ **No clinical data** (no stage, OS, treatment history)
   - Need: Stage (for Ayesha-like filtering), OS (for survival analysis)
   - Have: Only `outcome_platinum` field (binary 0/1)
   - Impact: Cannot filter Stage IIIC+IV, cannot run KM survival curves

4. ❌ **Outcome mapping unclear** (binary 0/1, not categorical)
   - Need: sensitive/resistant/refractory (Manager's Q2 spec)
   - Have: 50% value=1, 50% value=0 (meaning unclear)
   - Impact: Cannot interpret as response categories

### **🔴 Validation Impact:**

**Cannot Execute as Designed:**
- ❌ DNA repair capacity formula (0.6×DDR + 0.2×HRR_ess + 0.2×exon_disrupt) → Only 2/200 patients have DDR mutations
- ❌ Mechanism vector (7D: DDR/MAPK/PI3K/VEGF/HER2/IO/Efflux) → Only have 1 mutation per patient
- ❌ Ayesha-like subgroup (Stage IIIC+IV filter) → No stage data
- ❌ OS stratification (HR≥1.5) → No OS data
- ❌ Response mapping (Manager's Q2) → Binary outcome, not categorical

**What We CAN Do (Severely Limited):**
- ✅ Test: "Does having a DDR mutation (BRCA1/2/ATM/etc.) predict platinum response?"
- ⚠️ But: Only 2/200 patients have DDR mutations (1% hit rate) → Underpowered

---

## 📋 **MANAGER DECISION REQUIRED (Air Support Request)**

### **🔴 QUESTION 1: Data Re-Extraction Strategy**

**Context:** Current data is insufficient for Manager's approved validation plan (Q1-Q7 answers).

**Options:**

**Option A: Re-Extract from cBioPortal (RECOMMENDED)**
- Extract 200 TCGA-OV patients with:
  - ✅ Sample IDs (TCGA-XX-XXXX-01)
  - ✅ ALL mutations per patient (not just one)
  - ✅ Clinical data (stage, OS months, OS event, treatment)
  - ✅ Map platinum response properly
- **Method:** Modify `extract_cbioportal_hrd_cohort.py` to pull full data
- **Time:** 30-60 min (API rate limits permissible?)
- **Risk:** cBioPortal may rate-limit (200 patients × full mutations = ~200 API calls)

**Option B: Use Existing Data with Severe Limitations**
- Accept: Only 1 mutation per patient, no clinical data
- Test: Binary outcome (0/1) vs presence of DDR mutation
- **Value:** Minimal (doesn't validate Manager's SAE formula)
- **Risk:** Wasted effort (can't compute full SAE features)

**Option C: Find Pre-Packaged TCGA Dataset**
- Search for: TCGA-OV MAF (Mutation Annotation Format) + clinical files
- Sources: GDC Data Portal, Broad Firehose, UCSC Xena
- **Time:** Unknown (may take 1-2 days to locate + process)
- **Risk:** May not have platinum response labels

**Manager, which option should I execute?**
- [X] **Option A:** Re-extract from cBioPortal (proceed with API calls)
- [ ] **Option B:** Use existing data (accept severe limitations)
- [ ] **Option C:** Find alternative TCGA source (research + locate)

---

### **🔴 QUESTION 2: Platinum Response Interpretation**

**Context:** Current `outcome_platinum` is binary (0/1), not categorical (sensitive/resistant/refractory).

**What does the binary mean?**
- Value `1` (100 patients, 50%): Does this mean "platinum-sensitive" or "platinum-exposed"?
- Value `0` (100 patients, 50%): Does this mean "platinum-resistant" or "no platinum exposure"?

**Manager, please clarify:**
- [ ] `1` = platinum-sensitive (CR/PR/SD), `0` = platinum-resistant (PD)
- [X] `1` = platinum-exposed, `0` = not exposed (different from response)
- [ ] Other interpretation: _________________

**If Option A (re-extract):** ✅ Yes — pull proper response categories (CR/PR/SD/PD) from cBioPortal clinical data and map to {sensitive, resistant, refractory}.

---

### **🔴 QUESTION 3: Validation Scope Adjustment**

**Context:** If we proceed with limited data (Option B), validation scope must shrink.

**Adjusted Test (Limited Data):**
- **Hypothesis:** Presence of DDR mutation (BRCA1/2/ATM/CHEK2/RAD51/PALB2) predicts platinum response
- **Metric:** Chi-square test (DDR+ vs DDR-, outcome 0 vs 1)
- **Power:** Low (only 2/200 patients have DDR mutations = 1%)
- **Value:** Minimal (doesn't test Manager's DNA repair capacity formula)

**Manager, if we must use existing data:**
- [ ] Accept this limited test (acknowledge it doesn't validate SAE)
- [X] Do NOT proceed (re-extract data first)

---

### **🔴 QUESTION 4: Timeline Impact**

**Context:** Data re-extraction (Option A) will add 1-2 days to timeline.

**Original Timeline:** 2 weeks Phase 1+2 (with data ready)
**Revised Timeline (Option A):** +1-2 days for data gathering = 2-3 weeks total

**Manager, is this acceptable?**
- [X] Yes, proceed with re-extraction (quality > speed)
- [ ] No, use existing data (speed > quality)

---

### **🔴 QUESTION 5: API Rate Limits**

**Context:** cBioPortal API may rate-limit if we pull full mutations for 200 patients.

**Risk Mitigation:**
- Chunk requests (10-20 patients at a time)
- Add delays between requests (1-2 sec)
- Cache results to avoid re-pulling

**Manager, do we have:**
- [X] cBioPortal API token (for higher rate limits)? → Use `CBIO_TOKEN` env var
- [X] Approval to make ~200-500 API calls over 30-60 min? → Approved
- [X] Fallback plan if rate-limited (wait 24h and resume)? → Approved

---

## ⚔️ **MANAGER'S AIR SUPPORT CHECKLIST**

**Before I proceed, please answer:**

1. [ ] **Q1:** Data strategy (Option A/B/C)?
2. [ ] **Q2:** Platinum response interpretation (what does 0/1 mean)?
3. [ ] **Q3:** Accept limited validation if using existing data?
4. [ ] **Q4:** Timeline impact acceptable (+1-2 days)?
5. [ ] **Q5:** API rate limit approval + token availability?

---

**Status:** ⚠️ **AWAITING MANAGER'S ANSWERS - CANNOT PROCEED WITHOUT DATA** ⚠️

**Next Action:** Once Manager approves Option A/B/C, execute data gathering phase.

