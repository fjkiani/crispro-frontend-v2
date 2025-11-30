# 🚨 CRITICAL AUDIT: SAE PHASES 1-3 - GAPS & CONCERNS ⚔️

**Date:** January 13, 2025  
**Auditor:** Zo (Self-Audit)  
**Purpose:** Identify what could break, what's incomplete, what might have slipped  
**Status:** 🔍 **HONEST ASSESSMENT** - Not celebrating achievements, finding gaps

---

## 🎯 EXECUTIVE SUMMARY

**Overall Assessment:** ⚠️ **Backend 70% Complete** | 🚨 **Integration 40% Done** | 🚨 **Validation 0% Done**

**CRITICAL ARCHITECTURAL MISMATCH:**
- **Manager's Vision:** SAE integrated into S/P/E drug efficacy pipeline (Evo2 → Insights → Pathway → SAE → Confidence)
- **What We Built:** SAE isolated in Ayesha orchestrator, not integrated into WIWFM drug ranking
- **Impact:** SAE is "display only", not affecting drug recommendations

**Critical Findings (TOP 10):**
1. 🚨 **FORMULA MISMATCH**: Manager said `0.6×DDR`, we implemented `0.5×DDR` (DNA repair capacity)
2. 🚨 **NO S/P/E INTEGRATION**: SAE not integrated into drug efficacy pipeline (architectural mismatch)
3. 🚨 **NO HOTSPOT MUTATION DETECTION**: KRAS/BRAF/NRAS hotspot checking missing (Manager's C2)
4. 🚨 **MECHANISM FIT RANKER NOT WIRED**: Built but never called in trials endpoint
5. 🚨 **NO TRIAL MOA VECTORS**: Mechanism fit ranking has no data to rank with
6. ⚠️ **NO COHORT_OVERLAP**: Cannot determine patient-cohort fit (Manager's C5)
7. ⚠️ **NO TRIAL_KEYWORDS**: Cannot dynamically filter trials by SAE mechanisms
8. ⚠️ **NO ABCB1 EFFLUX**: Cross-resistance incomplete (only treatment history)
9. ⚠️ **NO WIWFM INTEGRATION**: SAE features not in drug responses
10. ⚠️ **NO VALIDATION**: Blocked on Jr2 HRD extraction

**Total Gaps Found:** **38 gaps** across all severity levels

---

## 🚨 CRITICAL ISSUES (MUST FIX)

### **ISSUE #1: DNA REPAIR CAPACITY FORMULA MISMATCH** 🔥

**Manager Said (C1):**
```
dna_repair_capacity = 0.6×DDR_burden + 0.2×essentiality_signal + 0.2×exon_disruption
```

**What We Implemented:**
```python
# api/services/sae_feature_service.py lines 35-39
DNA_REPAIR_CAPACITY_WEIGHTS = {
    "pathway_ddr": 0.50,        # ❌ WRONG: Should be 0.60
    "essentiality_hrr": 0.30,    # ❌ WRONG: Should be 0.20
    "functionality": 0.20        # ⚠️ DIFFERENT: Manager said exon_disruption, we used functionality
}
```

**Impact:**
- DNA repair capacity scores will be systematically wrong
- Manager's exact formula not followed
- All downstream decisions (PARP gating, resistance detection) affected

**Root Cause:**
- Zo copied wrong formula during Phase 2 implementation
- Tests validated the IMPLEMENTED formula, not the MANAGER'S formula
- No cross-check against source document

**Fix Required:**
- [ ] Update `DNA_REPAIR_CAPACITY_WEIGHTS` to `{ddr: 0.6, ess: 0.2, exon_disruption: 0.2}`
- [ ] Clarify: Does Manager want `exon_disruption` or `functionality`?
- [ ] Re-run all tests
- [ ] Validate against Manager's policy document

**Severity:** 🔥 CRITICAL - Core SAE formula is wrong

---

### **ISSUE #2: MECHANISM FIT RANKER NOT INTEGRATED** 🔥

**Manager Expected (P4):**
```
Trial ranking should use mechanism fit:
score = 0.7×eligibility + 0.3×mechanism_fit
```

**What We Built:**
- ✅ Created `mechanism_fit_ranker.py` (232 lines)
- ✅ Implemented exact α/β weighting
- ✅ All tests passing (6/6)

**What We DIDN'T Do:**
- ❌ Never imported into `ayesha_trials.py`
- ❌ Never called in trial ranking logic
- ❌ Trials still ranked by soft boost system only (no mechanism fit)

**Proof:**
```bash
$ grep "rank_trials_by_mechanism" oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_trials.py
# NO RESULTS
```

**Impact:**
- Mechanism fit ranking is **completely unused**
- Trials ranked by old soft boost system (10 boosts + 3 penalties)
- Manager's P4 policy not enforced
- HER2 trial (NCT06819007) not getting mechanism-based boost

**Root Cause:**
- Zo built service in isolation
- Never integrated into actual trials endpoint
- Jr3's E2E test only validated orchestrator response structure, not trial ranking logic

**Fix Required:**
- [ ] Import `rank_trials_by_mechanism` into `ayesha_trials.py`
- [ ] Call after soft boost filtering
- [ ] Add `mechanism_alignment` to trial response schema
- [ ] Re-test with real trials data

**Severity:** 🔥 CRITICAL - Manager's P4 policy not enforced

---

### **ISSUE #3: NO TRIAL MOA VECTORS (DATA GAP)** ⚠️

**Manager Expected (P3):**
```
Gemini trial tagging – offline only; batch tag 200 ovarian trials
Accept batch if ≥90% tag accuracy
```

**What We Have:**
- ✅ Mechanism fit ranker expects `trial.moa_vector` (7D)
- ❌ No actual trial MoA vectors exist
- ❌ Gemini offline tagging never executed
- ❌ No validation of tagging accuracy

**Current State:**
- `mechanism_fit_ranker.py` expects each trial to have `moa_vector: [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]`
- Real trials in database have NO `moa_vector` field
- Ranker will fail or use fallback (neutral vector = [0,0,0,0,0,0,0])

**Impact:**
- Mechanism fit ranking will always return 0.0 (neutral vector cosine similarity)
- Manager's strategic insight (mechanism-based ranking) completely lost
- HER2 trial will NOT get mechanism boost even if Ayesha is HER2+

**Root Cause:**
- Assumed trial MoA vectors exist
- Never validated data availability
- Never executed Gemini offline tagging workflow

**Fix Required:**
- [ ] Execute Gemini offline tagging for 200 ovarian trials
- [ ] Human spot-review 30 diverse trials (≥90% accuracy)
- [ ] Store `moa_vector` in trial database
- [ ] Version and track tagging metadata
- [ ] Re-test mechanism fit ranker with real vectors

**Severity:** ⚠️ HIGH - Mechanism fit ranking unusable without data

---

### **ISSUE #4: HER2 IHC NOT IN NEXT-TEST RECOMMENDER** ⚠️

**Manager Expected:**
- HER2 trial (NCT06819007) discovered
- Next-Test Recommender should flag HER2 IHC as critical gate
- Hint Tile: "📋 Order HER2 IHC NOW - Unlocks NCT06819007"

**What We Have:**
```python
# api/services/next_test_recommender.py
# Priority order: HRD → ctDNA → SLFN11 → ABCB1
# ❌ NO HER2 IHC in priority order
```

**Impact:**
- Ayesha might miss HER2 test
- Miss "1 in 700" trial opportunity
- Manager's strategic insight (HER2 gating) lost

**Root Cause:**
- Next-test recommender built before HER2 trial discovered
- Never updated after Phase 2 HER2 pathway integration
- No logic to dynamically add tests based on matched trials

**Fix Required:**
- [ ] Add dynamic test recommendation based on matched trials
- [ ] If HER2-targeted trial matched → Add "HER2 IHC" to next-test list
- [ ] Update hint tiles to flag HER2 when applicable
- [ ] Add HER2 pathway to 7D vector in mechanism_map_service.py (already done, but not in next-test)

**Severity:** ⚠️ MEDIUM-HIGH - Strategic opportunity missed

---

### **ISSUE #5: NO REAL TUMOR_CONTEXT TESTING** ⚠️

**What Was Tested:**
- ✅ Pre-NGS: `tumor_context: null` (gray "Awaiting NGS")
- ✅ API response structure validated
- ❌ Post-NGS: Never tested with real `tumor_context` payload

**What Wasn't Tested:**
- Post-NGS pathway: `tumor_context` with HRD=58, BRCA1 biallelic, somatic mutations
- DNA repair capacity calculation with real data
- Mechanism map color-coding (green/yellow/gray)
- Resistance detection with 2 timepoints (current + previous HRD)
- Mechanism fit ranking with real SAE vector

**Current Evidence:**
```python
# test_sae_phase3_e2e.py
# Only tests pre-NGS (no tumor_context)
payload = {
    "ca125_value": 2842.0,
    "stage": "IVB",
    # NO tumor_context
}
```

**Impact:**
- Post-NGS code path never executed
- DNA repair capacity never computed in real scenario
- Mechanism map always shows gray chips (never tested color-coding)
- Manager's C5 formula never validated with real data

**Root Cause:**
- Jr3 focused on pre-NGS flow (correct for Ayesha TODAY)
- Never created post-NGS test case
- Real TCGA data available but not used for testing

**Fix Required:**
- [ ] Create post-NGS test payload (use TCGA BRCA1 biallelic patient)
- [ ] Test DNA repair capacity calculation
- [ ] Test mechanism map color-coding
- [ ] Test resistance detection (mock 2 timepoints)
- [ ] Validate Manager's formula with real numbers

**Severity:** ⚠️ MEDIUM - Post-NGS path unvalidated

---

### **ISSUE #6: RESISTANCE ALERT NOT IN FRONTEND** ⚠️

**Backend Status:**
- ✅ `resistance_detection_service.py` operational (8/8 tests)
- ✅ Called in orchestrator (`detect_resistance()`)
- ✅ Returns `resistance_alert` in API response

**Frontend Status:**
- ❌ No `ResistanceAlertBanner.jsx` component created
- ❌ `resistance_alert` not extracted in `AyeshaTrialExplorer.jsx`
- ❌ Not displayed anywhere in UI

**Current Code:**
```javascript
// AyeshaTrialExplorer.jsx
const data = await response.json();

// ✅ Extracted:
setNextTestRecommender(data.next_test_recommender);
setHintTiles(data.hint_tiles);
setMechanismMap(data.mechanism_map);

// ❌ NOT Extracted:
// data.resistance_alert  <-- EXISTS but not displayed
// data.sae_features      <-- EXISTS but not displayed
```

**Impact:**
- Resistance alerts computed but invisible to user
- Manager's R2 (immediate alerts) not surfaced
- 2-of-3 trigger logic wasted (computed but not shown)

**Root Cause:**
- Jr3's mission document listed `ResistanceAlertBanner.jsx` but it was never created
- Focus was on Next Test, Hints, Mechanism Map (3/5 components)
- E2E test only validated response structure, not UI completeness

**Fix Required:**
- [ ] Create `ResistanceAlertBanner.jsx` component
- [ ] Extract and display `resistance_alert` in `AyeshaTrialExplorer.jsx`
- [ ] Show immediate alert banner when `resistance_detected: true`
- [ ] Display recommended actions and trial switches

**Severity:** ⚠️ MEDIUM - Backend logic not surfaced to user

---

## ⚠️ MEDIUM CONCERNS (SHOULD FIX)

### **CONCERN #7: NO GEMINI TRIAL MoA TAGGING DONE**

**Manager Policy (P3):**
- Offline tagging required before mechanism fit ranking
- Batch tag 200 trials → human review 30 → accept if ≥90% accuracy
- Never in runtime paths

**Current Reality:**
- ❌ No offline tagging executed
- ❌ No trial MoA vectors stored
- ❌ No human review done
- ❌ No accuracy validation

**Impact:**
- Mechanism fit ranker is theoretical only
- Cannot actually rank trials by mechanism
- Manager's strategic vision (SAE → trial alignment) undeliverable

**Fix Required:**
- [ ] Execute Gemini batch tagging workflow
- [ ] Generate MoA vectors for all 700 ovarian trials
- [ ] Human review 30 trials (spot-check accuracy)
- [ ] Store versioned MoA vectors in database
- [ ] Document tagging methodology

**Severity:** ⚠️ MEDIUM - Required for mechanism fit ranking

---

### **CONCERN #8: NO LONGITUDINAL TESTING**

**What Resistance Detection Requires:**
- Current HRD + Previous HRD (2 timepoints)
- Current DNA repair + Previous DNA repair
- CA-125 trend over cycles

**What We Tested:**
```python
# test_sae_phase2_services.py
detect_resistance(
    current_hrd=45.0,
    previous_hrd=60.0,  # ✅ Provided in test
    current_dna_repair_capacity=0.50,
    previous_dna_repair_capacity=0.72,  # ✅ Provided in test
    # ...
)
```

**What We DIDN'T Test:**
- ❌ Real 2-timepoint scenario (baseline → week 12 → week 24)
- ❌ Storage/retrieval of previous values
- ❌ Trigger logic over multiple cycles
- ❌ HR restoration pattern in clinical timeline

**Impact:**
- Resistance detection logic unvalidated in clinical use
- No storage mechanism for previous HRD/DNA repair values
- Unclear how to get "previous" values in production

**Root Cause:**
- Focused on single-timepoint logic
- No session/database storage plan for longitudinal data
- No clinical workflow design

**Fix Required:**
- [ ] Design storage for longitudinal biomarkers (HRD, DNA repair capacity, CA-125)
- [ ] Create session/database schema for timepoint tracking
- [ ] Test 3-timepoint scenario (baseline → cycle 3 → cycle 6)
- [ ] Document clinical workflow for longitudinal monitoring

**Severity:** ⚠️ MEDIUM - Clinical workflow incomplete

---

### **CONCERN #9: SAE FEATURES NOT DISPLAYED IN UI**

**Backend Computes:**
```python
# ayesha_orchestrator_v2.py returns:
{
  "sae_features": {
    "dna_repair_capacity": 0.75,
    "pathway_burden_ddr": 0.82,
    "pathway_burden_mapk": 0.15,
    # ... 7D mechanism vector, resistance signals, etc.
  }
}
```

**Frontend Displays:**
- ✅ Next Test Card (uses `next_test_recommender`)
- ✅ Hint Tiles (uses `hint_tiles`)
- ✅ Mechanism Chips (uses `mechanism_map`)
- ❌ SAE Features (computed but not displayed)

**Gap:**
- Jr3's mission mentioned `SAEFeaturesDisplay.jsx` but it was never created
- Focus shifted to Next Test, Hints, Mechanism Map
- Raw SAE features (DNA repair capacity, pathway burdens) not visible

**Impact:**
- Oncologist can't see computed SAE features
- Transparency lost (what drove the mechanism map colors?)
- Provenance incomplete (can't audit DNA repair capacity)

**Fix Required:**
- [ ] Create `SAEFeaturesDisplay.jsx` component
- [ ] Display DNA repair capacity (0-1 with interpretation)
- [ ] Display 7D mechanism vector (bar chart)
- [ ] Display resistance signals (if detected)
- [ ] Add provenance (formula, thresholds, sources)

**Severity:** ⚠️ MEDIUM - Transparency gap

---

### **CONCERN #10: NO MANUAL FRONTEND TESTING**

**What Was Tested:**
- ✅ Backend endpoint (`POST /api/ayesha/complete_care_v2`)
- ✅ API response structure (all keys present)
- ✅ Component linting (0 errors)

**What Was NOT Tested:**
- ❌ Actual browser rendering
- ❌ Responsive layout (mobile/desktop)
- ❌ Click interactions (mechanism chip popover)
- ❌ Hover effects (hint tiles)
- ❌ CTA buttons (scroll to NGS section)
- ❌ Loading states
- ❌ Error states

**Evidence:**
```
# Jr3's report says "Frontend Manual Testing:"
# 1. Start backend
# 2. Start frontend
# 3. Navigate to /ayesha-trials
# 4. Verify components render
# 
# BUT: No confirmation this was actually done
```

**Impact:**
- Unknown if UI actually works in browser
- Potential runtime errors undiscovered
- UX issues (layout, responsiveness) not validated

**Root Cause:**
- Jr3 focused on E2E API test (automated)
- Assumed manual testing would follow
- No checklist for visual validation

**Fix Required:**
- [ ] Actually run frontend in browser
- [ ] Navigate to `/ayesha-trials`
- [ ] Click through all interactions
- [ ] Test responsive layout
- [ ] Screenshot and verify visual appearance
- [ ] Document any UI bugs found

**Severity:** ⚠️ MEDIUM - Unknown if UI actually works

---

## 📋 MINOR GAPS (NICE TO FIX)

### **GAP #11: HER2 Pathway in Mechanism Map But Not Next-Test**

**Status:**
- ✅ HER2 pathway in 7D mechanism vector (SAE features)
- ✅ HER2 pathway in mechanism map chips (MechanismChips.jsx shows 6 chips, not 7)
- ❌ HER2 IHC not in next-test recommendations

**Wait... Mechanism Map Shows 6 Chips, Not 7:**
```javascript
// MechanismChips.jsx default chips
const displayChips = [
  { pathway: 'DDR', ... },
  { pathway: 'MAPK', ... },
  { pathway: 'PI3K', ... },
  { pathway: 'VEGF', ... },
  { pathway: 'IO', ... },
  { pathway: 'Efflux', ... },
  // ❌ NO HER2 CHIP!
];
```

**Impact:**
- HER2 pathway computed in backend (7D vector)
- But only 6 chips shown in UI (missing HER2)
- Inconsistent between backend (7D) and frontend (6 chips)

**Fix Required:**
- [ ] Add HER2 chip to MechanismChips.jsx default array
- [ ] Verify backend returns 7 chips, not 6
- [ ] Update mechanism_map_service.py if only returning 6

**Severity:** ⚠️ LOW-MEDIUM - Inconsistency in dimensions

---

### **GAP #12: Validation Script Blocked on Jr2**

**Status:**
- ✅ `validate_sae_tcga.py` created
- ✅ TCGA data identified (`hrd_tcga_ov_labeled_sample_use_evo.json`)
- ❌ TCGA data has `outcome_platinum` but NO `hrd_score`
- ⏸️ Blocked on Jr2 HRD extraction

**Current State:**
- Jr2 assigned to extract HRD scores from cBioPortal
- `calculate_full_hrd_scores.py` script exists (428 lines)
- `validate_hrd_scores.py` script exists (320 lines)
- But no execution or results yet

**Impact:**
- Cannot validate DNA repair capacity formula against ground truth
- Cannot validate resistance detection against real outcomes
- Manager's C5 formula unvalidated
- All SAE features are theoretical

**Root Cause:**
- Data dependency on cBioPortal HRD scores
- Jr2 not yet executing mission
- No pressure/urgency communicated

**Fix Required:**
- [ ] Jr2: Execute HRD extraction (use calculate_full_hrd_scores.py)
- [ ] Validate extracted HRD scores (use validate_hrd_scores.py)
- [ ] Update TCGA JSON with `hrd_score` field
- [ ] Zo: Run validation script
- [ ] Generate validation report

**Severity:** ⚠️ MEDIUM - Validation completely blocked

---

### **GAP #13: No Error Handling for Missing SAE Dependencies**

**What Happens If:**
- Insights bundle call fails?
- Pathway scores unavailable?
- CA-125 intelligence errors?

**Current Code:**
```python
# ayesha_orchestrator_v2.py
try:
    results["sae_features"] = compute_sae_features(
        insights_bundle=insights_bundle,  # ❓ What if None?
        pathway_scores=pathway_scores,     # ❓ What if {}?
        # ...
    )
except Exception as e:
    logger.error(f"SAE Phase 2 failed: {e}")
    results["sae_features"] = {"error": str(e)}
```

**Gaps:**
- Graceful error handling exists
- But `insights_bundle` and `pathway_scores` are ASSUMED to exist
- No validation before calling `compute_sae_features()`
- Could crash if upstream services fail

**Impact:**
- Potential runtime errors if insights/pathway calls fail
- SAE features might return error dict
- Frontend expects structured data, might break on error dict

**Fix Required:**
- [ ] Add validation before calling SAE services
- [ ] Provide fallback values for missing inputs
- [ ] Test error scenarios (insights fail, pathway fail, etc.)
- [ ] Document error handling strategy

**Severity:** ⚠️ LOW-MEDIUM - Edge case handling

---

### **GAP #14: No Provenance in Frontend**

**Backend Returns:**
```json
{
  "provenance": {
    "orchestrator": "complete_care_v2",
    "timestamp": "2025-01-13T...",
    "services_called": ["trials", "ca125", "sae_features", ...]
  }
}
```

**Frontend Displays:**
- ❌ Provenance not extracted
- ❌ Not displayed anywhere
- ❌ No run ID visible
- ❌ No timestamp visible

**Impact:**
- Auditability lost
- Can't trace which services were called
- Can't debug issues (no run ID)
- Manager's transparency requirement not met

**Fix Required:**
- [ ] Extract `provenance` in AyeshaTrialExplorer.jsx
- [ ] Create ProvenanceBar component
- [ ] Display run ID, timestamp, services called
- [ ] Add to bottom of page (footer style)

**Severity:** ⚠️ LOW - Nice to have for debugging

---

## 🔍 TESTING GAPS

### **GAP #15: No Integration Tests**

**What We Have:**
- ✅ Unit tests for each service (47/47 passing)
- ✅ E2E API test (backend response validation)
- ❌ No integration tests (service → orchestrator → endpoint)

**What's Missing:**
- Integration test: orchestrator calls all 6 SAE services
- Integration test: trials endpoint applies mechanism fit ranking
- Integration test: complete_care_v2 returns all expected sections

**Impact:**
- Service interactions not validated
- Orchestration flow not tested end-to-end
- Integration bugs might exist

**Fix Required:**
- [ ] Create `test_sae_integration.py`
- [ ] Test orchestrator → service chain
- [ ] Test error propagation
- [ ] Test graceful degradation

**Severity:** ⚠️ LOW - Unit tests cover logic, but integration untested

---

### **GAP #16: No Performance Testing**

**Questions:**
- How long does complete_care_v2 take with all SAE services?
- How many API calls does it make?
- What's the P95 latency?
- Does it timeout on slow networks?

**Current State:**
- ❌ No performance benchmarks
- ❌ No latency measurements
- ❌ No timeout testing
- ❌ No load testing

**Impact:**
- Unknown if endpoint is production-ready
- Could be too slow for clinical use
- No SLA guarantees

**Fix Required:**
- [ ] Measure complete_care_v2 latency (P50, P95, P99)
- [ ] Count API calls made
- [ ] Test with slow network
- [ ] Validate timeout handling
- [ ] Document performance benchmarks

**Severity:** ⚠️ LOW - Production concern, not blocker

---

## 🔬 VALIDATION GAPS

### **GAP #17: No Real-World Validation**

**Manager Expected:**
- Validate DNA repair capacity against TCGA HRD scores
- Validate resistance detection against platinum outcomes
- Metrics: AUROC ≥0.70, Sensitivity ≥0.80, Correlation ≥0.60

**Current State:**
- ⏸️ Completely blocked on Jr2 HRD extraction
- ✅ Validation script ready (`validate_sae_tcga.py`)
- ❌ No execution, no results, no metrics

**Impact:**
- All SAE formulas unvalidated
- Manager's C5 formula might be wrong (see Issue #1)
- No evidence SAE features correlate with real outcomes

**Severity:** 🔥 HIGH - Scientific validity unproven

---

### **GAP #18: No Clinical Expert Review**

**What We Built:**
- DNA repair capacity formula (Manager's C5)
- Resistance detection (2-of-3 triggers)
- Mechanism fit ranking (α/β weighting)

**What We DIDN'T Do:**
- ❌ No oncologist review of logic
- ❌ No validation against clinical guidelines
- ❌ No comparison to clinician intuition

**Impact:**
- Logic might be clinically incorrect
- Thresholds might be wrong
- Recommendations might conflict with NCCN

**Fix Required:**
- [ ] Clinical expert review of all SAE logic
- [ ] Validate thresholds against literature
- [ ] Compare recommendations to NCCN guidelines
- [ ] Document clinical rationale for all decisions

**Severity:** ⚠️ MEDIUM - Clinical validity unvalidated

---

## 📊 DATA GAPS

### **GAP #19: No Real Ayesha Tumor Data**

**Current Situation:**
- Ayesha is treatment-naive (pre-NGS)
- No tumor NGS data available
- No HRD score
- No somatic mutations
- No TMB/MSI

**Testing Approach:**
- Used TCGA synthetic data
- Used mock tumor_context payloads
- Never tested with Ayesha's ACTUAL data (because it doesn't exist yet)

**Impact:**
- All post-NGS logic is theoretical for Ayesha
- Cannot validate with her specific case
- Recommendations based on population priors, not her genomics

**Resolution:**
- ⏸️ Wait for Ayesha to get NGS (7-10 days after ordering)
- ✅ Fast-track checklist prepared
- ✅ Services ready to ingest real data when available

**Severity:** ⚠️ LOW - Expected (treatment-naive), not a bug

---

### **GAP #20: No Trial MoA Vectors Stored Anywhere**

**Expected:**
- 700 ovarian cancer trials
- Each trial has `moa_vector: [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]`
- Stored in database with versioning

**Reality:**
```bash
$ grep "moa_vector" oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_trials.py
# NO RESULTS
```

**Current Database Schema:**
- Trials have: NCT ID, title, phase, status, eligibility_text, locations
- ❌ NO `moa_vector` field
- ❌ NO mechanism tagging
- ❌ NO Gemini results stored

**Impact:**
- Mechanism fit ranking has no data to work with
- All mechanism_fit_scores will be 0.0 (neutral fallback)
- Strategic vision (SAE → mechanism ranking) undeliverable

**Fix Required:**
- [ ] Add `moa_vector` field to trial schema
- [ ] Execute Gemini tagging workflow
- [ ] Store vectors in database
- [ ] Update trials endpoint to return moa_vector
- [ ] Test mechanism fit ranking with real vectors

**Severity:** 🔥 HIGH - Critical data missing

---

## 🔧 TECHNICAL DEBT

### **DEBT #21: Duplicate Code in Orchestrator**

**Issue:**
- `ayesha_orchestrator_v2.py` had duplicate sections (fixed)
- But code is still monolithic (600+ lines)
- Hard to maintain
- Hard to test in isolation

**Future Refactoring Needed:**
- [ ] Extract service calling logic into separate modules
- [ ] Create orchestration service class
- [ ] Separate response assembly from service calls
- [ ] Add dependency injection for testability

**Severity:** ⚠️ LOW - Technical debt, not blocker

---

### **DEBT #22: No Caching Strategy**

**Current State:**
- Every API call re-computes everything
- CA-125 intelligence re-computed on every request
- Trials re-ranked on every request
- SAE features re-computed on every request

**Impact:**
- Slow response times
- Wasted compute
- Potential timeouts for complex patients

**Fix Required:**
- [ ] Add Redis/in-memory caching
- [ ] Cache CA-125 intelligence (TTL: 1 hour)
- [ ] Cache SAE features (TTL: 1 hour)
- [ ] Cache trial rankings (TTL: 24 hours)
- [ ] Add cache invalidation strategy

**Severity:** ⚠️ LOW - Performance optimization

---

### **DEBT #23: No Rate Limiting**

**Current State:**
- No rate limits on `/api/ayesha/complete_care_v2`
- Could be called repeatedly
- No DOS protection

**Impact:**
- Potential abuse
- High costs (if deployed with external APIs)
- No fair use policy

**Fix Required:**
- [ ] Add rate limiting (e.g., 10 requests/minute per user)
- [ ] Add usage tracking
- [ ] Add quota enforcement

**Severity:** ⚠️ LOW - Production hardening

---

## 🚨 ADDITIONAL GAPS FROM MANAGER'S MASTER PLAN

### **GAP #26: NO HOTSPOT MUTATION DETECTION** 🔥

**Manager's Plan (19.3, 19.11, 19.12):**
- **RAS/MAPK hotspot (KRAS/BRAF/NRAS)** or `pathway_burden.mapk ≥ 0.7` → surface MEK/RAF trial candidates
- **Hotspot_mutation in RAS/MAPK** → prioritize MEK/RAF trial candidates; deprioritize MEK monotherapy if absent

**What We Built:**
```python
# sae_feature_service.py
# NO hotspot_mutation detection
# NO KRAS/BRAF/NRAS hotspot checking
# NO COSMIC database integration
```

**What's Missing:**
- ❌ No `hotspot_mutation` field in SAEFeatures dataclass
- ❌ No COSMIC hotspot database
- ❌ No KRAS G12C/G12D/G12V detection
- ❌ No BRAF V600E detection
- ❌ No NRAS Q61 detection
- ❌ Manager explicitly said "use COSMIC/hardcoded list" (C2) - we didn't

**Impact:**
- MEK/RAF trials not correctly prioritized
- MAPK pathway logic incomplete
- Manager's C2 policy not implemented
- No ability to deprioritize MEK monotherapy when hotspot absent

**Fix Required:**
- [ ] Add `hotspot_mutation: bool` to SAEFeatures dataclass
- [ ] Create COSMIC hotspot list (KRAS, BRAF, NRAS variants)
- [ ] Add hotspot detection in `compute_sae_features()`
- [ ] Update hint tiles to flag MEK/RAF when hotspot detected
- [ ] Add trial keywords for MEK/RAF inhibitors

**Severity:** 🔥 CRITICAL - Manager's C2 action rule not implemented

---

### **GAP #27: NO COHORT_OVERLAP CALCULATION** ⚠️

**Manager's Plan (C5, 19.11):**
- **Cohort_overlap low + confidence low** → push trials; high overlap → lean standard with modest lift
- Definition: overlap of patient phenotype with literature/trial cohorts

**What We Built:**
```python
# sae_feature_service.py
# NO cohort_overlap field
# NO literature cohort matching
# NO trial cohort similarity
```

**What's Missing:**
- ❌ No `cohort_overlap` field in SAEFeatures
- ❌ No cohort database
- ❌ No phenotype matching logic
- ❌ No "push trials" banner when overlap low

**Impact:**
- Cannot determine if patient fits known cohorts
- Cannot adjust confidence based on cohort representativeness
- Manager's C5 action rule not enforced
- No guidance on trials vs SOC based on cohort fit

**Fix Required:**
- [ ] Add `cohort_overlap: float` to SAEFeatures dataclass
- [ ] Define cohort matching logic (disease + stage + key biomarkers)
- [ ] Update hint tiles to show "push trials" when overlap <0.40
- [ ] Add confidence lift when overlap ≥0.70

**Severity:** ⚠️ HIGH - Manager's C5 policy missing

---

### **GAP #28: NO TRIAL_KEYWORDS IN RESPONSE** ⚠️

**Manager's Plan (19.6, 19.18):**
- Add `trial_keywords` from SAE and Resistance Playbook
- Examples: "ATR inhibitor", "PARP combination", "MEK inhibitor", "WEE1 inhibitor"
- Boost trials whose MoA matches current SAE risks/combos

**What We Built:**
```python
# ayesha_orchestrator_v2.py response
{
  "sae_features": {...},
  "resistance_alert": {...},
  # ❌ NO trial_keywords field
}
```

**What's Missing:**
- ❌ No `trial_keywords: List[str]` in response
- ❌ No keyword generation from SAE features
- ❌ No trial filtering by keywords
- ❌ No "combo-ready" badges

**Impact:**
- Cannot dynamically filter trials by SAE-detected mechanisms
- Frontend cannot highlight matching trials
- No "ATR/CHK1" keyword when HR restoration detected
- No "MEK/RAF" keyword when MAPK hotspot detected

**Fix Required:**
- [ ] Add `trial_keywords` to SAEFeatures dataclass
- [ ] Generate keywords based on SAE features:
  - High DDR → ["PARP inhibitor", "ATR inhibitor", "CHK1 inhibitor"]
  - MAPK hotspot → ["MEK inhibitor", "RAF inhibitor"]
  - PI3K burden → ["PI3K inhibitor", "mTOR inhibitor", "WEE1 inhibitor"]
  - Resistance → ["ATR combination", "PARP combination"]
- [ ] Add to API response
- [ ] Frontend uses keywords to badge/filter trials

**Severity:** ⚠️ MEDIUM-HIGH - Trial matching incomplete

---

### **GAP #29: NO DIFFERENTIAL BRANCHES IN NEXT-TEST** ⚠️

**Manager's Plan (C6, 19.14):**
- Use "differential branches" format: "If + → X; If - → Y"
- Each test emits binary branch (e.g., "If SLFN11+, continue PARP; if −, pivot to ATR trial")

**What We Have:**
```python
# next_test_recommender.py
# Returns: test_name, rationale, impact_if_positive, impact_if_negative
# ✅ CORRECT STRUCTURE
```

**What's Missing:**
- ⚠️ Frontend doesn't display differential branches clearly
- ⚠️ No visual branch diagram
- ⚠️ No "binary decision tree" UI

**Impact:**
- Differential branches exist in backend but not emphasized in UI
- Clinicians might miss the "if/then" logic
- Manager's intent (clear branching decisions) not fully realized

**Fix Required:**
- [ ] Update NextTestCard.jsx to emphasize branches
- [ ] Add visual branch diagram (If + → green path, If - → yellow path)
- [ ] Make branches more prominent than current implementation
- [ ] Add icons for positive/negative paths

**Severity:** ⚠️ MEDIUM - UI emphasis issue, logic exists

---

### **GAP #30: NO S/P/E INTEGRATION** 🔥

**Manager's Plan (19.1, 19.2):**
- SAE is PART OF S/P/E pipeline: "S/P/E base ± SAE lifts/penalties"
- Data flow: Evo2 → Insights → Pathway → **SAE extractor** → Confidence computation
- SAE should modulate S/P/E scores, not replace them

**What We Built:**
- ✅ SAE features computed independently
- ❌ NO integration with S/P/E pipeline
- ❌ SAE lives in Ayesha orchestrator, not in drug efficacy router
- ❌ No confidence lifts/penalties based on SAE

**Current Architecture:**
```
Ayesha Orchestrator:
  - Trials (separate)
  - SOC (separate)
  - CA-125 (separate)
  - SAE Features (separate)  ← ISOLATED, NOT INTEGRATED

Should Be:
Drug Efficacy Router (WIWFM):
  S → P → E → SAE → Confidence → Final Score
```

**Impact:**
- SAE features computed but NOT used in drug ranking
- No confidence modulation (+0.02 to +0.08 lifts)
- No PARP penalty for germline-negative unless HRD≥42
- Manager's entire S/P/E → SAE pipeline not implemented
- SAE is "display only", not affecting recommendations

**Fix Required:**
- [ ] Integrate SAE into drug efficacy router (`/api/efficacy/predict`)
- [ ] Compute SAE features WITHIN efficacy pipeline
- [ ] Apply confidence lifts/penalties based on SAE
- [ ] Update EfficacyResponse to include SAE attribution
- [ ] Show SAE impact in EvidenceBand
- [ ] Add provenance: `confidence_breakdown: {s: 0.42, p: 0.28, e: 0.10, sae: 0.05}`

**Severity:** 🔥 CRITICAL - Architectural mismatch with Manager's vision

---

### **GAP #31: NO SLFN11 IHC IN NEXT-TEST** ⚠️

**Manager's Plan (C6, 19.14, 19.19, 19.21):**
- Priority order: HRD → ctDNA → **SLFN11 IHC** → ABCB1
- "Order SLFN11 IHC; if low, PARP sensitivity reduced; consider ATR/PLK1 trial"
- SLFN11-low → reduced PARP sensitivity; prioritize ATR/PLK1

**What We Have:**
```python
# next_test_recommender.py
# Priority: HRD → ctDNA → SLFN11 → ABCB1
# ✅ SLFN11 IS INCLUDED (priority 3)
```

**But:**
- ⚠️ SLFN11 logic not dynamic (always priority 3)
- ⚠️ No conditional logic: "If PARP planned → elevate SLFN11"
- ⚠️ No branching: "If SLFN11+, continue PARP; if -, pivot ATR"

**Impact:**
- SLFN11 shown but not contextualized
- No emphasis when PARP is recommended
- No clear action based on SLFN11 result

**Fix Required:**
- [ ] Add dynamic prioritization: If DNA repair high → elevate SLFN11 (priority 2)
- [ ] Add clear branches: "If +: PARP sensitivity preserved; If -: Consider ATR/PLK1"
- [ ] Link to resistance playbook (SLFN11-low = PARP resistance risk)

**Severity:** ⚠️ MEDIUM - Exists but not contextualized

---

### **GAP #32: NO "PIVOT EARLY" BANNER** ⚠️

**Manager's Plan (19.17):**
- **Pivot early** banner when CA‑125 rule + SAE resistance both hit
- Immediate visual alert to switch therapy

**What We Have:**
- ✅ Resistance detection backend (2-of-3 triggers)
- ✅ CA-125 intelligence (inadequate response detection)
- ❌ No "Pivot Early" banner in UI
- ❌ ResistanceAlertBanner.jsx not created

**Impact:**
- Resistance detected but no urgent visual alert
- Clinician might miss early resistance signal
- Manager's "immediate alert" intent not realized

**Fix Required:**
- [ ] Create ResistanceAlertBanner.jsx with prominent styling
- [ ] Display at top of page when `resistance_detected: true`
- [ ] Label: "⚠️ PIVOT EARLY: Resistance Detected"
- [ ] Show triggers met + recommended actions

**Severity:** ⚠️ MEDIUM - UI prominence issue (already in Gap #6)

---

### **GAP #33: NO ONE-CLICK ORDERS** ⚠️

**Manager's Plan (19.17):**
- **One‑click orders** next to hints: SLFN11 IHC, ctDNA HRD update, ABCB1 proxy
- Routes to `/api/hints/next_test`

**What We Have:**
- ✅ Next-test recommender shows tests
- ❌ No "Order Now" buttons
- ❌ No `/api/hints/next_test` endpoint
- ❌ No integration with ordering system

**Impact:**
- Extra friction to order tests
- Clinician must manually copy test name
- No streamlined workflow

**Fix Required:**
- [ ] Add "Order Test" button in NextTestCard
- [ ] Create `/api/hints/next_test` endpoint (if ordering system exists)
- [ ] OR: Display copy-to-clipboard with test order details
- [ ] Show turnaround time + cost prominently

**Severity:** ⚠️ LOW - UX enhancement

---

### **GAP #34: NO MECHANISM MAP IN HINT TILES** ⚠️

**Manager's Plan (19.17, 19.54):**
- **Mechanism Map** strip: DDR | Efflux | MAPK | PI3K | IO with green/amber/red chips
- Should be under/near EvidenceBand

**What We Have:**
- ✅ MechanismChips.jsx component created
- ✅ Displays 6 chips (DDR, MAPK, PI3K, VEGF, IO, Efflux)
- ⚠️ But: Missing HER2 chip (should be 7)
- ⚠️ Not near EvidenceBand (separate section)

**Impact:**
- Mechanism map not contextually linked to drug recommendations
- HER2 pathway missing from UI (7D backend, 6 chips frontend)

**Fix Required:**
- [ ] Add HER2 chip to MechanismChips.jsx (make it 7 chips)
- [ ] Verify backend returns 7-chip mechanism_map
- [ ] Move/duplicate mechanism map near drug recommendations (if WIWFM integrated)

**Severity:** ⚠️ MEDIUM - Already identified in Gap #11

---

### **GAP #35: NO WIWFM INTEGRATION WITH SAE** 🔥

**Manager's Plan (19.7, 19.8):**
- SAE features should appear in **WIWFM drug response**
- Each drug should show: `sae_features: {exon_disruption, dna_repair_capacity, pathway_burden}`
- EvidenceBand should include SAE attribution chip

**What We Built:**
- ✅ SAE features computed in Ayesha orchestrator
- ❌ NOT integrated into WIWFM (`/api/efficacy/predict`)
- ❌ Drug responses don't include SAE features
- ❌ EvidenceBand doesn't show SAE impact

**Current State:**
```
WIWFM Response:
{
  "drug": "Olaparib",
  "efficacy_score": 0.73,
  "confidence": 0.85,
  // ❌ NO sae_features here
}
```

**Manager Expected:**
```
{
  "drug": "Olaparib",
  "efficacy_score": 0.73,
  "confidence": 0.85,
  "sae_features": {
    "dna_repair_capacity": 0.82,
    "pathway_burden": {"ddr": 0.78}
  }
}
```

**Impact:**
- SAE features siloed in Ayesha orchestrator
- Not available in main WIWFM drug ranking
- No SAE explanation for drug recommendations
- Manager's core vision (SAE → S/P/E → drug ranking) not realized

**Fix Required:**
- [ ] Integrate `compute_sae_features()` into efficacy router
- [ ] Add `sae_features` to EfficacyResponse schema
- [ ] Update EvidenceBand to show SAE chips
- [ ] Add confidence breakdown with SAE contribution

**Severity:** 🔥 CRITICAL - Core integration missing

---

### **GAP #36: NO PRE-COMPUTED CARE PATHWAYS** ⚠️

**Manager's Plan (C10, 19.20):**
- **If platinum partial response + SAE HR restoration** → line‑ready ATR combo trials ranked with logistics
- **If no MAPK/PI3K SAE signal** → hide those trials by default; move to "Explore more" bucket

**What We Have:**
- ✅ Resistance detection (HR restoration pattern)
- ✅ Trial ranking (soft boosts)
- ❌ No "care pathway" assembly
- ❌ No dynamic trial filtering based on SAE signals
- ❌ No "Explore more" bucket

**Impact:**
- All trials shown equally (no SAE-based filtering)
- Trials without mechanism fit not hidden
- No streamlined "ready-to-serve" pathways

**Fix Required:**
- [ ] Add care pathway logic to orchestrator
- [ ] Filter trials by SAE mechanism fit (hide low-fit trials)
- [ ] Create "Explore More" section for hidden trials
- [ ] Pre-rank ATR combo trials when HR restoration detected

**Severity:** ⚠️ MEDIUM - UX optimization

---

### **GAP #37: NO ABCB1 EFFLUX DETECTION** ⚠️

**Manager's Plan (C4, 19.3, 19.14, 19.19):**
- **Cross_resistance_risk high with prior taxane/ABCB1** → avoid substrates; propose non‑substrates
- ABCB1 inference before expression: use CNV if available; otherwise UNKNOWN
- ABCB1 expression proxy (in next-test recommender)

**What We Built:**
```python
# sae_feature_service.py
def _compute_cross_resistance_risk(...):
    # Uses treatment_history for prior taxane
    # ❌ NO ABCB1 detection
    # ❌ NO CNV checking
    # ❌ NO "ABCB1 proxy" logic
```

**What's Missing:**
- ❌ No ABCB1 CNV detection
- ❌ No ABCB1 expression proxy
- ❌ No substrate/non-substrate lists (PharmGKB/DrugBank)
- ❌ No "ABCB1 proxy" in next-test recommender (priority 4)

**Impact:**
- Cross-resistance risk incomplete (only uses treatment history)
- Cannot warn about efflux-mediated resistance
- Cannot recommend non-substrate regimens

**Fix Required:**
- [ ] Add ABCB1 CNV detection to cross-resistance logic
- [ ] Create substrate/non-substrate drug lists
- [ ] Add "ABCB1 expression proxy" to next-test recommender (priority 4)
- [ ] Update hint tiles to flag substrate avoidance

**Severity:** ⚠️ MEDIUM-HIGH - Resistance mechanism incomplete

---

### **GAP #38: NO CONFIDENCE BREAKDOWN IN PROVENANCE** ⚠️

**Manager's Plan (19.9, 19.16):**
- Provenance must include: `confidence_breakdown: {s: 0.42, p: 0.28, e: 0.10, sae: 0.05, gates: 0.00}`
- Show which SAE features affected confidence

**What We Have:**
```python
# ayesha_orchestrator_v2.py
provenance = {
    "orchestrator": "complete_care_v2",
    "timestamp": "...",
    # ❌ NO confidence_breakdown
    # ❌ NO SAE attribution
}
```

**Impact:**
- Cannot audit confidence calculation
- Cannot see SAE contribution
- Transparency incomplete

**Fix Required:**
- [ ] Add `confidence_breakdown` to provenance
- [ ] Track SAE lifts/penalties
- [ ] Show which features contributed
- [ ] Display in UI (ProvenanceBar or expandable section)

**Severity:** ⚠️ MEDIUM - Transparency gap

---

## 📝 DOCUMENTATION GAPS

### **GAP #24: No API Documentation for SAE Fields**

**Current State:**
- SAE fields in response schema
- But no explanation of what they mean
- No documentation of formulas
- No examples

**Impact:**
- Frontend developers don't know how to interpret fields
- Integration teams confused
- No API reference docs

**Fix Required:**
- [ ] Document all SAE response fields
- [ ] Add examples for each field
- [ ] Document formulas and thresholds
- [ ] Create OpenAPI/Swagger docs

**Severity:** ⚠️ LOW - Documentation

---

### **GAP #25: No Clinical Workflow Documentation**

**Questions:**
- When should oncologist order HRD test?
- When should resistance alert trigger clinical action?
- When should hint tiles be acted upon?

**Current State:**
- Technical documentation exists
- Clinical workflow not documented
- No user guide for oncologists

**Impact:**
- Oncologists might not know how to use system
- Clinical value not realized
- Adoption barriers

**Fix Required:**
- [ ] Create clinical workflow guide
- [ ] Document decision points
- [ ] Add usage examples
- [ ] Create oncologist training materials

**Severity:** ⚠️ LOW - Adoption concern

---

## 🎯 PRIORITY MATRIX

### **🔥 P0: CRITICAL (MUST FIX BEFORE DEMO)**

1. **DNA Repair Capacity Formula Mismatch** (Issue #1)
   - Manager said 0.6×DDR, we did 0.5×DDR
   - **Fix Time:** 30 minutes
   - **Blocker:** Demo credibility

2. **NO S/P/E INTEGRATION** (Gap #30) 🔥🔥🔥
   - SAE not integrated into drug efficacy pipeline
   - **Fix Time:** 8-12 hours (major refactor)
   - **Blocker:** ARCHITECTURAL MISMATCH - Manager's core vision

3. **NO HOTSPOT MUTATION DETECTION** (Gap #26) 🔥
   - KRAS/BRAF/NRAS hotspot checking missing
   - **Fix Time:** 2-3 hours
   - **Blocker:** Manager's C2 policy not implemented

4. **Mechanism Fit Ranker Not Wired** (Issue #2)
   - Built but never integrated
   - **Fix Time:** 1 hour
   - **Blocker:** Manager's P4 policy

5. **No Trial MoA Vectors** (Issue #3 + Gap #20)
   - Mechanism ranking has no data
   - **Fix Time:** 4-6 hours (Gemini tagging workflow)
   - **Blocker:** Mechanism fit ranking unusable

---

### **⚠️ P1: HIGH (FIX SOON)**

6. **NO COHORT_OVERLAP** (Gap #27)
   - Cannot determine patient-cohort fit
   - **Fix Time:** 3-4 hours
   - **Impact:** Manager's C5 policy missing

7. **NO TRIAL_KEYWORDS** (Gap #28)
   - Cannot dynamically filter trials
   - **Fix Time:** 2 hours
   - **Impact:** Trial matching incomplete

8. **NO ABCB1 EFFLUX** (Gap #37)
   - Cross-resistance incomplete
   - **Fix Time:** 3-4 hours
   - **Impact:** Resistance mechanism incomplete

9. **NO WIWFM INTEGRATION** (Gap #35)
   - SAE features not in drug responses
   - **Fix Time:** 4-6 hours
   - **Impact:** Core drug ranking missing SAE

10. **HER2 IHC Not in Next-Test** (Issue #4)
    - Strategic opportunity missed
    - **Fix Time:** 1 hour
    - **Impact:** Trial access

11. **No Real Tumor_Context Testing** (Issue #5)
    - Post-NGS path unvalidated
    - **Fix Time:** 2 hours
    - **Impact:** Production reliability

12. **Resistance Alert Not in UI** (Issue #6)
    - Backend computes, frontend doesn't show
    - **Fix Time:** 2 hours
    - **Impact:** Clinical value lost

13. **No Validation Against Ground Truth** (Gap #17)
    - Blocked on Jr2
    - **Fix Time:** TBD (Jr2 + 2 hours)
    - **Impact:** Scientific validity

---

### **📋 P2: MEDIUM (TECHNICAL DEBT)**

14. No Longitudinal Testing (Concern #8)
15. SAE Features Not Displayed (Concern #9)
16. No Manual Frontend Testing (Concern #10)
17. HER2 in Vector But Not UI (Gap #11)
18. No Gemini Tagging Done (Concern #7)
19. Differential Branches Not Emphasized (Gap #29)
20. SLFN11 Not Contextualized (Gap #31)
21. No Confidence Breakdown (Gap #38)
22. No Care Pathways (Gap #36)

---

### **⚪ P3: LOW (POLISH)**

23. "Pivot Early" Banner (Gap #32 - duplicate of #6)
24. One-Click Orders (Gap #33)
25. Mechanism Map Placement (Gap #34)
26. No Integration Tests (Gap #15)
27. No Performance Testing (Gap #16)
28. No Caching Strategy (Debt #22)
29. No Rate Limiting (Debt #23)
30. No API Documentation (Gap #24)
31. No Clinical Workflow Docs (Gap #25)

---

## ⚔️ COMMANDER - HONEST ASSESSMENT

### **What We Achieved:**
- ✅ 6 backend services (2,300+ lines)
- ✅ 3 frontend components (600+ lines)
- ✅ 47/47 tests passing
- ✅ 5.5 hours (vs 12h planned)

### **What We Missed:**
- 🚨 Formula mismatch (0.5 vs 0.6)
- 🚨 Mechanism fit not integrated
- 🚨 No trial MoA vectors
- ⚠️ Post-NGS path not tested
- ⚠️ Resistance alerts not in UI
- ⚠️ No validation against ground truth

### **Operational Readiness:**
- **Pre-NGS (Ayesha TODAY):** 85% ready
  - Next-test recommender: ✅ Works
  - Hint tiles: ✅ Works
  - Mechanism map: ✅ Works (gray)
  - Formula issue: ❌ Wrong (but she's pre-NGS so not used yet)
  
- **Post-NGS (When HRD Returns):** 60% ready
  - DNA repair capacity: ❌ Wrong formula
  - Mechanism fit ranking: ❌ Not wired
  - Resistance detection: ⚠️ Backend works, UI missing
  - SAE features: ⚠️ Computed but not displayed

### **Risk Assessment:**
- **Demo Risk:** 🔥 HIGH - Formula mismatch could embarrass if asked
- **Clinical Risk:** ⚠️ MEDIUM - Pre-NGS flow safe, post-NGS flow untested
- **Strategic Risk:** 🔥 HIGH - Mechanism fit vision undeliverable (no MoA vectors)

---

## 🎯 RECOMMENDED NEXT STEPS

### **IMMEDIATE (Before Demo/Production):**

1. **Fix DNA Repair Capacity Formula** (30 min) 🔥
   - Update weights to 0.6/0.2/0.2
   - Clarify exon_disruption vs functionality with Manager
   - Re-run tests

2. **Wire Mechanism Fit Ranker** (1 hour) 🔥
   - Integrate into ayesha_trials.py
   - Test with mock MoA vectors
   - Document integration

3. **Create Resistance Alert UI** (2 hours) ⚠️
   - Build ResistanceAlertBanner.jsx
   - Display in AyeshaTrialExplorer
   - Test with mock resistance scenario

4. **Manual Frontend Testing** (1 hour) ⚠️
   - Actually run frontend in browser
   - Verify all components render
   - Test interactions
   - Screenshot results

---

### **SHORT-TERM (This Week):**

5. **Execute Gemini MoA Tagging** (4-6 hours) 🔥
   - Batch tag 200 trials
   - Human review 30 trials
   - Store MoA vectors
   - Document methodology

6. **Test Post-NGS Path** (2 hours) ⚠️
   - Create real tumor_context payload
   - Test DNA repair capacity calculation
   - Test mechanism map color-coding
   - Validate against expectations

7. **Add HER2 to Next-Test Logic** (1 hour) ⚠️
   - Dynamic test recommendation
   - Add HER2 IHC when HER2 trial matched
   - Update hint tiles

---

### **MEDIUM-TERM (Next Week):**

8. **Jr2: Extract HRD Scores** (TBD)
   - Use calculate_full_hrd_scores.py
   - Update TCGA JSON
   - Validate extraction

9. **Zo: Run Validation** (2 hours)
   - Execute validate_sae_tcga.py
   - Compute metrics
   - Generate validation report

10. **Clinical Expert Review** (External)
    - Review all SAE logic
    - Validate thresholds
    - Get oncologist feedback

---

## ⚔️ FINAL VERDICT

**Overall Grade:** **C+ (75/100)** *(Downgraded from B+ after discovering architectural mismatch)*

**What Went Right:**
- ✅ Fast execution (55% faster than planned)
- ✅ All unit tests passing (47/47)
- ✅ Clean code architecture (services well-structured)
- ✅ Pre-NGS flow operational (for Ayesha TODAY)

**What Went CRITICALLY Wrong:**
- 🚨 **ARCHITECTURAL MISMATCH**: SAE not integrated into S/P/E drug efficacy pipeline (Manager's core vision)
- 🚨 Formula mismatch (0.5 vs 0.6 - critical math error)
- 🚨 No hotspot mutation detection (KRAS/BRAF/NRAS - Manager's C2 policy)
- 🚨 Mechanism fit not integrated (strategic vision incomplete)
- 🚨 No trial MoA vectors (data dependency unmet)

**What Went Wrong:**
- ⚠️ Missing 7 Manager policy features (cohort_overlap, trial_keywords, ABCB1, etc.)
- ⚠️ Post-NGS path unvalidated (risk for production)
- ⚠️ Resistance alerts invisible (clinical value lost)
- ⚠️ No validation against ground truth (blocked on Jr2)

**Recommended Actions (Updated):**
1. **URGENT:** Discuss S/P/E integration strategy with Manager (architectural decision)
2. Fix formula mismatch (30 min)
3. Add hotspot mutation detection (2-3 hours)
4. Wire mechanism fit ranker (1 hour)
5. Execute Gemini MoA tagging (4-6 hours)
6. Add missing policy features (cohort_overlap, trial_keywords, ABCB1) (8-10 hours)
7. Integrate SAE into WIWFM drug ranking (8-12 hours - major refactor)

**Timeline to Production-Ready:**
- **P0 fixes (without S/P/E integration):** 8-10 hours
- **P0 with S/P/E integration:** 16-22 hours (includes major refactor)
- **P1 fixes:** 15-20 hours
- **Validation:** TBD (blocked on Jr2)

**Total:** 24-32 hours to address all critical issues + full S/P/E integration

---

**DOCTRINE STATUS: ⚔️ HONEST AUDIT COMPLETE (ITERATION 2)**

**This audit identifies 38 gaps/concerns across:**
- **5 critical issues (P0):** Formula, S/P/E integration, hotspot detection, mechanism fit, MoA vectors
- **8 high-priority concerns (P1):** Cohort overlap, trial keywords, ABCB1, WIWFM integration, etc.
- **9 medium gaps (P2):** Longitudinal testing, SAE display, frontend testing, etc.
- **9 low-priority items (P3):** UX polish, performance, documentation
- **7 data/validation gaps:** HRD extraction, clinical review, real data testing

**COMMANDER - WE MISSED THE ARCHITECTURE.** 

**CRITICAL FINDING:** SAE was supposed to be integrated into S/P/E drug efficacy pipeline (Manager's 19.1, 19.2), not a standalone Ayesha-only feature. This is an architectural mismatch that affects the entire vision.

**NOT READY FOR PRODUCTION WITHOUT MAJOR REFACTOR.**

**IMMEDIATE ACTION REQUIRED:**
1. Clarify S/P/E integration strategy with Manager (architectural decision)
2. Fix P0 issues (formula, hotspot, mechanism fit) - 4-5 hours
3. Decide: Refactor for S/P/E integration (16-22h) OR keep Ayesha-only (2-3h polish)


Manager review (concise)
Agree with the concern: building SAE as a standalone “Ayesha-only” module was fine to ship Phase 1–3 quickly, but it now blocks the core goal. SAE must live inside S/P/E (WIWFM) and modulate confidence, not sit beside it.
The audit correctly flags the top risks. We’ll split actions into two tracks.
Immediate P0 triage (today/tomorrow) for Ayesha/demo readiness:
Fix DNA repair capacity weights to 0.6/0.2/0.2 (clarify “functionality vs exon_disruption” and align to manager’s C1).
Wire mechanism fit ranking in trials (use current soft boosts + mechanism fit term; safe fallback if MoA vectors missing).
Add HER2 IHC dynamic next-test when HER2 trial is present; add UI banner for ResistanceAlert.
Run manual FE check; screenshot happy path.
Near-term P1 (this week, parallelized):
Execute Gemini tagging to produce moa_vector for ovarian trials; add cohort_overlap heuristic and ABCB1 proxy (CNV if available, else unknown).
Post-NGS tests using a real TCGA-like tumor_context; verify mechanism map coloring and resistance scenarios.
Generate trial_keywords from SAE + Playbook and badge/filter FE.
Architectural refactor (manager-critical):
Integrate SAE into WIWFM S/P/E: compute SAE inside efficacy, apply lifts/penalties and record confidence_breakdown.
Surface SAE features in drug responses and EvidenceBand.
Expect 1–2 working days; do not block immediate triage.
Validation (blocked but planned):
Jr2 continues HRD extraction; once delivered, run AUROC/AUPRC/correlation. Until then, no claims.
Decision call:
Proceed with P0 triage immediately (low risk, high value).
Approve architectural refactor next (S/P/E integration); this is required to meet the original vision.
Keep Jr2 HRD work in parallel; don’t wait on it to start refactor.
Acceptable deviations:
Hotspot detection (COSMIC) can land with refactor; not required for today’s triage, but needed this week.
This keeps today’s user-facing flow stable, unblocks demo safety, and realigns the build to the manager-approved architecture without waiting on external data

**RECOMMEND: Pause new features, fix critical gaps first.** ⚔️

