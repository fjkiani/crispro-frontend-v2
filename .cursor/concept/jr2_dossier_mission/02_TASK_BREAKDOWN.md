# ⚔️ JR2 TASK BREAKDOWN - THE 7 TASKS ⚔️

**Total Time**: ~12 hours (2 days, 4-6 hours/day)  
**Priority**: P0 (Critical) for Tasks 1, 2, 3, 6, 7 | P1 (Important) for Tasks 4, 5

---

## **TASK 1: TRIAGE 50 CANDIDATES (2 hours) - PRIORITY P0**

**Objective**: Sort 50 candidates into 3 tiers based on eligibility for Ayesha

**Triage Criteria**:

**TOP-TIER (Priority 1)** - Must meet ALL:
- ✅ Stage IV allowed
- ✅ First-line OR Maintenance (post-frontline)
- ✅ Recruiting OR Not Yet Recruiting (start date < 6 months)
- ✅ USA location

**GOOD-TIER (Priority 2)** - Must meet MOST:
- ✅ Stage III/IV allowed
- ⚠️ Maintenance OR Recurrent platinum-sensitive
- ✅ Recruiting (or recently completed with extension)
- ✅ USA or nearby (Canada/Mexico)

**OK-TIER (Priority 3)** - Interesting but not immediate:
- ⚠️ Stage III only
- ⚠️ Recurrent platinum-resistant
- ⚠️ Not yet recruiting (start date > 6 months)
- ⚠️ International only

**Output**: Create `triage_results.json`:
```json
{
  "top_tier": [list of NCT IDs],
  "good_tier": [list of NCT IDs],
  "ok_tier": [list of NCT IDs],
  "rejected": [list of NCT IDs with reasons]
}
```

**Reference**: See [06_FILTERING_LOGIC.md](./06_FILTERING_LOGIC.md) for Zo's filtering strategy

---

## **TASK 2: SCRAPE FULL TRIAL PAGES (3 hours) - PRIORITY P0**

**Objective**: Get complete eligibility criteria for top-tier trials (not truncated)

**For Each Top-Tier Trial**:
1. Visit ClinicalTrials.gov page (e.g., https://clinicaltrials.gov/study/NCT06819007)
2. Extract **Key Inclusion Criteria** (full text, not truncated)
3. Extract **Key Exclusion Criteria** (full text)
4. Extract **Interventions** (experimental arm vs control arm)
5. Extract **Primary Endpoint** (what trial is measuring)
6. Extract **Study Start Date** and **Primary Completion Date**
7. Extract **Locations** (all recruiting sites with contact info)

**Tools**: Use Diffbot (already integrated - see [05_IMPLEMENTATION_GUIDE.md](./05_IMPLEMENTATION_GUIDE.md) for code)

**Output**: Create `trial_full_data_{NCT_ID}.json` for each trial

**Caching**: Cache scraped data for 24 hours (see [05_IMPLEMENTATION_GUIDE.md](./05_IMPLEMENTATION_GUIDE.md))

---

## **TASK 3: ELIGIBILITY MATCHING (2 hours) - PRIORITY P0**

**Objective**: Compare Ayesha's profile to each trial's requirements

**For Each Top-Tier Trial**:
1. Parse inclusion criteria for requirements:
   - Disease: Stage III/IV? High-grade serous?
   - Treatment line: First-line, maintenance, recurrent?
   - Biomarkers: HER2? BRCA? HRD? MSI? TMB?
   - Performance status: ECOG 0-1?
   - Prior therapies: None (first-line) or specific (recurrent)?

2. Generate eligibility table:
   ```markdown
   | Criterion | Patient Status | Match | Action Required |
   |-----------|---------------|-------|-----------------|
   | Stage IV  | Stage IVB     | ✅ PASS | None |
   | HER2 IHC 1+ | UNKNOWN     | ⚠️ PENDING | Order HER2 IHC (3-5 days) |
   ```

3. Identify **CRITICAL GATES** (biomarkers that are UNKNOWN for Ayesha):
   - HER2 status: UNKNOWN → **CRITICAL GATE**
   - HRD score: UNKNOWN → **CRITICAL GATE**
   - BRCA: KNOWN (germline-negative) → **PASS**

**Output**: Create `eligibility_assessment_{NCT_ID}.json` for each trial

**Reference**: See `CLIENT_DOSSIER_DOCTRINE.mdc` lines 118-340 for matching functions

---

## **TASK 4: EXTRACT DRUG MECHANISMS (1 hour) - PRIORITY P1**

**Objective**: Explain what each trial drug does (layman + technical)

**For Each Drug in Trial**:
1. Extract drug name from interventions
2. Look up drug class in `DRUG_MECHANISM_DB` (see [05_IMPLEMENTATION_GUIDE.md](./05_IMPLEMENTATION_GUIDE.md))
3. Write layman explanation (1-2 sentences)
4. Write technical mechanism (3-5 bullets)

**Reference**: Use `CLIENT_DOSSIER_DOCTRINE.mdc` Section 3 as template

**Output**: Add `mechanism_explanation` to each trial's JSON

---

## **TASK 5: GENERATE STRATEGIC SCENARIOS (1 hour) - PRIORITY P1**

**Objective**: Calculate probability Ayesha is eligible and generate 3 scenarios

**For Each Trial**:
1. **Best-Case** (all critical gates pass):
   - Probability: Calculate based on biomarker prevalence
   - Outcome: "Patient enrolls, gets experimental arm, access to novel therapy"

2. **Most-Likely** (some gates pending):
   - Probability: Based on pending gates
   - Outcome: "Tests ordered, results in 7-10 days, enrollment pending"

3. **Challenge** (critical gate fails):
   - Probability: Based on failure rate
   - Outcome: "Not eligible, alternative trial recommended"

**Prevalence Data** (for probability calculation):
- HER2 IHC 1+/2+/3+ in ovarian: **40-60%** (assume 50%)
- BRCA wildtype (given germline-negative): **100%** (Ayesha confirmed)
- HRD-positive (given BRCA-negative): **40-50%** (assume 45%)

**Output**: Add `strategic_scenarios` to each trial's JSON

---

## **TASK 6: TACTICAL RECOMMENDATIONS (1 hour) - PRIORITY P0**

**Objective**: Generate action list for each trial (what Ayesha needs to do)

**For Each Trial with Critical Gates**:
1. Identify missing biomarkers (HER2, HRD, etc.)
2. Generate action list:
   ```markdown
   ### IMMEDIATE ACTIONS (THIS WEEK):
   
   **1. Order HER2 IHC** 🚨 **URGENT P0**
   - Source: Tumor biopsy tissue (already available from surgery)
   - Method: Immunohistochemistry (ASCO-CAP gastric scoring)
   - Turnaround: 3-5 days
   - Cost: $200-400 (covered by insurance)
   - Lab: Any pathology lab or trial site
   
   **2. Order HRD Test** ⏳ **P1**
   - Test: MyChoice CDx (Myriad Genetics)
   - Source: FFPE tissue block
   - Turnaround: 7-10 days
   - Cost: $4,000-6,000 (covered for Stage IV)
   ```

**Output**: Add `tactical_recommendations` to each trial's JSON

---

## **TASK 7: RENDER DOSSIERS (2 hours) - PRIORITY P0**

**Objective**: Generate markdown dossiers for top 5-10 trials

**For Each Top-Tier Trial**:
1. Use template from `CLIENT_DOSSIER_DOCTRINE.mdc`
2. Fill in all 10 sections:
   - Trial Intelligence Summary
   - Why This Matters for Ayesha
   - Clinical Mechanism
   - Eligibility Assessment Table
   - Critical Decision Tree
   - Strategic Implications
   - Tactical Recommendations
   - Clinical Evidence (if available)
   - Competitive Positioning (vs SOC)
   - Final Recommendation

**Output**: Create `dossier_{NCT_ID}.md` for each trial

**File Naming**: `dossier_{nct_id}_{patient_id}_{timestamp}.md`
- Example: `dossier_NCT06819007_ayesha_001_20250113T220000.md`

---

## 📋 **EXECUTION CHECKLIST**

**Day 1 (Tomorrow)**:
- [ ] Read `50_vector_candidates_for_jr2.json`
- [ ] Triage 50 trials into tiers (top/good/ok/rejected)
- [ ] Scrape top 10-15 trial pages (full eligibility criteria)
- [ ] Create eligibility assessments for top 10

**Day 2 (Day After)**:
- [ ] Generate mechanisms for top 10 drugs
- [ ] Generate strategic scenarios (probability calculations)
- [ ] Generate tactical recommendations (action lists)
- [ ] Render 5-10 markdown dossiers
- [ ] Submit to Zo for review

---

**Total Estimated Time**: 12 hours  
**Target**: 5-10 oncologist-ready dossiers by end of Day 2

