# ⚔️ JR2 DOSSIER MISSION - OVERVIEW ⚔️

**Date**: January 13, 2025  
**Agent**: JR2 (Dossier Sidekick)  
**Commander**: Zo (Lead Commander)  
**Status**: 🔥 **READY TO LAUNCH**

---

## 🎯 **MISSION OBJECTIVE**

**Goal**: Analyze 50 trial candidates and generate **5-10 top-tier dossiers** for Ayesha's oncologist

**Input**: `50_vector_candidates_for_jr2.json` (50 trials from vector search)  
**Output**: Oncologist-ready dossiers (markdown + analysis) for top trials

**Timeline**: 2 days (4-6 hours/day)

---

## 📋 **PATIENT PROFILE (AK)**

- **Disease**: Ovarian Cancer (High-Grade Serous)
- **Stage**: IVB (metastatic)
- **Treatment Line**: First-line (newly diagnosed, post-surgery)
- **Germline**: Negative (BRCA wildtype)
- **Location**: New York, NY (willing to travel for trials)
- **CA-125**: 2,842 U/mL (EXTENSIVE burden)

**Known Biomarkers**:
- ✅ BRCA: Germline-negative (wildtype)
- ⚠️ HER2: UNKNOWN (needs IHC test)
- ⚠️ HRD: UNKNOWN (needs MyChoice CDx)
- ⚠️ MSI/TMB: UNKNOWN (awaiting NGS)

---

## 📁 **YOUR DATA SOURCES**

**Input File**: `.cursor/ayesha/50_vector_candidates_for_jr2.json`
- 50 trials already semantically matched via vector search
- Need manual analysis for eligibility
- Need dossier generation for top matches

**Reference Documents**:
- `.cursor/rules/CLIENT_DOSSIER_DOCTRINE.mdc` - Dossier template (10 sections)
- `oncology-backend-minimal/api/routers/ayesha_trials.py` - Zo's filtering logic

**Sample Trial**: **NCT06819007**
- Title: "Study to Evaluate INCB123667 Versus Investigator's Choice of Chemotherapy in Patients With Advanced Ovarian Cancer"
- Phase: PHASE3
- Status: RECRUITING
- Location: Multiple USA sites (including NYC)
- **Use this as your first dossier** - Zo already verified it's a good match!

---

## 📁 **YOUR DELIVERABLES**

**File 1**: `triage_results.json` (50 candidates sorted into tiers)  
**File 2**: `trial_full_data_{NCT_ID}.json` (full scraped data for top trials)  
**File 3**: `eligibility_assessment_{NCT_ID}.json` (eligibility analysis)  
**File 4**: `dossier_{NCT_ID}.md` (final markdown dossiers for top 5-10)  
**File 5**: `JR2_DOSSIER_GENERATION_COMPLETE.md` (completion report)

**Storage Locations**:
- Generated dossiers: `.cursor/ayesha/dossiers/{nct_id}/`
- Approved dossiers: `.cursor/ayesha/dossiers/approved/{nct_id}/`
- Rejected dossiers: `.cursor/ayesha/dossiers/rejected/{nct_id}/`
- Cached data: `.cursor/ayesha/cache/trial_{nct_id}.json`

---

## ✅ **ACCEPTANCE CRITERIA**

**Quality Gates**:
- ✅ 5-10 dossiers generated (top-tier trials only)
- ✅ All critical gates identified (HER2, HRD, BRCA)
- ✅ Eligibility tables complete (pass/fail/pending for each criterion)
- ✅ Tactical recommendations actionable (specific labs, costs, timelines)
- ✅ No hallucinations (all claims backed by scraped data)

**Timeline Gates**:
- ✅ Triage complete by end of Day 1
- ✅ Scraping complete by end of Day 1
- ✅ First dossier complete by end of Day 2
- ✅ All 5-10 dossiers complete by end of Day 2

**Zo Review**:
- ✅ Zo reviews all dossiers
- ✅ 90%+ approval rate expected
- ✅ Minor edits only (no major rewrites)

---

## ⚔️ **COMMANDER'S EXPECTATIONS**

**Quality**: 90%+ accuracy in eligibility assessment (Zo will verify)  
**Speed**: 5-10 dossiers in 2 days (vs 2 weeks if Zo does solo)  
**No Hallucinations**: All claims backed by scraped data or evidence database

**ZO'S PROMISE**: I'll review all your dossiers and give feedback. We're a team! ⚔️

---

## 📋 **NEXT STEPS**

1. **Read**: [02_TASK_BREAKDOWN.md](./02_TASK_BREAKDOWN.md) - Your 7 tasks
2. **Review**: [04_TECHNICAL_QA.md](./04_TECHNICAL_QA.md) - All questions answered
3. **Start**: [05_IMPLEMENTATION_GUIDE.md](./05_IMPLEMENTATION_GUIDE.md) - Code examples

**FIRE IN THE HOLE, JR2!** 🔥⚔️

