# ⚔️ JR2 PREPARATION CHECKLIST - GET READY FOR ZO'S CANDIDATES ⚔️

**Date**: January 13, 2025 - 22:30 EST  
**Status**: 🔥 **PREPARING** - Zo is seeding, JR2 preparing pipeline  
**Expected Candidates**: 100 trials (tier-tagged) by midnight

---

## 🎯 **UNDERSTANDING ZO'S AUTONOMOUS WORK**

### **What Zo Is Doing RIGHT NOW**:
- ✅ **Iteration 1** (22:30-23:30): Seeding 755 → 1000 trials (SQLite + AstraDB)
- ⏳ **Iteration 2** (23:30-23:50): Uploading to AstraDB with vectors
- ⏳ **Iteration 3** (23:50-23:55): Vector search → 100 candidates
- ⏳ **Iteration 4** (23:55-00:00): Testing 5 filter strategies
- ⏳ **Iteration 5** (00:00-00:15): Quality analysis + tier tagging

### **What Zo Will Deliver**:
- ✅ `100_vector_candidates_for_jr2_FULL.json` - 100 trials with tier tags (TOP/GOOD/OK)
- ✅ `filtering_strategy_guide.md` - How to replicate Zo's filtering
- ✅ `quality_analysis_report.md` - Which trials to prioritize
- ✅ Updated `00_ZO_JR2_SYNC.json` - Status updates

### **Zo's Recommendation**: **Option 5 - Multi-Tier Strategy**
- **Top-Tier**: 10-15 trials (Stage IV, first-line, recruiting, USA)
- **Good-Tier**: 10-15 trials (maintenance, upcoming, platinum-sensitive)
- **OK-Tier**: 10-20 trials (interesting but conditional)
- **Total**: 30-50 trials for dossier generation

---

## 📋 **JR2 PREPARATION TASKS (DO NOW)**

### **✅ TASK 1: Review All Documentation** (30 min)
- [ ] Read [01_MISSION_OVERVIEW.md](./01_MISSION_OVERVIEW.md) - Understand mission
- [ ] Read [02_TASK_BREAKDOWN.md](./02_TASK_BREAKDOWN.md) - Your 7 tasks
- [ ] Read [04_TECHNICAL_QA.md](./04_TECHNICAL_QA.md) - All questions answered
- [ ] Read [05_IMPLEMENTATION_GUIDE.md](./05_IMPLEMENTATION_GUIDE.md) - Code examples
- [ ] Read [06_FILTERING_LOGIC.md](./06_FILTERING_LOGIC.md) - Replicate Zo's filters
- [ ] Read `CLIENT_DOSSIER_DOCTRINE.mdc` - Dossier template (10 sections)

**Status**: 🔄 **IN PROGRESS** - Reviewing now

---

### **✅ TASK 2: Set Up Project Structure** (15 min)
- [ ] Create folder: `oncology-backend-minimal/api/services/client_dossier/`
- [ ] Create folder: `oncology-backend-minimal/api/resources/drug_mechanism/`
- [ ] Create folder: `.cursor/ayesha/dossiers/` (for output)
- [ ] Create folder: `.cursor/ayesha/cache/` (for scraped data)
- [ ] Verify BeautifulSoup is installed: `pip install beautifulsoup4 lxml`

**Files to Create**:
- `api/services/client_dossier/trial_scraper.py` - BeautifulSoup scraper
- `api/services/client_dossier/eligibility_matcher.py` - Matching logic
- `api/services/client_dossier/dossier_generator.py` - Main orchestrator
- `api/resources/drug_mechanism/drug_mechanism_db.json` - 20 drugs

---

### **✅ TASK 3: Build Multi-Tier Filtering Logic** (45 min)
- [ ] Implement `filter_50_candidates()` from [06_FILTERING_LOGIC.md](./06_FILTERING_LOGIC.md)
- [ ] Add tier classification (TOP/GOOD/OK/REJECTED)
- [ ] Test with sample data (NCT06819007)
- [ ] Verify tier counts match Zo's expectations

**Expected Output**:
```json
{
  "top_tier": [10-15 trials],
  "good_tier": [10-15 trials],
  "ok_tier": [10-20 trials],
  "rejected": [remaining trials]
}
```

---

### **✅ TASK 4: Build Trial Scraper** (30 min)
- [ ] Implement `scrape_trial_page(nct_id)` using **Diffbot** (already integrated!)
- [ ] Use existing endpoint: `POST /api/evidence/extract` or call Diffbot API directly
- [ ] Extract: inclusion/exclusion, interventions, dates, locations from Diffbot HTML
- [ ] Add caching (24hr TTL) - see [05_IMPLEMENTATION_GUIDE.md](./05_IMPLEMENTATION_GUIDE.md)
- [ ] Test with NCT06819007 (Zo's verified trial)

**Existing Diffbot Integration**:
- File: `api/routers/evidence/extraction.py`
- Endpoint: `POST /api/evidence/extract` with `{"url": "https://clinicaltrials.gov/study/NCT06819007"}`
- Config: `DIFFBOT_TOKEN` from environment

**Test Command**:
```python
from api.services.client_dossier.trial_scraper import scrape_trial_page
data = await scrape_trial_page("NCT06819007")
assert 'inclusion_criteria_full' in data
assert 'study_start_date' in data
assert 'full_html' in data  # From Diffbot
```

---

### **✅ TASK 5: Build Eligibility Matcher** (45 min)
- [ ] Implement `assess_disease_match()` - Stage matching
- [ ] Implement `assess_treatment_line_match()` - Treatment line
- [ ] Implement `assess_biomarker_match()` - HER2, BRCA, HRD gates
- [ ] Generate eligibility table (pass/fail/pending)

**Reference**: `CLIENT_DOSSIER_DOCTRINE.mdc` lines 118-340

**Test**: Compare Ayesha's profile to NCT06819007 requirements

---

### **✅ TASK 6: Create Drug Mechanism Database** (30 min)
- [ ] Create `drug_mechanism_db.json` with 20 drugs
- [ ] Structure: `{drug_name: {class, target, moa, layman, breakthrough}}`
- [ ] Include: T-DXd, Olaparib, Bevacizumab, Pembrolizumab, etc.
- [ ] Add fallback logic (3-tier: DB → Evidence Service → Flag for Zo)

**Reference**: [05_IMPLEMENTATION_GUIDE.md](./05_IMPLEMENTATION_GUIDE.md) for complete list

---

### **✅ TASK 7: Prepare for Tier-Based Generation** (15 min)
- [ ] Understand tier priorities:
  - **Top-Tier**: Generate ALL dossiers (10-15 trials)
  - **Good-Tier**: Generate top 5-10 dossiers (highest match scores)
  - **OK-Tier**: Generate only if Commander requests (conditional)
- [ ] Set up dossier output structure:
  - `dossiers/top_tier/{nct_id}/`
  - `dossiers/good_tier/{nct_id}/`
  - `dossiers/ok_tier/{nct_id}/`

---

## 🎯 **READY FOR ZO'S DELIVERABLES**

### **What to Check at Midnight (00:00 EST)**:
1. ✅ `100_vector_candidates_for_jr2_FULL.json` - 100 trials with tier tags
2. ✅ `01_ZO_ITERATION_LOG.md` - Updated with search results
3. ✅ `00_ZO_JR2_SYNC.json` - Status updates
4. ✅ `filtering_strategy_guide.md` - Zo's filtering recommendations
5. ✅ `quality_analysis_report.md` - Which trials to prioritize

### **What to Do When Candidates Arrive**:
1. ✅ Load `100_vector_candidates_for_jr2_FULL.json`
2. ✅ Run multi-tier filtering (verify Zo's tier tags)
3. ✅ Prioritize top-tier trials (10-15 dossiers)
4. ✅ Start scraping top-tier trials (full eligibility criteria)
5. ✅ Generate first dossier (NCT06819007 as test)

---

## 📊 **EXPECTED WORKFLOW (AFTER MIDNIGHT)**

### **Phase 1: Filter & Triage** (1 hour)
- Load 100 candidates
- Run multi-tier filtering
- Verify tier classification
- Select top 10-15 for dossier generation

### **Phase 2: Scrape & Assess** (2 hours)
- Scrape top 10-15 trial pages
- Generate eligibility assessments
- Identify critical gates (HER2, HRD, BRCA)
- Cache all scraped data

### **Phase 3: Generate Dossiers** (3 hours)
- Generate 10-15 top-tier dossiers
- Use CLIENT_DOSSIER_DOCTRINE.mdc template
- All 10 sections per dossier
- Submit to Zo for review

**Total Time**: ~6 hours (can be done in 2 days)

---

## 🔥 **IMMEDIATE ACTIONS (DO NOW)**

### **Priority P0 (Before Midnight)**:
1. ✅ Review all modularized documents
2. ✅ Set up project structure
3. ✅ Build trial scraper (test with NCT06819007)
4. ✅ Create drug mechanism database (20 drugs)
5. ✅ Build filtering logic (multi-tier)

### **Priority P1 (After Candidates Arrive)**:
1. ⏳ Load and filter 100 candidates
2. ⏳ Scrape top-tier trials
3. ⏳ Generate eligibility assessments
4. ⏳ Generate first 5 dossiers

---

## ⚔️ **ZO'S NOTES FOR JR2**

**While I'm Seeding** (22:30-00:30):
- ✅ Build your folder structure
- ✅ Review CLIENT_DOSSIER_DOCTRINE.mdc
- ✅ Prepare Diffbot scraper (already integrated!)
- ✅ Set up drug mechanism database
- ✅ Build filtering logic

**When I'm Done** (by midnight):
- ✅ Check `100_vector_candidates_for_jr2_FULL.json`
- ✅ Check iteration log (see filtering strategy)
- ✅ Start your pipeline (filter → scrape → assess → generate)

**We're a Team**: I find gold, you refine it into diamonds. **LET'S GO!** ⚔️

---

## 📋 **PREPARATION STATUS**

**Completed**:
- ✅ Modularized documentation (10 files created)
- ✅ Master index with all references
- ✅ Technical Q&A (all questions answered)
- ✅ Implementation guide (code examples)

**In Progress**:
- 🔄 Reviewing Zo's strategic options
- 🔄 Understanding multi-tier strategy
- 🔄 Preparing filtering logic

**Next Steps**:
- ⏳ Set up project structure
- ⏳ Build trial scraper
- ⏳ Create drug mechanism database
- ⏳ Wait for Zo's candidates (midnight)

---

**STATUS**: 🔥 **PREPARING** - Ready to consume Zo's candidates when they arrive!

**Last Updated**: January 13, 2025 - 22:45 EST  
**Next Check**: Midnight (when Zo exports candidates)

