# ⚔️ AGENT JR - PRE-EXECUTION CLARIFICATIONS ⚔️

**Purpose**: Answer common questions before starting execution  
**Created By**: Zo  
**Date**: January 12, 2025

---

## ✅ **VERIFIED - READY TO PROCEED**

### **1. HybridTrialSearchService Structure** ✅
**Location**: `api/services/hybrid_trial_search.py`

**Method to Call**:
```python
from api.services.hybrid_trial_search import HybridTrialSearchService

service = HybridTrialSearchService()
trials = await service.search_optimized(
    query="ovarian cancer first-line",
    patient_context={
        "condition": "ovarian_cancer_high_grade_serous",
        "disease_category": "ovarian_cancer",
        "location_state": "NY"
    },
    germline_status="negative",  # ✅ Already supported!
    tumor_context={},  # Empty for now (NGS pending)
    top_k=50  # Get 50 candidates for filtering
)
```

**Returns**: `List[Dict[str, Any]]` - Each trial dict has:
- `nct_id` (str)
- `title` (str)
- `status` (str) - "Recruiting", "Active, not recruiting", etc.
- `phase` (str) - "Phase 1", "Phase 2", "Phase 3"
- `description` (str) - Full trial description
- `eligibility_text` (str) - Eligibility criteria
- `interventions` (List[str]) - Drug names
- `locations` (List[Dict]) - Site locations
- `optimization_score` (float) - Graph-based score
- `sporadic_filtering_applied` (bool) - True if germline filtering used

**✅ CONFIRMED**: Service already supports `germline_status="negative"` filtering!

---

### **2. Router Registration Pattern** ✅
**Location**: `api/main.py`

**Pattern** (see lines 32-44):
```python
from .routers import ayesha_trials as ayesha_trials_router

# Register router
app.include_router(ayesha_trials_router.router)
```

**✅ CONFIRMED**: Simple pattern - just import and include!

---

### **3. Frontend Routing Pattern** ✅
**Location**: `src/App.jsx`

**Pattern** (see lines 56, 139-152):
```jsx
import AyeshaTrialExplorer from './pages/AyeshaTrialExplorer';

// In Routes:
<Route path="/ayesha-trials" element={<AyeshaTrialExplorer />} />
```

**Sidebar Integration**: Add to `constants/index.js` navigation array

**✅ CONFIRMED**: Standard React Router pattern!

---

### **4. AstraDB Seeding** ✅ **FOUND**

**Status**: Seeding scripts exist!

**Available Scripts**:
1. `scripts/seed_astradb_from_sqlite.py` - Seeds from SQLite database
2. `scripts/seed_trials_standalone.py` - Standalone seeding
3. `scripts/seed_trials_simple.py` - Simple seeding

**Recommended Approach**:
1. **FIRST**: Check if AstraDB already has ≥200 ovarian trials
   ```python
   # Quick check via ClinicalTrialSearchService
   service = ClinicalTrialSearchService()
   results = await service.search_trials(
       query="ovarian cancer",
       disease_category="ovarian_cancer",
       top_k=200
   )
   # If len(results['data']['found_trials']) >= 200 → Skip seeding
   ```

2. **IF NEEDED**: Use existing seeding script
   ```bash
   cd oncology-coPilot/oncology-backend-minimal
   venv/bin/python scripts/seed_astradb_from_sqlite.py --disease ovarian --count 200
   ```

**✅ CONFIRMED**: Seeding scripts exist - check first, seed if needed!

---

### **5. Trial Data Structure** ✅

**From HybridTrialSearchService, each trial Dict has**:
```python
{
    "nct_id": "NCT12345678",
    "title": "Trial Title",
    "status": "Recruiting",
    "phase": "Phase 2",
    "description": "Full description text...",
    "eligibility_text": "Inclusion criteria...",
    "interventions": ["Drug1", "Drug2"],
    "locations": [
        {
            "facility": "Memorial Sloan Kettering",
            "city": "New York",
            "state": "NY",
            "zip": "10065"
        }
    ],
    "optimization_score": 0.85,
    "sporadic_filtering_applied": True
}
```

**✅ CONFIRMED**: Use `.get()` with defaults for all fields!

---

### **6. Frontend Component Styling** ✅

**Reference Components**:
- `src/components/research/ResultsDisplay.jsx` - Trial card patterns
- `src/components/sporadic/SporadicProvenanceCard.jsx` - Card styling
- `src/components/ayesha/DrugRankingPanel.jsx` - Ayesha-specific styling

**MUI Components to Use**:
- `Card`, `CardContent` - Container
- `Typography` - Text
- `Chip` - Badges/tags
- `LinearProgress` - Match score bar
- `Alert` - Warnings/errors
- `Box` - Layout

**✅ CONFIRMED**: Follow existing patterns!

---

## ✅ **ZO'S ANSWERS TO JR'S QUESTIONS**

### **Q1: AstraDB Seeding** ✅ ANSWERED
**Answer**: Use existing `scripts/seed_trials_simple.py`. Jr previously seeded 30 trials successfully. If needed, run:
```bash
cd oncology-coPilot/oncology-backend-minimal
venv/bin/python scripts/seed_trials_simple.py --disease ovarian --count 200
```
**BUT**: According to `ZO_ASTRADB_SEEDING_STATUS.md`, Jr already completed this! Check count first:
```python
# In your code:
results = await service.search_trials(query="ovarian cancer", top_k=200)
# If len(results) >= 200 → Skip seeding
```

### **Q2: Trial Contact Info** ✅ ANSWERED
**Answer**: **Leave blank for now**. ClinicalTrials.gov data structure doesn't always include contact_name/phone/email in structured fields. Instead:
1. Display `locations[*].facility` and `locations[*].city, state` clearly
2. Add "View on ClinicalTrials.gov" link using `nct_id`
3. In dossier export (Markdown), include note: "Contact info available on ClinicalTrials.gov: https://clinicaltrials.gov/study/{nct_id}"

### **Q3: Location Distance Calculation** ✅ ANSWERED
**Answer**: **Hardcode 'NYC metro' filter** for V1. Manager's decision (line 514 of plan):
- Distance filtering: `≤50 miles` from NYC
- Backend: Filter `locations` array where `state == "NY" OR state == "NJ" OR state == "CT"`
- Frontend: Display "📍 NYC Metro" badge for trials with NY/NJ/CT sites
- V2 (future): Add proper geopy distance calculation

### **Q4: Conditional NGS Features** ✅ ANSWERED
**Answer**: **Show with "Awaiting NGS" warning**. Manager's decision (lines 719-722):
- UI displays: Grayed WIWFM panel with "🔒 Unlock with NGS" banner
- Text: "Personalized drug predictions available when tumor NGS completes"
- NGS Fast-Track Checklist displayed prominently (ctDNA, HRD, IHC)
- Eligibility checklist STILL WORKS (uses clinical criteria only)
- CA-125 tracker STILL WORKS (uses CA-125 value only)

---

## 🔧 **ADDITIONAL CLARIFICATIONS FROM MANAGER**

### **A1: CA-125 Intelligence Thresholds** ✅
**From Manager** (lines 391-395):
```python
# Burden classes:
if ca125 < 100: burden_class = "MINIMAL"
elif ca125 < 500: burden_class = "MODERATE"
elif ca125 < 1000: burden_class = "SIGNIFICANT"
else: burden_class = "EXTENSIVE"  # Ayesha is here (2842)

# Forecast:
cycle3_expected_drop = "≥70%"  # Chemo-sensitive
cycle6_expected_drop = "≥90%"  # Complete response target
target = "<35 U/mL"  # CR threshold

# Resistance signal:
if on_therapy_rise OR drop_less_than_50_percent_by_cycle_3:
    flag = "⚠️ RESISTANCE SIGNAL - Consider alternative therapy"
```

### **A2: Confidence Gates Formula** ✅
**From Manager** (lines 427-434):
```python
# Formula: confidence = max(gates) with cap 1.0
# NOTE: This is DIFFERENT from efficacy confidence (which uses weighted S/P/E)
# Trial matching uses gate-based system (highest gate wins)
gates = []
if soc_aligned_nccn: gates.append(0.95)
if frontline_trial_eligibility >= 0.80: gates.append(0.90)
# NYC proximity and CA-125 monitoring are DISPLAY badges (+0.05 each), NOT stacked

confidence = min(max(gates), 1.0)  # Cap at 1.0
```

**⚠️ CRITICAL**: This is **trial matching confidence**, NOT drug efficacy confidence. Different formula for different use case!

### **A3: Hard/Soft Criteria Scoring** ✅
**From Manager** (lines 1049-1058):
```python
# Hard criteria (MUST pass): Stage, Treatment line, Major exclusions
# Soft criteria (% match): ECOG, Age, Distance, Biomarkers, Organ function

if all_hard_pass:
    if soft_percent >= 0.80: eligibility_gate = 0.90
    elif soft_percent >= 0.60: eligibility_gate = 0.85  # Yellow notice
    else: eligibility_gate = 0.75  # Yellow notice
else:
    # Any hard fail → trial EXCLUDED (red)
    eligibility_gate = 0.0
```

### **A4: Gemini for Eligibility Parsing** ✅
**From Manager** (lines 987-992):
```python
# Use Gemini (free tier) for eligibility parsing
# Caching: AstraDB - store in `structured_criteria` field
# Processing: Offline pre-processing (batch 200 trials → human review → cache)
# Runtime: NEVER call Gemini during trial search (serve cached only)
```

**⚠️ CRITICAL**: **DO NOT** call Gemini API at runtime! Only use cached `structured_criteria` from AstraDB. If missing, use text-based keyword matching instead.

---

## 🎯 **EXECUTION STRATEGY**

### **Start Here** (First 30 min):
1. ✅ **Verify AstraDB**: Check if ≥200 ovarian trials exist
2. ✅ **Read HybridTrialSearchService**: Understand exact return structure
3. ✅ **Read ResultsDisplay.jsx**: Understand frontend patterns
4. ✅ **Create first module**: Start with schemas (no dependencies)

### **If Blocked**:
- Create `AGENT_JR_TRIALS_QUESTIONS.md` immediately
- Don't wait - ask Zo right away
- Include: What you tried, error message, what you need

### **Progress Updates**:
- Every 2 hours: Update `AGENT_JR_TRIALS_PROGRESS.md`
- Mark steps as ✅/🔄/❌
- Note any deviations from plan

---

## ✅ **READY TO START**

**All Critical Questions Answered**: ✅  
**Code References Verified**: ✅  
**Patterns Confirmed**: ✅  
**AstraDB Seeding**: ⚠️ Check first, ask if needed

**Agent Jr - You're cleared to proceed!** ⚔️

**Start with**: Step 1 (Schemas) - 30 min, no dependencies!

---

**Last Updated**: January 12, 2025  
**By**: Zo

