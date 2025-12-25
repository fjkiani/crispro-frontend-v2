# Frontend Mechanism Fit Display Audit

**Date:** January 28, 2025  
**Purpose:** Audit frontend components for mechanism fit display capabilities  
**Focus:** Gap 2 from DELIVERABLES_GAPS_ANALYSIS.md - Mechanism Fit Validation

---

## 🎯 Executive Summary

**Status:** ⚠️ **PARTIAL IMPLEMENTATION** - Some components support mechanism fit, but not all features documented

**What Exists:**
- ✅ `ClinicalTrialMatchingSection.jsx` - **FULL SUPPORT** (combined_score, mechanism_fit_score, sorting)
- ✅ `TrialMatchesCard.jsx` - **PARTIAL SUPPORT** (shows mechanism_fit_score, but no breakdown)
- ⚠️ `TrialMatchCard.jsx` - **NO SUPPORT** (only shows match_score, no mechanism fit)
- ⚠️ `AyeshaTrialExplorer.jsx` - **NO SUPPORT** (uses TrialMatchCard which doesn't show mechanism fit)

**What's Missing:**
- ❌ Mechanism alignment breakdown display (per-pathway alignment)
- ❌ "Low mechanism fit" warning badge
- ❌ Mechanism-aligned badge/chip
- ❌ Mechanism vector visualization
- ❌ Full integration in all trial display components

---

## 📊 Component-by-Component Audit

### ✅ 1. ClinicalTrialMatchingSection.jsx (Clinical Dossier)

**Location:** `oncology-coPilot/oncology-frontend/src/components/ClinicalDossier/components/ClinicalTrialMatchingSection.jsx`

**Status:** ✅ **FULL SUPPORT** - Best implementation

**What It Supports:**
- ✅ Combined score display (0.7×eligibility + 0.3×mechanism_fit)
- ✅ Mechanism fit score display
- ✅ Eligibility score display
- ✅ Sorting by combined_score, eligibility, or mechanism_fit
- ✅ Filtering by phase
- ✅ Score color coding (success/warning/default)
- ✅ Formula explanation: "0.7 × Eligibility + 0.3 × Mechanism Fit"

**What's Missing:**
- ❌ Mechanism alignment breakdown (per-pathway: DDR, MAPK, PI3K, etc.)
- ❌ "Low mechanism fit" warning badge (when mechanism_fit < 0.50)
- ❌ "Mechanism-aligned" badge (when mechanism_fit ≥ 0.50)
- ❌ Mechanism vector visualization
- ❌ Tooltip explaining mechanism fit meaning

**Code Evidence:**
```jsx
// Lines 45-56: Combined score calculation
const calculateCombinedScore = (trial) => {
  if (trial.combined_score !== undefined) {
    return trial.combined_score;
  }
  const eligibility = trial.eligibility_score || trial.match_score || 0;
  const mechanismFit = trial.mechanism_fit_score || 0;
  return 0.7 * eligibility + 0.3 * mechanismFit;
};

// Lines 286-303: Mechanism fit display
<Grid item xs={12} sm={4}>
  <Box>
    <Typography variant="caption" color="text.secondary">
      Mechanism Fit
    </Typography>
    <Typography variant="body2" sx={{ fontWeight: 600 }}>
      {(mechanismFitScore * 100).toFixed(0)}%
    </Typography>
    <LinearProgress
      variant="determinate"
      value={mechanismFitScore * 100}
      color={getScoreColor(mechanismFitScore)}
    />
  </Box>
</Grid>
```

**Recommendation:** ✅ **KEEP** - This is the best implementation. Add missing features (breakdown, badges).

---

### ⚠️ 2. TrialMatchesCard.jsx (Orchestrator)

**Location:** `oncology-coPilot/oncology-frontend/src/components/orchestrator/Analysis/TrialMatchesCard.jsx`

**Status:** ⚠️ **PARTIAL SUPPORT** - Shows mechanism_fit_score but no breakdown

**What It Supports:**
- ✅ Mechanism fit score display (if present)
- ✅ Combined score display (if present)
- ✅ Score color coding

**What's Missing:**
- ❌ Mechanism alignment breakdown
- ❌ "Low mechanism fit" warning
- ❌ "Mechanism-aligned" badge
- ❌ Formula explanation
- ❌ Mechanism vector visualization

**Code Evidence:**
```jsx
// Lines 161-175: Mechanism fit display (conditional)
{trial.mechanism_fit_score !== undefined && (
  <Box>
    <Typography variant="caption">Mechanism Fit</Typography>
    <LinearProgress
      variant="determinate"
      value={trial.mechanism_fit_score * 100}
      color={getScoreColor(trial.mechanism_fit_score)}
    />
    <Typography variant="caption">
      {(trial.mechanism_fit_score * 100).toFixed(0)}%
    </Typography>
  </Box>
)}
```

**Recommendation:** ⚠️ **ENHANCE** - Add breakdown, badges, and formula explanation.

---

### ❌ 3. TrialMatchCard.jsx (Ayesha Trial Explorer)

**Location:** `oncology-coPilot/oncology-frontend/src/components/trials/TrialMatchCard.jsx`

**Status:** ❌ **NO SUPPORT** - Only shows match_score, no mechanism fit

**What It Supports:**
- ✅ Match score display
- ✅ Reasoning display (why_eligible, why_good_fit)
- ✅ Eligibility checklist
- ✅ Evidence tier and enrollment likelihood

**What's Missing:**
- ❌ Combined score display
- ❌ Mechanism fit score display
- ❌ Mechanism alignment breakdown
- ❌ "Low mechanism fit" warning
- ❌ "Mechanism-aligned" badge
- ❌ Formula explanation

**Code Evidence:**
```jsx
// Lines 36-49: Only match_score, no mechanism fit
const {
  nct_id,
  title,
  phase,
  status,
  interventions,
  locations,
  match_score,  // ❌ Only match_score, no mechanism_fit_score
  reasoning,
  source_url,
} = trial;
```

**Recommendation:** ❌ **NEEDS UPDATE** - Add mechanism fit support to match ClinicalTrialMatchingSection.jsx.

---

### ❌ 4. AyeshaTrialExplorer.jsx (Page)

**Location:** `oncology-coPilot/oncology-frontend/src/pages/AyeshaTrialExplorer.jsx`

**Status:** ❌ **NO SUPPORT** - Uses TrialMatchCard which doesn't support mechanism fit

**What It Does:**
- ✅ Loads trials from `/api/ayesha/complete_care_v2`
- ✅ Displays trials using `TrialMatchCard` component
- ✅ Shows other sections (CA-125, SOC, Next Test, etc.)

**What's Missing:**
- ❌ Mechanism fit display (inherits from TrialMatchCard limitation)
- ❌ Mechanism alignment breakdown
- ❌ "Low mechanism fit" warning
- ❌ "Mechanism-aligned" badge

**Code Evidence:**
```jsx
// Lines 193-199: Uses TrialMatchCard (no mechanism fit support)
{trials.map((trial, index) => (
  <TrialMatchCard
    key={trial.nct_id || index}
    trial={trial}
    rank={index + 1}
  />
))}
```

**Recommendation:** ❌ **NEEDS UPDATE** - Either update TrialMatchCard or switch to ClinicalTrialMatchingSection.

---

### ⚠️ 5. TrialsListCard.jsx (Clinical Genomics Command Center)

**Location:** `oncology-coPilot/oncology-frontend/src/components/ClinicalGenomicsCommandCenter/cards/TrialsListCard.jsx`

**Status:** ⚠️ **UNKNOWN** - Need to check implementation

**Recommendation:** ⚠️ **AUDIT NEEDED** - Check if it supports mechanism fit.

---

## 📋 Missing Features (From Documentation)

### From ADVANCED_CARE_PLAN_MECHANISM_TRIAL_MATCHING.md:

**Expected Features:**
1. ✅ Combined score display (0.7×eligibility + 0.3×mechanism_fit) - **PARTIAL** (ClinicalTrialMatchingSection has it)
2. ❌ Mechanism alignment breakdown (per-pathway: DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux) - **MISSING**
3. ❌ "Low mechanism fit" warning badge (when mechanism_fit < 0.50) - **MISSING**
4. ❌ "Mechanism-aligned" badge (when mechanism_fit ≥ 0.50) - **MISSING**
5. ❌ Mechanism vector visualization - **MISSING**
6. ❌ Tooltip explaining mechanism fit meaning - **MISSING**

### From mechanism_trial_matching_contribution.mdc:

**Expected Display:**
- ✅ Combined score (0.7×eligibility + 0.3×mechanism_fit) - **PARTIAL**
- ❌ Mechanism fit score (0.92 avg for DDR-high patients) - **PARTIAL** (score shown, but no context)
- ❌ Mechanism alignment breakdown - **MISSING**
- ❌ Shortlist compression indicator (50+ → 5-12 trials) - **MISSING**

---

## 🎯 Gap Analysis: What Needs to Be Built

### Priority 1: Enhance Existing Components (2-3 hours)

**1.1: Update TrialMatchCard.jsx** ⚠️ **HIGH PRIORITY**
- Add `mechanism_fit_score` display
- Add `combined_score` display
- Add mechanism alignment breakdown (if `trial.mechanism_alignment` exists)
- Add "Low mechanism fit" warning badge
- Add "Mechanism-aligned" badge
- Add formula explanation tooltip

**1.2: Enhance TrialMatchesCard.jsx** ⚠️ **MEDIUM PRIORITY**
- Add mechanism alignment breakdown
- Add "Low mechanism fit" warning
- Add "Mechanism-aligned" badge
- Add formula explanation

**1.3: Enhance ClinicalTrialMatchingSection.jsx** ⚠️ **MEDIUM PRIORITY**
- Add mechanism alignment breakdown display
- Add "Low mechanism fit" warning badge
- Add "Mechanism-aligned" badge
- Add mechanism vector visualization
- Add tooltip explaining mechanism fit

### Priority 2: Create New Components (3-4 hours)

**2.1: MechanismAlignmentBreakdown.jsx** ⚠️ **HIGH PRIORITY**
- Display per-pathway alignment (DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux)
- Show pathway scores with color coding
- Tooltip explaining each pathway

**2.2: MechanismFitBadge.jsx** ⚠️ **MEDIUM PRIORITY**
- "Mechanism-aligned" badge (when mechanism_fit ≥ 0.50)
- "Low mechanism fit" warning badge (when mechanism_fit < 0.50)
- Color coding (green/yellow/red)

**2.3: MechanismVectorVisualization.jsx** ⚠️ **LOW PRIORITY**
- Visual representation of 7D mechanism vector
- Bar chart or radar chart
- Pathway labels

---

## 📊 Component Usage Map

| Component | Used In | Mechanism Fit Support | Status |
|-----------|---------|----------------------|--------|
| `ClinicalTrialMatchingSection.jsx` | Clinical Dossier | ✅ Full | ✅ **BEST** |
| `TrialMatchesCard.jsx` | Orchestrator | ⚠️ Partial | ⚠️ **ENHANCE** |
| `TrialMatchCard.jsx` | AyeshaTrialExplorer | ❌ None | ❌ **UPDATE** |
| `TrialsListCard.jsx` | Clinical Genomics Command Center | ⚠️ Unknown | ⚠️ **AUDIT** |

---

## 🎯 Recommended Action Plan

### Week 1: Enhance Existing Components (2-3 hours)

**Day 1: Update TrialMatchCard.jsx** (1-2 hours)
- Add mechanism_fit_score and combined_score display
- Add mechanism alignment breakdown (if available)
- Add "Low mechanism fit" warning badge
- Add "Mechanism-aligned" badge
- Test with AyeshaTrialExplorer

**Day 2: Enhance TrialMatchesCard.jsx** (1 hour)
- Add mechanism alignment breakdown
- Add badges
- Add formula explanation

**Day 3: Enhance ClinicalTrialMatchingSection.jsx** (1 hour)
- Add mechanism alignment breakdown
- Add badges
- Add tooltip

### Week 2: Create New Components (3-4 hours)

**Day 1-2: Create MechanismAlignmentBreakdown.jsx** (2 hours)
- Display per-pathway alignment
- Color coding
- Tooltips

**Day 3: Create MechanismFitBadge.jsx** (1 hour)
- Badge components
- Color coding

**Day 4: Create MechanismVectorVisualization.jsx** (1 hour)
- Visual representation
- Chart component

---

## ✅ Success Criteria

**Component Complete When:**
1. ✅ All trial display components show mechanism_fit_score
2. ✅ All trial display components show combined_score
3. ✅ Mechanism alignment breakdown displayed (when available)
4. ✅ "Low mechanism fit" warning shown (when mechanism_fit < 0.50)
5. ✅ "Mechanism-aligned" badge shown (when mechanism_fit ≥ 0.50)
6. ✅ Formula explanation visible (0.7×eligibility + 0.3×mechanism_fit)
7. ✅ All components tested with real API responses

**Demo Ready When:**
1. ✅ Navigate to AyeshaTrialExplorer
2. ✅ See mechanism fit scores in trial cards
3. ✅ See mechanism alignment breakdown
4. ✅ See "Mechanism-aligned" badges on high-fit trials
5. ✅ See "Low mechanism fit" warnings on low-fit trials
6. ✅ Sort by mechanism fit works

---

## 📝 Notes

**Backend Response Format:**
Based on testing, backend returns:
```json
{
  "trials": [
    {
      "nct_id": "...",
      "title": "...",
      "mechanism_fit_score": 0.92,
      "combined_score": 0.87,
      "mechanism_alignment": {
        "DDR": 0.95,
        "MAPK": 0.10
      },
      "low_mechanism_fit_warning": false,
      "mechanism_boost_applied": true
    }
  ]
}
```

**Frontend Should Display:**
- ✅ `mechanism_fit_score` - Score (0-1)
- ✅ `combined_score` - Combined score (0-1)
- ✅ `mechanism_alignment` - Per-pathway breakdown (object)
- ✅ `low_mechanism_fit_warning` - Boolean flag
- ✅ `mechanism_boost_applied` - Boolean flag

---

*Document Author: Zo*  
*Last Updated: January 28, 2025*  
*Status: ⚠️ PARTIAL IMPLEMENTATION - Needs Enhancement*


