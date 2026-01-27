# 🎨 AYESHA FRONTEND PRODUCTION AUDIT - Component Status & Deliverables

**Last Updated:** January 26, 2026  
**Synced With:** `AYESHA_ENGINE_DEEP_AUDIT.md`  
**Purpose:** Track frontend component status and define what needs to be wired

---

## 🚨 CRITICAL FINDING: BACKEND CONFIG ISSUE

**From Engine Audit:** MBD4 is NOT in the DDR gene lists in `ddr_config.py`!

**Impact on Frontend:** If we wire DDR components without fixing config:
- `DDRStatusCard` will show "Unknown" instead of "DDR Defective"
- `DDRTreatmentEligibility` will show "ELIGIBILITY UNKNOWN"

**Solution:** Either fix backend config first, OR use Synthetic Lethality result (which already detects MBD4+TP53)

---

## 📋 EXECUTIVE SUMMARY

### Component Status Overview

| Status | Count | Description |
|--------|-------|-------------|
| ✅ Working + Displayed | 12 | On Ayesha's page, showing data |
| ⚠️ Working + Awaiting Data | 4 | Built but needs user input or backend fix |
| ❌ Built + Wrong Location | 6 | Exist on separate pages |
| 🔴 Missing | 2 | Need to be created |

---

## 🎯 AYESHA'S MAIN PAGE: `AyeshaTrialExplorer.jsx`

**Location:** `oncology-frontend/src/pages/AyeshaTrialExplorer.jsx` (852 lines)

### Tab Structure

| Tab | Name | Content Status |
|-----|------|----------------|
| 0 | Overview | ✅ Mechanism + SOC + Hints (missing DDR status) |
| 1 | Trials | ✅ Working - Holistic score + PGx |
| 2 | Treatment | ⚠️ SOC + WIWFM ("awaiting NGS") + Timing ("treatment-naive") |
| 3 | Monitoring | 🔴 CA-125 blocked (null value) |
| 4 | Resistance | ✅ Playbook + Prophet |
| 5 | Synthetic Lethality | ✅ MBD4+TP53 detected ← **BEST SOURCE FOR PARP ELIGIBILITY** |

---

## ✅ COMPONENTS WORKING + DISPLAYED

### Tab 0: Overview

| Component | File | Status | Notes |
|-----------|------|--------|-------|
| `MechanismVectorVisualization` | `ayesha/MechanismVectorVisualization.jsx` | ✅ Working | Shows DDR=0.88 |
| `PathwayDisruptionCard` | `ayesha/PathwayDisruptionCard.jsx` | ✅ Working | From SL result |
| `EssentialPathwaysCard` | `ayesha/EssentialPathwaysCard.jsx` | ✅ Working | From SL result |
| `SOCRecommendationCard` | `ayesha/SOCRecommendationCard.jsx` | ✅ Working | Carboplatin+Paclitaxel+Bevacizumab |
| `NextTestCard` | `ayesha/NextTestCard.jsx` | ✅ Working | |
| `HintTilesPanel` | `ayesha/HintTilesPanel.jsx` | ✅ Working | |
| `MechanismChips` | `ayesha/MechanismChips.jsx` | ✅ Working | |

### Tab 5: Synthetic Lethality ← **KEY TAB**

| Component | File | Status | Notes |
|-----------|------|--------|-------|
| `SyntheticLethalityCard` | `orchestrator/Analysis/SyntheticLethalityCard.jsx` | ✅ Working | **Detects MBD4+TP53!** |
| `SLOpportunityBanner` | `ayesha/SLOpportunityBanner.jsx` | ✅ Working | Shows PARP recommended |
| `SLDrugRecommendations` | `ayesha/SLDrugRecommendations.jsx` | ✅ Working | |

---

## ❌ COMPONENTS ON WRONG PAGE (DDR)

### Located on `/ddr-status` page only

| Component | File | What It Shows | Will Show For Ayesha |
|-----------|------|---------------|---------------------|
| `DDRStatusCard` | `ddr/DDRStatusCard.jsx` | DDR_bin status | ⚠️ "Unknown" (MBD4 not in config) |
| `DDRTreatmentEligibility` | `ddr/DDRTreatmentEligibility.jsx` | PARP eligible badge | ⚠️ "Unknown" (MBD4 not in config) |
| `DDRFeatureBreakdown` | `ddr/DDRFeatureBreakdown.jsx` | Feature breakdown | Shows nothing |
| `DDRRecommendationsPanel` | `ddr/DDRRecommendationsPanel.jsx` | Recommendations | Shows nothing |

**Export Verified:** `components/ddr/index.js` exports all components correctly

---

## 🎯 REVISED DELIVERABLES

### D0: Fix MBD4 in Backend Config (30 min) - **BEFORE FRONTEND**

**Why:** Without this, DDR components show "Unknown"

**File:** `api/services/resistance/config/ddr_config.py`

**Change:** Add MBD4 to extended_ddr_genes for ovary:
```python
"extended_ddr_genes": [
    "ATM", "ATR", "CHEK1", "CHEK2",
    "FANCA", "FANCD2", "RAD50", "MRE11", "NBN", "POLQ",
    "MBD4"  # Add this
],
```

**Alternative:** Skip DDR components, use SL banner instead (already works!)

---

### D1: Add DDR Components to Tab 0 (4-6 hours)

**Prerequisites:** D0 complete OR decision to use SL as primary source

**Implementation (if using DDR endpoint):**

```jsx
// In AyeshaTrialExplorer.jsx

// Add imports
import { DDRStatusCard, DDRTreatmentEligibility } from '../components/ddr';
import { useDDRStatus } from '../hooks/useDDRStatus';

// Inside component, after other hooks
const { ddrStatus, loading: ddrLoading, calculateDDRStatus } = useDDRStatus();

// Add useEffect to compute DDR on load
useEffect(() => {
  const mutations = (AYESHA_11_17_25_PROFILE.germline?.mutations || []).map(m => ({
    gene_symbol: m.gene,
    variant_classification: m.classification || 'pathogenic'
  }));

  if (mutations.length > 0 && !ddrStatus && !ddrLoading) {
    calculateDDRStatus({
      patient_id: AYESHA_11_17_25_PROFILE.patient?.patient_id || 'AK',
      disease_site: 'ovary',
      tumor_subtype: 'HGSOC',
      mutations
    });
  }
}, []);

// In Tab 0 (Overview), add after Mechanism Intelligence section:
<Box mb={3}>
  <Typography variant="h5" gutterBottom>🧬 DDR Status</Typography>
  <Grid container spacing={3}>
    <Grid item xs={12} md={6}>
      <DDRStatusCard ddrStatus={ddrStatus} />
    </Grid>
    <Grid item xs={12} md={6}>
      <DDRTreatmentEligibility 
        ddrStatus={ddrStatus} 
        onViewTrials={() => setActiveTab(1)} 
      />
    </Grid>
  </Grid>
</Box>
```

**Alternative (if NOT fixing config):**
Use SL banner which already shows PARP eligible from MBD4+TP53 detection.

---

### D1-ALT: Promote SL Banner to Tab 0 (2 hours)

**Why:** SL already detects MBD4+TP53, no backend change needed

**Implementation:**
```jsx
// SLOpportunityBanner already exists and is used on page
// Just move it from just above tabs to Tab 0 Overview section
// OR make it more prominent

// In Tab 0 (Overview), add prominently:
{slResult?.synthetic_lethality_detected && (
  <Box mb={3}>
    <SLOpportunityBanner
      slDetected={slResult.synthetic_lethality_detected}
      suggestedTherapy={slResult.suggested_therapy}
      doubleHitDescription={slResult.double_hit_description}
      confidence={slResult.recommended_drugs?.[0]?.confidence}
      onViewDetails={() => setActiveTab(5)}
    />
  </Box>
)}
```

---

### D2: Create Treatment Options Summary (6-8 hours)

**Data Sources Available:**
- `socRecommendation` - From API ✅
- `trials` - Array from API ✅
- `slResult.recommended_drugs` - Drug recommendations ✅
- `wiwfm` - Shows "awaiting_ngs" ⚠️

---

### D3: Create CA-125 Entry Form (4 hours)

**Why:** CA-125 Intelligence blocked by null value

**CA125Tracker Props (verified):**
```jsx
<CA125Tracker
  current_value={number | null}
  burden_class={string}
  forecast={object}
  resistance_rule={object}
  monitoring_strategy={object}
/>
```

---

## 📊 HOOKS AVAILABLE

| Hook | File | Endpoint | Status |
|------|------|----------|--------|
| `useDDRStatus` | `hooks/useDDRStatus.js` | `/api/resistance/ddr-status` | ⚠️ Returns "unknown" for MBD4 |
| `useSyntheticLethality` | `hooks/useSyntheticLethality.js` | SL endpoint | ✅ Detects MBD4+TP53 |
| `useTimingChemoFeatures` | `hooks/useTimingChemoFeatures.js` | `/api/resistance/timing-chemo-features` | ✅ (treatment-naive shows alert) |

---

## 📁 COMPONENT IMPORT PATHS (Verified)

```javascript
// DDR Components
import { 
  DDRStatusCard, 
  DDRTreatmentEligibility,
  DDRFeatureBreakdown,
  DDRRecommendationsPanel,
  HRDPanel,
  DDRMutationSummary
} from '../components/ddr';  // ✅ Exports verified

// Existing Ayesha Components
import {
  MechanismVectorVisualization,
  SLOpportunityBanner,
  PathwayDisruptionCard,
  EssentialPathwaysCard,
  SLDrugRecommendations,
} from '../components/ayesha';  // ✅ Already imported

// Hooks
import { useDDRStatus } from '../hooks/useDDRStatus';  // ✅ Verified
import { useSyntheticLethality } from '../hooks/useSyntheticLethality';  // ✅ Already used
```

---

## 📊 EFFORT SUMMARY - REVISED

| Deliverable | Hours | Priority | Blocked By |
|-------------|-------|----------|------------|
| D0: Fix MBD4 in config | 0.5h | **P0** | Nothing |
| D1: DDR on Ayesha's page | 4-6h | **P0** | D0 |
| D1-ALT: Promote SL Banner | 2h | **P0** | Nothing |
| D2: Treatment Options | 6-8h | **P1** | Nothing |
| D3: CA-125 Entry | 4h | **P1** | Nothing |
| **TOTAL** | **16-20h** | | |

---

## ✅ DECISION MATRIX

### Option A: Full DDR Integration
1. Fix MBD4 in ddr_config.py (30 min)
2. Wire DDR components to Tab 0 (4h)
3. Test end-to-end (1h)
**Total: 5.5 hours**
**Result: DDRStatusCard shows "DDR Defective", PARP Eligible badge**

### Option B: Use SL as Primary
1. Promote SLOpportunityBanner to Tab 0 (2h)
2. Add "PARP ELIGIBLE" messaging based on slResult (1h)
**Total: 3 hours**
**Result: SL banner shows PARP recommended based on MBD4+TP53**

### Recommendation: **Option A** (Fix config + DDR components)
- More comprehensive (DDR score, features, recommendations)
- Aligns with clinical expectations (DDR deficiency status)
- Option B is good fallback if time-constrained

---

## ✅ WHAT'S ALREADY EXCELLENT

1. **Synthetic Lethality Tab** - Detects MBD4+TP53 ✅
2. **Trials Display** - Holistic scores, PGx safety ✅
3. **SOC Recommendation** - NCCN-aligned ✅
4. **Mechanism Vector** - Shows DDR=0.88 ✅
5. **Resistance Playbook** - 7 combos, 6 switches ✅

---

**SYNC STATUS:** Aligned with `AYESHA_ENGINE_DEEP_AUDIT.md`  
**Last Updated:** January 26, 2026
