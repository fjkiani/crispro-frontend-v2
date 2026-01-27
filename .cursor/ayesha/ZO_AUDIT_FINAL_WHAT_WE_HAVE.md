# 🔬 ZO FINAL AUDIT: WHAT WE HAVE vs WHAT'S DISPLAYED

**Date:** January 26, 2026  
**Auditor:** Zo  
**Mission:** Identify built engines NOT being utilized for Ayesha

---

## 🎯 KEY FINDING: ENGINES ARE BUILT, JUST NOT IN THE RIGHT PLACE

### The Problem Is NOT "Missing Engines"

The problem is **fragmented display** - engines exist but are:
1. On separate pages (not on Ayesha's dashboard)
2. Blocked by data that exists elsewhere in the profile
3. Not mapped to show Ayesha's ACTUAL available treatment options

---

## ✅ ENGINES BUILT AND WORKING

| Engine | Backend | Frontend | Location | Ayesha's Page? |
|--------|---------|----------|----------|----------------|
| **Timing & Chemosensitivity** | ✅ `timing_chemo_features.py` | ✅ `TimingFeaturesCard.jsx` | AyeshaTrialExplorer Tab 2 | ✅ YES (but shows "treatment-naive") |
| **Resistance Playbook** | ✅ `resistance_playbook_service.py` | ✅ `ResistancePlaybook.jsx` | AyeshaTrialExplorer Tab 4 | ✅ YES |
| **Resistance Prophet** | ✅ `resistance_prophet_service.py` | ✅ Inline in Tab 4 | AyeshaTrialExplorer Tab 4 | ✅ YES |
| **Synthetic Lethality** | ✅ `synthetic_lethality/` | ✅ `SyntheticLethalityCard.jsx` | AyeshaTrialExplorer Tab 5 | ✅ YES |
| **DDR_bin Classification** | ✅ `ddr_bin_engine.py` | ✅ `DDRStatusPage.jsx` | **SEPARATE PAGE** `/ddr-status` | ❌ **NOT ON AYESHA** |
| **Sporadic Gates** | ✅ `sporadic_gates.py` | ✅ `SporadicProvenanceCard.jsx` | DrugRankingPanel | ✅ YES (when WIWFM shows drugs) |
| **Holistic Score** | ✅ `holistic_score_service.py` | ✅ `TrialMatchCard.jsx` | AyeshaTrialExplorer Tab 1 | ✅ YES |
| **Mechanism Vector** | ✅ In orchestrator | ✅ `MechanismVectorVisualization.jsx` | AyeshaTrialExplorer Tab 0 | ✅ YES |
| **CA-125 Intelligence** | ✅ `ca125_intelligence.py` | ✅ `CA125Tracker.jsx` | AyeshaTrialExplorer Tab 3 | ⚠️ Shows empty (ca125=null) |
| **SAE Features** | ✅ `sae_service.py` | ✅ `AyeshaSAEFeaturesCard.jsx` | AyeshaTrialExplorer Tab 0,4 | ⚠️ Shows "awaiting_ngs" |

---

## ❌ WHAT'S NOT BEING DISPLAYED FOR AYESHA

### Issue 1: DDR Status Engine NOT on Ayesha's Page

**Engine exists:** `DDRStatusPage.jsx` - Full DDR classification page
- `DDRStatusCard.jsx` - Shows DDR_bin status
- `DDRFeatureBreakdown.jsx` - Shows feature breakdown
- `DDRTreatmentEligibility.jsx` - Shows "PARP ELIGIBLE" badge
- `DDRRecommendationsPanel.jsx` - Shows recommendations

**Not utilized because:** It's on a separate page (`/ddr-status`), not integrated into `AyeshaTrialExplorer.jsx`

**FIX:** Import and display DDR components on Tab 0 (Overview)

---

### Issue 2: Available Treatment Options NOT Computed

**Document exists:** `AYESHA_AVAILABLE_OPTIONS.md` - Lists all options
- SOC: Carboplatin + Paclitaxel + Bevacizumab
- PARP Trials: Olaparib, Niraparib, Rucaparib
- DDR Combinations: PARP+ATR, PARP+WEE1
- Expected Holistic Scores

**Not displayed because:** No component renders this document's computed options

**FIX:** Create `TreatmentOptionsCard.jsx` that computes and displays these with confidence

---

### Issue 3: Drug Efficacy Shows "Awaiting NGS" But Has Germline Data

**What we have:**
```javascript
germline: {
  mutations: [{
    gene: "MBD4",
    variant: "c.1293delA",            // Full genomic coordinates!
    protein_change: "p.K431Nfs*54",   // Protein change!
    classification: "pathogenic"
  }]
}
```

**What happens now:** WIWFM shows "Awaiting NGS" because `tumor_context.somatic_mutations` is incomplete

**What should happen:** WIWFM should use germline MBD4 to:
1. Compute DDR pathway score (0.88)
2. Show PARP drugs with confidence (germline-based)
3. Flag as "Germline-Only Mode" with "Add NGS to unlock full ranking"

**FIX:** Backend already has `sporadic_gates` - just need to ensure germline passes through

---

### Issue 4: MBD4 Significance NOT Explained

**What we have:**
- MBD4 homozygous pathogenic (c.1293delA)
- p53 IHC mutant type

**What should be explained:**
- MBD4 = Base Excision Repair (BER) pathway → DDR deficiency
- Homozygous = Both copies broken
- DDR_defective → PARP eligible
- MBD4 + TP53 = Synthetic lethality opportunity

**Current display:** Just shows "Germline: Positive (MBD4)" chip

**FIX:** Create `GermlineMutationCard.jsx` that explains significance

---

## 🎯 THE REAL GAPS (Priority Order)

### P0: Add DDR Components to Ayesha's Page (4 hours)

Just import existing DDR components from `/pages/DDRStatusPage.jsx` into `AyeshaTrialExplorer.jsx` Overview tab:

```jsx
// In AyeshaTrialExplorer.jsx Tab 0 (Overview)
import { DDRStatusCard, DDRTreatmentEligibility } from '../components/ddr';
import { useDDRStatus } from '../hooks/useDDRStatus';

// Use the hook with Ayesha's profile
const { ddrStatus, loading: ddrLoading, calculateDDRStatus } = useDDRStatus();

useEffect(() => {
  // Compute DDR status from Ayesha's mutations
  calculateDDRStatus({
    mutations: [
      ...AYESHA_11_17_25_PROFILE.germline?.mutations || [],
      ...AYESHA_11_17_25_PROFILE.tumor_context?.somatic_mutations || []
    ]
  });
}, []);

// Then display:
<Grid item xs={12} md={6}>
  <DDRStatusCard ddrStatus={ddrStatus} />
</Grid>
<Grid item xs={12} md={6}>
  <DDRTreatmentEligibility ddrStatus={ddrStatus} />
</Grid>
```

### P1: Wire Treatment Options to Engine Output (6 hours)

Create `TreatmentOptionsCard.jsx` that:
1. Takes trials + SOC + drug efficacy as input
2. Groups by pathway alignment (DDR-targeted, Platinum, IO, etc.)
3. Shows expected holistic scores
4. Indicates immediate availability (SOC=immediate, trials=enrollment needed)

### P2: Fix WIWFM to Use Germline When NGS Missing (8 hours)

Modify `drug_efficacy_service.py` to:
1. Check if `tumor_context` is incomplete
2. If germline has DDR mutations (MBD4, BRCA, etc.), use that for DDR pathway
3. Return drugs with "Germline-Only Mode" provenance
4. Flag which additional tests unlock full ranking

### P3: Add Mutation Significance Card (4 hours)

Create `MutationSignificanceCard.jsx` that:
1. Takes a mutation (germline or somatic)
2. Looks up pathway impact (DDR, MAPK, PI3K, etc.)
3. Explains what it means (DDR_defective, synthetic lethality, etc.)
4. Shows treatment implications

---

## 📊 CURRENT vs TARGET STATE

### Current Experience

```
User sees Ayesha's page:
┌───────────────────────────────────────────────────────────┐
│ Profile: Stage IVB, Germline: Positive (MBD4)            │
│                                                           │
│ [Overview] [Trials] [Treatment] [Monitoring] [Resistance] │
│                                                           │
│ 🧬 Mechanism Intelligence                                 │
│   - Mechanism vector visualization (DDR=0.88)             │
│   - Pathway disruption card                               │
│                                                           │
│ 💊 Drug Efficacy: "Awaiting NGS"  ← BLOCKED               │
│                                                           │
│ 🔬 SAE Features: "awaiting_ngs"   ← BLOCKED               │
│                                                           │
└───────────────────────────────────────────────────────────┘
```

### Target Experience (After Fixes)

```
User sees Ayesha's page:
┌───────────────────────────────────────────────────────────┐
│ Profile: Stage IVB, DDR_DEFECTIVE ← NEW                  │
│                                                           │
│ [Overview] [Trials] [Treatment] [Monitoring] [Resistance] │
│                                                           │
│ 🧬 DDR Status: DDR_DEFECTIVE (from MBD4+TP53)  ← NEW     │
│   ├─ PARP ELIGIBLE badge                                  │
│   ├─ Synthetic Lethality: MBD4+TP53 → PARP vulnerable    │
│   └─ Recommendations panel                                │
│                                                           │
│ 🎯 Treatment Options (computed from engines):  ← NEW     │
│   ┌─ SOC: Carboplatin+Paclitaxel+Bevacizumab (immediate) │
│   ├─ PARP Trials: 0.90-0.95 holistic score               │
│   ├─ DDR Combos: PARP+ATR (0.80-0.85)                    │
│   └─ Lower: Non-DDR trials (0.45-0.60)                   │
│                                                           │
│ 💊 Drug Efficacy (Germline Mode):  ← NEW                 │
│   Based on MBD4 germline → DDR=0.88                       │
│   "Add NGS to unlock full WIWFM ranking"                  │
│                                                           │
└───────────────────────────────────────────────────────────┘
```

---

## 🚀 ACTION PLAN (22 hours total)

| Priority | Task | Hours | Components to Use |
|----------|------|-------|-------------------|
| **P0** | Add DDR components to Ayesha Overview | 4h | `DDRStatusCard`, `DDRTreatmentEligibility` (existing) |
| **P0** | Auto-compute DDR status from mutations | 2h | `useDDRStatus` hook (existing) |
| **P1** | Create TreatmentOptionsCard | 6h | New component, uses trial + SOC + WIWFM data |
| **P2** | Fix WIWFM germline fallback | 6h | Backend `drug_efficacy_service.py` |
| **P3** | MutationSignificanceCard | 4h | New component |

**Total: 22 hours** - All using EXISTING engines!

---

## 📝 SUMMARY

**The engines exist.** We don't need to reinvent anything.

**What's missing:**
1. DDR components not on Ayesha's main page (separate page `/ddr-status`)
2. No unified "Treatment Options" view that computes from all engines
3. WIWFM doesn't fall back to germline when NGS missing
4. Mutation significance not explained

**Fix strategy:** Wire existing components to existing engines, add 2 new cards.

---

**AUDIT COMPLETE** 💀
