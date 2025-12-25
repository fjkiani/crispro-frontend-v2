# 🧬 Synthetic Lethality Frontend Enhancement Plan

**Date:** January 28, 2025  
**Status:** ✅ IMPLEMENTED  
**Goal:** Create a production-ready Synthetic Lethality Analyzer for clinicians

---

## ✅ IMPLEMENTATION COMPLETE

### Files Created:

| File | Description |
|------|-------------|
| `hooks/useSyntheticLethality.js` | Main analysis hook with API orchestration |
| `components/EssentialityScoreCard.jsx` | Visual gene essentiality display |
| `components/PathwayDependencyDiagram.jsx` | Broken → Essential → Drug flow diagram |
| `components/TherapyRecommendationList.jsx` | Ranked drug recommendations |
| `components/MutationInputForm.jsx` | Multi-mutation input with gene selector |
| `components/ClinicalDossierModal.jsx` | Export-ready clinical report |
| `SyntheticLethalityAnalyzer.jsx` | Main page component |
| `index.js` | Module exports |

### Route Added:
- **URL:** `/synthetic-lethality`
- **Component:** `<SyntheticLethalityAnalyzer />`

### Features:
- ✅ Multi-gene mutation input
- ✅ Disease context selection (ovarian, breast, etc.)
- ✅ Real-time essentiality scoring
- ✅ Pathway dependency visualization
- ✅ Ranked drug recommendations
- ✅ Clinical dossier generation (copy/download/print)
- ✅ "Load Example" button for Ayesha's MBD4+TP53 case

---

---

## 📊 Current State Assessment

### ✅ What We Already Have

| Component | Status | Location |
|-----------|--------|----------|
| **Backend Endpoint** | ✅ Working | `/api/guidance/synthetic_lethality` |
| **Demo Component** | ✅ Exists | `SyntheticLethalityDetective.jsx` |
| **Essentiality Hook** | ✅ Working | `useEssentiality` in `useInsights.js` |
| **CoPilot Integration** | ✅ Configured | Intent routing in `Q2CRouter/intents.js` |

### Current `SyntheticLethalityDetective.jsx` Features:
- ✅ Step-by-step analysis visualization
- ✅ Calls backend endpoint
- ✅ Shows damage report & essentiality scores
- ✅ Displays therapy recommendation
- ❌ Hardcoded to BRCA1 C61G (no user input)
- ❌ Single gene only (not multi-gene)
- ❌ No disease selector
- ❌ No clinical dossier output

### Current Backend `/api/guidance/synthetic_lethality` Returns:
```json
{
  "suggested_therapy": "platinum",
  "damage_report": [
    {"variant": {...}, "vep": {...}, "functionality": {...}}
  ],
  "essentiality_report": [
    {"gene": "BRCA1", "result": {"essentiality_score": 0.85, "flags": {...}}}
  ],
  "guidance": {...}
}
```

---

## 🎯 Vision: What We Want

### For Clinicians (Ayesha's Oncologist):

```
┌─────────────────────────────────────────────────────────────────┐
│  🔬 SYNTHETIC LETHALITY ANALYZER                                │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  PATIENT PROFILE                                                │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ Disease: [Ovarian Cancer ▼]  Stage: [IVB]                │  │
│  │ Subtype: [High-Grade Serous ▼]                           │  │
│  └──────────────────────────────────────────────────────────┘  │
│                                                                 │
│  GENETIC MUTATIONS                                              │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ + Add Mutation                                            │  │
│  │ ┌────────────────────────────────────────────────────┐   │  │
│  │ │ MBD4 | c.1239delA (frameshift) | Germline | ✅     │   │  │
│  │ └────────────────────────────────────────────────────┘   │  │
│  │ ┌────────────────────────────────────────────────────┐   │  │
│  │ │ TP53 | p.R175H (missense)     | Somatic  | ✅     │   │  │
│  │ └────────────────────────────────────────────────────┘   │  │
│  └──────────────────────────────────────────────────────────┘  │
│                                                                 │
│  [🔬 Analyze Synthetic Lethality]                              │
│                                                                 │
├─────────────────────────────────────────────────────────────────┤
│  ANALYSIS RESULTS                                               │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ Step 1: Damage Assessment      ✅ Complete               │  │
│  │ Step 2: Pathway Mapping        ✅ Complete               │  │
│  │ Step 3: Essentiality Scoring   ✅ Complete               │  │
│  │ Step 4: Synthetic Lethality    ✅ Complete               │  │
│  │ Step 5: Drug Recommendations   ✅ Complete               │  │
│  └──────────────────────────────────────────────────────────┘  │
│                                                                 │
│  ESSENTIALITY SCORES                                            │
│  ┌─────────────────────┐  ┌─────────────────────┐              │
│  │ MBD4               │  │ TP53               │              │
│  │ Score: 0.80        │  │ Score: 0.75        │              │
│  │ ████████░░ HIGH    │  │ ███████░░░ HIGH    │              │
│  │ Frameshift → LoF   │  │ Hotspot → Inactive │              │
│  │ BER: NON-FUNCTIONAL│  │ Checkpoint: BYPASS │              │
│  └─────────────────────┘  └─────────────────────┘              │
│                                                                 │
│  SYNTHETIC LETHALITY OPPORTUNITIES                              │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │  🎯 DOUBLE-HIT EFFECT DETECTED                           │  │
│  │                                                          │  │
│  │  MBD4 Loss ──┐                                          │  │
│  │              ├──► HR Pathway Essential ──► PARP Target  │  │
│  │  TP53 Loss ──┤                                          │  │
│  │              └──► ATR/CHK1 Essential ──► ATR Target     │  │
│  └──────────────────────────────────────────────────────────┘  │
│                                                                 │
│  RECOMMENDED THERAPIES                                          │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ #1 OLAPARIB (PARP)     Confidence: 89%  [FDA ✓]        │  │
│  │ #2 NIRAPARIB (PARP)    Confidence: 87%  [FDA ✓]        │  │
│  │ #3 CERALASERTIB (ATR)  Confidence: 72%  [Clinical]     │  │
│  └──────────────────────────────────────────────────────────┘  │
│                                                                 │
│  [📄 Generate Clinical Dossier]  [📊 Export PDF]               │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## 🏗️ Implementation Plan

### Option A: Enhance Existing Component (Quick Win)
**Timeline:** 1-2 days  
**Effort:** Low

Modify `SyntheticLethalityDetective.jsx` to:
1. Add variant input form (gene, hgvs_p, consequence)
2. Add disease selector
3. Support multiple mutations
4. Keep existing step visualization

**Pros:** Quick, reuses existing UI  
**Cons:** Demo aesthetic, not production-ready

---

### Option B: Create New Production Component (Recommended)
**Timeline:** 3-5 days  
**Effort:** Medium

Create `SyntheticLethalityAnalyzer.jsx`:

#### Components to Create:

| Component | Purpose |
|-----------|---------|
| `SyntheticLethalityAnalyzer.jsx` | Main page container |
| `PatientProfileForm.jsx` | Disease, stage, subtype inputs |
| `MutationInputPanel.jsx` | Multi-mutation input (reuse from MyelomaDigitalTwin) |
| `EssentialityScoreCard.jsx` | Display gene essentiality with visual gauge |
| `PathwayDependencyDiagram.jsx` | Visual of broken vs essential pathways |
| `SyntheticLethalityCard.jsx` | Show double-hit effect explanation |
| `TherapyRecommendationList.jsx` | Ranked drugs with confidence |
| `ClinicalDossierModal.jsx` | Generate the formatted dossier |

#### File Structure:
```
src/components/SyntheticLethality/
├── SyntheticLethalityAnalyzer.jsx    # Main page
├── components/
│   ├── PatientProfileForm.jsx
│   ├── MutationInputPanel.jsx
│   ├── EssentialityScoreCard.jsx
│   ├── PathwayDependencyDiagram.jsx
│   ├── SyntheticLethalityCard.jsx
│   ├── TherapyRecommendationList.jsx
│   └── ClinicalDossierModal.jsx
├── hooks/
│   └── useSyntheticLethality.js      # Orchestrates API calls
└── index.js
```

---

### Option C: Integrate into MyelomaDigitalTwin (Full Integration)
**Timeline:** 2-3 days  
**Effort:** Medium

Add "Synthetic Lethality" tab to existing `MyelomaDigitalTwin.jsx`:
- Reuse existing variant input
- Add essentiality analysis section
- Show synthetic lethality results inline

**Pros:** Consistent UX, no new page  
**Cons:** MyelomaDigitalTwin is already complex

---

## 📋 Recommended Approach: Option B

### Phase 1: Core Infrastructure (Day 1)

1. **Create hook `useSyntheticLethality.js`:**
```javascript
export function useSyntheticLethality({ disease, mutations }) {
  const [loading, setLoading] = useState(false);
  const [results, setResults] = useState(null);
  const [error, setError] = useState(null);

  const analyze = useCallback(async () => {
    setLoading(true);
    try {
      // Call /api/guidance/synthetic_lethality
      const response = await fetch(`${API_ROOT}/api/guidance/synthetic_lethality`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ disease, mutations })
      });
      const data = await response.json();
      setResults(data);
    } catch (err) {
      setError(err.message);
    } finally {
      setLoading(false);
    }
  }, [disease, mutations]);

  return { analyze, loading, results, error };
}
```

2. **Create main page `SyntheticLethalityAnalyzer.jsx`:**
   - Patient profile form
   - Mutation input (reuse `VariantInputList`)
   - Analysis trigger button

### Phase 2: Results Display (Day 2)

3. **Create `EssentialityScoreCard.jsx`:**
```jsx
const EssentialityScoreCard = ({ gene, score, flags, rationale, pathwayImpact }) => (
  <Card>
    <Typography variant="h6">{gene}</Typography>
    <LinearProgress value={score * 100} color={score >= 0.7 ? 'error' : 'warning'} />
    <Typography>Score: {score.toFixed(2)}</Typography>
    <Chip label={flags.frameshift ? 'Frameshift' : flags.hotspot ? 'Hotspot' : 'Variant'} />
    <Typography variant="body2">{pathwayImpact}</Typography>
  </Card>
);
```

4. **Create `PathwayDependencyDiagram.jsx`:**
   - Visual showing broken pathways
   - Arrows to essential backup pathways
   - Connection to drug targets

5. **Create `TherapyRecommendationList.jsx`:**
   - Ranked drug list with confidence
   - FDA approval status
   - Mechanism of action

### Phase 3: Clinical Output (Day 3)

6. **Create `ClinicalDossierModal.jsx`:**
   - Generate formatted clinical dossier
   - Export as PDF/Markdown
   - Based on template from `AYESHA_CLINICAL_DOSSIER_ESSENTIALITY.md`

### Phase 4: Integration (Day 4-5)

7. **Add route to app:**
```jsx
<Route path="/synthetic-lethality" element={<SyntheticLethalityAnalyzer />} />
```

8. **Add to navigation menu**

9. **Connect to CoPilot for conversational access**

---

## 🎨 UI Design Guidelines

### Color Scheme (Consistent with ClinicalDossier):
- **High Essentiality (≥0.7):** Red/Error
- **Moderate (0.5-0.7):** Orange/Warning
- **Low (<0.5):** Green/Success

### Visual Elements:
- Gauge charts for essentiality scores
- Step-by-step progress (like existing Detective)
- Pathway diagrams (SVG or MUI icons)
- Card-based layout for recommendations

### Accessibility:
- Clear labels on all interactive elements
- Color-blind friendly indicators
- Screen reader support

---

## 📡 API Enhancements (Optional)

### Current Endpoint Works, But Could Add:

1. **Enhanced response format:**
```json
{
  "suggested_therapies": [
    {"drug": "Olaparib", "target": "PARP", "confidence": 0.89, "evidence_tier": "I"},
    {"drug": "Ceralasertib", "target": "ATR", "confidence": 0.72, "evidence_tier": "II"}
  ],
  "pathway_analysis": {
    "broken_pathways": ["BER", "G1/S Checkpoint"],
    "essential_pathways": ["HR", "ATR/CHK1"],
    "double_hit_detected": true
  },
  "synthetic_lethality_explanation": "MBD4 loss + TP53 loss creates dependency on HR and ATR pathways..."
}
```

2. **Clinical dossier generation endpoint:**
```
POST /api/dossier/generate_synthetic_lethality_report
```

---

## ✅ Success Criteria

| Criteria | Measurement |
|----------|-------------|
| **User Input** | Can input any gene/variant combination |
| **Multi-Gene** | Supports 2+ mutations simultaneously |
| **Disease Context** | Can select disease type and stage |
| **Essentiality Display** | Shows scores with visual gauges |
| **Pathway Visualization** | Clear diagram of broken vs essential |
| **Drug Recommendations** | Ranked list with confidence |
| **Clinical Dossier** | Exportable formatted report |
| **Doctor Usable** | Oncologist can use without technical help |

---

## 🚀 Next Steps

1. **Choose approach:** Option B recommended (production component)
2. **Create folder structure:** `src/components/SyntheticLethality/`
3. **Start with hook:** `useSyntheticLethality.js`
4. **Build incrementally:** One component at a time
5. **Test with Ayesha case:** MBD4 + TP53

---

## 📁 Related Files

| File | Purpose |
|------|---------|
| `SyntheticLethalityDetective.jsx` | Existing demo (reference) |
| `MyelomaDigitalTwin.jsx` | Pattern for variant input |
| `useInsights.js` | Essentiality hook pattern |
| `AYESHA_CLINICAL_DOSSIER_ESSENTIALITY.md` | Dossier template |
| `guidance.py` | Backend endpoint |

---

**Ready to proceed?** Let me know which option you prefer, and I'll start building!

