# METASTASIS INTERCEPTION - FRONTEND COMPLETION PLAN

**Goal:** Create demo-ready metastasis dashboard with dedicated route  
**Timeline:** 4-6 hours  
**Priority:** P1 (publication complete, demo needed for partners)

---

## CURRENT STATE

### **✅ What Exists (Backend)**
- All endpoints operational (100% test coverage)
- `/api/metastasis/assess` - 8-step cascade risk assessment
- `/api/metastasis_interception/intercept` - CRISPR weapon design
- Real ClinVar data (14 pathogenic variants)
- Publication figures, tables, datasets

### **⚠️ What Exists (Frontend)**
- `MetastasisReport.jsx` (966 lines) - 8-step risk visualization
- `MetastasisInterceptionPanel.jsx` (966 lines) - weapon design panel
- Hooks: `useMetastasisAssess.js`, `useMetastasisInterception.js`
- **Problem:** Buried in VUS Explorer (not exposed, hard to demo)

### **❌ What's Missing**
- Dedicated `/metastasis` route
- Standalone dashboard (patient variant input → assessment → weapon design)
- Quick demo flow (not nested 5 layers deep)

---

## SOLUTION: METASTASIS DASHBOARD

### **New Route:** `/metastasis`

**Component:** `MetastasisDashboard.jsx`

**Layout:**
```
┌──────────────────────────────────────────────────┐
│ METASTASIS INTERCEPTION PLATFORM (RUO)           │
└──────────────────────────────────────────────────┘

┌────────────────────────────────────────┐
│ STEP 1: PATIENT VARIANT INPUT          │
│ ┌────────────────────────────────────┐ │
│ │ Gene: [BRAF ▼]                     │ │
│ │ Variant: [V600E]                   │ │
│ │ OR paste coords: chr7:140753336 T>A│ │
│ │                                    │ │
│ │ [Analyze Metastatic Risk →]        │ │
│ └────────────────────────────────────┘ │
└────────────────────────────────────────┘

┌────────────────────────────────────────┐
│ STEP 2: 8-STEP CASCADE ASSESSMENT      │
│ ┌────────────────────────────────────┐ │
│ │ [MetastasisReport component]       │ │
│ │ - Risk bars per step               │ │
│ │ - Driver genes table               │ │
│ │ - Provenance                       │ │
│ └────────────────────────────────────┘ │
└────────────────────────────────────────┘

┌────────────────────────────────────────┐
│ STEP 3: CRISPR WEAPON DESIGN           │
│ ┌────────────────────────────────────┐ │
│ │ Select highest-risk step:          │ │
│ │ [Angiogenesis] [EMT] [Invasion] .. │ │
│ │                                    │ │
│ │ [MetastasisInterceptionPanel]      │ │
│ │ - Target lock scores               │ │
│ │ - Ranked guide candidates          │ │
│ │ - Assassin scores                  │ │
│ └────────────────────────────────────┘ │
└────────────────────────────────────────┘
```

---

## IMPLEMENTATION PLAN

### **Phase 1: Create Dashboard Component (2 hours)**

**File:** `oncology-coPilot/oncology-frontend/src/pages/MetastasisDashboard.jsx`

**Structure:**
```jsx
import React, { useState } from 'react';
import MetastasisReport from '../components/metastasis/MetastasisReport.jsx';
import MetastasisInterceptionPanel from '../components/metastasis/MetastasisInterceptionPanel.jsx';
import { useMetastasisAssess } from '../hooks/useMetastasis.js';
import { useMetastasisInterception } from '../hooks/useMetastasisInterception.js';
import RUOLabel from '../components/common/RUOLabel.jsx';

const KNOWN_VARIANTS = [
  { label: 'BRAF V600E', gene: 'BRAF', hgvsp: 'V600E', chrom: '7', pos: 140753336, ref: 'T', alt: 'A' },
  { label: 'KRAS G12D', gene: 'KRAS', hgvsp: 'G12D', chrom: '12', pos: 25398284, ref: 'C', alt: 'T' },
  // ... (14 real ClinVar variants)
];

const MISSION_STEPS = [
  'primary_growth', 'local_invasion', 'intravasation',
  'survival_in_circulation', 'extravasation',
  'micrometastasis_formation', 'angiogenesis', 'metastatic_colonization'
];

export default function MetastasisDashboard() {
  const [selectedVariant, setSelectedVariant] = useState(null);
  const [selectedMissionStep, setSelectedMissionStep] = useState(null);

  // Assessment data
  const assessmentData = useMetastasisAssess({
    mutations: selectedVariant ? [selectedVariant] : [],
    disease: 'PanCancer',
    options: { profile: 'baseline' }
  });

  // Interception data (only fetch when step selected)
  const interceptionData = useMetastasisInterception({
    missionStep: selectedMissionStep,
    mutations: selectedVariant ? [selectedVariant] : [],
    options: { profile: 'baseline' }
  }, !!selectedMissionStep);

  return (
    <div className="min-h-screen bg-slate-900 p-6">
      <RUOLabel className="mb-6" />
      
      {/* Header */}
      <div className="mb-8">
        <h1 className="text-4xl font-bold text-cyan-300 mb-2">
          Metastasis Interception Platform
        </h1>
        <p className="text-slate-400">
          Assess metastatic cascade risk → Design CRISPR interventions
        </p>
      </div>

      {/* Step 1: Variant Input */}
      <div className="bg-slate-800 rounded-lg p-6 mb-6">
        <h2 className="text-2xl font-bold text-white mb-4">Step 1: Select Variant</h2>
        <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
          {KNOWN_VARIANTS.map(v => (
            <button
              key={v.label}
              onClick={() => setSelectedVariant(v)}
              className={`p-4 rounded transition-colors ${
                selectedVariant?.label === v.label
                  ? 'bg-cyan-600 text-white'
                  : 'bg-slate-700 text-slate-300 hover:bg-slate-600'
              }`}
            >
              {v.label}
            </button>
          ))}
        </div>
      </div>

      {/* Step 2: Assessment (only show when variant selected) */}
      {selectedVariant && (
        <div className="bg-slate-800 rounded-lg p-6 mb-6">
          <h2 className="text-2xl font-bold text-white mb-4">
            Step 2: Cascade Assessment
          </h2>
          <MetastasisReport
            data={assessmentData.data}
            loading={assessmentData.loading}
            error={assessmentData.error}
          />
        </div>
      )}

      {/* Step 3: Weapon Design (only show when variant selected) */}
      {selectedVariant && (
        <div className="bg-slate-800 rounded-lg p-6">
          <h2 className="text-2xl font-bold text-white mb-4">
            Step 3: Design CRISPR Weapons
          </h2>
          <p className="text-slate-400 mb-4">
            Select a high-risk cascade step to design targeted interventions:
          </p>
          <div className="flex flex-wrap gap-2 mb-6">
            {MISSION_STEPS.map(step => (
              <button
                key={step}
                onClick={() => setSelectedMissionStep(step === selectedMissionStep ? null : step)}
                className={`px-4 py-2 rounded transition-colors ${
                  selectedMissionStep === step
                    ? 'bg-cyan-600 text-white'
                    : 'bg-slate-700 text-slate-300 hover:bg-slate-600'
                }`}
              >
                {step.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase())}
              </button>
            ))}
          </div>

          {selectedMissionStep && (
            <MetastasisInterceptionPanel
              data={interceptionData.data}
              loading={interceptionData.loading}
              error={interceptionData.error}
            />
          )}
        </div>
      )}
    </div>
  );
}
```

---

### **Phase 2: Add Route (10 minutes)**

**File:** `oncology-coPilot/oncology-frontend/src/App.jsx`

**Add:**
```jsx
import MetastasisDashboard from './pages/MetastasisDashboard.jsx';

// In <Routes>:
<Route path="/metastasis" element={<MetastasisDashboard />} />
```

---

### **Phase 3: Add Nav Link (5 minutes)**

**File:** Navigation component (wherever main nav is)

**Add:**
```jsx
<Link to="/metastasis" className="...">
  Metastasis Interception
</Link>
```

---

### **Phase 4: Polish Components (1-2 hours)**

**Minor Updates Needed:**

1. **`MetastasisReport.jsx`:**
   - ✅ Already complete (displays 8-step bars, drivers table, provenance)
   - Add: Quick "Select for Interception" buttons on high-risk steps

2. **`MetastasisInterceptionPanel.jsx`:**
   - ✅ Already complete (target lock scores, guide candidates, assassin scores)
   - Add: "Download Designs" button (export CSV)

3. **Loading States:**
   - Ensure spinners show during API calls
   - Add skeleton loaders for better UX

4. **Error Handling:**
   - Toast notifications for API failures
   - Retry buttons

---

## ACCEPTANCE CRITERIA

### **Demo Flow (2 minutes):**
1. Navigate to `/metastasis`
2. Click "BRAF V600E"
3. See 8-step cascade assessment (risk bars)
4. Click highest-risk step (e.g., "Angiogenesis")
5. See ranked targets (VEGFA top)
6. See ranked guide candidates with assassin scores
7. Click "Download Designs" → get CSV

### **Technical:**
- All API calls succeed (no 404s, no CORS errors)
- Loading states visible during fetch
- Provenance bar shows run_id, profile, methods
- RUO label prominent

---

## ALTERNATIVE: MINIMAL DEMO ROUTE (1 hour)

**If time-constrained, create simpler version:**

**File:** `oncology-coPilot/oncology-frontend/src/pages/MetastasisDemo.jsx`

**Content:**
```jsx
// Hardcoded BRAF V600E example
// Show pre-fetched assessment + pre-fetched interception results
// No live API calls (static demo)
// Just display components with mock data
```

**Benefits:**
- Fastest to ship (1 hour vs 4-6 hours)
- No API stability concerns
- Perfect for screenshots/videos

**Downsides:**
- Not interactive
- Can't test with different variants

---

## RECOMMENDATION

**Go with Full Dashboard (4-6 hours):**
- Backend is 100% stable (21/21 tests passing)
- Components already exist (just need wiring)
- Interactive demo more impressive than static
- Reusable for future demos/trials

**Timeline:**
- Phase 1 (Dashboard): 2 hours
- Phase 2 (Route): 10 minutes
- Phase 3 (Nav): 5 minutes
- Phase 4 (Polish): 1-2 hours
- **Total:** 4-6 hours

**Next Action:**
1. Create `MetastasisDashboard.jsx` with 3-step flow
2. Wire to existing components
3. Add route + nav link
4. Test demo flow end-to-end
5. Polish loading/error states

---

**STATUS:** 🟡 **READY TO IMPLEMENT**  
**BLOCKER:** None (all dependencies complete)  
**PRIORITY:** P1 (publication done, demo needed)


