# Frontend Integration Guide

## Component Overview

We've created 3 new components for treatment line integration:

1. **TreatmentHistoryForm.jsx** - Collects patient treatment line context
2. **TreatmentLineProvenance.jsx** - Displays treatment line calculations
3. **SAETreatmentLineChips.jsx** - Shows SAE treatment line features

---

## Integration Steps

### Step 1: Copy Components to Frontend

Copy the 3 component files to your frontend components directory:

```bash
# From: .cursor/ayesha/treatment_lines/frontend/components/
# To: oncology-coPilot/oncology-frontend/src/components/

cp .cursor/ayesha/treatment_lines/frontend/components/TreatmentHistoryForm.jsx \
   oncology-coPilot/oncology-frontend/src/components/

cp .cursor/ayesha/treatment_lines/frontend/components/TreatmentLineProvenance.jsx \
   oncology-coPilot/oncology-frontend/src/components/

cp .cursor/ayesha/treatment_lines/frontend/components/SAETreatmentLineChips.jsx \
   oncology-coPilot/oncology-frontend/src/components/
```

---

### Step 2: Wire Treatment History Form

**Option A: Add to VUS Explorer (Mutation Explorer)**

In `oncology-frontend/src/components/vus/AnalysisResults.jsx` or similar:

```jsx
import TreatmentHistoryForm from '../TreatmentHistoryForm';

// Inside your component:
const [treatmentHistory, setTreatmentHistory] = useState(null);

// Handler
const handleTreatmentHistorySubmit = (history) => {
    setTreatmentHistory(history);
    // Re-run analysis with treatment history
    fetchAnalysis({ ...currentRequest, treatment_history: history });
};

// Render
<TreatmentHistoryForm 
    onSubmit={handleTreatmentHistorySubmit}
    disease={disease} // e.g., "ovarian_cancer"
    defaultLine={1}
    defaultPriorTherapies={[]}
/>
```

**Option B: Add to Myeloma Digital Twin**

In `oncology-frontend/src/components/MyelomaDigitalTwin.jsx`:

```jsx
import TreatmentHistoryForm from './TreatmentHistoryForm';

// Add to state
const [treatmentHistory, setTreatmentHistory] = useState(null);

// Update efficacy prediction call
const predictEfficacy = async () => {
    const payload = {
        mutations: variants,
        disease: "multiple_myeloma",
        model_id: "evo2_1b",
        options: { adaptive: true },
        treatment_history: treatmentHistory  // ← Add this
    };
    
    const response = await fetch(`${API_BASE}/api/efficacy/predict`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload)
    });
    
    // ... handle response
};

// Render form before prediction button
<TreatmentHistoryForm 
    onSubmit={setTreatmentHistory}
    disease="multiple_myeloma"
/>
```

---

### Step 3: Display Treatment Line Provenance

In your efficacy results display (e.g., `EfficacyModal.jsx` or per-drug cards):

```jsx
import TreatmentLineProvenance from './TreatmentLineProvenance';

// For each drug result:
{drug.treatment_line_provenance && (
    <TreatmentLineProvenance 
        provenance={drug.treatment_line_provenance} 
    />
)}
```

**Example full drug card:**

```jsx
<div className="drug-card">
    <h3>{drug.drug_name}</h3>
    
    <div className="scores">
        <span>Efficacy: {(drug.efficacy_score * 100).toFixed(0)}%</span>
        <span>Confidence: {(drug.confidence * 100).toFixed(0)}%</span>
    </div>
    
    {/* Treatment Line Analysis */}
    {drug.treatment_line_provenance && (
        <TreatmentLineProvenance 
            provenance={drug.treatment_line_provenance} 
        />
    )}
    
    {/* Other drug info */}
    <div className="rationale">
        {drug.rationale.map((r, i) => <p key={i}>{r.text}</p>)}
    </div>
</div>
```

---

### Step 4: Add SAE Treatment Line Chips

If you have SAE features displayed, integrate the treatment line chips:

```jsx
import SAETreatmentLineChips from './SAETreatmentLineChips';

// In your SAE features section:
<div className="sae-features">
    {/* Existing SAE chips (exon, hotspot, etc.) */}
    <SAEFeatureChips features={saeFeatures} />
    
    {/* Treatment line chips */}
    <SAETreatmentLineChips saeFeatures={saeFeatures} />
</div>
```

**Note**: The `SAETreatmentLineChips` component automatically filters for treatment line features (`line_appropriateness`, `cross_resistance_risk`, `sequencing_fitness`), so you can pass the full SAE features array.

---

### Step 5: Update API Call

Ensure your efficacy prediction API call includes `treatment_history`:

```javascript
const efficacyRequest = {
    mutations: [
        {
            gene: "BRCA1",
            hgvs_p: "p.Gln356Ter",
            chrom: "17",
            pos: 43094464,
            ref: "C",
            alt: "T"
        }
    ],
    disease: "ovarian_cancer",
    model_id: "evo2_1b",
    options: {
        adaptive: true,
        ensemble: true
    },
    treatment_history: {  // ← Add this
        current_line: 2,
        prior_therapies: ["carboplatin", "paclitaxel"]
    }
};

const response = await fetch('http://localhost:8000/api/efficacy/predict', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(efficacyRequest)
});

const data = await response.json();

// Each drug in data.drugs will have treatment_line_provenance:
data.drugs.forEach(drug => {
    console.log(drug.treatment_line_provenance);
    // {
    //     current_line: 2,
    //     prior_therapies: ["carboplatin", "paclitaxel"],
    //     line_appropriateness: 1.0,
    //     cross_resistance_risk: 0.4,
    //     sequencing_fitness: 0.6,
    //     nccn_category: "1",
    //     confidence_penalty: 0.08,
    //     rationale: "Reduced by 8.0% due to cross-resistance risk"
    // }
});
```

---

## Example: Complete Integration in Ayesha's VUS Explorer

```jsx
import React, { useState } from 'react';
import TreatmentHistoryForm from './TreatmentHistoryForm';
import TreatmentLineProvenance from './TreatmentLineProvenance';
import SAETreatmentLineChips from './SAETreatmentLineChips';

const AyeshaVUSExplorer = () => {
    const [treatmentHistory, setTreatmentHistory] = useState(null);
    const [efficacyResults, setEfficacyResults] = useState(null);
    const [saeFeatures, setSaeFeatures] = useState([]);

    const runAnalysis = async () => {
        const response = await fetch('http://localhost:8000/api/efficacy/predict', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                mutations: [
                    {
                        gene: "BRCA1",
                        hgvs_p: "p.Gln356Ter",
                        chrom: "17",
                        pos: 43094464,
                        ref: "C",
                        alt: "T"
                    }
                ],
                disease: "ovarian_cancer",
                model_id: "evo2_1b",
                treatment_history: treatmentHistory
            })
        });
        
        const data = await response.json();
        setEfficacyResults(data);
        
        // If SAE features are included in response
        if (data.sae_features) {
            setSaeFeatures(data.sae_features);
        }
    };

    return (
        <div className="ayesha-vus-explorer">
            <h2>Ayesha's VUS Explorer</h2>
            
            {/* Treatment History Form */}
            <TreatmentHistoryForm 
                onSubmit={(history) => {
                    setTreatmentHistory(history);
                    runAnalysis();
                }}
                disease="ovarian_cancer"
                defaultLine={2}
                defaultPriorTherapies={["carboplatin", "paclitaxel"]}
            />
            
            {/* SAE Features with Treatment Line Chips */}
            {saeFeatures.length > 0 && (
                <div className="sae-section">
                    <h3>SAE Features</h3>
                    <SAETreatmentLineChips saeFeatures={saeFeatures} />
                </div>
            )}
            
            {/* Efficacy Results */}
            {efficacyResults && (
                <div className="efficacy-results">
                    <h3>Drug Rankings</h3>
                    {efficacyResults.drugs.map((drug, index) => (
                        <div key={index} className="drug-card">
                            <h4>{drug.drug_name}</h4>
                            <p>Confidence: {(drug.confidence * 100).toFixed(0)}%</p>
                            
                            {/* Treatment Line Provenance */}
                            {drug.treatment_line_provenance && (
                                <TreatmentLineProvenance 
                                    provenance={drug.treatment_line_provenance} 
                                />
                            )}
                        </div>
                    ))}
                </div>
            )}
        </div>
    );
};

export default AyeshaVUSExplorer;
```

---

## Testing Checklist

- [ ] Treatment history form renders correctly
- [ ] Form validates input (line > 1 requires prior therapies)
- [ ] Quick-add therapy suggestions work
- [ ] Form submission triggers API call with `treatment_history`
- [ ] Backend returns `treatment_line_provenance` for each drug
- [ ] Provenance display shows scores, NCCN category, confidence penalty
- [ ] SAE treatment line chips render with correct colors
- [ ] Confidence penalty is visible in drug rankings
- [ ] Rationale explains confidence adjustment

---

## Styling Notes

All components use inline styles for portability. To integrate with your theme:

1. Replace inline styles with your CSS classes
2. Update color schemes to match your design system
3. Adjust spacing/padding to match existing components

---

## Backend Compatibility

These components expect the backend to return:

```typescript
{
    drugs: [
        {
            drug_name: string,
            efficacy_score: number,
            confidence: number,
            treatment_line_provenance: {
                current_line: number,
                prior_therapies: string[],
                line_appropriateness: number,
                cross_resistance_risk: number,
                sequencing_fitness: number,
                nccn_category: string,
                confidence_penalty: number,
                rationale: string
            }
        }
    ],
    sae_features: [
        {
            id: string,
            name: string,
            activation: number,
            impact: "positive" | "negative" | "neutral",
            explanation: string
        }
    ]
}
```

This schema is already implemented in Phase 3 (backend integration complete).

---

## Questions?

See:
- Phase 3 completion report: `.cursor/ayesha/treatment_lines/docs/PHASE3_COMPLETION.md`
- Backend integration test: `.cursor/ayesha/treatment_lines/backend/tests/test_phase3_integration.py`
- Execution plan: `.cursor/ayesha/treatment_lines/EXECUTION_PLAN.md`


