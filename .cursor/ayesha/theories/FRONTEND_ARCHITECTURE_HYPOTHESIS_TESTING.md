# 🎯 Universal Hypothesis Testing Frontend Architecture

**Mission:** Build a frontend that can test **ANY** hypothesis, claim, or theory with the same ease as testing a single compound.

**Current Status:** Single-compound testing works ✅  
**Target:** Batch testing, comparative analysis, theory validation ✅

---

## 🏗️ **ARCHITECTURE OVERVIEW**

### **Core Principle: Test Anything, Anytime**

The frontend should support:
1. **Single Item Testing** - Test one food/compound/drug (existing)
2. **Batch Testing** - Test multiple items simultaneously
3. **Theory Validation** - Test entire theories (like "10 Cancer-Fighting Foods")
4. **Comparative Analysis** - Compare results side-by-side
5. **Mechanism Exploration** - Filter and explore by biological mechanisms

---

## 📊 **COMPONENT ARCHITECTURE**

### **Layer 1: Core Testing Engine**

```
UniversalHypothesisTester
├── InputLayer (Single/Batch/Theory)
├── ProcessingLayer (Orchestration)
└── OutputLayer (Results/Comparison)
```

#### **1. InputLayer Components:**

**a) SingleTestInput** (Existing - `DynamicFoodValidator.jsx`)
- Single compound input
- Disease context
- Patient medications
- ✅ **Already built**

**b) BatchTestInput** (NEW - `BatchFoodValidator.jsx`)
- Multi-compound input (paste list, upload CSV, select from library)
- Shared disease context
- Batch processing status
- Progress tracking

**c) TheoryTestInput** (NEW - `TheoryValidator.jsx`)
- Theory definition (JSON/YAML)
- List of claims to test
- Expected outcomes
- Comparative framework

#### **2. ProcessingLayer Components:**

**a) TestOrchestrator** (NEW - `HypothesisTestOrchestrator.jsx`)
- Manages parallel/sequential API calls
- Error handling and retry logic
- Progress tracking
- Result aggregation

**b) BatchProcessor** (NEW - `BatchProcessor.jsx`)
- Processes multiple items efficiently
- Rate limiting and queuing
- Partial result handling
- Resume capability

#### **3. OutputLayer Components:**

**a) SingleResultDisplay** (Existing - `DynamicFoodValidator.jsx` results)
- PercentileBar
- EvidenceQualityChips
- MechanismPanel
- ✅ **Already built**

**b) BatchResultsDisplay** (NEW - `BatchResultsTable.jsx`)
- Sortable/filterable table
- Score comparison
- Mechanism grouping
- Export functionality

**c) ComparativeAnalysis** (NEW - `ComparativeAnalysisPanel.jsx`)
- Side-by-side comparison
- Radar charts
- Mechanism heatmaps
- Ranking visualization

**d) TheoryValidationReport** (NEW - `TheoryValidationReport.jsx`)
- Overall theory score
- Claim-by-claim breakdown
- Evidence synthesis
- Recommendations

---

## 🎨 **FRONTEND ORGANIZATION**

### **File Structure:**

```
oncology-frontend/src/
├── pages/
│   ├── DynamicFoodValidator.jsx          ✅ Single compound (existing)
│   ├── BatchFoodValidator.jsx            🆕 Batch testing (NEW)
│   ├── TheoryValidator.jsx               🆕 Theory validation (NEW)
│   └── MechanismExplorer.jsx             🆕 Mechanism-based filtering (NEW)
│
├── components/
│   ├── food/                              ✅ Existing components
│   │   ├── PercentileBar.jsx
│   │   ├── EvidenceQualityChips.jsx
│   │   └── MechanismPanel.jsx
│   │
│   ├── batch/                             🆕 Batch testing components (NEW)
│   │   ├── BatchTestInput.jsx
│   │   ├── BatchResultsTable.jsx
│   │   ├── BatchProgressTracker.jsx
│   │   └── BatchExportButton.jsx
│   │
│   ├── theory/                            🆕 Theory validation components (NEW)
│   │   ├── TheoryDefinitionForm.jsx
│   │   ├── TheoryValidationReport.jsx
│   │   ├── ClaimBreakdown.jsx
│   │   └── TheoryScoreCard.jsx
│   │
│   └── comparison/                        🆕 Comparative analysis (NEW)
│       ├── ComparativeAnalysisPanel.jsx
│       ├── MechanismHeatmap.jsx
│       ├── RadarChart.jsx
│       └── RankingVisualization.jsx
│
├── services/
│   ├── hypothesisTestService.js          🆕 Unified API service (NEW)
│   ├── batchProcessor.js                  🆕 Batch orchestration (NEW)
│   └── theoryValidator.js                 🆕 Theory validation logic (NEW)
│
└── hooks/
    ├── useFoodValidation.js               ✅ Single compound (existing)
    ├── useBatchValidation.js              🆕 Batch processing (NEW)
    └── useTheoryValidation.js             🆕 Theory validation (NEW)
```

---

## 🔌 **BACKEND INTEGRATION**

### **API Endpoints (Existing):**

```javascript
// Single compound testing
POST /api/hypothesis/validate_food_dynamic

// Batch testing (NEW - needs backend)
POST /api/hypothesis/validate_food_batch
{
  "compounds": ["Green Tea", "Broccoli", "Papaya", ...],
  "disease_context": {...},
  "options": {
    "parallel": true,
    "max_concurrent": 5
  }
}

// Theory validation (NEW - needs backend)
POST /api/hypothesis/validate_theory
{
  "theory_name": "10 Cancer-Fighting Foods",
  "claims": [
    {"compound": "Green Tea", "mechanism": "anti-angiogenic", ...},
    {"compound": "Broccoli", "mechanism": "immune-boosting", ...},
    ...
  ],
  "disease_context": {...}
}
```

### **Frontend Service Layer:**

```javascript
// services/hypothesisTestService.js
export class HypothesisTestService {
  // Single test
  async testSingle(compound, diseaseContext) {
    return await fetch('/api/hypothesis/validate_food_dynamic', {...});
  }
  
  // Batch test
  async testBatch(compounds, diseaseContext, options) {
    return await fetch('/api/hypothesis/validate_food_batch', {...});
  }
  
  // Theory validation
  async validateTheory(theoryDefinition) {
    return await fetch('/api/hypothesis/validate_theory', {...});
  }
}
```

---

## 🎯 **USE CASE: "10 CANCER-FIGHTING FOODS"**

### **Scenario 1: Batch Test All 10 Foods**

**User Flow:**
1. User navigates to "Batch Food Validator"
2. Pastes list: "Green Tea, Broccoli, Papaya, Purple Potatoes, Pomegranates, ..."
3. Selects disease: "Ovarian Cancer (HGS)"
4. Clicks "Test All"
5. Frontend shows:
   - Progress bar (10/10 complete)
   - Real-time results as each completes
   - Final ranked table

**Implementation:**
```jsx
// pages/BatchFoodValidator.jsx
import { useState } from 'react';
import BatchTestInput from '../components/batch/BatchTestInput';
import BatchResultsTable from '../components/batch/BatchResultsTable';
import { useBatchValidation } from '../hooks/useBatchValidation';

export default function BatchFoodValidator() {
  const [compounds, setCompounds] = useState([]);
  const { results, loading, progress, testBatch } = useBatchValidation();
  
  const handleTest = async () => {
    await testBatch(compounds, {
      disease: 'ovarian_cancer_hgs',
      // ... other context
    });
  };
  
  return (
    <Box>
      <BatchTestInput 
        onCompoundsChange={setCompounds}
        onTest={handleTest}
      />
      <BatchProgressTracker progress={progress} />
      <BatchResultsTable results={results} />
    </Box>
  );
}
```

### **Scenario 2: Theory Validation**

**User Flow:**
1. User navigates to "Theory Validator"
2. Loads theory definition (JSON/YAML) or fills form
3. Clicks "Validate Theory"
4. Frontend shows:
   - Overall theory score
   - Claim-by-claim breakdown
   - Evidence synthesis
   - Recommendations

**Implementation:**
```jsx
// pages/TheoryValidator.jsx
import TheoryDefinitionForm from '../components/theory/TheoryDefinitionForm';
import TheoryValidationReport from '../components/theory/TheoryValidationReport';
import { useTheoryValidation } from '../hooks/useTheoryValidation';

export default function TheoryValidator() {
  const { result, loading, validateTheory } = useTheoryValidation();
  
  const handleValidate = async (theoryDefinition) => {
    await validateTheory(theoryDefinition);
  };
  
  return (
    <Box>
      <TheoryDefinitionForm onSubmit={handleValidate} />
      <TheoryValidationReport result={result} />
    </Box>
  );
}
```

### **Scenario 3: Mechanism Explorer**

**User Flow:**
1. User navigates to "Mechanism Explorer"
2. Filters by mechanism: "Anti-angiogenic"
3. Frontend shows all foods with anti-angiogenic properties
4. User can compare scores, evidence, pathways

**Implementation:**
```jsx
// pages/MechanismExplorer.jsx
import MechanismFilter from '../components/comparison/MechanismFilter';
import MechanismHeatmap from '../components/comparison/MechanismHeatmap';
import ComparativeAnalysisPanel from '../components/comparison/ComparativeAnalysisPanel';

export default function MechanismExplorer() {
  const [selectedMechanism, setSelectedMechanism] = useState(null);
  const [results, setResults] = useState([]);
  
  // Filter results by mechanism
  const filteredResults = results.filter(r => 
    r.mechanisms.includes(selectedMechanism)
  );
  
  return (
    <Box>
      <MechanismFilter 
        onSelect={setSelectedMechanism}
        mechanisms={['anti-angiogenic', 'immune-boosting', 'anti-inflammatory']}
      />
      <MechanismHeatmap results={filteredResults} />
      <ComparativeAnalysisPanel results={filteredResults} />
    </Box>
  );
}
```

---

## 🎨 **UI/UX DESIGN PRINCIPLES**

### **1. Consistency Across All Testing Modes**
- Same visual components (PercentileBar, EvidenceQualityChips, MechanismPanel)
- Same color coding (green=high, yellow=moderate, red=low)
- Same data structure (S/P/E scores, evidence, mechanisms)

### **2. Progressive Disclosure**
- **Single Test:** Full detail view
- **Batch Test:** Summary table → click for detail
- **Theory Test:** Overall score → claim breakdown → evidence

### **3. Real-Time Feedback**
- Progress bars for batch processing
- Live updates as results complete
- Error handling with retry options

### **4. Export & Sharing**
- Export results to CSV/JSON
- Shareable links for theory validations
- PDF reports for documentation

---

## 🔄 **DATA FLOW**

### **Single Test Flow:**
```
User Input → SingleTestInput → API Call → SingleResultDisplay
```

### **Batch Test Flow:**
```
User Input → BatchTestInput → BatchProcessor → 
  Parallel API Calls → Result Aggregation → BatchResultsTable
```

### **Theory Validation Flow:**
```
Theory Definition → TheoryDefinitionForm → 
  BatchProcessor (for each claim) → 
  Result Aggregation → TheoryValidationReport
```

---

## 📋 **IMPLEMENTATION PRIORITY**

### **Phase 1: Batch Testing (P0 - 1 week)**
- [ ] Create `BatchFoodValidator.jsx` page
- [ ] Create `BatchTestInput.jsx` component
- [ ] Create `BatchResultsTable.jsx` component
- [ ] Create `useBatchValidation.js` hook
- [ ] Create backend `/api/hypothesis/validate_food_batch` endpoint
- [ ] Test with "10 Cancer-Fighting Foods" list

### **Phase 2: Comparative Analysis (P1 - 3 days)**
- [ ] Create `ComparativeAnalysisPanel.jsx`
- [ ] Create `MechanismHeatmap.jsx`
- [ ] Create `RadarChart.jsx` (using recharts or chart.js)
- [ ] Integrate into `BatchFoodValidator.jsx`

### **Phase 3: Theory Validation (P1 - 1 week)**
- [ ] Create `TheoryValidator.jsx` page
- [ ] Create `TheoryDefinitionForm.jsx`
- [ ] Create `TheoryValidationReport.jsx`
- [ ] Create backend `/api/hypothesis/validate_theory` endpoint
- [ ] Test with "10 Cancer-Fighting Foods" theory

### **Phase 4: Mechanism Explorer (P2 - 3 days)**
- [ ] Create `MechanismExplorer.jsx` page
- [ ] Create `MechanismFilter.jsx` component
- [ ] Integrate filtering with batch results

---

## 🎯 **SUCCESS CRITERIA**

### **User Can:**
1. ✅ Test single compound (existing)
2. 🆕 Test 10 foods simultaneously (batch)
3. 🆕 Validate entire theory (10 foods as one theory)
4. 🆕 Compare results side-by-side
5. 🆕 Filter by mechanism type
6. 🆕 Export results for documentation

### **Frontend Provides:**
1. ✅ Consistent UI across all testing modes
2. 🆕 Real-time progress tracking
3. 🆕 Error handling and retry
4. 🆕 Export functionality (CSV/JSON/PDF)
5. 🆕 Responsive design (mobile-friendly)

---

## 🚀 **NEXT STEPS**

1. **Immediate:** Implement Phase 1 (Batch Testing) for "10 Cancer-Fighting Foods"
2. **Short-term:** Add comparative analysis and mechanism filtering
3. **Long-term:** Theory validation framework for any hypothesis

---

**DOCTRINE STATUS:** Ready for Implementation  
**LAST UPDATED:** 2025-01-XX



