# ⚔️ Platinum Response Frontend Integration - Audit & Implementation Plan

**Date:** January 13, 2025  
**Mission:** Build frontend capability to query and visualize platinum response data (like Streamlit dashboard)  
**Status:** 🔍 **AUDIT COMPLETE** - Ready for implementation

---

## 📊 CURRENT STATE AUDIT

### ✅ **What We Have (Backend)**

1. **Platinum Response Data (Jr2's Extraction)**
   - **Location:** `data/validation/tcga_ov_platinum_response_labels.json`
   - **Structure:**
     ```json
     {
       "metadata": {
         "n_patients": 469,
         "response_distribution": { "sensitive": 396, "resistant": 31, "refractory": 42 },
         "source_distribution": {...}
       },
       "patients": [
         {
           "tcga_sample_id": "TCGA-XX-XXXX",
           "patient_id": "...",
           "platinum_response": "sensitive|resistant|refractory",
           "response_value": "CR|PR|SD|PD|...",
           "source": "gdc_xml|pybioportal|broad_firehose|cbio",
           "extraction_date": "..."
         }
       ]
     }
     ```
   - **Overlap with Zo's Dataset:** 161 patients (exceeds N=40 threshold)

2. **Existing Backend Endpoint**
   - **Endpoint:** `POST /api/datasets/extract_and_benchmark`
   - **Current Capability:**
     - Extracts mutations from cBioPortal (TCGA-OV)
     - Labels with `outcome_platinum` (binary: 1=sensitive, 0=resistant)
     - Returns `{ rows, count, metrics, profile }`
   - **Limitation:** Does NOT query Jr2's platinum response JSON directly

3. **Zo's Mutation Dataset**
   - **Location:** `data/validation/tcga_ov_full_validation_dataset.json`
   - **Structure:** Contains mutations, sample IDs, clinical data
   - **Overlap:** 161 patients match Jr2's platinum response data

### ✅ **What We Have (Frontend)**

1. **CohortLab.jsx Component**
   - **Location:** `oncology-coPilot/oncology-frontend/src/components/research/CohortLab.jsx`
   - **Current Features:**
     - Study ID input (default: `ov_tcga`)
     - Genes input (comma-separated)
     - Limit, Mode, Profile selectors
     - Calls `/api/datasets/extract_and_benchmark`
     - Displays simple metrics (count, positives, prevalence, AUPRC proxy)
     - Shows basic table (sample, gene, HGVS, coord, platinum)
   - **Limitations:**
     - No rich visualizations (pie charts, bar charts, overlap analysis)
     - No platinum response distribution charts
     - No validation status display
     - No patient explorer with filters
     - No overlap analysis with Zo's dataset

2. **Research Portal Page**
   - **Location:** `oncology-coPilot/oncology-frontend/src/pages/ResearchPortal/ResearchPortal.jsx`
   - **Status:** Contains CohortLab component (basic integration)

### ❌ **What's Missing (Gaps)**

1. **Backend Gaps:**
   - ❌ No endpoint to query Jr2's platinum response JSON directly
   - ❌ No endpoint to merge Jr2's data with Zo's mutation dataset
   - ❌ No endpoint to compute overlap statistics
   - ❌ No endpoint to return response distribution charts data

2. **Frontend Gaps:**
   - ❌ No visualization components (pie charts, bar charts, Venn diagrams)
   - ❌ No hero metrics cards (like Streamlit dashboard)
   - ❌ No response distribution charts
   - ❌ No overlap analysis visualization
   - ❌ No patient explorer with search/filter
   - ❌ No validation status component
   - ❌ No integration with Zo's mutation dataset

3. **Data Flow Gaps:**
   - ❌ Frontend cannot query platinum response data directly
   - ❌ No way to merge platinum response + mutation data in frontend
   - ❌ No way to display overlap statistics

---

## 🎯 TARGET STATE (What We Want)

### **User Experience Goal**

Users should be able to:
1. **Query Platinum Response Data** - Select study, view platinum response distribution
2. **View Rich Visualizations** - Pie charts, bar charts, response distribution
3. **Explore Overlap** - See how Jr2's data overlaps with Zo's mutation dataset
4. **Validate Readiness** - Check if sample size meets statistical thresholds
5. **Browse Patients** - Search, filter, and explore individual patient data
6. **Export Data** - Download merged datasets for analysis

### **Visual Components Needed (Mirror Streamlit Dashboard)**

1. **Hero Metrics Cards** - 4 key stats (Jr2 patients, Zo patients, Overlap, Patients with mutations)
2. **Response Distribution Charts** - Pie chart + bar chart of sensitive/resistant/refractory
3. **Overlap Analysis** - Venn diagram + match rate metrics
4. **Validation Status** - Gauge chart + readiness assessment
5. **Patient Table** - Searchable, filterable table with merged data
6. **Raw Data Display** - JSON viewer for merged patient data

---

## 🔧 IMPLEMENTATION PLAN

### **Phase 1: Backend API Enhancement (2-3 hours)**

#### **Task 1.1: Create Platinum Response Query Endpoint**
- **File:** `oncology-coPilot/oncology-backend-minimal/api/routers/datasets.py`
- **New Endpoint:** `POST /api/datasets/platinum_response/query`
- **Input:**
  ```json
  {
    "study_id": "tcga_ov",
    "filters": {
      "response_type": ["sensitive", "resistant", "refractory"],  // optional
      "source": ["gdc_xml", "pybioportal", "broad_firehose", "cbio"]  // optional
    }
  }
  ```
- **Output:**
  ```json
  {
    "metadata": {
      "n_patients": 469,
      "response_distribution": {...},
      "source_distribution": {...}
    },
    "patients": [...],
    "statistics": {
      "sensitive_count": 396,
      "resistant_count": 31,
      "refractory_count": 42
    }
  }
  ```
- **Implementation:**
  - Load `data/validation/tcga_ov_platinum_response_labels.json`
  - Apply filters (response_type, source)
  - Return filtered data + statistics

#### **Task 1.2: Create Overlap Analysis Endpoint**
- **File:** `oncology-coPilot/oncology-backend-minimal/api/routers/datasets.py`
- **New Endpoint:** `POST /api/datasets/platinum_response/overlap`
- **Input:**
  ```json
  {
    "merge_with_zo": true  // optional, default true
  }
  ```
- **Output:**
  ```json
  {
    "overlap_stats": {
      "jr2_total": 469,
      "zo_total": 212,
      "overlap_count": 161,
      "match_rate_jr2": 34.3,
      "match_rate_zo": 75.9
    },
    "merged_patients": [...],  // 161 patients with both platinum response + mutations
    "validation_readiness": {
      "is_ready": true,
      "overlap_count": 161,
      "min_required": 40,
      "status": "READY_FOR_VALIDATION"
    }
  }
  ```
- **Implementation:**
  - Load both JSON files
  - Match by `tcga_sample_id` (Jr2) ↔ `sample_id` (Zo)
  - Compute overlap statistics
  - Return merged patient list + validation readiness

#### **Task 1.3: Enhance Existing Endpoint**
- **File:** `oncology-coPilot/oncology-backend-minimal/api/routers/datasets.py`
- **Enhancement:** Add `platinum_response_data` option to `extract_and_benchmark`
- **New Input Option:**
  ```json
  {
    "mode": "both",
    "study_id": "tcga_ov",
    "include_platinum_response": true,  // NEW
    "merge_with_zo": true  // NEW
  }
  ```
- **Enhanced Output:**
  ```json
  {
    "rows": [...],
    "metrics": {...},
    "platinum_response": {  // NEW
      "distribution": {...},
      "overlap_stats": {...},
      "merged_patients": [...]
    }
  }
  ```

### **Phase 2: Frontend Component Development (4-6 hours)**

#### **Task 2.1: Create Visualization Components**

**2.1.1: HeroMetricsCards Component**
- **File:** `oncology-coPilot/oncology-frontend/src/components/research/platinum/HeroMetricsCards.jsx`
- **Props:** `{ stats: { jr2_total, zo_total, overlap_count, patients_with_mutations } }`
- **Features:**
  - 4 metric cards (styled like Streamlit)
  - Color-coded values
  - Responsive grid layout

**2.1.2: ResponseDistributionCharts Component**
- **File:** `oncology-coPilot/oncology-frontend/src/components/research/platinum/ResponseDistributionCharts.jsx`
- **Props:** `{ distribution: { sensitive, resistant, refractory }, rawResponses?: [...] }`
- **Features:**
  - Pie chart (using Recharts or Plotly)
  - Bar chart
  - Raw response value breakdown table
  - Color-coded by response type

**2.1.3: OverlapAnalysis Component**
- **File:** `oncology-coPilot/oncology-frontend/src/components/research/platinum/OverlapAnalysis.jsx`
- **Props:** `{ overlapData: { jr2_total, zo_total, overlap_count, match_rates } }`
- **Features:**
  - 3 metric cards (Jr2, Overlap, Zo)
  - Venn diagram (using `react-venn-diagram` or similar)
  - Match rate breakdown
  - Fallback text if Venn library not available

**2.1.4: ValidationStatus Component**
- **File:** `oncology-coPilot/oncology-frontend/src/components/research/platinum/ValidationStatus.jsx`
- **Props:** `{ stats: { overlap_count, min_required }, isReady: boolean }`
- **Features:**
  - Status badge (READY / INSUFFICIENT)
  - Gauge chart (using Recharts)
  - Manager's criteria review
  - Success/warning messages

**2.1.5: PatientTable Component**
- **File:** `oncology-coPilot/oncology-frontend/src/components/research/platinum/PatientTable.jsx`
- **Props:** `{ patients: [...], mergedData?: [...] }`
- **Features:**
  - Searchable table (by sample ID, patient ID)
  - Filterable (by response type, mutation count)
  - Sortable columns
  - Pagination (show 100 at a time)
  - Export to CSV button

#### **Task 2.2: Create Main Platinum Response Page**

**2.2.1: PlatinumResponsePage Component**
- **File:** `oncology-coPilot/oncology-frontend/src/pages/ResearchPortal/PlatinumResponsePage.jsx`
- **Features:**
  - Tab navigation (Overview, Patients, Overlap, Validation)
  - Query form (study selection, filters)
  - Integration of all visualization components
  - Loading states
  - Error handling
  - Export functionality

**2.2.2: Update Research Portal**
- **File:** `oncology-coPilot/oncology-frontend/src/pages/ResearchPortal/ResearchPortal.jsx`
- **Enhancement:** Add "Platinum Response" tab alongside existing tabs
- **Integration:** Route to `PlatinumResponsePage` component

#### **Task 2.3: Enhance CohortLab Component**

**2.3.1: Add Platinum Response Option**
- **File:** `oncology-coPilot/oncology-frontend/src/components/research/CohortLab.jsx`
- **Enhancement:**
  - Add checkbox: "Include Platinum Response Data"
  - Add checkbox: "Merge with Zo's Mutation Dataset"
  - When checked, call new endpoints and display enhanced visualizations
  - Show response distribution charts when platinum data available

### **Phase 3: Data Integration & Testing (2-3 hours)**

#### **Task 3.1: Data Flow Integration**
- Test backend endpoints with real data
- Verify JSON loading and filtering
- Test overlap calculation accuracy
- Validate merged patient data structure

#### **Task 3.2: Frontend Integration Testing**
- Test all visualization components render correctly
- Test data flow from backend → frontend
- Test filtering and search functionality
- Test export functionality

#### **Task 3.3: E2E Testing**
- Test complete user flow: Query → Visualize → Export
- Test with different filter combinations
- Test error handling (missing data, API failures)
- Test responsive design (mobile/desktop)

---

## 📋 DETAILED TASK BREAKDOWN

### **Backend Tasks (Priority Order)**

| Task | File | Time | Status |
|------|------|------|--------|
| 1.1: Platinum Response Query Endpoint | `api/routers/datasets.py` | 1h | ⏸️ Pending |
| 1.2: Overlap Analysis Endpoint | `api/routers/datasets.py` | 1h | ⏸️ Pending |
| 1.3: Enhance Existing Endpoint | `api/routers/datasets.py` | 30min | ⏸️ Pending |

### **Frontend Tasks (Priority Order)**

| Task | File | Time | Status |
|------|------|------|--------|
| 2.1.1: HeroMetricsCards | `components/research/platinum/HeroMetricsCards.jsx` | 30min | ⏸️ Pending |
| 2.1.2: ResponseDistributionCharts | `components/research/platinum/ResponseDistributionCharts.jsx` | 1h | ⏸️ Pending |
| 2.1.3: OverlapAnalysis | `components/research/platinum/OverlapAnalysis.jsx` | 1h | ⏸️ Pending |
| 2.1.4: ValidationStatus | `components/research/platinum/ValidationStatus.jsx` | 1h | ⏸️ Pending |
| 2.1.5: PatientTable | `components/research/platinum/PatientTable.jsx` | 1.5h | ⏸️ Pending |
| 2.2.1: PlatinumResponsePage | `pages/ResearchPortal/PlatinumResponsePage.jsx` | 1.5h | ⏸️ Pending |
| 2.2.2: Update Research Portal | `pages/ResearchPortal/ResearchPortal.jsx` | 30min | ⏸️ Pending |
| 2.3.1: Enhance CohortLab | `components/research/CohortLab.jsx` | 1h | ⏸️ Pending |

### **Testing Tasks**

| Task | File | Time | Status |
|------|------|------|--------|
| 3.1: Data Flow Integration | Test scripts | 1h | ⏸️ Pending |
| 3.2: Frontend Integration Testing | Test scripts | 1h | ⏸️ Pending |
| 3.3: E2E Testing | Test scripts | 1h | ⏸️ Pending |

**Total Estimated Time:** 12-15 hours

---

## 🎨 VISUAL DESIGN SPECIFICATIONS

### **Color Scheme (Match Streamlit Dashboard)**
- **Primary:** `#FF4B4B` (Red)
- **Secondary:** `#00C0F2` (Blue)
- **Sensitive:** `#28a745` (Green)
- **Resistant:** `#dc3545` (Red)
- **Refractory:** `#ffc107` (Yellow/Orange)
- **Unknown:** `#6c757d` (Grey)
- **Background:** `#0e1117` (Dark)
- **Text:** `#FAFAFA` (Light)

### **Component Styling**
- Use Material-UI (MUI) components for consistency
- Dark theme to match existing frontend
- Responsive grid layouts (mobile/desktop)
- Loading spinners for async operations
- Error banners for failures

### **Chart Libraries**
- **Option 1:** Recharts (lightweight, React-native)
- **Option 2:** Plotly.js (more features, heavier)
- **Recommendation:** Recharts for simplicity, Plotly for advanced features

---

## 🔌 API CONTRACTS (Detailed)

### **Endpoint 1: Query Platinum Response Data**

```typescript
POST /api/datasets/platinum_response/query

Request:
{
  study_id?: string;  // default: "tcga_ov"
  filters?: {
    response_type?: ("sensitive" | "resistant" | "refractory")[];
    source?: ("gdc_xml" | "pybioportal" | "broad_firehose" | "cbio")[];
  };
}

Response:
{
  metadata: {
    n_patients: number;
    response_distribution: {
      sensitive: number;
      resistant: number;
      refractory: number;
    };
    source_distribution: Record<string, number>;
    extraction_date: string;
  };
  patients: Array<{
    tcga_sample_id: string;
    patient_id: string;
    platinum_response: "sensitive" | "resistant" | "refractory";
    response_value: string;
    source: string;
    extraction_date: string;
  }>;
  statistics: {
    sensitive_count: number;
    resistant_count: number;
    refractory_count: number;
  };
}
```

### **Endpoint 2: Overlap Analysis**

```typescript
POST /api/datasets/platinum_response/overlap

Request:
{
  merge_with_zo?: boolean;  // default: true
}

Response:
{
  overlap_stats: {
    jr2_total: number;
    zo_total: number;
    overlap_count: number;
    jr2_only: number;
    zo_only: number;
    match_rate_jr2: number;  // percentage
    match_rate_zo: number;   // percentage
  };
  merged_patients: Array<{
    // Jr2 fields
    tcga_sample_id: string;
    platinum_response: string;
    response_value: string;
    // Zo fields
    sample_id: string;
    mutations: Array<{...}>;
    // ... other merged fields
  }>;
  validation_readiness: {
    is_ready: boolean;
    overlap_count: number;
    min_required: number;
    status: "READY_FOR_VALIDATION" | "INSUFFICIENT_SAMPLE_SIZE";
  };
}
```

### **Endpoint 3: Enhanced Extract & Benchmark**

```typescript
POST /api/datasets/extract_and_benchmark

Request (Enhanced):
{
  mode: "extract_only" | "run_only" | "both";
  study_id: string;
  genes?: string[];
  limit?: number;
  profile?: "baseline" | "richer_s" | "fusion";
  include_platinum_response?: boolean;  // NEW
  merge_with_zo?: boolean;  // NEW
}

Response (Enhanced):
{
  rows?: Array<{...}>;
  count?: number;
  metrics?: {...};
  profile: string;
  platinum_response?: {  // NEW (when include_platinum_response=true)
    distribution: {...};
    overlap_stats: {...};
    merged_patients: Array<{...}>;
  };
}
```

---

## 📊 COMPARISON: STREAMLIT vs FRONTEND

### **Streamlit Dashboard Features**
✅ Hero metrics (4 cards)  
✅ Response distribution (pie + bar charts)  
✅ Overlap analysis (Venn diagram + metrics)  
✅ Validation status (gauge + readiness)  
✅ Patient table (searchable, filterable)  
✅ Raw data display (JSON viewer)  
✅ Multi-page navigation (Overview, Patients, Overlap, Validation)

### **Frontend Target Features**
✅ Hero metrics (4 cards) - **MIRROR**  
✅ Response distribution (pie + bar charts) - **MIRROR**  
✅ Overlap analysis (Venn diagram + metrics) - **MIRROR**  
✅ Validation status (gauge + readiness) - **MIRROR**  
✅ Patient table (searchable, filterable) - **MIRROR**  
✅ Raw data display (JSON viewer) - **MIRROR**  
✅ Tab navigation (Overview, Patients, Overlap, Validation) - **MIRROR**  
✅ Integration with existing CohortLab - **ENHANCEMENT**  
✅ Export to CSV/JSON - **ENHANCEMENT**  
✅ Real-time API queries (not file-based) - **ENHANCEMENT**

---

## 🚨 CRITICAL GAPS FROM `cohort_context_concept.mdc`

### **What the Concept Document Should Include (But Doesn't)**

The `cohort_context_concept.mdc` file is mostly empty (just metadata). Based on the search results and existing code, here's what it SHOULD include:

1. **Study Selection UI**
   - Dropdown/autocomplete for available studies
   - Study metadata display (n patients, date, source)

2. **Metrics Display**
   - Prevalence calculations
   - Response rates
   - Outcome statistics (OS, PFS)

3. **Artifacts Management**
   - Download extracted datasets
   - Download benchmark results
   - Download merged datasets

4. **Mapping & Coverage**
   - Gene coverage indicators
   - Variant coverage indicators
   - Sample overlap visualization

5. **Cohort Overlays**
   - Overlay cohort context on variant analysis
   - Show prevalence in cohort
   - Show response rates for variants

### **What We're Building (Beyond Concept)**

Our implementation goes BEYOND the concept document by:
- ✅ Querying platinum response data directly (not just mutations)
- ✅ Merging with Zo's mutation dataset
- ✅ Visualizing overlap statistics
- ✅ Validating statistical readiness
- ✅ Rich interactive visualizations (not just tables)

---

## ✅ ACCEPTANCE CRITERIA

### **Backend**
- [ ] `POST /api/datasets/platinum_response/query` returns filtered platinum response data
- [ ] `POST /api/datasets/platinum_response/overlap` returns accurate overlap statistics
- [ ] Enhanced `extract_and_benchmark` includes platinum response when requested
- [ ] All endpoints handle missing data gracefully
- [ ] All endpoints return proper error messages

### **Frontend**
- [ ] Hero metrics cards display correctly
- [ ] Response distribution charts render (pie + bar)
- [ ] Overlap analysis shows Venn diagram + metrics
- [ ] Validation status shows gauge + readiness assessment
- [ ] Patient table is searchable and filterable
- [ ] Tab navigation works (Overview, Patients, Overlap, Validation)
- [ ] Export functionality works (CSV/JSON)
- [ ] Loading states display during API calls
- [ ] Error handling displays user-friendly messages
- [ ] Responsive design works (mobile/desktop)

### **Integration**
- [ ] Frontend queries backend endpoints successfully
- [ ] Data flows correctly from backend → frontend
- [ ] Visualizations update when filters change
- [ ] Export downloads correct data format
- [ ] E2E user flow works end-to-end

---

## 🎯 NEXT STEPS

### **Immediate (Start Now)**
1. ✅ **AUDIT COMPLETE** - This document
2. ⏭️ **Phase 1.1** - Create platinum response query endpoint (1h)
3. ⏭️ **Phase 1.2** - Create overlap analysis endpoint (1h)

### **Short-Term (This Week)**
4. ⏭️ **Phase 1.3** - Enhance existing endpoint (30min)
5. ⏭️ **Phase 2.1** - Create visualization components (4-5h)
6. ⏭️ **Phase 2.2** - Create main page (1.5h)
7. ⏭️ **Phase 2.3** - Enhance CohortLab (1h)

### **Testing (After Implementation)**
8. ⏭️ **Phase 3** - Data integration & testing (3h)

---

## 📝 NOTES

- **Data Location:** Ensure `data/validation/tcga_ov_platinum_response_labels.json` is accessible from backend
- **Chart Libraries:** Install `recharts` or `plotly.js` for frontend visualizations
- **Venn Diagram:** Install `react-venn-diagram` or use `@upsetjs/react` for overlap visualization
- **Export:** Use `papaparse` for CSV export, native JSON for JSON export
- **Styling:** Use MUI components for consistency with existing frontend

---

**DOCTRINE STATUS: ACTIVE** ⚔️  
**READY FOR IMPLEMENTATION** ✅


