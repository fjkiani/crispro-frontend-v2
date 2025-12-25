# Universal Pages Audit & Deliverables Plan

**Date:** January 28, 2025  
**Purpose:** Complete audit of Universal frontend pages and comprehensive deliverables plan  
**Status:** 🔴 **AUDIT COMPLETE** - Deliverables defined  
**Location:** `.cursor/MOAT/UNIVERSAL_PAGES_AUDIT_AND_DELIVERABLES.md`

---

## 🎯 Executive Summary

**Audited Components:**
1. `UniversalCompleteCare.jsx` - Complete care plan orchestration
2. `UniversalTrialIntelligence.jsx` - Trial filtering and dossier generation
3. `UniversalDossierBrowser.jsx` - Dossier browsing and management
4. `UniversalDossierDetail.jsx` - Individual dossier detail view

**Key Findings:**
- ✅ **Working:** Basic functionality, API integration, component structure
- ⚠️ **Incomplete:** Missing components (SAE Features, Resistance Playbook), no orchestrator integration
- 🔴 **Gaps:** No testing, limited error handling, missing real-time updates, no orchestrator state management

**Priority Deliverables:**
1. **Complete Missing Components** (2-3 days)
2. **Orchestrator Integration** (3-4 days)
3. **Testing & Validation** (2-3 days)
4. **Enhancements** (2-3 days)

**Total Estimated Time:** 9-13 days

---

## 📋 Component Audit

### 1. UniversalCompleteCare.jsx

#### **Current State**

**Purpose:** Universal complete care plan page that orchestrates all MOAT capabilities via `/api/complete_care/v2`

**What Works:**
- ✅ Patient profile loading (default AK profile for demo)
- ✅ API integration with `/api/complete_care/v2`
- ✅ Component structure (DrugRankingPanel, BiomarkerCard, ResistanceCard, etc.)
- ✅ Export JSON functionality
- ✅ Provenance modal
- ✅ Summary stats display
- ✅ Grid layout for analysis components

**What's Missing/Incomplete:**
- ⚠️ **TODO:** Resistance Playbook component (line 499) - Shows alert placeholder
- ⚠️ **TODO:** SAE Features component (line 519) - Shows alert placeholder
- ⚠️ **Missing:** Orchestrator integration (`/api/orchestrate/full` not used)
- ⚠️ **Missing:** Real-time pipeline status updates
- ⚠️ **Missing:** Error recovery/retry mechanisms
- ⚠️ **Missing:** Loading states for individual components
- ⚠️ **Missing:** Patient profile upload (file upload capability)
- ⚠️ **Missing:** State persistence (refresh loses data)

**API Integration:**
- ✅ Uses: `/api/complete_care/v2` (POST)
- ❌ Not using: `/api/orchestrate/full` (orchestrator endpoint)
- ❌ Not using: `/api/orchestrate/status/{patient_id}` (status polling)
- ❌ Not using: `/api/patients/{patient_id}` (state retrieval)

**Component Dependencies:**
- ✅ `DrugRankingPanel` - Working
- ✅ `BiomarkerCard` - Working
- ✅ `ResistanceCard` - Working
- ✅ `ToxicityRiskCard` - Working
- ✅ `TrialMatchesCard` - Working
- ✅ `SOCRecommendationCard` - Working
- ✅ `NextTestCard` - Working
- ✅ `HintTilesPanel` - Working
- ✅ `MechanismChips` - Working
- ⚠️ `ResistancePlaybook` - **MISSING** (TODO)
- ⚠️ `SAEFeatures` - **MISSING** (TODO)

**Data Flow:**
```
User → Generate Plan Button
  → POST /api/complete_care/v2
  → Response: Complete care plan with all sections
  → Display in grid layout
```

**Issues:**
1. **No Orchestrator Integration:** Uses legacy `/api/complete_care/v2` instead of orchestrator
2. **No State Management:** Doesn't track pipeline phases or agent execution
3. **No Real-time Updates:** Can't show progress of individual agents
4. **Missing Components:** Resistance Playbook and SAE Features show placeholders
5. **No File Upload:** Can't upload VCF/PDF/MAF files for data extraction
6. **No Error Recovery:** Limited error handling, no retry logic

---

### 2. UniversalTrialIntelligence.jsx

#### **Current State**

**Purpose:** Complete interface for filtering trials and generating dossiers for any patient

**What Works:**
- ✅ Tab-based workflow (Patient Profile → Filter Trials → Generate Dossiers → Autonomous Flow)
- ✅ Patient profile form integration
- ✅ Trial filtering via `/api/dossiers/intelligence/filter`
- ✅ Dossier generation via `/api/dossiers/intelligence/generate`
- ✅ Batch dossier generation
- ✅ Autonomous flow (`/api/trials/agent/generate-dossiers`)
- ✅ JSON trial candidate input
- ✅ Dossier summary cards display

**What's Missing/Incomplete:**
- ⚠️ **Missing:** Integration with orchestrator trial matching agent
- ⚠️ **Missing:** Mechanism fit display in filtered results
- ⚠️ **Missing:** Real-time progress for dossier generation
- ⚠️ **Missing:** Error handling for partial failures
- ⚠️ **Missing:** Trial search integration (relies on external search)
- ⚠️ **Missing:** Patient profile persistence
- ⚠️ **Missing:** Dossier quality indicators

**API Integration:**
- ✅ Uses: `/api/dossiers/intelligence/filter` (POST)
- ✅ Uses: `/api/dossiers/intelligence/generate` (POST)
- ✅ Uses: `/api/dossiers/intelligence/batch-generate` (POST)
- ✅ Uses: `/api/trials/agent/generate-dossiers` (POST)
- ❌ Not using: `/api/trials/agent/search` (mechanism fit ranking)
- ❌ Not using: `/api/orchestrate/full` (orchestrator pipeline)

**Component Dependencies:**
- ✅ `PatientProfileForm` - Working
- ✅ `UniversalDossierSummaryCard` - Working
- ⚠️ Missing: Trial search component integration
- ⚠️ Missing: Mechanism fit display components

**Data Flow:**
```
Tab 1: Patient Profile
  → Save profile
  → Store in state

Tab 2: Filter Trials
  → Paste JSON candidates
  → POST /api/dossiers/intelligence/filter
  → Display filtered results (top_tier, good_tier)

Tab 3: Generate Dossiers
  → Select NCT IDs
  → POST /api/dossiers/intelligence/batch-generate
  → Display generated dossiers

Tab 4: Autonomous Flow
  → POST /api/trials/agent/generate-dossiers
  → Auto-search → Filter → Generate
```

**Issues:**
1. **No Orchestrator Integration:** Doesn't use orchestrator trial matching agent
2. **No Mechanism Fit Display:** Filtered results don't show mechanism fit scores
3. **Manual Trial Input:** Requires manual JSON paste instead of integrated search
4. **No Progress Tracking:** Can't see dossier generation progress
5. **Limited Error Handling:** Basic error display, no recovery options

---

### 3. UniversalDossierBrowser.jsx

#### **Current State**

**Purpose:** Browse and manage trial intelligence dossiers for any patient

**What Works:**
- ✅ Patient ID input
- ✅ Dossier listing via `/api/dossiers/intelligence/list/{patient_id}`
- ✅ Search functionality (NCT ID, keywords)
- ✅ Tier filtering (ALL, TOP_TIER, GOOD_TIER)
- ✅ Dossier summary cards
- ✅ Patient profile form integration
- ✅ Export dossier functionality

**What's Missing/Incomplete:**
- ⚠️ **Missing:** Real-time updates (new dossiers don't appear automatically)
- ⚠️ **Missing:** Dossier status indicators (generating, complete, failed)
- ⚠️ **Missing:** Bulk operations (delete, export multiple)
- ⚠️ **Missing:** Sorting options (by date, tier, match score)
- ⚠️ **Missing:** Pagination for large dossier lists
- ⚠️ **Missing:** Dossier comparison view

**API Integration:**
- ✅ Uses: `/api/dossiers/intelligence/list/{patient_id}` (GET)
- ✅ Uses: `/api/dossiers/intelligence/{patient_id}/{nct_id}` (GET) - via export
- ❌ Not using: WebSocket/SSE for real-time updates
- ❌ Not using: Orchestrator status endpoints

**Component Dependencies:**
- ✅ `UniversalDossierSummaryCard` - Working
- ✅ `PatientProfileForm` - Working
- ⚠️ Missing: Dossier status indicators
- ⚠️ Missing: Bulk action components

**Data Flow:**
```
User enters Patient ID
  → GET /api/dossiers/intelligence/list/{patient_id}
  → Display dossier list
  → Apply filters (tier, search)
  → Display filtered results
```

**Issues:**
1. **No Real-time Updates:** Dossiers generated elsewhere don't appear automatically
2. **No Status Tracking:** Can't see if dossier generation is in progress
3. **Limited Filtering:** Only tier and search, no date/sort options
4. **No Pagination:** All dossiers loaded at once (performance issue for large lists)

---

### 4. UniversalDossierDetail.jsx

#### **Current State**

**Purpose:** Display full dossier markdown for any patient

**What Works:**
- ✅ Route parameters (patientId, nct_id)
- ✅ Dossier loading via `/api/dossiers/intelligence/{patient_id}/{nct_id}`
- ✅ Markdown rendering (ReactMarkdown with remarkGfm)
- ✅ Export markdown functionality
- ✅ Back navigation
- ✅ Metadata display (tier, match score)

**What's Missing/Incomplete:**
- ⚠️ **Missing:** Interactive elements (expandable sections, tabs)
- ⚠️ **Missing:** Dossier versioning/history
- ⚠️ **Missing:** Comparison with other dossiers
- ⚠️ **Missing:** Print-friendly view
- ⚠️ **Missing:** Share functionality
- ⚠️ **Missing:** Annotations/comments

**API Integration:**
- ✅ Uses: `/api/dossiers/intelligence/{patient_id}/{nct_id}` (GET)
- ❌ Not using: Dossier versioning endpoints
- ❌ Not using: Dossier comparison endpoints

**Component Dependencies:**
- ✅ `ReactMarkdown` - Working
- ⚠️ Missing: Interactive markdown components
- ⚠️ Missing: Dossier comparison components

**Data Flow:**
```
Route: /dossiers/{patientId}/{nct_id}
  → GET /api/dossiers/intelligence/{patientId}/{nct_id}
  → Render markdown
  → Display metadata
```

**Issues:**
1. **Static Display:** Markdown is read-only, no interactivity
2. **No Versioning:** Can't see dossier history or updates
3. **Limited Navigation:** Only back button, no related dossiers
4. **No Comparison:** Can't compare multiple dossiers side-by-side

---

## 🔗 Orchestrator Integration Analysis

### **Current State vs Orchestrator Capabilities**

| Feature | Current Pages | Orchestrator Provides | Gap |
|---------|--------------|----------------------|-----|
| **Pipeline Execution** | `/api/complete_care/v2` | `/api/orchestrate/full` | ❌ Not integrated |
| **Status Tracking** | None | `/api/orchestrate/status/{patient_id}` | ❌ Not integrated |
| **State Management** | Local React state | `PatientState` object | ❌ Not integrated |
| **Real-time Updates** | None | WebSocket/SSE (planned) | ❌ Not integrated |
| **Agent Execution** | Single API call | 7-phase pipeline | ❌ Not integrated |
| **Error Recovery** | Basic | Agent-level error handling | ❌ Not integrated |
| **File Upload** | None | VCF/PDF/MAF parsing | ❌ Not integrated |

### **What Orchestrator Provides (Not Used)**

1. **7-Phase Pipeline:**
   - Phase 1: Data Extraction (VCF/PDF/MAF)
   - Phase 2: Parallel Analysis (Biomarker, Resistance, Nutrition)
   - Phase 3: Drug Ranking (S/P/E)
   - Phase 4: Trial Matching (Mechanism Fit)
   - Phase 5: Care Plan Generation
   - Phase 6: Monitoring Setup
   - Phase 7: Complete

2. **State Management:**
   - `PatientState` object tracks all agent outputs
   - Full audit trail of state changes
   - History tracking for debugging

3. **Agent Coordination:**
   - Parallel execution where dependencies allow
   - Sequential execution for dependent agents
   - Error handling and recovery
   - Execution timing and performance tracking

4. **API Endpoints:**
   - `POST /api/orchestrate/full` - Run complete pipeline
   - `GET /api/orchestrate/status/{patient_id}` - Get status
   - `GET /api/patients/{patient_id}` - Get full state
   - `GET /api/patients/{patient_id}/care-plan` - Get care plan

---

## 📊 Deliverables Plan

### **Phase 1: Complete Missing Components** (2-3 days)

#### **Deliverable 1.1: Resistance Playbook Component** 🔴 **HIGH PRIORITY**

**Objective:** Create `ResistancePlaybook.jsx` component to display resistance playbook data

**Requirements:**
- Display alternative treatment options if resistance detected
- Show resistance pathways and mechanisms
- Display drug recommendations for resistance scenarios
- Include evidence and citations

**Implementation:**
```jsx
// Location: oncology-frontend/src/components/orchestrator/Analysis/ResistancePlaybook.jsx

// Props:
interface ResistancePlaybookProps {
  resistance_playbook: {
    resistance_detected: boolean;
    pathways: Array<{
      pathway: string;
      resistance_mechanism: string;
      alternative_drugs: Array<{
        drug_name: string;
        rationale: string;
        evidence_tier: string;
      }>;
    }>;
    recommendations: Array<{
      action: string;
      priority: 'high' | 'medium' | 'low';
      rationale: string;
    }>;
  };
}
```

**Acceptance Criteria:**
- ✅ Component displays resistance playbook data
- ✅ Shows pathways and mechanisms
- ✅ Lists alternative drug recommendations
- ✅ Includes evidence tiers and citations
- ✅ Integrated into `UniversalCompleteCare.jsx`

**Estimated Time:** 4-6 hours

---

#### **Deliverable 1.2: SAE Features Component** 🔴 **HIGH PRIORITY**

**Objective:** Create `SAEFeatures.jsx` component to display SAE features and diagnostics

**Requirements:**
- Display DNA repair capacity
- Show pathway burden scores (DDR, MAPK, PI3K, etc.)
- Display DDR_bin score (if TRUE SAE available)
- Show SAE provenance (TRUE SAE vs PROXY SAE)
- Include validation metrics

**Implementation:**
```jsx
// Location: oncology-frontend/src/components/orchestrator/Analysis/SAEFeatures.jsx

// Props:
interface SAEFeaturesProps {
  sae_features: {
    provenance: {
      sae: 'true_sae' | 'proxy';
      sae_version?: string;
      mapping_version?: string;
    };
    sae_diagnostics: {
      ddr_bin_score?: number;
      ddr_sae_score?: number;
      mapk_sae_score?: number;
      // ... other pathway scores
    };
    dna_repair_capacity?: number;
    pathway_burden?: {
      ddr: number;
      mapk: number;
      pi3k: number;
      // ... other pathways
    };
  };
}
```

**Acceptance Criteria:**
- ✅ Component displays SAE features and diagnostics
- ✅ Shows DNA repair capacity gauge
- ✅ Displays pathway burden scores
- ✅ Shows DDR_bin score with TRUE SAE indicator
- ✅ Includes validation metrics (AUROC, etc.)
- ✅ Integrated into `UniversalCompleteCare.jsx`

**Estimated Time:** 6-8 hours

---

### **Phase 2: Orchestrator Integration** (3-4 days)

#### **Deliverable 2.1: Replace `/api/complete_care/v2` with Orchestrator** 🔴 **CRITICAL**

**Objective:** Integrate orchestrator pipeline into `UniversalCompleteCare.jsx`

**Requirements:**
- Replace `/api/complete_care/v2` with `/api/orchestrate/full`
- Add file upload capability (VCF/PDF/MAF)
- Implement status polling (`/api/orchestrate/status/{patient_id}`)
- Display pipeline phase progress
- Show agent execution status

**Implementation Steps:**

1. **Add File Upload:**
```jsx
// Add to UniversalCompleteCare.jsx
import { PatientUpload } from '../components/orchestrator/Patient/PatientUpload';

// Add file upload section before "Generate Plan" button
<PatientUpload
  onUploadComplete={(patientId) => {
    setPatientProfile({ ...patientProfile, patient_id: patientId });
    handleGeneratePlan();
  }}
  patientId={patientProfile?.patient_id}
/>
```

2. **Replace API Call:**
```jsx
// OLD:
const response = await fetch(`${API_ROOT}/api/complete_care/v2`, {
  method: 'POST',
  body: JSON.stringify({ patient_profile: patientProfile, ... })
});

// NEW:
const response = await fetch(`${API_ROOT}/api/orchestrate/full`, {
  method: 'POST',
  body: JSON.stringify({
    patient_id: patientProfile.patient_id,
    disease: patientProfile.disease?.type,
    mutations: patientProfile.tumor_context?.somatic_mutations || [],
    // ... other fields
  })
});
```

3. **Add Status Polling:**
```jsx
// Add status polling hook
const [pipelineStatus, setPipelineStatus] = useState(null);

useEffect(() => {
  if (patientProfile?.patient_id && !result) {
    const interval = setInterval(async () => {
      const statusResponse = await fetch(
        `${API_ROOT}/api/orchestrate/status/${patientProfile.patient_id}`
      );
      const status = await statusResponse.json();
      setPipelineStatus(status);
      
      if (status.phase === 'COMPLETE') {
        // Load full state
        const stateResponse = await fetch(
          `${API_ROOT}/api/patients/${patientProfile.patient_id}`
        );
        const state = await stateResponse.json();
        setResult(state);
        clearInterval(interval);
      }
    }, 2000); // Poll every 2 seconds
    
    return () => clearInterval(interval);
  }
}, [patientProfile, result]);
```

4. **Display Pipeline Status:**
```jsx
// Add pipeline status component
import { PipelineStatusCard } from '../components/orchestrator/Dashboard/PipelineStatusCard';

{pipelineStatus && (
  <PipelineStatusCard status={pipelineStatus} />
)}
```

**Acceptance Criteria:**
- ✅ File upload works (VCF/PDF/MAF)
- ✅ Orchestrator pipeline executes successfully
- ✅ Status polling shows real-time progress
- ✅ All agent outputs displayed correctly
- ✅ Error handling for failed agents
- ✅ Backward compatible with existing patient profiles

**Estimated Time:** 8-10 hours

---

#### **Deliverable 2.2: Integrate Orchestrator Trial Matching** 🟡 **HIGH PRIORITY**

**Objective:** Use orchestrator trial matching agent in `UniversalTrialIntelligence.jsx`

**Requirements:**
- Use `/api/orchestrate/full` for trial matching phase
- Display mechanism fit scores in filtered results
- Show mechanism alignment breakdown
- Integrate with orchestrator state management

**Implementation:**
```jsx
// In UniversalTrialIntelligence.jsx

// Replace manual trial filtering with orchestrator
const handleOrchestratorTrialMatching = async () => {
  // Run orchestrator pipeline (trial matching phase only)
  const response = await fetch(`${API_ROOT}/api/orchestrate/full`, {
    method: 'POST',
    body: JSON.stringify({
      patient_id: patientProfile.patient_id,
      disease: patientProfile.disease?.type,
      mutations: patientProfile.tumor_context?.somatic_mutations || [],
      skip_agents: ['01', '02', '03', '04', '06', '07', '08'], // Only run trial matching
    })
  });
  
  const state = await response.json();
  // Extract trials from state.trial_matching
  setFilteredResults({
    top_tier: state.trial_matching?.trials?.filter(t => t.combined_score >= 0.8) || [],
    good_tier: state.trial_matching?.trials?.filter(t => t.combined_score >= 0.6 && t.combined_score < 0.8) || [],
  });
};
```

**Acceptance Criteria:**
- ✅ Orchestrator trial matching agent used
- ✅ Mechanism fit scores displayed
- ✅ Mechanism alignment breakdown shown
- ✅ Combined scores (0.7×eligibility + 0.3×mechanism_fit) displayed
- ✅ TRUE SAE indicators shown (if available)

**Estimated Time:** 6-8 hours

---

### **Phase 3: Testing & Validation** (2-3 days)

#### **Deliverable 3.1: Component Testing** 🟡 **HIGH PRIORITY**

**Objective:** Create comprehensive tests for all Universal pages

**Test Coverage:**

1. **UniversalCompleteCare.jsx:**
   - ✅ Renders correctly with default AK profile
   - ✅ Generates care plan successfully
   - ✅ Displays all components (biomarker, resistance, trials, etc.)
   - ✅ Handles API errors gracefully
   - ✅ Exports JSON correctly
   - ✅ Shows provenance modal
   - ✅ File upload works (VCF/PDF/MAF)
   - ✅ Orchestrator integration works
   - ✅ Status polling updates correctly

2. **UniversalTrialIntelligence.jsx:**
   - ✅ Patient profile form works
   - ✅ Trial filtering works
   - ✅ Dossier generation works
   - ✅ Batch generation works
   - ✅ Autonomous flow works
   - ✅ Error handling works

3. **UniversalDossierBrowser.jsx:**
   - ✅ Dossier listing works
   - ✅ Search functionality works
   - ✅ Tier filtering works
   - ✅ Export works

4. **UniversalDossierDetail.jsx:**
   - ✅ Dossier loading works
   - ✅ Markdown rendering works
   - ✅ Export works
   - ✅ Navigation works

**Test Framework:**
- React Testing Library
- Jest
- MSW (Mock Service Worker) for API mocking

**Estimated Time:** 12-16 hours

---

#### **Deliverable 3.2: End-to-End Testing** 🟡 **HIGH PRIORITY**

**Objective:** Test complete workflows from patient upload to care plan generation

**Test Scenarios:**

1. **Complete Care Plan Workflow:**
   - Upload VCF file → Extract mutations → Run pipeline → Display care plan
   - Verify all agents execute correctly
   - Verify all components display data
   - Verify export functionality

2. **Trial Intelligence Workflow:**
   - Create patient profile → Search trials → Filter → Generate dossiers
   - Verify mechanism fit scores
   - Verify dossier quality
   - Verify batch generation

3. **Dossier Management Workflow:**
   - Generate dossier → Browse → View detail → Export
   - Verify search and filtering
   - Verify navigation

**Estimated Time:** 8-10 hours

---

### **Phase 4: Enhancements** (2-3 days)

#### **Deliverable 4.1: Real-time Updates** 🟢 **MEDIUM PRIORITY**

**Objective:** Add real-time updates for dossier generation and pipeline status

**Requirements:**
- WebSocket/SSE connection for real-time updates
- Auto-refresh dossier list when new dossiers generated
- Real-time pipeline status updates
- Progress indicators for long-running operations

**Implementation:**
```jsx
// Add WebSocket hook
const useDossierUpdates = (patientId) => {
  const [dossiers, setDossiers] = useState([]);
  
  useEffect(() => {
    if (!patientId) return;
    
    const ws = new WebSocket(`ws://localhost:8000/ws/dossiers/${patientId}`);
    
    ws.onmessage = (event) => {
      const data = JSON.parse(event.data);
      if (data.type === 'dossier_generated') {
        setDossiers(prev => [...prev, data.dossier]);
      }
    };
    
    return () => ws.close();
  }, [patientId]);
  
  return dossiers;
};
```

**Estimated Time:** 6-8 hours

---

#### **Deliverable 4.2: Enhanced Error Handling** 🟢 **MEDIUM PRIORITY**

**Objective:** Improve error handling and recovery mechanisms

**Requirements:**
- Retry logic for failed API calls
- Partial failure handling (some agents succeed, others fail)
- User-friendly error messages
- Error recovery suggestions

**Estimated Time:** 4-6 hours

---

#### **Deliverable 4.3: Performance Optimizations** 🟢 **MEDIUM PRIORITY**

**Objective:** Optimize performance for large datasets

**Requirements:**
- Pagination for dossier lists
- Virtual scrolling for large trial lists
- Lazy loading for components
- Memoization for expensive computations

**Estimated Time:** 6-8 hours

---

#### **Deliverable 4.4: Enhanced Dossier Detail View** 🟢 **LOW PRIORITY**

**Objective:** Add interactive features to dossier detail view

**Requirements:**
- Expandable sections
- Tabbed navigation for different sections
- Print-friendly view
- Share functionality
- Comparison view (multiple dossiers)

**Estimated Time:** 8-10 hours

---

## 📊 Priority Matrix

| Deliverable | Priority | Time | Dependencies | Impact |
|-------------|----------|------|--------------|--------|
| **1.1: Resistance Playbook Component** | 🔴 CRITICAL | 4-6h | None | High - Blocks care plan display |
| **1.2: SAE Features Component** | 🔴 CRITICAL | 6-8h | None | High - Blocks care plan display |
| **2.1: Orchestrator Integration** | 🔴 CRITICAL | 8-10h | 1.1, 1.2 | Very High - Core functionality |
| **2.2: Orchestrator Trial Matching** | 🟡 HIGH | 6-8h | 2.1 | High - Better trial matching |
| **3.1: Component Testing** | 🟡 HIGH | 12-16h | 1.1, 1.2, 2.1 | High - Quality assurance |
| **3.2: End-to-End Testing** | 🟡 HIGH | 8-10h | 3.1 | High - Quality assurance |
| **4.1: Real-time Updates** | 🟢 MEDIUM | 6-8h | 2.1 | Medium - UX improvement |
| **4.2: Enhanced Error Handling** | 🟢 MEDIUM | 4-6h | 2.1 | Medium - Reliability |
| **4.3: Performance Optimizations** | 🟢 MEDIUM | 6-8h | 2.1 | Medium - Scalability |
| **4.4: Enhanced Dossier Detail** | 🟢 LOW | 8-10h | None | Low - Nice to have |

---

## 🎯 Success Criteria

### **Phase 1 Complete When:**
- ✅ Resistance Playbook component displays data correctly
- ✅ SAE Features component displays data correctly
- ✅ Both components integrated into UniversalCompleteCare
- ✅ No placeholder alerts visible

### **Phase 2 Complete When:**
- ✅ UniversalCompleteCare uses orchestrator pipeline
- ✅ File upload works (VCF/PDF/MAF)
- ✅ Status polling shows real-time progress
- ✅ All agent outputs displayed correctly
- ✅ UniversalTrialIntelligence uses orchestrator trial matching
- ✅ Mechanism fit scores displayed

### **Phase 3 Complete When:**
- ✅ All components have >80% test coverage
- ✅ End-to-end tests pass
- ✅ No critical bugs
- ✅ Performance benchmarks met

### **Phase 4 Complete When:**
- ✅ Real-time updates working
- ✅ Enhanced error handling implemented
- ✅ Performance optimizations complete
- ✅ Enhanced dossier detail view available

---

## 🔗 Related Documents

- **Orchestration Scope:** `.cursor/MOAT/ORCHESTRATION_SCOPE_SYNTHESIS.md`
- **Orchestrator Implementation:** `api/services/orchestrator/orchestrator.py`
- **Component Library:** `oncology-frontend/src/components/orchestrator/`
- **API Documentation:** `api/routers/orchestrate.py`

---

## 📝 Implementation Notes

### **Key Decisions:**

1. **Orchestrator Integration Strategy:**
   - Replace `/api/complete_care/v2` with `/api/orchestrate/full`
   - Maintain backward compatibility with existing patient profiles
   - Add file upload capability for data extraction

2. **State Management:**
   - Use orchestrator `PatientState` for single source of truth
   - Poll status endpoint for real-time updates
   - Cache state locally for offline access

3. **Component Architecture:**
   - Keep existing component structure
   - Add new components to `orchestrator/Analysis/` directory
   - Follow existing patterns (DrugRankingCard, BiomarkerCard, etc.)

4. **Testing Strategy:**
   - Unit tests for each component
   - Integration tests for API calls
   - End-to-end tests for complete workflows
   - Mock orchestrator responses for testing

---

**Document Status:** ✅ **AUDIT COMPLETE** - Deliverables defined  
**Last Updated:** January 28, 2025  
**Next Review:** After Phase 1 completion

---

## 🔗 Integration with Orchestration Plan

**This deliverable is part of Module 12 (UI Dashboard) in the MOAT Orchestration System.**

**See:** `.cursor/MOAT/ORCHESTRATION_SCOPE_SYNTHESIS.md` for complete orchestration plan

**Status in Orchestration Plan:**
- **Module 12 (UI Dashboard):** ⏳ IN PROGRESS - 30% complete
- **Priority:** 🟡 HIGH (moved from MEDIUM - critical for user adoption)
- **Estimated Time:** 9-13 days (4 phases)
- **Dependencies:** None (can proceed in parallel with other modules)

**Integration Points:**
- Phase 2 requires orchestrator endpoints (`/api/orchestrate/full`, `/api/orchestrate/status/{patient_id}`)
- Phase 2 requires file upload capability (VCF/PDF/MAF parsing)
- All phases can proceed independently of other orchestration modules

