# 🚀 Frontend MOAT Capabilities - Production Implementation Plan

**Date:** January 2025  
**Status:** 🔴 **IN PROGRESS**  
**Based on:** `.cursor/MOAT/FRONTEND_MOAT_CAPABILITIES_AUDIT.md`

---

## 🎯 Implementation Strategy

### **Phase 1: Orchestrator Integration** 🔴 **CRITICAL** (Week 1)
**Goal:** Replace legacy endpoint with orchestrator pipeline

### **Phase 2: Missing Components** 🔴 **HIGH** (Week 2)
**Goal:** Complete Resistance Playbook and SAE Features integration

### **Phase 3: Trial Matching Enhancement** 🟡 **MEDIUM** (Week 3)
**Goal:** Add mechanism fit scores and TRUE SAE indicators

### **Phase 4: Testing & Polish** 🟢 **ONGOING** (Week 4)
**Goal:** Comprehensive testing and UX improvements

---

## 📋 Detailed Implementation Tasks

### **Task 1.1: Update UniversalCompleteCare to Use Orchestrator** 🔴 **CRITICAL**

**File:** `oncology-coPilot/oncology-frontend/src/pages/UniversalCompleteCare.jsx`

**Changes:**
1. Replace `/api/complete_care/v2` with `/api/orchestrate/full`
2. Update request format to match `OrchestratePipelineRequest` schema
3. Update response handling to match `OrchestratePipelineResponse` schema
4. Map orchestrator state to existing component props

**Request Format:**
```typescript
{
  patient_id?: string,
  disease: string,
  mutations: Array<{
    gene: string,
    hgvs_p?: string,
    consequence?: string,
    zygosity?: string
  }>,
  treatment_line?: number,
  prior_therapies?: string[],
  current_regimen?: string,
  current_drug_class?: string,
  clinical_data?: {...},
  run_async?: boolean,
  skip_agents?: string[]
}
```

**Response Format:**
```typescript
{
  job_id: string,
  patient_id: string,
  status: string,
  phase: string,
  progress: number,
  alerts: Array<{...}>,
  biomarker_profile?: {...},
  resistance_prediction?: {...},
  drug_ranking?: Array<{...}>,
  trial_matches?: Array<{...}>,
  care_plan?: {...},
  // ... other agent outputs
}
```

**Estimated Time:** 4-6 hours

---

### **Task 1.2: Add File Upload Component** 🔴 **CRITICAL**

**File:** `oncology-coPilot/oncology-frontend/src/pages/UniversalCompleteCare.jsx`

**Changes:**
1. Import `PatientUpload` component (already exists)
2. Add file upload section before "Generate Plan" button
3. Handle file upload completion
4. Pass file to orchestrator API

**Implementation:**
```jsx
import { PatientUpload } from '../components/orchestrator/Patient/PatientUpload';

// Add before Generate Plan button
<PatientUpload
  onUploadComplete={(patientId) => {
    setPatientProfile({ ...patientProfile, patient_id: patientId });
    handleGeneratePlan();
  }}
  patientId={patientProfile?.patient_id}
/>
```

**Estimated Time:** 2-3 hours

---

### **Task 1.3: Create Pipeline Status Polling Hook** 🔴 **CRITICAL**

**File:** `oncology-coPilot/oncology-frontend/src/hooks/usePipelineStatus.js` (NEW)

**Purpose:** Poll `/api/orchestrate/status/{patient_id}` for real-time progress

**Implementation:**
```jsx
import { useState, useEffect, useRef } from 'react';

export const usePipelineStatus = (patientId, enabled = true) => {
  const [status, setStatus] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const intervalRef = useRef(null);

  useEffect(() => {
    if (!patientId || !enabled) return;

    const pollStatus = async () => {
      try {
        setLoading(true);
        const response = await fetch(`${API_ROOT}/api/orchestrate/status/${patientId}`);
        if (response.ok) {
          const data = await response.json();
          setStatus(data);
          setError(null);
          
          // Stop polling if complete or error
          if (data.phase === 'complete' || data.phase === 'error') {
            if (intervalRef.current) {
              clearInterval(intervalRef.current);
            }
          }
        }
      } catch (err) {
        setError(err.message);
      } finally {
        setLoading(false);
      }
    };

    // Poll immediately, then every 2 seconds
    pollStatus();
    intervalRef.current = setInterval(pollStatus, 2000);

    return () => {
      if (intervalRef.current) {
        clearInterval(intervalRef.current);
      }
    };
  }, [patientId, enabled]);

  return { status, loading, error };
};
```

**Estimated Time:** 2-3 hours

---

### **Task 1.4: Create PipelineStatusCard Component** 🔴 **CRITICAL**

**File:** `oncology-coPilot/oncology-frontend/src/components/orchestrator/Dashboard/PipelineStatusCard.jsx` (NEW)

**Purpose:** Display pipeline phase progress, agent execution status, alerts

**Implementation:**
```jsx
import React from 'react';
import {
  Card,
  CardContent,
  Typography,
  LinearProgress,
  Box,
  Chip,
  Alert,
  List,
  ListItem,
  ListItemText
} from '@mui/material';
import CheckCircleIcon from '@mui/icons-material/CheckCircle';
import ErrorIcon from '@mui/icons-material/Error';
import HourglassEmptyIcon from '@mui/icons-material/HourglassEmpty';

const PHASE_LABELS = {
  initialized: 'Initialized',
  extracting: 'Extracting Data',
  analyzing: 'Analyzing',
  ranking: 'Ranking Drugs',
  matching: 'Matching Trials',
  planning: 'Generating Care Plan',
  monitoring: 'Setting Up Monitoring',
  complete: 'Complete',
  error: 'Error'
};

export const PipelineStatusCard = ({ status }) => {
  if (!status) return null;

  const { phase, progress, alerts = [], agent_executions = [] } = status;

  return (
    <Card sx={{ mb: 3 }}>
      <CardContent>
        <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 2 }}>
          <Typography variant="h6">
            Pipeline Status
          </Typography>
          <Chip 
            label={PHASE_LABELS[phase] || phase}
            color={phase === 'complete' ? 'success' : phase === 'error' ? 'error' : 'primary'}
          />
        </Box>

        <Box sx={{ mb: 2 }}>
          <Box sx={{ display: 'flex', justifyContent: 'space-between', mb: 1 }}>
            <Typography variant="body2" color="text.secondary">
              Progress
            </Typography>
            <Typography variant="body2" fontWeight="bold">
              {Math.round(progress * 100)}%
            </Typography>
          </Box>
          <LinearProgress 
            variant="determinate" 
            value={progress * 100} 
            sx={{ height: 8, borderRadius: 4 }}
          />
        </Box>

        {alerts.length > 0 && (
          <Box sx={{ mb: 2 }}>
            {alerts.slice(-3).map((alert, idx) => (
              <Alert 
                key={idx} 
                severity={alert.severity || 'info'} 
                sx={{ mb: 1 }}
              >
                {alert.message}
              </Alert>
            ))}
          </Box>
        )}

        {agent_executions.length > 0 && (
          <Box>
            <Typography variant="subtitle2" gutterBottom>
              Agent Execution Status
            </Typography>
            <List dense>
              {agent_executions.slice(-5).map((exec, idx) => (
                <ListItem key={idx}>
                  <ListItemText
                    primary={exec.agent_name || exec.agent_id}
                    secondary={exec.status}
                  />
                  {exec.status === 'completed' && <CheckCircleIcon color="success" />}
                  {exec.status === 'failed' && <ErrorIcon color="error" />}
                  {exec.status === 'running' && <HourglassEmptyIcon color="primary" />}
                </ListItem>
              ))}
            </List>
          </Box>
        )}
      </CardContent>
    </Card>
  );
};
```

**Estimated Time:** 3-4 hours

---

### **Task 2.1: Verify ResistancePlaybook Integration** 🟡 **HIGH**

**File:** `oncology-coPilot/oncology-frontend/src/pages/UniversalCompleteCare.jsx`

**Changes:**
1. Check if `ResistancePlaybook` component receives data from orchestrator
2. Map orchestrator `resistance_prediction` to component props
3. Remove placeholder alert if data is available

**Current State:**
- Component exists at `components/ayesha/ResistancePlaybook.jsx`
- Expects `resistance_playbook` prop
- Currently shows placeholder alert

**Action:**
- Map `result.resistance_prediction` to `resistance_playbook` format
- Or update component to accept orchestrator format

**Estimated Time:** 2-3 hours

---

### **Task 2.2: Create/Adapt SAEFeatures Component** 🟡 **HIGH**

**File:** `oncology-coPilot/oncology-frontend/src/components/orchestrator/Analysis/SAEFeatures.jsx` (NEW or adapt existing)

**Changes:**
1. Check if existing `SAEFeaturesCard` can be adapted
2. Create new component if needed for orchestrator format
3. Display:
   - DNA repair capacity
   - Pathway burden scores (DDR, MAPK, PI3K, etc.)
   - DDR_bin score (if TRUE SAE available)
   - SAE provenance (TRUE SAE vs PROXY SAE)
   - Validation metrics

**Data Source:**
- From orchestrator: `state.sae_features` or `result.sae_features`
- Format: See `.cursor/SUPABASE_SCHEMA_UPDATES_MOAT.sql` for schema

**Estimated Time:** 4-6 hours

---

### **Task 3.1: Integrate Orchestrator Trial Matching** 🟡 **MEDIUM**

**File:** `oncology-coPilot/oncology-frontend/src/pages/UniversalTrialIntelligence.jsx`

**Changes:**
1. Add option to use orchestrator trial matching
2. Call `/api/orchestrate/full` with `skip_agents` to only run trial matching
3. Extract trials from orchestrator response
4. Display mechanism fit scores

**Estimated Time:** 4-6 hours

---

### **Task 3.2: Add Mechanism Fit Display** 🟡 **MEDIUM**

**Files:** 
- `oncology-coPilot/oncology-frontend/src/components/orchestrator/Analysis/TrialMatchesCard.jsx`
- `oncology-coPilot/oncology-frontend/src/pages/UniversalTrialIntelligence.jsx`

**Changes:**
1. Display `mechanism_fit_score` in trial cards
2. Show mechanism alignment breakdown (per-pathway)
3. Add TRUE SAE indicators (if available)
4. Show combined score formula (0.7×eligibility + 0.3×mechanism_fit)

**Estimated Time:** 3-4 hours

---

## 🔧 Technical Implementation Details

### **API Endpoint Mapping**

| Legacy Endpoint | New Orchestrator Endpoint | Notes |
|----------------|---------------------------|-------|
| `POST /api/complete_care/v2` | `POST /api/orchestrate/full` | Main pipeline execution |
| N/A | `GET /api/orchestrate/status/{patient_id}` | Status polling |
| N/A | `GET /api/patients/{patient_id}` | Full state retrieval |

### **Data Mapping**

**From Legacy Response to Orchestrator Response:**

```javascript
// Legacy format
{
  biomarker_intelligence: {...},
  resistance_prediction: {...},
  wiwfm: {...},
  trials: {...}
}

// Orchestrator format
{
  biomarker_profile: {...},
  resistance_prediction: {...},
  drug_ranking: [...],
  trial_matches: [...]
}
```

**Mapping Function:**
```javascript
const mapOrchestratorToLegacy = (orchestratorResponse) => {
  return {
    biomarker_intelligence: orchestratorResponse.biomarker_profile,
    resistance_prediction: orchestratorResponse.resistance_prediction,
    wiwfm: {
      drugs: orchestratorResponse.drug_ranking || [],
      evidence_tier: 'Supported' // Default
    },
    trials: {
      trials: orchestratorResponse.trial_matches || []
    },
    sae_features: orchestratorResponse.sae_features,
    resistance_playbook: orchestratorResponse.resistance_playbook,
    care_plan: orchestratorResponse.care_plan
  };
};
```

---

## 📝 File Structure

```
oncology-coPilot/oncology-frontend/src/
├── pages/
│   ├── UniversalCompleteCare.jsx          [MODIFY] - Orchestrator integration
│   └── UniversalTrialIntelligence.jsx      [MODIFY] - Trial matching enhancement
├── components/
│   ├── orchestrator/
│   │   ├── Analysis/
│   │   │   ├── SAEFeatures.jsx            [CREATE] - SAE features display
│   │   │   └── TrialMatchesCard.jsx      [MODIFY] - Add mechanism fit
│   │   ├── Dashboard/
│   │   │   └── PipelineStatusCard.jsx     [CREATE] - Pipeline status display
│   │   └── Patient/
│   │       └── PatientUpload.jsx          [EXISTS] - File upload
│   └── ayesha/
│       └── ResistancePlaybook.jsx         [VERIFY] - Check integration
└── hooks/
    └── usePipelineStatus.js               [CREATE] - Status polling hook
```

---

## ✅ Success Criteria

### **Phase 1 Complete When:**
- ✅ UniversalCompleteCare uses `/api/orchestrate/full`
- ✅ File upload works (VCF/PDF/MAF)
- ✅ Status polling shows real-time progress
- ✅ PipelineStatusCard displays correctly
- ✅ All agent outputs displayed

### **Phase 2 Complete When:**
- ✅ ResistancePlaybook displays data (no placeholder)
- ✅ SAEFeatures component displays data (no placeholder)
- ✅ Both integrated into UniversalCompleteCare

### **Phase 3 Complete When:**
- ✅ UniversalTrialIntelligence uses orchestrator trial matching
- ✅ Mechanism fit scores displayed
- ✅ TRUE SAE indicators shown

---

## 🚨 Risk Mitigation

### **Backward Compatibility**
- Keep legacy endpoint as fallback
- Add feature flag to toggle between endpoints
- Gradual migration path

### **Error Handling**
- Graceful degradation if orchestrator fails
- Clear error messages for users
- Retry logic for transient failures

### **Performance**
- Async pipeline execution option
- Status polling with exponential backoff
- Cache orchestrator responses

---

## 📊 Progress Tracking

| Task | Status | Time Spent | Blockers |
|------|--------|------------|----------|
| 1.1: Update UniversalCompleteCare | 🔴 In Progress | 0h | None |
| 1.2: Add File Upload | ⏳ Pending | 0h | None |
| 1.3: Create Status Hook | ⏳ Pending | 0h | None |
| 1.4: Create Status Card | ⏳ Pending | 0h | None |
| 2.1: Verify ResistancePlaybook | ⏳ Pending | 0h | None |
| 2.2: Create SAEFeatures | ⏳ Pending | 0h | None |
| 3.1: Integrate Trial Matching | ⏳ Pending | 0h | None |
| 3.2: Add Mechanism Fit | ⏳ Pending | 0h | None |

**Total Estimated Time:** 25-35 hours  
**Current Progress:** 0% (0/8 tasks)

---

*This plan will be updated as implementation progresses.*

