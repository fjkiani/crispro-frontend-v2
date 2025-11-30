# ⚔️ SPORADIC CANCER - AGENT ARCHITECTURE & WORKFLOWS

**Status**: ✅ **COMPLETE**  
**Duration**: 3-4 hours  
**Created**: January 14, 2025  
**Cycle**: SC-I5 (Cycle 9)

---

## **9.1 AGENT SYSTEM ARCHITECTURE**

### **9.1.1 Agent Router** (`routers/agents.py`)

**Purpose**: API endpoints for agent management, execution, and results

**Endpoints**:
- `POST /api/agents`: Create new agent
- `GET /api/agents`: List all agents for current user
- `GET /api/agents/{agent_id}`: Get agent details
- `PUT /api/agents/{agent_id}`: Update agent
- `DELETE /api/agents/{agent_id}`: Delete agent
- `POST /api/agents/{agent_id}/execute`: Execute agent
- `GET /api/agents/{agent_id}/results`: Get agent execution results

**Authentication**: Requires user authentication (via `get_optional_user` middleware)

---

### **9.1.2 Agent Types**

#### **1. PubMed Sentinel** (`pubmed_sentinel`):
- **Purpose**: Monitors PubMed for new papers
- **Execution**: Scheduled (hourly, daily, weekly, monthly)
- **Output**: Alerts for new publications matching criteria

#### **2. Trial Scout** (`trial_scout`):
- **Purpose**: Searches for new clinical trials
- **Execution**: Scheduled or on-demand
- **Output**: Trial matches based on patient profile

#### **3. VUS Vigil** (planned):
- **Purpose**: Resolve Variants of Uncertain Significance through continuous literature monitoring
- **Execution**: Scheduled
- **Output**: VUS resolution updates

---

### **9.1.3 Agent Services**

#### **Agent Manager** (`services/agent_manager.py`):
- **Purpose**: CRUD operations for agents
- **Operations**:
  - Create agent
  - Get user agents
  - Update agent
  - Delete agent
  - Get agent details

#### **Agent Executor** (`services/agent_executor.py`):
- **Purpose**: Agent execution logic
- **Operations**:
  - Execute agent
  - Get execution results
  - Handle agent-specific logic

#### **Agent Scheduler** (`services/agent_scheduler.py`):
- **Purpose**: Scheduled agent execution
- **Operations**:
  - Schedule agent runs
  - Manage execution frequency
  - Track execution history

---

## **9.2 AGENT WORKFLOWS FOR AYESHA**

### **9.2.1 Agent 3 - E2E Testing Mission**

**Status**: ⏸️ **ASSIGNED - PENDING EXECUTION**  
**Timeline**: 4-6 hours  
**Priority**: P1 (High value)

**Mission Objectives**:
1. Validate complete Ayesha demo workflow end-to-end
2. Create provider report template (Markdown + PDF)
3. Validate demo data quality
4. Test edge cases and error handling

**Task Breakdown**:
1. **E2E Smoke Testing** (2-3 hours):
   - Level 0 Flow: Quick Intake → Efficacy with PARP penalty
   - Level 2 Flow: NGS Upload → Re-analyze with PARP rescue
   - Clinical Trials: Verify filtering + badges
   - Edge Case Testing: Germline positive, TMB-high, missing data

2. **Provider Report Template** (2-3 hours):
   - Design Report Template (Markdown)
   - Create Report Generator Backend (`/api/reports/provider`)
   - Create Frontend Export Button

3. **Demo Data Validation** (30 min):
   - Verify test data quality
   - Cross-validate with literature

4. **Documentation** (30 min):
   - Create `E2E_SMOKE_TEST_RESULTS.md`

---

### **9.2.2 Agent Jr - Mission 4 Completion**

**Status**: ✅ **100% COMPLETE**  
**Mission**: Wire WIWFM (HypothesisValidator.jsx) to SporadicContext

**Deliverables**:
1. **BiomarkerSummaryWidget Component** ✅
   - Displays TMB, HRD, MSI status with color-coded chips
   - Shows data level (L0/L1/L2) and completeness score
   - Shows germline status

2. **HypothesisValidator.jsx Transformation** ✅
   - Removed ToolRunner wrapper
   - Added `useSporadic()` hook integration
   - Mutation input form with parsing
   - Efficacy API call with `getEfficacyPayload()` injection
   - Drug results display with cards
   - SporadicProvenanceCard below each drug

3. **Backend Router Update** ✅
   - Extract `germline_status` from request
   - Extract `tumor_context` from request
   - Pass both to `EfficacyRequest` model

---

### **9.2.3 Agent JR1 - Trial Seeding Mission**

**Status**: 🔄 **PRE-FLIGHT CHECK**  
**Mission**: Seed 200 trials for GTM (Go-To-Market)

**Tasks**:
1. **Schema Verification**:
   - Verify 28-column schema
   - Add missing base columns if needed

2. **Parsing Logic**:
   - PI name extraction (overallOfficials plural)
   - Contact information parsing
   - Primary endpoint extraction
   - Mechanism tagging
   - Biomarker extraction

3. **Database Seeding**:
   - SQLite seeding (28 columns)
   - AstraDB seeding (with vector embeddings)
   - Vector field at root level (not nested)

**Status**: ✅ **COMPLETE** - 525 trials seeded successfully

---

### **9.2.4 Agent Jr - Integration Testing Mission**

**Status**: ✅ **100% COMPLETE**  
**Mission**: Frontend integration testing for sporadic cancer workflow

**Deliverables**:
1. **Frontend Wired to Backend** ✅
   - SporadicContext integration
   - API calls with sporadic fields
   - Error handling and loading states

2. **Bugs Fixed** ✅
   - SOC bug fixed (duplicate code removed)
   - Loading states added
   - Error states added

3. **Provenance Cards** ✅
   - SporadicProvenanceCard displays
   - BiomarkerSummaryWidget displays
   - Complete audit trail

---

## **9.3 AGENT INTEGRATION WITH PLATFORM**

### **9.3.1 Agent → Platform Integration Points**

1. **Agent → Efficacy Orchestrator**:
   - Agents can trigger efficacy predictions
   - Agents can pass sporadic context (germline_status, tumor_context)

2. **Agent → Trial Matching**:
   - Agents can trigger trial searches
   - Agents can filter by sporadic status

3. **Agent → Evidence RAG**:
   - Agents can query literature
   - Agents can add variants to knowledge base

4. **Agent → Food Validator**:
   - Agents can validate compounds
   - Agents can pass treatment line context

---

### **9.3.2 Agent Execution Patterns**

#### **Scheduled Execution**:
```python
# Agent runs on schedule (hourly, daily, weekly, monthly)
agent = {
    "agent_type": "pubmed_sentinel",
    "run_frequency": "daily",
    "config": {
        "query": "ovarian cancer AND PARP inhibitor",
        "alert_threshold": 0.8
    }
}
```

#### **On-Demand Execution**:
```python
# Agent runs on-demand via API
POST /api/agents/{agent_id}/execute
{
    "trigger": "manual",
    "context": {
        "patient_id": "ayesha",
        "germline_status": "negative"
    }
}
```

---

### **9.3.3 Agent Results Storage**

**Storage**: Supabase (`public.agent_results` table)

**Schema**:
- `agent_id`: Foreign key to agents
- `execution_id`: UUID for this execution
- `status`: "pending" | "running" | "completed" | "failed"
- `results`: JSON blob with agent-specific results
- `created_at`: Timestamp
- `updated_at`: Timestamp

---

## **9.4 AGENT MISSIONS FOR AYESHA**

### **9.4.1 Mission 1: Trial Seeding (Agent JR1)**

**Objective**: Seed 200+ ovarian cancer trials for GTM

**Status**: ✅ **COMPLETE** - 525 trials seeded

**Deliverables**:
- SQLite database with 28 columns
- AstraDB collection with vector embeddings
- PI names extracted (overallOfficials plural fix)
- Mechanism tagging complete
- Biomarker extraction complete

---

### **9.4.2 Mission 2: WIWFM Integration (Agent Jr)**

**Objective**: Wire WIWFM to SporadicContext

**Status**: ✅ **COMPLETE**

**Deliverables**:
- BiomarkerSummaryWidget component
- HypothesisValidator.jsx transformation
- Backend router update (germline_status + tumor_context)

---

### **9.4.3 Mission 3: Integration Testing (Agent Jr)**

**Objective**: Frontend integration testing

**Status**: ✅ **COMPLETE**

**Deliverables**:
- Frontend wired to backend
- Bugs fixed (SOC duplicate code)
- Loading/error states added
- Provenance cards displayed

---

### **9.4.4 Mission 4: E2E Testing (Agent 3)**

**Objective**: Complete workflow validation

**Status**: ⏸️ **PENDING**

**Deliverables**:
- E2E smoke test results
- Provider report template
- Demo data validation
- Documentation

---

## **9.5 AGENT ARCHITECTURE SUMMARY**

### **9.5.1 Agent System Components**:

1. **Agent Router**: API endpoints for agent management
2. **Agent Manager**: CRUD operations for agents
3. **Agent Executor**: Agent execution logic
4. **Agent Scheduler**: Scheduled agent execution
5. **Agent Results**: Results storage in Supabase

### **9.5.2 Agent Types**:

1. **PubMed Sentinel**: Literature monitoring
2. **Trial Scout**: Trial matching
3. **VUS Vigil**: VUS resolution (planned)

### **9.5.3 Agent Workflows**:

1. **Scheduled**: Agents run on schedule (hourly, daily, weekly, monthly)
2. **On-Demand**: Agents run on-demand via API
3. **Event-Driven**: Agents triggered by events (planned)

### **9.5.4 Integration Points**:

1. **Agent → Efficacy Orchestrator**: Trigger efficacy predictions
2. **Agent → Trial Matching**: Trigger trial searches
3. **Agent → Evidence RAG**: Query literature
4. **Agent → Food Validator**: Validate compounds

---

**Status**: ✅ **CYCLE 9 COMPLETE** - Agent Architecture & Workflows  
**Next**: Cycle 10 (SC-I6) - Synthesis & Alignment

---



