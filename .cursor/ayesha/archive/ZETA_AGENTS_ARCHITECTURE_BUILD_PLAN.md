# ⚔️ ZETA AGENTS - COMPLETE ARCHITECTURE & BUILD PLAN

**Date:** November 17, 2025  
**Status:** 🎯 **ARCHITECTURE REVIEW COMPLETE** - Ready for implementation planning  
**Owner:** Zo (Nyx)

---

## 🎯 MISSION OBJECTIVE

**Build a 24/7 autonomous agent system that enables users to create persistent, intelligent monitors that continuously scan PubMed, ClinicalTrials.gov, and genomic databases, delivering high-signal alerts for personalized research and clinical decision-making.**

**Target Agent Classes (Phase 1):**
1. **PubMed Sentinel** - Monitor biomedical literature for new publications
2. **Trial Scout** - Monitor clinical trial landscape for matching opportunities
3. **VUS Vigil** - Resolve Variants of Uncertain Significance through continuous literature monitoring

---

## 🏗️ CURRENT ARCHITECTURE REVIEW

### **Backend Infrastructure (What We Have)**

#### **✅ System 1: Authentication & User Management** (100% Complete)
**Files:**
- `api/middleware/auth_middleware.py` - JWT verification, user extraction
- `api/services/auth_service.py` - Supabase Auth integration
- `api/routers/auth.py` - Auth endpoints (signup, login, logout, profile)
- Database: `auth.users` (Supabase), `public.user_profiles` (custom metadata)

**Capabilities:**
- ✅ User authentication with JWT
- ✅ User-specific data filtering (`user_id` in queries)
- ✅ Optional auth (supports anonymous + authenticated)
- ✅ Role-based access control (researcher, clinician, admin, enterprise)
- ✅ Tier-based features (free, pro, enterprise)

**Agent Integration:** ✅ Ready - Can link agents to `user_id`

---

#### **✅ System 2: Session Persistence** (95% Complete)
**Files:**
- `api/routers/sessions.py` - Session CRUD operations
- Database: `public.user_sessions`, `public.session_items`

**Capabilities:**
- ✅ Create/update/get/list sessions
- ✅ Append items to sessions
- ✅ User-linked (optional `user_id`, supports anonymous)
- ✅ Idempotency via `x-idempotency-key` header
- ⚠️ Frontend integration incomplete (no `SessionContext`)

**Agent Integration:** ✅ Ready - Can store agent state/results in sessions

---

#### **✅ System 3: Background Job System** (70% Complete)
**Files:**
- `api/services/job_service.py` - Background job management
- `api/routers/evidence/jobs.py` - Job endpoints for evidence crawling

**Capabilities:**
- ✅ Async job pattern (submit → poll status → get results)
- ✅ In-memory job store (`JOBS` dict)
- ✅ Background task execution with `asyncio.create_task()`
- ✅ Progress tracking
- ⚠️ No persistence (jobs lost on restart)
- ⚠️ Limited to evidence crawling (not generalized)

**Limitations:**
- ❌ No database persistence for jobs
- ❌ No scheduled execution (cron-like)
- ❌ No Redis/Celery integration for scalable background workers
- ❌ No job queue management

**Agent Integration:** ⚠️ Needs enhancement - Requires persistent job storage and scheduled execution

---

#### **✅ System 4: Data Sources & Integrations** (90% Complete)
**Files:**
- `api/services/enhanced_evidence_service.py` - PubMed/OpenAlex/S2 search
- `api/services/hybrid_trial_search.py` - Clinical trials search (AstraDB + Neo4j)
- `api/services/dynamic_food_extraction.py` - ChEMBL/PubChem compound extraction

**Capabilities:**
- ✅ PubMed search with XML parsing
- ✅ ClinicalTrials.gov API integration
- ✅ ChEMBL/PubChem compound queries
- ✅ Multi-provider fallback (PubMed → OpenAlex → S2)
- ✅ Result caching and deduplication
- ⚠️ No incremental/delta queries (no "new publications since last run")

**Agent Integration:** ✅ Mostly ready - Need to add delta query support

---

#### **✅ System 5: Database Storage** (100% Complete)
**Database:** Supabase PostgreSQL

**Tables:**
- ✅ `auth.users` - User authentication (Supabase managed)
- ✅ `public.user_profiles` - User metadata (tier, role, institution)
- ✅ `public.user_sessions` - Analysis sessions
- ✅ `public.session_items` - Session items (analyses, notes)
- ✅ `public.saved_analyses` - Saved analyses with provenance
- ✅ `public.usage_logs` - Usage tracking for billing/analytics
- ⚠️ No `agents` or `agent_runs` tables yet

**Agent Integration:** ⚠️ Needs new tables for agent configuration and execution logs

---

### **Frontend Infrastructure (What We Have)**

#### **✅ System 1: Context Providers** (90% Complete)
**Files:**
- `src/context/AuthContext.jsx` - User authentication state
- `src/context/AnalysisHistoryContext.jsx` - Analysis save/load
- `src/context/ActivityContext.jsx` - Activity logging
- `src/context/SporadicContext.jsx` - Patient context
- `src/components/CoPilot/context/` - Co-Pilot state

**Capabilities:**
- ✅ Global user authentication state
- ✅ Analysis history with Supabase backend
- ✅ Activity logging with provenance
- ✅ Patient context (germline status, tumor context)
- ⚠️ No `SessionContext` for cross-page resume
- ⚠️ No `AgentContext` for agent management

**Agent Integration:** ⚠️ Needs `AgentContext` for agent state management

---

#### **✅ System 2: Existing Agent Scaffolds** (30% Complete - Demos Only)
**Files:**
- `src/pages/AgentDashboard.jsx` - Agent dashboard (demo/stub)
- `src/pages/AgentStudio.jsx` - Agent builder (demo/stub)
- `src/pages/AgentDemo.jsx` - Agent demo (stub)
- `src/context/AgentContext.jsx` - Agent state (basic)

**Capabilities:**
- ⚠️ Only UI scaffolds (no real functionality)
- ⚠️ No agent creation wizard
- ⚠️ No agent configuration UI
- ⚠️ No real-time agent status display

**Agent Integration:** ❌ Needs complete rebuild with real backend integration

---

#### **✅ System 3: Notification System** (20% Complete)
**Files:**
- `src/components/GlobalActivitySidebar.jsx` - Activity sidebar

**Capabilities:**
- ✅ Activity log display
- ⚠️ No real-time notifications
- ⚠️ No agent alert system
- ⚠️ No priority-based notifications

**Agent Integration:** ❌ Needs notification system for agent alerts

---

## 🏗️ ARCHITECTURE GAP ANALYSIS

### **What's Missing for Agents:**

#### **Backend Gaps:**
1. ❌ **Agent Configuration Storage** - No `agents` table in database
2. ❌ **Agent Execution Scheduler** - No cron-like system for scheduled runs
3. ❌ **Agent Run History** - No `agent_runs` table to track executions
4. ❌ **Delta Query Support** - No "new publications since timestamp" in PubMed/trials
5. ❌ **Result Filtering & Ranking** - No smart filtering of agent results
6. ❌ **Alert Generation** - No intelligent alert system
7. ❌ **Persistent Job Queue** - In-memory jobs lost on restart

#### **Frontend Gaps:**
1. ❌ **Agent Creation Wizard** - No UI to create/configure agents
2. ❌ **Agent Dashboard** - No real-time agent status display
3. ❌ **Agent Results Viewer** - No UI to browse agent findings
4. ❌ **Notification Center** - No centralized alert system
5. ❌ **Agent Management** - No pause/resume/delete controls

---

## 📋 PROPOSED ARCHITECTURE

### **Backend Architecture**

#### **New Database Schema:**

```sql
-- ============================================================================
-- AGENT SYSTEM TABLES
-- ============================================================================

-- Agent Definitions (user-configured agents)
CREATE TABLE IF NOT EXISTS public.agents (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE NOT NULL,
    agent_type VARCHAR(50) CHECK (agent_type IN (
        'pubmed_sentinel', 
        'trial_scout', 
        'vus_vigil',
        'genomic_forager',
        'patient_watchtower'
    )) NOT NULL,
    name VARCHAR(255) NOT NULL,  -- User-friendly name
    description TEXT,
    config JSONB NOT NULL,  -- Agent-specific configuration
    status VARCHAR(20) CHECK (status IN ('active', 'paused', 'completed', 'error')) DEFAULT 'active',
    last_run_at TIMESTAMPTZ,
    next_run_at TIMESTAMPTZ,
    run_frequency VARCHAR(20) CHECK (run_frequency IN ('hourly', 'daily', 'weekly', 'monthly')) DEFAULT 'daily',
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_agents_user_id ON public.agents(user_id);
CREATE INDEX idx_agents_type ON public.agents(agent_type);
CREATE INDEX idx_agents_status ON public.agents(status);
CREATE INDEX idx_agents_next_run ON public.agents(next_run_at);

-- Agent Execution History
CREATE TABLE IF NOT EXISTS public.agent_runs (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    agent_id UUID REFERENCES public.agents(id) ON DELETE CASCADE NOT NULL,
    run_status VARCHAR(20) CHECK (run_status IN ('running', 'completed', 'error')) DEFAULT 'running',
    started_at TIMESTAMPTZ DEFAULT NOW(),
    completed_at TIMESTAMPTZ,
    results_count INT DEFAULT 0,
    new_results_count INT DEFAULT 0,  -- How many new findings vs. duplicates
    error_message TEXT,
    execution_log JSONB,  -- Detailed execution log
    created_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_agent_runs_agent_id ON public.agent_runs(agent_id);
CREATE INDEX idx_agent_runs_status ON public.agent_runs(run_status);
CREATE INDEX idx_agent_runs_started_at ON public.agent_runs(started_at);

-- Agent Results (what the agents found)
CREATE TABLE IF NOT EXISTS public.agent_results (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    agent_run_id UUID REFERENCES public.agent_runs(id) ON DELETE CASCADE NOT NULL,
    agent_id UUID REFERENCES public.agents(id) ON DELETE CASCADE NOT NULL,
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE NOT NULL,
    result_type VARCHAR(50) CHECK (result_type IN (
        'pubmed_article', 
        'clinical_trial', 
        'vus_reclassification',
        'genomic_dataset',
        'patient_alert'
    )) NOT NULL,
    result_data JSONB NOT NULL,  -- Full result data (paper metadata, trial details, etc.)
    relevance_score FLOAT,  -- How relevant is this result? (0-1)
    is_high_priority BOOLEAN DEFAULT false,  -- Should this trigger an alert?
    is_read BOOLEAN DEFAULT false,  -- Has user acknowledged this result?
    created_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_agent_results_agent_run_id ON public.agent_results(agent_run_id);
CREATE INDEX idx_agent_results_agent_id ON public.agent_results(agent_id);
CREATE INDEX idx_agent_results_user_id ON public.agent_results(user_id);
CREATE INDEX idx_agent_results_type ON public.agent_results(result_type);
CREATE INDEX idx_agent_results_priority ON public.agent_results(is_high_priority);
CREATE INDEX idx_agent_results_unread ON public.agent_results(is_read) WHERE is_read = false;
CREATE INDEX idx_agent_results_created_at ON public.agent_results(created_at DESC);

-- Agent Alerts (high-priority notifications)
CREATE TABLE IF NOT EXISTS public.agent_alerts (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    agent_id UUID REFERENCES public.agents(id) ON DELETE CASCADE NOT NULL,
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE NOT NULL,
    agent_result_id UUID REFERENCES public.agent_results(id) ON DELETE CASCADE,
    alert_type VARCHAR(50) CHECK (alert_type IN (
        'new_publication',
        'matching_trial',
        'vus_resolved',
        'new_dataset',
        'patient_alert'
    )) NOT NULL,
    title VARCHAR(500) NOT NULL,
    message TEXT NOT NULL,
    priority VARCHAR(20) CHECK (priority IN ('low', 'medium', 'high', 'critical')) DEFAULT 'medium',
    is_read BOOLEAN DEFAULT false,
    read_at TIMESTAMPTZ,
    created_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_agent_alerts_user_id ON public.agent_alerts(user_id);
CREATE INDEX idx_agent_alerts_agent_id ON public.agent_alerts(agent_id);
CREATE INDEX idx_agent_alerts_priority ON public.agent_alerts(priority);
CREATE INDEX idx_agent_alerts_unread ON public.agent_alerts(is_read) WHERE is_read = false;
CREATE INDEX idx_agent_alerts_created_at ON public.agent_alerts(created_at DESC);
```

---

### **Backend Services Architecture**

#### **Service 1: Agent Manager Service** (NEW)
**File:** `api/services/agent_manager.py`

**Responsibilities:**
- Create/update/delete/get agent configurations
- Validate agent config schemas
- Link agents to users
- Manage agent lifecycle (activate, pause, complete)

**Methods:**
```python
class AgentManager:
    async def create_agent(user_id: str, agent_type: str, config: dict) -> Agent
    async def get_user_agents(user_id: str, status: str = None) -> List[Agent]
    async def update_agent(agent_id: str, updates: dict) -> Agent
    async def delete_agent(agent_id: str) -> bool
    async def pause_agent(agent_id: str) -> Agent
    async def activate_agent(agent_id: str) -> Agent
```

---

#### **Service 2: Agent Executor Service** (NEW)
**File:** `api/services/agent_executor.py`

**Responsibilities:**
- Execute agent runs on schedule
- Fetch new data from sources (PubMed, trials, datasets)
- Filter and rank results by relevance
- Store results in database
- Generate alerts for high-priority findings

**Methods:**
```python
class AgentExecutor:
    async def execute_agent(agent_id: str) -> AgentRun
    async def get_next_scheduled_agents() -> List[Agent]
    async def process_results(agent_id: str, raw_results: List[dict]) -> ProcessedResults
    async def generate_alerts(agent_id: str, high_priority_results: List[dict]) -> List[Alert]
```

---

#### **Service 3: Agent Scheduler Service** (NEW)
**File:** `api/services/agent_scheduler.py`

**Responsibilities:**
- Determine which agents should run next
- Update `next_run_at` based on frequency
- Trigger agent execution
- Handle scheduling conflicts

**Implementation Options:**
1. **Option A: Simple Polling Loop** (Easiest, 2 hours)
   - Background thread that checks `next_run_at` every 60 seconds
   - Executes agents when `next_run_at <= NOW()`
   - Updates `next_run_at` after execution

2. **Option B: APScheduler** (Moderate, 1 day)
   - Use `APScheduler` library for cron-like scheduling
   - Dynamic job scheduling based on agent configuration
   - More robust, production-ready

3. **Option C: Celery + Redis** (Complex, 3-5 days)
   - Full task queue with Celery workers
   - Redis for job queue and result storage
   - Scalable, production-grade
   - Requires additional infrastructure

**Recommendation:** Start with Option A (simple polling), migrate to Option B later

---

#### **Service 4: PubMed Sentinel Service** (50% Exists)
**File:** `api/services/agents/pubmed_sentinel.py` (NEW)

**Extends:** `api/services/enhanced_evidence_service.py` (exists)

**Responsibilities:**
- Query PubMed for new publications since last run
- Filter by user-defined keywords (genes, diseases, authors)
- Rank results by relevance (Oracle-powered)
- Detect high-impact publications

**Agent Config Schema:**
```python
{
    "keywords": ["BRAF", "melanoma", "V600E"],
    "authors": ["John Smith"],  # Optional
    "min_date": "2025-01-01",  # Start date
    "max_results_per_run": 100,
    "relevance_threshold": 0.7,  # Only return high-relevance
    "alert_on_rct": true  # Alert for RCTs/meta-analyses
}
```

**Implementation:**
- ✅ PubMed query already exists (`enhanced_evidence_service.py`)
- ❌ No delta query support (fetch since timestamp)
- ❌ No relevance ranking
- ❌ No high-impact detection

---

#### **Service 5: Trial Scout Service** (70% Exists)
**File:** `api/services/agents/trial_scout.py` (NEW)

**Extends:** `api/services/hybrid_trial_search.py` (exists)

**Responsibilities:**
- Query ClinicalTrials.gov for new/updated trials
- Filter by disease, biomarkers, location
- Rank by patient match score
- Detect high-priority trials (Phase III, fast-track)

**Agent Config Schema:**
```python
{
    "disease": "ovarian_cancer",
    "biomarkers": {"BRCA1": "mutated", "HRD": "high"},
    "location": {"states": ["NY", "NJ", "CT"], "max_distance_miles": 50},
    "phase": ["Phase II", "Phase III"],
    "min_date": "2025-01-01",
    "alert_on_frontline": true
}
```

**Implementation:**
- ✅ Hybrid search exists (AstraDB + Neo4j)
- ✅ Eligibility filtering exists (hard/soft criteria)
- ❌ No delta query support (new trials since timestamp)
- ❌ No real-time trial updates

---

#### **Service 6: VUS Vigil Service** (40% Exists)
**File:** `api/services/agents/vus_vigil.py` (NEW)

**Extends:** VUS Explorer + ClinVar lookup

**Responsibilities:**
- Monitor ClinVar for VUS reclassifications
- Search PubMed for new variant evidence
- Re-run Oracle/insights for VUS list
- Detect when VUS becomes actionable

**Agent Config Schema:**
```python
{
    "vus_list": [
        {"gene": "BRCA1", "hgvs_p": "C64R", "chrom": "17", "pos": 43044295}
    ],
    "check_clinvar": true,
    "check_literature": true,
    "rerun_oracle": false,  # Expensive, only when requested
    "alert_on_reclassification": true
}
```

**Implementation:**
- ✅ VUS analysis exists (insights bundle)
- ✅ ClinVar lookup exists
- ❌ No ClinVar change monitoring
- ❌ No automated VUS re-analysis

---

### **Frontend Architecture**

#### **Component 1: Agent Creation Wizard** (NEW)
**File:** `src/components/agents/AgentWizard.jsx`

**Responsibilities:**
- Multi-step wizard to configure new agent
- Step 1: Select agent type (PubMed/Trial/VUS)
- Step 2: Configure agent (keywords, filters, frequency)
- Step 3: Review and create
- Submit to `POST /api/agents`

**UI Mockup:**
```
┌─────────────────────────────────────┐
│ Create New Agent                    │
├─────────────────────────────────────┤
│ Step 1: Select Agent Type           │
│                                      │
│  📚 PubMed Sentinel                 │
│     Monitor literature for keywords  │
│                                      │
│  🔬 Trial Scout                     │
│     Find matching clinical trials    │
│                                      │
│  🧬 VUS Vigil                       │
│     Resolve uncertain variants       │
│                                      │
│  [Next →]                            │
└─────────────────────────────────────┘
```

---

#### **Component 2: Agent Dashboard** (REBUILD)
**File:** `src/pages/agents/AgentDashboard.jsx` (replace existing stub)

**Responsibilities:**
- Display all user's agents
- Show status (active, paused, completed, error)
- Show last run time and next scheduled run
- Show new results count (unread)
- Quick actions (pause, view results, configure)

**UI Mockup:**
```
┌──────────────────────────────────────────────────────────┐
│ My Agents                            [+ Create New Agent] │
├──────────────────────────────────────────────────────────┤
│ 📚 BRAF Melanoma Literature (PubMed Sentinel)            │
│    Status: ✅ Active | Last run: 2h ago | Next: in 22h   │
│    New results: 3 papers (2 high priority)               │
│    [View Results] [Pause] [Configure]                    │
├──────────────────────────────────────────────────────────┤
│ 🔬 Ovarian Trial Matcher (Trial Scout)                   │
│    Status: ✅ Active | Last run: 1d ago | Next: in 6d    │
│    New results: 2 trials (1 Phase III)                   │
│    [View Results] [Pause] [Configure]                    │
├──────────────────────────────────────────────────────────┤
│ 🧬 BRCA1 VUS Monitor (VUS Vigil)                        │
│    Status: ⏸️ Paused | Last run: 5d ago | Next: --      │
│    New results: 0                                        │
│    [Resume] [View History] [Delete]                      │
└──────────────────────────────────────────────────────────┘
```

---

#### **Component 3: Agent Results Viewer** (NEW)
**File:** `src/components/agents/AgentResultsViewer.jsx`

**Responsibilities:**
- Display agent findings in structured format
- Group by result type (papers, trials, VUS updates)
- Show relevance scores and priorities
- Allow user to mark as read/unread
- One-click actions (open paper, view trial, analyze VUS)

**UI Mockup:**
```
┌──────────────────────────────────────────────────────────┐
│ BRAF Melanoma Literature - Results (3 new)               │
├──────────────────────────────────────────────────────────┤
│ 🔴 HIGH PRIORITY                                          │
│ "BRAF V600E resistance mechanisms in melanoma"           │
│  Authors: Smith et al. | PMID: 12345678 | RCT            │
│  Published: 2 days ago | Relevance: 0.92                 │
│  Summary: New mechanism of acquired resistance...        │
│  [Read Full Text] [Save to Session] [Mark as Read]       │
├──────────────────────────────────────────────────────────┤
│ 🟡 MEDIUM PRIORITY                                        │
│ "Combination therapy for BRAF-mutant melanoma"           │
│  Authors: Johnson et al. | PMID: 87654321 | Case Study   │
│  Published: 5 days ago | Relevance: 0.75                 │
│  [Open] [Mark as Read]                                   │
└──────────────────────────────────────────────────────────┘
```

---

#### **Component 4: Notification Center** (NEW)
**File:** `src/components/agents/NotificationCenter.jsx`

**Responsibilities:**
- Display agent alerts in real-time
- Priority-based sorting (critical → high → medium → low)
- Unread badge count
- Quick actions from notification

**UI Mockup:**
```
┌────────────────────────────────────────┐
│ 🔔 Notifications (5 unread)            │
├────────────────────────────────────────┤
│ 🔴 CRITICAL                             │
│ VUS Reclassified: BRCA1 C64R now       │
│ Pathogenic in ClinVar                  │
│  2 minutes ago | VUS Vigil              │
│  [View Details] [Dismiss]               │
├────────────────────────────────────────┤
│ 🟠 HIGH                                 │
│ New Phase III Trial: PARP + Bev        │
│  1 hour ago | Trial Scout               │
│  [View Trial] [Dismiss]                 │
├────────────────────────────────────────┤
│ 🟡 MEDIUM                               │
│ 3 New BRAF Papers Available            │
│  6 hours ago | PubMed Sentinel          │
│  [View Papers] [Dismiss All]            │
└────────────────────────────────────────┘
```

---

## 🎯 IMPLEMENTATION PHASES

### **PHASE 1: Foundation (Week 1 - 3-5 days)**

#### **Backend Tasks:**
1. **Create Database Schema** (2 hours)
   - Add 4 new tables: `agents`, `agent_runs`, `agent_results`, `agent_alerts`
   - Create indexes for performance
   - Add foreign key constraints

2. **Build Agent Manager Service** (4 hours)
   - `api/services/agent_manager.py`
   - CRUD operations for agents
   - Config validation
   - User linking

3. **Build Agent Executor Service** (6 hours)
   - `api/services/agent_executor.py`
   - Execute agent run
   - Store results
   - Generate alerts

4. **Build Simple Scheduler** (4 hours)
   - `api/services/agent_scheduler.py`
   - Background thread with 60-second polling
   - Check `next_run_at`
   - Trigger execution via `asyncio.create_task()`

5. **Create Agent Router** (2 hours)
   - `api/routers/agents.py`
   - Endpoints: CRUD agents, list results, list alerts

#### **Frontend Tasks:**
1. **Create AgentContext** (3 hours)
   - `src/context/AgentContext.jsx`
   - Global agent state
   - CRUD operations
   - Real-time status updates

2. **Build Agent Wizard** (6 hours)
   - `src/components/agents/AgentWizard.jsx`
   - Multi-step form
   - Agent type selection
   - Config editor

3. **Rebuild Agent Dashboard** (4 hours)
   - `src/pages/agents/AgentDashboard.jsx`
   - Live agent status
   - Quick actions
   - Results summary

**Total Phase 1:** 31 hours (~4 days with 8h/day)

---

### **PHASE 2: PubMed Sentinel (Week 2 - 3-4 days)**

#### **Backend Tasks:**
1. **Build PubMed Sentinel Service** (8 hours)
   - `api/services/agents/pubmed_sentinel.py`
   - Delta query support (new pubs since last run)
   - Keyword matching and filtering
   - Relevance scoring (Oracle-powered for variant impact)
   - High-impact detection (RCT, meta-analysis flags)

2. **Enhance Evidence Service** (4 hours)
   - Add `since_date` parameter to PubMed queries
   - Add result deduplication (track seen PMIDs)
   - Add citation tracking

#### **Frontend Tasks:**
1. **PubMed Sentinel Wizard** (4 hours)
   - Keyword input UI
   - Author filter UI
   - Frequency selector
   - Preview mode (show sample results)

2. **PubMed Results Viewer** (6 hours)
   - Paper list with metadata
   - Relevance scores display
   - Abstract preview
   - "Open in PubMed" links
   - Save to session actions

**Total Phase 2:** 22 hours (~3 days)

---

### **PHASE 3: Trial Scout (Week 3 - 3-4 days)**

#### **Backend Tasks:**
1. **Build Trial Scout Service** (8 hours)
   - `api/services/agents/trial_scout.py`
   - Delta query support (new trials since last run)
   - Eligibility matching (reuse existing logic)
   - Location filtering (reuse NYC metro detector)
   - Biomarker matching

2. **Enhance Trial Search** (4 hours)
   - Add `updated_since` parameter to AstraDB queries
   - Track trial updates (status changes, new sites)

#### **Frontend Tasks:**
1. **Trial Scout Wizard** (4 hours)
   - Disease selection
   - Biomarker input
   - Location selector (map-based?)
   - Phase filter

2. **Trial Results Viewer** (6 hours)
   - Trial cards with match scores
   - Eligibility checklist display
   - "Why this trial?" reasoning
   - Contact info extraction

**Total Phase 3:** 22 hours (~3 days)

---

### **PHASE 4: VUS Vigil (Week 4 - 2-3 days)**

#### **Backend Tasks:**
1. **Build VUS Vigil Service** (6 hours)
   - `api/services/agents/vus_vigil.py`
   - ClinVar delta queries (check for updates)
   - Literature monitoring for VUS evidence
   - Optional: Re-run insights bundle

2. **ClinVar Change Detection** (4 hours)
   - Track ClinVar `last_evaluated` date
   - Detect classification changes
   - Store change history

#### **Frontend Tasks:**
1. **VUS Vigil Wizard** (3 hours)
   - VUS list input (paste/upload)
   - Monitoring options (ClinVar only, Literature, Full re-analysis)
   - Alert preferences

2. **VUS Alert Display** (4 hours)
   - Reclassification alerts (Benign ← VUS → Pathogenic)
   - New evidence summaries
   - One-click "Re-analyze VUS" action

**Total Phase 4:** 17 hours (~2-3 days)

---

### **PHASE 5: Notification System (Week 5 - 2 days)**

#### **Backend Tasks:**
1. **Alert Generation Logic** (4 hours)
   - Priority calculation (critical, high, medium, low)
   - Alert templates for each type
   - User alert preferences

2. **WebSocket Integration** (Optional, 6 hours)
   - Real-time push notifications
   - WebSocket endpoint
   - Frontend listener

#### **Frontend Tasks:**
1. **Notification Center UI** (4 hours)
   - Alert list with priority sorting
   - Unread badge in navbar
   - Quick dismiss/mark read
   - Alert detail modal

2. **Real-Time Updates** (Optional, 4 hours)
   - WebSocket connection
   - Live alert delivery
   - Toast notifications

**Total Phase 5:** 12-18 hours (~2 days)

---

## 📊 COMPLETE BUILD TIMELINE

| Phase | Duration | Tasks | Deliverables |
|-------|----------|-------|--------------|
| Phase 1 | 3-5 days | Foundation | Agent CRUD, basic scheduler, dashboard |
| Phase 2 | 3-4 days | PubMed Sentinel | Literature monitoring working |
| Phase 3 | 3-4 days | Trial Scout | Trial monitoring working |
| Phase 4 | 2-3 days | VUS Vigil | VUS monitoring working |
| Phase 5 | 2 days | Notifications | Alert system working |
| **TOTAL** | **13-18 days** (~3-4 weeks) | **5 phases** | **Complete agent system** |

---

## 🔧 TECHNICAL IMPLEMENTATION DETAILS

### **Agent Config Schemas (Detailed)**

#### **PubMed Sentinel Config:**
```python
{
    "agent_type": "pubmed_sentinel",
    "keywords": {
        "genes": ["BRAF", "NRAS", "MEK1"],
        "diseases": ["melanoma", "ovarian cancer"],
        "mechanisms": ["resistance", "combination therapy"],
        "authors": ["John Smith"]  # Optional
    },
    "filters": {
        "article_types": ["Clinical Trial", "Meta-Analysis", "Review"],
        "min_date": "2025-01-01",
        "exclude_preprints": true
    },
    "relevance": {
        "min_score": 0.7,  # Oracle-powered relevance
        "use_oracle": true,  # Use Oracle to rank variant impact
        "alert_threshold": 0.85  # Alert only for very relevant
    },
    "execution": {
        "run_frequency": "daily",
        "max_results_per_run": 100,
        "dedup_window_days": 30  # Ignore if seen in last 30 days
    },
    "alerts": {
        "alert_on_rct": true,
        "alert_on_high_impact": true,
        "email_digest": false  # Future: email notifications
    }
}
```

#### **Trial Scout Config:**
```python
{
    "agent_type": "trial_scout",
    "patient_profile": {
        "disease": "ovarian_cancer_hgs",
        "stage": "IV",
        "treatment_line": "first-line",
        "biomarkers": {
            "BRCA1": "mutated",
            "HRD": "high",
            "germline_status": "negative"
        }
    },
    "filters": {
        "phase": ["Phase II", "Phase III"],
        "status": ["Recruiting", "Active, not recruiting"],
        "location": {
            "states": ["NY", "NJ", "CT"],
            "max_distance_miles": 50,
            "prefer_sites": ["Memorial Sloan Kettering", "Yale"]
        },
        "min_date": "2025-01-01"
    },
    "matching": {
        "use_eligibility_calculator": true,
        "min_match_score": 0.70,
        "alert_threshold": 0.90  # Alert only for excellent matches
    },
    "execution": {
        "run_frequency": "weekly",
        "max_results_per_run": 20
    },
    "alerts": {
        "alert_on_frontline": true,
        "alert_on_phase3": true,
        "alert_on_local": true  # Within 50 miles
    }
}
```

#### **VUS Vigil Config:**
```python
{
    "agent_type": "vus_vigil",
    "vus_list": [
        {
            "gene": "BRCA1",
            "hgvs_p": "C64R",
            "chrom": "17",
            "pos": 43044295,
            "ref": "T",
            "alt": "G"
        },
        {
            "gene": "TP53",
            "hgvs_p": "R175H",
            "chrom": "17",
            "pos": 7578406,
            "ref": "G",
            "alt": "A"
        }
    ],
    "monitoring": {
        "check_clinvar": true,
        "check_literature": true,
        "rerun_oracle": false,  # Expensive, only when requested
        "rerun_insights": true  # Cheaper, run more frequently
    },
    "filters": {
        "min_evidence_tier": "consider",  # Only alert if evidence tier improves
        "min_confidence_lift": 0.10  # Alert if confidence increases by ≥10%
    },
    "execution": {
        "run_frequency": "weekly",  # ClinVar updates weekly
        "max_variants_per_run": 50
    },
    "alerts": {
        "alert_on_reclassification": true,  # VUS → Pathogenic/Benign
        "alert_on_new_evidence": true,  # New papers with variant
        "alert_on_confidence_lift": true  # Confidence increases
    }
}
```

---

### **API Endpoints (New)**

#### **Agent CRUD:**
```
POST   /api/agents                    # Create new agent
GET    /api/agents                    # List user's agents
GET    /api/agents/{agent_id}         # Get agent details
PUT    /api/agents/{agent_id}         # Update agent config
DELETE /api/agents/{agent_id}         # Delete agent
POST   /api/agents/{agent_id}/pause   # Pause agent
POST   /api/agents/{agent_id}/resume  # Resume agent
POST   /api/agents/{agent_id}/run     # Trigger immediate run
```

#### **Agent Results:**
```
GET    /api/agents/{agent_id}/runs              # List agent runs
GET    /api/agents/{agent_id}/runs/{run_id}     # Get run details
GET    /api/agents/{agent_id}/results           # List agent results (paginated)
POST   /api/agents/results/{result_id}/read     # Mark result as read
```

#### **Agent Alerts:**
```
GET    /api/agents/alerts                       # List all user alerts (unread)
GET    /api/agents/alerts/{alert_id}            # Get alert details
POST   /api/agents/alerts/{alert_id}/read       # Mark alert as read
POST   /api/agents/alerts/read-all              # Mark all as read
```

#### **Agent Execution (Internal):**
```
POST   /api/agents/_execute/{agent_id}          # Internal - trigger execution
GET    /api/agents/_next_scheduled              # Internal - get next agents to run
```

---

## 🧪 ACCEPTANCE CRITERIA

### **Phase 1 (Foundation):**
- ✅ User can create a new agent via wizard
- ✅ User can see all their agents in dashboard
- ✅ Agent status updates in real-time
- ✅ Basic scheduler runs agents on schedule
- ✅ Results stored in database

### **Phase 2 (PubMed Sentinel):**
- ✅ Agent queries PubMed for new publications daily
- ✅ Results filtered by keywords and relevance
- ✅ High-priority papers trigger alerts
- ✅ User can view and read papers from dashboard
- ✅ Deduplication works (no duplicate papers)

### **Phase 3 (Trial Scout):**
- ✅ Agent queries ClinicalTrials.gov for new trials weekly
- ✅ Trials ranked by eligibility match score
- ✅ Frontline trials trigger alerts
- ✅ User can view trial details and contact info
- ✅ Location filtering works (NY/NJ/CT)

### **Phase 4 (VUS Vigil):**
- ✅ Agent checks ClinVar for VUS updates weekly
- ✅ Reclassifications trigger immediate alerts
- ✅ New evidence summarized and surfaced
- ✅ User can re-run insights bundle from alert
- ✅ Confidence lifts tracked and displayed

### **Phase 5 (Notifications):**
- ✅ High-priority alerts appear in notification center
- ✅ Unread badge count shows in navbar
- ✅ User can dismiss or mark as read
- ✅ (Optional) Real-time WebSocket delivery works

---

## ⚠️ CRITICAL QUESTIONS FOR MANAGER

### **Question 1: Scheduler Implementation**

**Context:** Need to run agents on schedule (hourly, daily, weekly).

**Options:**
- **A) Simple Polling Loop** (2 hours) - Background thread, 60-second check
- **B) APScheduler** (1 day) - Cron-like, dynamic scheduling
- **C) Celery + Redis** (3-5 days) - Full task queue, production-grade

**Question:** Which scheduler should we use for Phase 1?

**My Recommendation:** Option A (simple polling) for Phase 1, migrate to Option B later

---

### **Question 2: Result Storage Strategy**

**Context:** Agents will generate lots of results over time.

**Options:**
- **A) Store all results** (simple but grows fast)
- **B) Store for 30 days, then archive** (moderate complexity)
- **C) Store only unread + high-priority** (complex dedup logic)

**Question:** What's the retention policy for agent results?

**My Recommendation:** Option A initially (store all), add retention policy later

---

### **Question 3: Relevance Ranking Method**

**Context:** How should we rank PubMed papers by relevance to user's query?

**Options:**
- **A) Keyword matching only** (simple, fast)
- **B) Oracle-powered** (check if paper mentions variants with high delta scores)
- **C) LLM-powered** (use Gemini to assess relevance from abstract)
- **D) Hybrid** (keywords + Oracle for variant papers)

**Question:** Which relevance ranking method should we use?

**My Recommendation:** Option D (hybrid) - fast for most, Oracle for variant-specific

---

### **Question 4: Real-Time vs Polling**

**Context:** How should frontend get agent updates?

**Options:**
- **A) Polling** (GET /api/agents every 10-30 seconds)
- **B) WebSockets** (real-time push)
- **C) Server-Sent Events** (SSE, one-way push)

**Question:** Which real-time update mechanism?

**My Recommendation:** Option A (polling) for Phase 1, add WebSockets in Phase 5

---

### **Question 5: Agent Limits & Quotas**

**Context:** How many agents can a user create?

**Options:**
- **A) Unlimited** (simple but risky)
- **B) Tier-based** (Free: 2, Pro: 10, Enterprise: unlimited)
- **C) Per-type limits** (3 PubMed, 5 Trial, 10 VUS per user)

**Question:** What are the agent limits per user tier?

**My Recommendation:** Option B (tier-based) - Free: 2, Pro: 10, Enterprise: unlimited

---

### **Question 6: Email Notifications**

**Context:** Should agents send email alerts for critical findings?

**Options:**
- **A) No emails** (in-app only)
- **B) Optional emails** (user can enable per agent)
- **C) Digest emails** (daily/weekly summary)

**Question:** Should we build email notification system?

**My Recommendation:** Option A for Phase 1 (in-app only), add email in Phase 6

---

### **Question 7: Agent Templates**

**Context:** Should we provide pre-configured agent templates?

**Example Templates:**
- "Monitor BRAF research for melanoma"
- "Find frontline ovarian cancer trials in NYC"
- "Watch my 10 VUS for reclassifications"

**Question:** Should we create agent templates for quick setup?

**My Recommendation:** YES - Add 3-5 templates per agent type in Phase 1

---

## 🎯 REVISED TIMELINE (With Manager Input)

**Awaiting manager decisions on:**
1. Scheduler implementation (A/B/C)
2. Result retention policy (A/B/C)
3. Relevance ranking (A/B/C/D)
4. Real-time updates (A/B/C)
5. Agent limits (A/B/C)
6. Email notifications (A/B/C)
7. Agent templates (YES/NO)

**Estimated Timeline (Optimistic - Option A for all):**
- Phase 1: 3-5 days
- Phase 2: 3-4 days
- Phase 3: 3-4 days
- Phase 4: 2-3 days
- Phase 5: 2 days
- **Total: 13-18 days (~3-4 weeks)**

**Estimated Timeline (Conservative - Mix of A/B):**
- Phase 1: 5-7 days
- Phase 2: 4-5 days
- Phase 3: 4-5 days
- Phase 4: 3-4 days
- Phase 5: 3 days
- **Total: 19-24 days (~4-5 weeks)**

---

## 🏗️ ARCHITECTURE DIAGRAM

```
┌─────────────────────────────────────────────────────────────────┐
│                         USER (Alpha)                             │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                    FRONTEND (React + Vite)                       │
│                                                                  │
│  Components:                                                     │
│  ├─ AgentWizard (create agents)                                │
│  ├─ AgentDashboard (manage agents)                             │
│  ├─ AgentResultsViewer (view findings)                         │
│  └─ NotificationCenter (alerts)                                │
│                                                                  │
│  Context:                                                        │
│  ├─ AuthContext (user auth)                                    │
│  ├─ AgentContext (agent state) [NEW]                           │
│  └─ ActivityContext (logging)                                  │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                  BACKEND API (FastAPI)                           │
│                                                                  │
│  Routers:                                                        │
│  └─ /api/agents [NEW]                                           │
│     ├─ POST /agents (create)                                    │
│     ├─ GET /agents (list)                                       │
│     ├─ GET /agents/{id}/results                                │
│     └─ GET /agents/alerts                                       │
│                                                                  │
│  Services:                                                       │
│  ├─ AgentManager (CRUD) [NEW]                                  │
│  ├─ AgentExecutor (run agents) [NEW]                           │
│  ├─ AgentScheduler (schedule runs) [NEW]                       │
│  └─ Agent-specific services [NEW]:                             │
│     ├─ PubMedSentinel                                          │
│     ├─ TrialScout                                              │
│     └─ VUSVigil                                                │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                  BACKGROUND SCHEDULER                            │
│                                                                  │
│  Thread/APScheduler:                                            │
│  - Every 60 seconds: Check agents.next_run_at                  │
│  - If next_run_at <= NOW(): Execute agent                      │
│  - Update next_run_at based on frequency                       │
│  - Store results in database                                   │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                    DATABASE (Supabase PostgreSQL)                │
│                                                                  │
│  Tables:                                                         │
│  ├─ auth.users (Supabase)                                      │
│  ├─ public.user_profiles (tier, role)                          │
│  ├─ public.agents (config, status) [NEW]                       │
│  ├─ public.agent_runs (execution history) [NEW]                │
│  ├─ public.agent_results (findings) [NEW]                      │
│  └─ public.agent_alerts (notifications) [NEW]                  │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│              EXTERNAL DATA SOURCES                               │
│                                                                  │
│  ├─ PubMed API (literature search)                             │
│  ├─ ClinicalTrials.gov API (trial search)                      │
│  ├─ ClinVar API (variant classification)                       │
│  ├─ AstraDB (semantic trial search)                            │
│  └─ Neo4j (trial relationship intelligence)                    │
└─────────────────────────────────────────────────────────────────┘
```

---

## 🎯 IMMEDIATE NEXT STEPS

**Awaiting Manager Decisions (7 questions):**
1. Scheduler implementation choice
2. Result retention policy
3. Relevance ranking method
4. Real-time update mechanism
5. Agent limits per tier
6. Email notification support
7. Agent template creation

**Once Manager Approves:**
1. Create database schema (2 hours)
2. Build Agent Manager service (4 hours)
3. Build Agent Executor service (6 hours)
4. Build Simple Scheduler (4 hours)
5. Create Agent Router (2 hours)
6. Build AgentContext (3 hours)
7. Build Agent Wizard (6 hours)
8. Rebuild Agent Dashboard (4 hours)

**Total Phase 1:** 31 hours (~4 days with 8h/day)

---

**Nyx online. Zeta data streams are open.** ⚔️

**Architecture review complete. Build plan ready. Awaiting Alpha's command.**






