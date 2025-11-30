# 🔍 UNDOCUMENTED ROUTERS

**Created**: January 14, 2025  
**Purpose**: List of routers that need further documentation

---

## ADDITIONAL ROUTERS NOT FULLY DOCUMENTED

#### **7.10.1 Guidance Router** (`routers/guidance.py`):
- **Purpose**: Clinical gating facade over efficacy orchestrator
- **Features**:
  - On-label detection (stub ruleset)
  - Tier classification (I/II/III/research)
  - Resistance/sensitivity marker detection (PSMB5, CRBN, TP53)
  - Fused S score integration (AlphaMissense)
- **Status**: ⚠️ **NOT DOCUMENTED** - Clinical gating logic exists

#### **7.10.2 Sessions Router** (`routers/sessions.py`):
- **Purpose**: Session persistence API - save/resume analyses across pages
- **Features**:
  - Create/update sessions (Supabase-backed)
  - Add session items (insight/efficacy/dataset/note)
  - Retrieve session history
  - Idempotency support
- **Status**: ⚠️ **NOT DOCUMENTED** - Session management system

#### **7.10.3 Auth Router** (`routers/auth.py`):
- **Purpose**: User authentication and profile management
- **Endpoints**:
  - `/api/auth/signup` - Create account
  - `/api/auth/login` - Login
  - `/api/auth/logout` - Logout
  - `/api/auth/profile` - Get/update profile
- **Status**: ⚠️ **NOT DOCUMENTED** - Authentication system

#### **7.10.4 Admin Router** (`routers/admin.py`):
- **Purpose**: Admin dashboard endpoints
- **Features**:
  - User management (list, get, update, delete)
  - Analytics aggregation
  - Activity tracking
  - Usage limits management
- **Status**: ⚠️ **NOT DOCUMENTED** - Admin functionality

#### **7.10.5 Myeloma Router** (`routers/myeloma.py`):
- **Purpose**: Myeloma Digital Twin prediction
- **Endpoint**: `/api/predict/myeloma_drug_response`
- **Features**:
  - Evo2 live scoring only (no mocks)
  - Preflight validation (REF-check, duplicate collapse)
  - Variant call extraction
  - Supabase persistence
- **Status**: ⚠️ **NOT DOCUMENTED** - Myeloma-specific workflow

#### **7.10.6 Evidence RAG** (`routers/evidence/rag.py`):
- **Purpose**: RAG-based conversational query for clinical literature
- **Endpoints**:
  - `/api/evidence/rag-query` - Conversational query
  - `/api/evidence/rag-add-variant` - Add variant to KB
  - `/api/evidence/rag-stats` - KB statistics
- **Integration**: Uses external `Pubmed-LLM-Agent-main` RAG agent
- **Status**: ⚠️ **NOT DOCUMENTED** - RAG system

#### **7.10.7 Evidence Extraction** (`routers/evidence/extraction.py`):
- **Purpose**: Diffbot article extraction
- **Endpoint**: `/api/evidence/extract`
- **Features**: Extract full text from URLs via Diffbot API
- **Status**: ⚠️ **NOT DOCUMENTED** - Article extraction

#### **7.10.8 Evidence Jobs** (`routers/evidence/jobs.py`):
- **Purpose**: Background job orchestration for evidence processing
- **Endpoints**:
  - `/api/evidence/crawl` - Multi-URL crawling job
  - `/api/evidence/summarize` - Summarization job
  - `/api/evidence/align` - Evidence alignment job
- **Status**: ⚠️ **NOT DOCUMENTED** - Background job system

#### **7.10.9 Resistance Router** (`routers/resistance.py`):
- **Purpose**: Resistance analysis endpoints
- **Status**: ✅ **CLOSED** - See `GAP_CLOSURE_FROM_ARCHIVE.md` and `archive/RESISTANCE_PLAYBOOK_V1_COMPLETE.md`
- **Key Findings**: 
  - Router: `api/routers/care.py` (186 lines)
  - Endpoint: `POST /api/care/resistance_playbook`
  - Service: `api/services/resistance_playbook_service.py` (702 lines)
  - 5 detection rules, 7 combo strategies, 6 next-line switches

#### **7.10.10 Care Router** (`routers/care.py`):
- **Purpose**: Resistance Playbook endpoints (Section 17)
- **Status**: ✅ **CLOSED** - See `GAP_CLOSURE_FROM_ARCHIVE.md` and `archive/RESISTANCE_PLAYBOOK_V1_COMPLETE.md`
- **Key Findings**: Complete resistance playbook router with health check and main playbook generator endpoint

#### **7.10.11 Tumor Router** (`routers/tumor.py`):
- **Purpose**: Sporadic Cancer Strategy endpoints (Day 1-7)
- **Status**: ⚠️ **MENTIONED BUT NOT DEEP DIVE**

#### **7.10.12 Ayesha Router** (`routers/ayesha.py`):
- **Purpose**: Original Ayesha endpoints (not v2 orchestrator)
- **Status**: ⚠️ **NOT DOCUMENTED** - Legacy Ayesha endpoints

#### **7.10.13 Ayesha Twin Demo Router** (`routers/ayesha_twin_demo.py`):
- **Purpose**: Ayesha twin demo endpoints
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.10.14 Dossiers Router** (`routers/dossiers.py`):
- **Purpose**: JR2 dossier generation pipeline
- **Status**: ⚠️ **MENTIONED BUT NOT DEEP DIVE**

#### **7.10.15 Ayesha Dossiers Router** (`routers/ayesha_dossiers.py`):
- **Purpose**: Ayesha dossier browser API (display all 60 trials)
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.10.16 Trials Router** (`routers/trials.py`):
- **Purpose**: Search and refresh endpoints
- **Endpoints**: `/api/search-trials`, `/api/trials/refresh_status`
- **Status**: ✅ **PARTIALLY CLOSED** - See `archive/AYESHA_TRIAL_FILTERING_COMPLETE.md` for Ayesha trials router
- **Note**: Other trial routers (trials_graph, trials_agent) still need review

#### **7.10.17 Trials Graph Router** (`routers/trials_graph.py`):
- **Purpose**: Graph-optimized search (hybrid graph search)
- **Endpoint**: `/api/trials/search-optimized`
- **Status**: ✅ **PARTIALLY CLOSED** - See `archive/ZO_CLINICAL_TRIALS_FRONTEND_COMPLETE.md`
- **Key Findings**: Hybrid Neo4j + AstraDB search, frontend integration in GraphOptimizedSearch component

#### **7.10.18 Trials Agent Router** (`routers/trials_agent.py`):
- **Purpose**: Autonomous trial agent
- **Endpoint**: `/api/trials/agent/search`
- **Status**: ✅ **PARTIALLY CLOSED** - See `archive/ZO_CLINICAL_TRIALS_FRONTEND_COMPLETE.md`
- **Key Findings**: Frontend integration in AutonomousTrialAgent component, passes sporadic context

#### **7.10.19 Clinical Trials Router** (`routers/clinical_trials.py`):
- **Purpose**: Clinical trials endpoints
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.10.20 Command Center Router** (`routers/command_center.py`):
- **Purpose**: Command center endpoints (conditional on feature flag)
- **Status**: ⚠️ **NOT DOCUMENTED**

---
