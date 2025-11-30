# 🔍 UNDOCUMENTED SERVICES

**Created**: January 14, 2025  
**Purpose**: List of services that need further documentation

---

## ADDITIONAL SERVICES NOT FULLY DOCUMENTED

#### **7.11.1 Admin Service** (`services/admin_service.py`):
- **Purpose**: Admin operations (user management, analytics)
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.2 Auth Service** (`services/auth_service.py`):
- **Purpose**: Authentication operations (Supabase Auth integration)
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.3 Agent Services**:
- **Agent Executor** (`services/agent_executor.py`): Execute agent tasks
- **Agent Manager** (`services/agent_manager.py`): CRUD for agents
- **Agent Scheduler** (`services/agent_scheduler.py`): Scheduled execution
- **Status**: ✅ **CLOSED** - See `GAP_CLOSURE_FROM_ARCHIVE.md` and `archive/AYESHA_AGENT_MISSIONS_MASTER.md`
- **Key Findings**: 
  - Agent types: `pubmed_sentinel`, `trial_scout`, `genomic_forager`
  - Scheduler: 60s polling loop, executes up to 5 agents concurrently
  - Storage: Supabase-backed (agents, runs, results)
  - Tier limits: Free (3), Pro (10), Enterprise (unlimited)

#### **7.11.4 Autonomous Trial Agent** (`services/autonomous_trial_agent.py`):
- **Purpose**: Autonomous trial search agent
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.5 Cache Service** (`services/cache_service.py`):
- **Purpose**: Caching operations
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.6 Compound Services**:
- **Compound Alias Resolver** (`services/compound_alias_resolver.py`): Resolve drug aliases
- **Compound Calibration** (`services/compound_calibration.py`): Drug-specific calibration
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.7 Dietician Recommendations** (`services/dietician_recommendations.py`):
- **Purpose**: Dietary recommendations
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.8 Dynamic Food Extraction** (`services/dynamic_food_extraction.py`):
- **Purpose**: Extract food/supplement information
- **Status**: ✅ **CLOSED** - See `GAP_CLOSURE_FROM_ARCHIVE.md` and `archive/FOOD_VALIDATOR_E2E_INTEGRATION_COMPLETE.md`
- **Key Findings**: Full E2E integration complete, endpoint `/api/hypothesis/validate_food_dynamic`, LLM literature mining, SAE features

#### **7.11.9 Food S/P/E Integration** (`services/food_spe_integration.py`):
- **Purpose**: S/P/E framework for foods
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.10 Hotspot Detector** (`services/hotspot_detector.py`):
- **Purpose**: Detect mutation hotspots
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.11 LLM Literature Service** (`services/llm_literature_service.py`):
- **Purpose**: LLM-based literature synthesis
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.12 Mechanism Fit Ranker** (`services/mechanism_fit_ranker.py`):
- **Purpose**: Rank drugs by mechanism fit
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.13 Metastasis Service** (`services/metastasis_service.py`):
- **Purpose**: Metastasis cascade risk assessment (NOT interception)
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.14 NGS Fast Track** (`services/ngs_fast_track.py`):
- **Purpose**: Fast-track NGS recommendations
- **Status**: ✅ **CLOSED** - See `archive/ZO_BACKEND_COMPLETE_STATUS.md`
- **Key Findings**: 
  - ctDNA (Guardant360), HRD (MyChoice CDx), IHC panel recommendations
  - Parallel execution (~10 days total), cost estimates, ordering info
  - 300+ lines, fully operational

#### **7.11.15 Safety Validator** (`services/safety_validator.py`):
- **Purpose**: Safety validation logic
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.16 Therapeutic Prompt Builder** (`services/therapeutic_prompt_builder.py`):
- **Purpose**: Build prompts for therapeutic generation
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.17 Trial Intelligence Pipeline** (`services/trial_intelligence/`):
- **Purpose**: Multi-stage trial processing pipeline
- **Stages**: 6 stages (hard filters → trial type → location → eligibility → LLM analysis → dossier)
- **Status**: ✅ **CLOSED** - See `GAP_CLOSURE_FROM_ARCHIVE.md` and `archive/AYESHA_TRIAL_FILTERING_COMPLETE.md`
- **Key Findings**: 
  - Complete 6-stage pipeline documented (MatchOrchestrator, EligibilityFilters, ScoringEngine, ReasoningGenerator, CA125Intelligence, SOC Recommendation)
  - Backend: 7 modules, Frontend: 4 components, Full E2E testing

#### **7.11.18 Trial Refresh Service** (`services/trial_refresh/`):
- **Purpose**: Refresh trial data from ClinicalTrials.gov
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.19 Tumor Quick Intake** (`services/tumor_quick_intake.py`):
- **Purpose**: Generate initial tumor context using disease priors
- **Status**: ✅ **CLOSED** - See `archive/ZO_BACKEND_COMPLETE_STATUS.md`
- **Key Findings**: Level 0/1 intake endpoint, 15 cancers with TCGA priors, 25 test scenarios, 100% pass rate

#### **7.11.20 Logging Services** (`services/logging/`):
- **Efficacy Logger**: Log efficacy predictions
- **Evidence Logger**: Log evidence queries
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.11.21 Infrastructure Services**:
- **Supabase Service** (`services/supabase_service.py`): Supabase operations
- **Neo4j Connection** (`services/neo4j_connection.py`): Neo4j graph database
- **Neo4j Graph Loader** (`services/neo4j_graph_loader.py`): Load graph data
- **Database Connections** (`services/database_connections.py`): DB connection management
- **Enformer Client** (`services/enformer_client.py`): Chromatin accessibility client
- **Status**: ⚠️ **NOT DOCUMENTED**

---
