# 🔍 INTEGRATION PATTERNS NOT DOCUMENTED

**Created**: January 14, 2025  
**Purpose**: Integration patterns that need documentation

---

## INTEGRATION PATTERNS NOT DOCUMENTED

#### **7.13.1 Supabase Integration**:
- **Tables**: `user_sessions`, `mdt_runs`, `mdt_run_variants`, `job_results`, `agents`, `agent_runs`
- **Operations**: Select, insert, update, events
- **Status**: ✅ **CLOSED** - See `GAP_CLOSURE_FROM_ARCHIVE.md`
- **Key Findings**: 
  - Singleton pattern via `get_supabase_client()`
  - Used for agent storage, session management, job results
  - CRUD operations via Supabase Python client

#### **7.13.2 Neo4j Graph Integration**:
- **Graph Structure**: Trial relationships, eligibility graph
- **Queries**: Graph-optimized trial search
- **Status**: ✅ **PARTIALLY CLOSED** - See `GAP_CLOSURE_FROM_ARCHIVE.md`
- **Key Findings**: 
  - Connection: `9669e5f3.databases.neo4j.io` (Neo4j Cloud)
  - 910 relationships, 30 trials
  - Endpoint: `/api/trials/search-optimized` (hybrid graph search)
  - Graceful degradation: Connection failures don't crash app
  - **Note**: Graph query patterns still need documentation

#### **7.13.3 AstraDB Integration**:
- **Collections**: `clinical_trials_eligibility` (with vector search)
- **Operations**: Vector search, document storage
- **Status**: ✅ **CLOSED** - See `GAP_CLOSURE_FROM_ARCHIVE.md` and `archive/ZO_ASTRADB_SEEDING_STATUS.md`
- **Key Findings**: 
  - 768-dim vectors using `text-embedding-004` (no quota limits)
  - Vector field: `$vector` (root level, not nested)
  - Upsert: `find_one_and_update` with `upsert=True`
  - Seeding script: `scripts/seed_astradb_from_sqlite.py`
  - 30 trials successfully seeded

#### **7.13.4 SQLite Integration**:
- **Tables**: `clinical_trials` (28 columns)
- **Operations**: Trial storage, querying
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.13.5 Background Job System**:
- **Job Types**: Crawl, summarize, align
- **Job Store**: In-memory (JOBS dict) or Supabase
- **Status**: ⚠️ **NOT DOCUMENTED**

#### **7.13.6 Agent System Architecture**:
- **Agent Types**: `pubmed_sentinel`, `trial_scout`, `genomic_forager`
- **Execution**: Scheduled vs on-demand
- **Status**: ✅ **CLOSED** - See `GAP_CLOSURE_FROM_ARCHIVE.md` and `archive/AYESHA_AGENT_MISSIONS_MASTER.md`
- **Key Findings**: 
  - Scheduler: 60s polling loop, executes up to 5 agents concurrently
  - Storage: Supabase (agents, runs, results)
  - Manager: CRUD operations, tier-based limits
  - Executor: Agent-specific execution logic

---