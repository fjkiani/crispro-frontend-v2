# ⚔️ CLINICAL TRIALS SYSTEM - MASTER DOCUMENT

**Last Updated:** January 13, 2025  
**Status:** ✅ **100% MIGRATION COMPLETE** | ✅ **GRAPH OPTIMIZATION IMPLEMENTED** | ⚠️ **ASTRADB SEEDING REQUIRED**  
**Single Source of Truth:** This document consolidates all clinical trials migration, graph optimization, and integration documentation.  
**Location:** `.cursor/MOAT/CLINICAL_TRIALS/01_SYSTEM_ARCHITECTURE.md`  
**See Also:** [00_MISSION.mdc](00_MISSION.mdc) for mission overview, [02_IMPLEMENTATION_STATUS.mdc](02_IMPLEMENTATION_STATUS.mdc) for implementation status

---

## 📚 TABLE OF CONTENTS

1. [Executive Summary](#executive-summary)
2. [Migration History](#migration-history)
3. [Current Architecture](#current-architecture)
4. [Graph Optimization Strategy](#graph-optimization-strategy)
5. [Implementation Status](#implementation-status)
6. [Integration Points](#integration-points)
7. [Testing & Verification](#testing--verification)
8. [Known Issues & Next Steps](#known-issues--next-steps)
9. [Reference Documents](#reference-documents)

---

## 🎯 EXECUTIVE SUMMARY

### **Mission Complete: Clinical Trials System Migration**

**What Was Achieved:**
- ✅ **100% Migration Complete**: All clinical trials functionality moved from main backend to minimal backend
- ✅ **ChromaDB → AstraDB**: Eliminated local file dependencies, enabled serverless deployment
- ✅ **Graph Optimization**: Neo4j + AstraDB hybrid architecture for intelligent trial discovery
- ✅ **Agent 1 Integration**: Seeding scripts updated to use minimal backend database paths
- ✅ **Frontend Integration**: ResearchPortal.jsx fully wired to new endpoints

**Current Capabilities:**
- **Hybrid Search**: Semantic search (AstraDB) + graph optimization (Neo4j)
- **Autonomous Agent**: AI-driven trial matching from patient data
- **Graph Algorithms**: PageRank, centrality, community detection for trial ranking
- **Relationship Modeling**: Trial-Sponsor, Trial-Site, Trial-PI, Trial-Drug relationships

**Critical Gaps:**
- ⚠️ **AstraDB Not Seeded**: Frontend will return zero results until seeding script is run
- ⚠️ **Neo4j Connection**: Graceful degradation implemented, but graph features limited if unavailable
- ⚠️ **PI Extraction**: Some Principal Investigator names may be missing from graph

---

## 📖 MIGRATION HISTORY

### **Phase 1: ChromaDB → AstraDB Migration**

**Problem:** ChromaDB used local files, incompatible with serverless deployment.

**Solution:**
- Migrated vector storage to AstraDB (cloud-based, serverless-compatible)
- Standardized on 768-dim Google Embeddings (`models/embedding-001`)
- Created `seed_astradb_from_sqlite.py` script for data migration
- Updated all search services to use AstraDB DataAPIClient

**Files Created:**
- `api/services/database_connections.py` - Unified database connection manager
- `scripts/seed_astradb_from_sqlite.py` - Migration script

**Status:** ✅ **COMPLETE**

---

### **Phase 2: Agent 1 Relocation**

**Problem:** Agent 1 seeding scripts were tightly coupled to main backend structure.

**Solution:**
- Updated `SQLITE_DB_PATH` in `oncology-backend/scripts/agent_1_seeding/config.py`
- Points to minimal backend database: `../oncology-backend-minimal/data/clinical_trials.db`
- Deprecated ChromaDB usage in favor of AstraDB
- Agent 1 remains in main backend as data pipeline (not migrated)

**Status:** ✅ **COMPLETE**

---

### **Phase 3: Backend Architecture Consolidation**

**Problem:** Clinical trials functionality split between main backend (proxy) and minimal backend (service).

**Solution:**
- **Eliminated Proxy Dependencies**: Removed all `MAIN_BACKEND_URL` proxy code
- **Direct Service Usage**: `api/routers/trials.py` now uses `ClinicalTrialSearchService` directly
- **Self-Contained**: Zero external dependencies on main backend
- **Unified Router**: All trial endpoints in minimal backend

**Files Modified:**
- `api/routers/trials.py` - Removed proxy code, added direct service calls
- `api/main.py` - Router registration verified

**Status:** ✅ **COMPLETE**

---

## 🏗️ CURRENT ARCHITECTURE

### **Backend Services**

#### **1. Database Connections** (`api/services/database_connections.py`)
- **SQLite**: Local database for trial metadata (1000+ trials)
- **AstraDB**: Vector storage for semantic search (768-dim embeddings)
- **Neo4j**: Graph database for relationship optimization
- **Singleton Pattern**: Reusable connection manager
- **Graceful Degradation**: Continues if Neo4j unavailable

#### **2. Clinical Trial Search Service** (`api/services/clinical_trial_search_service.py`)
- **Semantic Search**: Google embeddings (768-dim) via AstraDB
- **Disease Category Filtering**: Automatic category mapping
- **State Filtering**: Location-based trial filtering
- **Similarity Scoring**: `$similarity` extraction from AstraDB results
- **Error Handling**: Graceful degradation on service failures

#### **3. Hybrid Trial Search Service** (`api/services/hybrid_trial_search.py`)
- **Hybrid Architecture**: AstraDB semantic search + Neo4j graph optimization
- **Graph Algorithms**: PageRank, centrality, community detection
- **Relationship Leveraging**: Trial-Sponsor, Trial-Site, Trial-PI relationships
- **Ranking**: Combines semantic similarity + graph optimization scores

#### **4. Autonomous Trial Agent** (`api/services/autonomous_trial_agent.py`)
- **AI-Driven Matching**: Generates search queries from patient data
- **Context Extraction**: Mutations, disease, biomarkers → search queries
- **Multi-Query Strategy**: 1-3 queries per patient profile
- **Deduplication**: Merges results from multiple queries
- **Provenance Tracking**: Full audit trail of query generation

### **API Endpoints**

#### **`POST /api/trials/search`** (Basic Search)
- **Input**: `{ query, patient_context, top_k }`
- **Output**: `{ success, data: { found_trials[], count } }`
- **Method**: Semantic search via AstraDB

#### **`POST /api/trials/search-optimized`** (Graph-Optimized)
- **Input**: `{ query, patient_context, top_k }`
- **Output**: `{ success, data: { found_trials[], optimization_method, count } }`
- **Method**: Hybrid search (AstraDB + Neo4j)

#### **`POST /api/trials/agent/search`** (Autonomous Agent)
- **Input**: `{ mutations[], disease, state, biomarkers[] }`
- **Output**: `{ success, data: { matched_trials[], queries_used[], patient_context, total_found } }`
- **Method**: AI-driven query generation + hybrid search

### **Frontend Integration**

#### **ResearchPortal.jsx** (`oncology-frontend/src/pages/ResearchPortal.jsx`)
- **Status:** ✅ **FULLY IMPLEMENTED & WIRED**
- **Three Tabs:**
  1. **Manual Search**: Direct query input
  2. **Graph-Optimized**: Hybrid search with graph ranking
  3. **Autonomous Agent**: AI-driven matching from patient data
- **API Integration:** All tabs call `/api/trials/*` endpoints
- **Display:** Trial cards with NCT ID, title, phase, status, location

**Critical Note:** Frontend will return zero results until AstraDB is seeded.

---

## 🧠 GRAPH OPTIMIZATION STRATEGY

### **Graph Schema**

#### **Node Types:**
- **Trial**: Clinical trial with metadata (NCT ID, title, phase, status)
- **Sponsor**: Trial sponsor organization
- **Site**: Clinical trial site (location, state)
- **Principal Investigator (PI)**: Lead researcher
- **Drug**: Therapeutic agent
- **Condition**: Disease/condition

#### **Relationship Types:**
- **SPONSORED_BY**: Trial → Sponsor
- **CONDUCTED_AT**: Trial → Site
- **LED_BY**: Trial → PI
- **TESTS**: Trial → Drug
- **TARGETS**: Trial → Condition

### **Graph Algorithms**

#### **1. PageRank**
- **Purpose**: Identify high-authority trials based on network structure
- **Use Case**: Boost trials with strong sponsor/site connections

#### **2. Centrality Measures**
- **Degree Centrality**: Trials with many relationships
- **Betweenness Centrality**: Trials bridging different communities
- **Use Case**: Identify influential trials in network

#### **3. Community Detection**
- **Purpose**: Group related trials (same sponsor, similar conditions)
- **Use Case**: Recommend trials from same research network

#### **4. Shortest Path**
- **Purpose**: Find connections between patient profile and trials
- **Use Case**: "Similar patients enrolled in these trials"

### **Hybrid Search Flow**

1. **AstraDB Semantic Search**: Find ~50 candidate trials via vector similarity
2. **Neo4j Graph Optimization**: Rank candidates using graph algorithms
3. **Score Combination**: `final_score = semantic_score * 0.6 + graph_score * 0.4`
4. **Top-K Selection**: Return top 10 trials with optimization scores

---

## ✅ IMPLEMENTATION STATUS

### **Component 1: Relationship Extraction** ✅
- **Status:** ✅ Complete
- **File:** `api/services/trials_graph.py`
- **Capabilities:**
  - Extracts Trial-Sponsor relationships from ClinicalTrials.gov data
  - Extracts Trial-Site relationships (location, state)
  - Extracts Trial-PI relationships (Principal Investigator names)
  - Extracts Trial-Drug relationships from intervention data
  - Extracts Trial-Condition relationships from condition data

### **Component 2: Neo4j Setup** ✅
- **Status:** ✅ Complete
- **File:** `api/services/neo4j_connection.py`
- **Capabilities:**
  - Connection management with graceful degradation
  - Database verification and health checks
  - Error handling for connection failures

### **Component 3: Graph Loader** ✅
- **Status:** ✅ Complete
- **File:** `scripts/load_trials_to_neo4j.py`
- **Capabilities:**
  - Batch loading of trials, sponsors, sites, PIs, drugs, conditions
  - Relationship creation (SPONSORED_BY, CONDUCTED_AT, LED_BY, TESTS, TARGETS)
  - Idempotent (safe to re-run)

### **Component 4: Hybrid Search Service** ✅
- **Status:** ✅ Complete
- **File:** `api/services/hybrid_trial_search.py`
- **Capabilities:**
  - AstraDB semantic search (primary)
  - Neo4j graph optimization (secondary)
  - Score combination and ranking
  - Optimization method tracking

### **Component 5: Autonomous Trial Agent** ✅
- **Status:** ✅ Complete
- **File:** `api/services/autonomous_trial_agent.py`
- **Capabilities:**
  - Patient context extraction (mutations, disease, biomarkers)
  - AI-driven query generation (Gemini API)
  - Multi-query execution with deduplication
  - Provenance tracking

### **Frontend Integration** ✅
- **Status:** ✅ Complete
- **File:** `oncology-frontend/src/pages/ResearchPortal.jsx`
- **Capabilities:**
  - Three-tab interface (Manual, Graph-Optimized, Autonomous)
  - API integration with all endpoints
  - Trial card display with metadata
  - Error handling and loading states

---

## 🔗 INTEGRATION POINTS

### **Ayesha Complete Care Workflow**

**Integration Strategy:** Hybrid approach - trials integrated into Ayesha's care page while maintaining standalone Research Portal.

**Where Trials Fit:**
1. **Initial Diagnosis**: Show relevant trials during first-line treatment planning
2. **Treatment Response**: Update trial matches based on CA-125 trends
3. **Resistance Detection**: Suggest trials for resistant patients
4. **Biomarker Updates**: Re-rank trials when NGS results arrive

**Implementation:**
- Ayesha Complete Care page calls `/api/ayesha/trials/search` (Ayesha-specific filtering)
- Research Portal calls `/api/trials/*` (general search)
- Both use same backend services (hybrid search, autonomous agent)

### **Sporadic Cancer Integration**

**Biomarker Boosts:**
- HRD ≥42: PARP inhibitor trial boost (1.3x)
- TMB ≥20: Immunotherapy trial boost (1.3x)
- MSI-H: Immunotherapy trial boost (1.3x)

**Germline Filtering:**
- Germline-negative: Filter out germline-only trials
- Germline-positive: Include all relevant trials

---

## 🧪 TESTING & VERIFICATION

### **End-to-End Test Status**

**Test 1: Hybrid Graph-Optimized Search**
- **Endpoint:** `POST /api/trials/search-optimized`
- **Status:** ⚠️ **REQUIRES BACKEND RESTART** (API key fix applied)
- **Expected:** HTTP 200, `found_trials[]` with `optimization_score`
- **Performance:** < 2 seconds response time

**Test 2: Autonomous Trial Agent**
- **Endpoint:** `POST /api/trials/agent/search`
- **Status:** ⚠️ **REQUIRES BACKEND RESTART** (API key fix applied)
- **Expected:** HTTP 200, `matched_trials[]` with `queries_used[]`
- **Performance:** < 3 seconds response time

### **Known Issues**

1. **API Key Configuration:**
   - Service expects `GOOGLE_API_KEY` but `.env` has `GEMINI_API_KEY`
   - **Fix Applied:** Updated `clinical_trial_search_service.py` to support both
   - **Action Required:** Restart backend server

2. **AstraDB Seeding:**
   - **Status:** ⚠️ **NOT SEEDED**
   - **Impact:** Frontend will return zero results
   - **Fix:** Run `scripts/seed_astradb_from_sqlite.py`
   - **Time:** ~16 minutes for 1000 trials

3. **Neo4j Connection:**
   - **Status:** ✅ **GRACEFUL DEGRADATION IMPLEMENTED**
   - **Impact:** Graph features limited if unavailable
   - **Behavior:** System continues with semantic search only

4. **PI Extraction:**
   - **Status:** ⚠️ **SOME PIs MISSING**
   - **Impact:** Some trials may not have PI relationships in graph
   - **Fix:** Improve extraction logic in `trials_graph.py`

---

## ⚠️ KNOWN ISSUES & NEXT STEPS

### **Critical Actions Required**

1. **Seed AstraDB** (P0 - Blocks Frontend)
   ```bash
   cd oncology-coPilot/oncology-backend-minimal
   python scripts/seed_astradb_from_sqlite.py
   ```
   - **Time:** ~16 minutes
   - **Result:** Frontend will return trial results

2. **Restart Backend** (P0 - Enables Testing)
   ```bash
   cd oncology-coPilot/oncology-backend-minimal
   python -m uvicorn api.main:app --reload --host 0.0.0.0 --port 8000
   ```
   - **Reason:** API key fix requires restart
   - **Result:** Autonomous agent will work

3. **Run End-to-End Tests** (P1 - Verification)
   ```bash
   cd oncology-coPilot/oncology-backend-minimal
   ./test_end_to_end.sh
   ```
   - **Expected:** Both test endpoints return HTTP 200
   - **Result:** Verify hybrid search and autonomous agent work

### **Enhancement Opportunities**

1. **Improve PI Extraction** (P2)
   - Enhance `trials_graph.py` to extract more PI names
   - Add fallback parsing for different PI name formats

2. **Graph Algorithm Tuning** (P2)
   - Experiment with different score combination weights
   - A/B test semantic vs graph ranking performance

3. **Caching Strategy** (P2)
   - Cache graph optimization results for common queries
   - Reduce Neo4j query load

4. **Frontend Enhancements** (P2)
   - Add trial comparison view
   - Add "Save Trial" functionality
   - Add trial detail modal with full eligibility criteria

---

## 📖 REFERENCE DOCUMENTS

### **Essential References (Kept for Detailed Context)**

1. **`GRAPH_CONQUEST_MASTER_PLAN.md`** (48K)
   - Complete implementation details for all 5 components
   - File locations, code snippets, testing procedures
   - **Use When:** Need detailed implementation reference

2. **`GRAPH_OPTIMIZATION_COMPLETE_STRATEGY.md`** (25K)
   - Complete graph schema design
   - Algorithm explanations and use cases
   - **Use When:** Need graph architecture details

3. **`INTEGRATION_STRATEGY.md`** (17K)
   - Ayesha workflow integration points
   - Sporadic cancer biomarker boosts
   - **Use When:** Need integration planning

### **Archived Documents (Content Consolidated Above)**

All other documents have been archived to `.cursor/rules/clinical_trials_agents/archive/`:
- Migration plans and audits
- Status reports
- Test reports
- Frontend status documents

See `archive/ARCHIVE_INDEX.md` for complete list.

---

## ⚔️ DOCTRINE STATUS: ACTIVE

**LAST UPDATED:** January 13, 2025  
**APPLIES TO:** All clinical trials functionality, graph optimization, and integration  
**ENFORCEMENT:** Mandatory reference for all clinical trials development

**This master document represents the complete consolidation of all clinical trials migration, graph optimization, and integration knowledge. Every strategy, implementation detail, and status update is preserved and organized for maximum clarity and actionability.**


