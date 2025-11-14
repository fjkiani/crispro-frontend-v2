# 🏛️ CLINICAL TRIALS BACKEND ARCHITECTURE - STRATEGIC PLAN

## **⚔️ MISSION: Eliminate Monolith, Optimize for Production**

**Goal:** Determine optimal backend architecture for clinical trials as we upgrade from monolith to modular services.

---

## **📊 CURRENT STATE AUDIT**

### **Backend 1: `oncology-backend` (Main Backend - Port 8001)**
**Purpose:** Full-featured agent system with orchestration

**Clinical Trials Components:**
- ✅ `ClinicalTrialAgent` (567 lines) - Full vector search + LLM assessment
- ✅ `/api/search-trials` endpoint
- ✅ Uses AstraDB + SQLite (via `DatabaseConnections`)
- ✅ SentenceTransformer (384-dim) - needs migration
- ⚠️ Heavy dependencies: AgentInterface, orchestrator, patient management

**Other Features:**
- 15+ AI agents (genomic analyst, consultation synthesizer, etc.)
- Patient management system
- Agent orchestrator
- Population-level analytics

### **Backend 2: `oncology-backend-minimal` (Minimal Backend - Port 8000)**
**Purpose:** Production-focused, Vercel-deployable, modular API

**Clinical Trials Components:**
- ⚠️ `/api/search-trials` - **PROXIES to main backend** (tight coupling!)
- ✅ `/api/trials/refresh_status` - Self-contained service
- ✅ `/api/clinical_trials/match` - Direct CT.gov API calls
- ✅ `trial_refresh` service (modular, reusable)
- ✅ No agent dependencies

**Other Features:**
- S/P/E framework (efficacy orchestrator)
- Evidence services (ClinVar, literature)
- Insights bundle (4-chip system)
- Knowledge Base (KB system)
- Modular, production-ready architecture

---

## **🚨 CURRENT PROBLEMS**

### **Problem 1: Tight Coupling via Proxy**
```python
# Minimal backend proxies to main backend:
MAIN_BACKEND_URL = os.getenv("MAIN_BACKEND_URL", "http://localhost:8001")
# This creates deployment dependency!
```

**Issues:**
- Minimal backend can't run standalone
- Deployment requires both backends running
- Single point of failure
- Network latency for every search

### **Problem 2: Duplicate Implementations**
- Two different CT.gov API implementations:
  - `clinical_trials.py` (minimal backend) - direct API
  - `ClinicalTrialAgent` (main backend) - vector search + LLM
- Confusing which endpoint to use
- Inconsistent response formats

### **Problem 3: Agent Dependency**
- `ClinicalTrialAgent` depends on `AgentInterface` base class
- Tightly coupled to main backend's agent system
- Can't run in minimal backend without bringing entire agent framework

---

## **🎯 RECOMMENDED ARCHITECTURE: CONSOLIDATE INTO MINIMAL BACKEND**

### **✅ OPTION A: MIGRATE CLINICAL TRIALS TO MINIMAL BACKEND** (RECOMMENDED)

**Strategy:** Move ClinicalTrialAgent functionality to minimal backend as a **service** (not agent), eliminating proxy dependency.

**Benefits:**
- ✅ **Self-contained**: Minimal backend runs standalone
- ✅ **Production-ready**: Vercel-deployable, no external dependencies
- ✅ **Consistent**: All clinical trials code in one place
- ✅ **Modular**: Fits existing services architecture (trial_refresh, evidence, etc.)
- ✅ **Scalable**: Cloud-native (AstraDB, no local files)

**Implementation:**
1. Extract ClinicalTrialAgent logic → `ClinicalTrialService` (no AgentInterface dependency)
2. Move to `oncology-backend-minimal/api/services/clinical_trials/`
3. Update `/api/search-trials` to use service directly (no proxy)
4. Migrate to Google Embeddings (768-dim) for consistency
5. Remove proxy code

**New Structure:**
```
oncology-backend-minimal/
├── api/
│   ├── routers/
│   │   ├── trials.py              # Search endpoint (no proxy)
│   │   └── clinical_trials.py     # Matching endpoint
│   └── services/
│       ├── clinical_trials/
│       │   ├── search_service.py  # Vector search + LLM (from agent)
│       │   ├── database.py        # DatabaseConnections adapter
│       │   └── models.py           # Request/response schemas
│       └── trial_refresh/          # Already exists ✅
```

---

### **❌ OPTION B: SEPARATE CLINICAL TRIALS BACKEND** (NOT RECOMMENDED)

**Why Not:**
- ❌ Over-engineering: Clinical trials is just one feature
- ❌ Deployment complexity: Need 3 backends (main, minimal, trials)
- ❌ Duplication: Would duplicate DatabaseConnections, configs
- ❌ Network overhead: Extra HTTP hop for every request
- ❌ Maintenance burden: 3 codebases to update

**When It Makes Sense:**
- If clinical trials becomes 50%+ of backend code
- If it needs separate scaling/performance requirements
- If you're building microservices architecture (probably premature)

---

### **❌ OPTION C: KEEP CURRENT PROXY ARCHITECTURE** (NOT RECOMMENDED)

**Why Not:**
- ❌ Tight coupling: Minimal backend depends on main backend
- ❌ Deployment complexity: Both must run
- ❌ Single point of failure
- ❌ Network latency: Extra HTTP hop
- ❌ Confusing: Two backends, unclear ownership

---

## **📋 MIGRATION PLAN: CLINICAL TRIALS → MINIMAL BACKEND**

### **Phase 1: Extract Service Logic** (3-4 hours)

**File:** `oncology-backend-minimal/api/services/clinical_trials/search_service.py`

**Extract from:** `oncology-backend/backend/agents/clinical_trial_agent.py`

**Key Changes:**
1. Remove `AgentInterface` dependency (just a regular class)
2. Replace `SentenceTransformer` → `Google Embeddings` (768-dim)
3. Adapt `DatabaseConnections` pattern (or create adapter)
4. Keep all business logic (vector search, LLM assessment, SQLite fetch)

**Code Structure:**
```python
class ClinicalTrialSearchService:
    """Clinical trial search service - no agent dependency."""
    
    def __init__(self):
        self.db_connections = self._init_database_connections()
        self.google_ef = self._init_google_embeddings()
        self.llm_client = self._init_gemini_client()
    
    async def search_trials(
        self,
        query: str,
        patient_context: Optional[Dict[str, Any]] = None,
        page_state: Optional[str] = None
    ) -> Dict[str, Any]:
        """Search trials with vector search + LLM assessment."""
        # 1. Build rich query text
        rich_query = self._build_query_text(query, patient_context)
        
        # 2. Generate embedding (Google API - 768-dim)
        query_embedding = self.google_ef([rich_query])[0]
        
        # 3. Vector search in AstraDB
        similar_trials = await self._vector_search(query_embedding, page_state)
        
        # 4. Fetch details from SQLite
        trial_details = await self._fetch_trial_details(similar_trials)
        
        # 5. LLM eligibility assessment
        assessed_trials = await self._assess_eligibility(trial_details, patient_context)
        
        return {
            "results": assessed_trials,
            "next_page_state": next_page_state
        }
```

---

### **Phase 2: Create Database Connections Adapter** (1 hour)

**File:** `oncology-backend-minimal/api/services/clinical_trials/database.py`

**Purpose:** Port `DatabaseConnections` helper to minimal backend (or create adapter)

**Options:**
- **Option A**: Copy `DatabaseConnections` to minimal backend
- **Option B**: Create adapter that uses same pattern
- **Option C**: Use existing AstraDB client directly (if already available)

**Recommended:** Option A (copy the helper - it's clean and reusable)

---

### **Phase 3: Update Router** (1 hour)

**File:** `oncology-backend-minimal/api/routers/trials.py`

**Change:** Remove proxy, use service directly

**Before (Proxy):**
```python
@router.post("/api/search-trials")
async def search_trials(request: TrialSearchRequest):
    # Proxy to main backend
    response = await client.post(f"{MAIN_BACKEND_URL}/api/search-trials", ...)
```

**After (Direct Service):**
```python
from ..services.clinical_trials.search_service import ClinicalTrialSearchService

search_service = ClinicalTrialSearchService()

@router.post("/api/search-trials")
async def search_trials(request: TrialSearchRequest):
    """Search trials using AstraDB vector search (self-contained)."""
    results = await search_service.search_trials(
        query=request.query,
        patient_context=request.patient_context,
        page_state=request.page_state
    )
    return {
        "success": True,
        "data": {
            "found_trials": results["results"],
            "next_page_state": results.get("next_page_state")
        }
    }
```

---

### **Phase 4: Migrate Embeddings** (1 hour)

**Update:** SentenceTransformer → Google Embeddings (768-dim)

**File:** `oncology-backend-minimal/api/services/clinical_trials/search_service.py`

**Before:**
```python
from sentence_transformers import SentenceTransformer
self.embedding_model = SentenceTransformer('all-MiniLM-L6-v2')  # 384-dim
query_embedding = self.embedding_model.encode(query).tolist()
```

**After:**
```python
from chromadb.utils import embedding_functions
self.google_ef = embedding_functions.GoogleGenerativeAiEmbeddingFunction(
    api_key=os.getenv("GOOGLE_API_KEY"),
    model_name="models/embedding-001"  # 768-dim
)
query_embedding = self.google_ef([query])[0]  # Returns 768-dim vector
```

**Also Update AstraDB Collection:**
- Re-create collection with 768 dimensions (if currently 384-dim)
- Re-seed trials with 768-dim embeddings (via Agent 1 migration)

---

### **Phase 5: Remove Proxy & Update Main Backend** (30 min)

**Changes:**
1. Remove `/api/search-trials` from `oncology-backend/main.py` (or mark deprecated)
2. Update `AgentOrchestrator` to note ClinicalTrialAgent moved
3. Update documentation

**Keep for Compatibility:**
- Can keep ClinicalTrialAgent in main backend for now (if other agents use it)
- Or remove entirely if no dependencies

---

### **Phase 6: Testing & Validation** (1 hour)

**Tests:**
1. Unit tests: Service methods with mocked databases
2. Integration tests: End-to-end search → results
3. Performance tests: Latency <500ms
4. Compatibility tests: Frontend still works

---

## **📊 COMPARISON: OPTIONS**

| Aspect | Option A: Migrate to Minimal | Option B: Separate Backend | Option C: Keep Proxy |
|--------|------------------------------|--------------------------|-------------------|
| **Deployment** | ✅ Single backend (self-contained) | ❌ 3 backends to deploy | ❌ 2 backends required |
| **Dependencies** | ✅ Clean (no agent framework) | ❌ Shared dependencies | ❌ Cross-backend dependency |
| **Latency** | ✅ Direct (no proxy hop) | ❌ Extra HTTP hop | ❌ Proxy latency |
| **Maintainability** | ✅ One codebase | ❌ 3 codebases | ⚠️ 2 codebases |
| **Scalability** | ✅ Auto-scales with minimal backend | ⚠️ Separate scaling | ⚠️ Both must scale |
| **Complexity** | ✅ Low | ❌ High | ⚠️ Medium |

---

## **🎯 FINAL RECOMMENDATION**

### **✅ MIGRATE CLINICAL TRIALS TO MINIMAL BACKEND**

**Rationale:**
1. **Production Focus**: Minimal backend is your Vercel-deployable, production-ready backend
2. **Modular Architecture**: Already has services pattern (trial_refresh, evidence, etc.)
3. **Self-Contained**: Can run standalone (no proxy dependencies)
4. **Consistency**: All clinical trials code in one place
5. **Future-Proof**: Easier to scale and maintain

**What Stays in Main Backend:**
- Agent orchestrator (if ClinicalTrialAgent used by other agents)
- Patient management (if trials tied to patient workflows)
- Population analytics (different use case)

**What Moves to Minimal Backend:**
- Clinical trial search (vector search + LLM)
- Trial refresh (already there ✅)
- Trial matching (already there ✅)
- Database connections (adapted)

---

## **📋 MIGRATION CHECKLIST**

### **Pre-Migration:**
- [ ] Audit ClinicalTrialAgent dependencies
- [ ] Identify all usages in main backend
- [ ] Verify DatabaseConnections can be ported
- [ ] Check if other agents depend on ClinicalTrialAgent

### **Migration Steps:**
- [ ] Create `api/services/clinical_trials/search_service.py`
- [ ] Port DatabaseConnections helper to minimal backend
- [ ] Migrate embeddings: SentenceTransformer → Google (768-dim)
- [ ] Update `/api/search-trials` router (remove proxy)
- [ ] Update AstraDB collection to 768-dim (if needed)
- [ ] Update tests
- [ ] Remove proxy code from trials.py

### **Post-Migration:**
- [ ] Mark ClinicalTrialAgent in main backend as deprecated (if kept)
- [ ] Update documentation
- [ ] Update frontend (if any changes needed)
- [ ] Performance validation
- [ ] Deployment test (Vercel/local)

---

## **🔧 IMPLEMENTATION DETAILS**

### **New Service Structure:**

```
api/services/clinical_trials/
├── __init__.py
├── search_service.py       # Main search logic (from agent)
├── database.py             # DatabaseConnections adapter
├── models.py               # Request/response schemas
└── README.md              # Service documentation
```

### **Dependencies to Add:**

```python
# requirements.txt additions:
chromadb>=0.4.0  # For Google Embeddings
google-generativeai>=0.3.0  # Already have
astrapy>=2.0.1  # Already have
```

### **Remove from Requirements:**
```python
# Remove if not used elsewhere:
sentence-transformers  # Replace with Google Embeddings
```

---

## **💰 COST/BENEFIT ANALYSIS**

### **Benefits:**
- ✅ **Zero deployment dependencies** (minimal backend standalone)
- ✅ **Faster queries** (no proxy latency)
- ✅ **Easier scaling** (one backend to scale)
- ✅ **Cleaner architecture** (modular services)
- ✅ **Better maintainability** (one codebase)

### **Costs:**
- ⚠️ **Migration time**: ~6-8 hours
- ⚠️ **Testing effort**: ~2 hours
- ⚠️ **Risk**: Need to validate search quality matches current

### **ROI:**
- **High** - Production deployment becomes much simpler
- **Medium-term**: Easier maintenance, faster development
- **Long-term**: Scalable architecture foundation

---

## **📝 NOTES**

### **Why Not Separate Backend:**
- Clinical trials is ~5-10% of total backend code (not large enough)
- Minimal backend already has modular services architecture
- Would create more complexity than value

### **Why Minimal Backend:**
- Already production-focused (Vercel deployment)
- Already has modular services pattern
- Already has trial_refresh service (natural fit)
- Clean separation from agent framework

### **Migration Complexity:**
- **Low-Medium**: Most logic is straightforward port
- Main challenge: DatabaseConnections adapter
- Embedding migration already planned (separate task)

---

**STATUS**: 📋 **PLAN COMPLETE** - Ready for execution  
**ESTIMATED TIME**: 6-8 hours total  
**PRIORITY**: P1 (Improves deployment, removes coupling)

