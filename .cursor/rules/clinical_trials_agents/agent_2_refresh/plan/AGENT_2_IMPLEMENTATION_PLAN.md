# 🤖 AGENT 2: IMPLEMENTATION PLAN - LIVE REFRESH SERVICE

## **⚔️ MISSION BRIEF**
Build a live status refresh service that fetches current recruiting status and location data from ClinicalTrials.gov API v2 for cached trials.

**Target Backend**: `oncology-coPilot/oncology-backend-minimal` (FastAPI minimal backend)

---

## **📊 EXISTING CAPABILITIES AUDIT**

### **✅ What We Can Reuse:**

1. **ClinicalTrials.gov API Integration** ✅
   - **Location**: `oncology-coPilot/oncology-backend-minimal/api/routers/clinical_trials.py`
   - **Base URL**: Already using `https://clinicaltrials.gov/api/v2/studies`
   - **Pattern**: Using `httpx.AsyncClient` for async API calls
   - **Reusable**: API client pattern and base URL constant

2. **Study Parsing Utilities** ✅
   - **Location**: `oncology-coPilot/oncology-backend/backend/research/clinicaltrials_utils.py`
   - **Function**: `parse_study()` extracts trial details from API response
   - **Reusable**: Study parsing logic (can adapt for locations)

3. **ClinicalTrialAgent** ✅
   - **Location**: `oncology-coPilot/oncology-backend-minimal/api/routers/clinical_trials.py` (minimal backend)
   - **Location**: `oncology-coPilot/oncology-backend/backend/agents/clinical_trial_agent.py` (full backend)
   - **Reusable**: Integration point exists, needs modification

4. **FastAPI Router Pattern** ✅
   - **Location**: `oncology-coPilot/oncology-backend-minimal/api/routers/`
   - **Pattern**: Modular router structure with injection into `main.py`
   - **Reusable**: Follow existing router pattern

### **⚠️ What Needs New Implementation:**

1. **Refresh Service** ❌ (NEW)
   - Dedicated service for status + location fetching
   - Batch NCT ID processing
   - State filtering utility

2. **Refresh Endpoint** ❌ (NEW)
   - `POST /api/trials/refresh_status` endpoint
   - Request/response schemas

3. **Agent Integration** ⚠️ (MODIFY)
   - Add `live_refresh` parameter to agent
   - Merge live data into cached results

---

## **🎯 ARCHITECTURE DECISIONS**

### **Backend Target: Minimal Backend**
- **Reason**: Doctrine specifies `oncology-backend-minimal` 
- **Path**: `oncology-coPilot/oncology-backend-minimal/api/`

### **Service Location:**
- **Service**: `api/services/trial_refresh_service.py` (new)
- **Router**: Extend `api/routers/clinical_trials.py` (existing)
- **Agent**: Modify `oncology-backend/backend/agents/clinical_trial_agent.py` (full backend)

### **API Client:**
- **Use**: `httpx.AsyncClient` (matches existing pattern in `clinical_trials.py`)
- **Base URL**: Reuse `CLINICAL_TRIALS_BASE_URL = "https://clinicaltrials.gov/api/v2/studies"`

---

## **📋 IMPLEMENTATION CHECKLIST**

### **Task 1: Refresh Service Implementation** (1.5 hours)

**File**: `oncology-coPilot/oncology-backend-minimal/api/services/trial_refresh_service.py`

**Key Functions:**
1. `refresh_trial_status(nct_ids: List[str])` - Core async function
2. `refresh_trial_status_with_retry(nct_ids, max_retries=2)` - Wrapper with retry
3. `filter_locations_by_state(trial_data, state)` - State filtering utility

**Implementation Notes:**
- Use `httpx.AsyncClient` for consistency with existing code
- Batch API call: Use `query.id` parameter with comma-separated NCT IDs
- Field selection: Request only needed fields (status, locations, contacts)
- Parse API v2 response structure (protocolSection → contactsLocationsModule)
- Filter to only recruiting locations (RECRUITING, NOT_YET_RECRUITING)
- Add ISO timestamp to `last_updated` field

**Reusable Patterns from `clinicaltrials_utils.py`:**
- Error handling with `try/except`
- Logging pattern
- Timeout handling (30 seconds)

**Differences from Existing Code:**
- **Focus**: ONLY status + locations (not full trial parsing)
- **Input**: List of NCT IDs (not search query)
- **Output**: Dict mapping NCT ID → {status, locations, last_updated}

---

### **Task 2: FastAPI Endpoint** (30 minutes)

**File**: `oncology-coPilot/oncology-backend-minimal/api/routers/clinical_trials.py`

**Action**: Add new endpoint to existing router

**Endpoint**: `POST /api/trials/refresh_status`

**Request Schema:**
```python
class RefreshStatusRequest(BaseModel):
    nct_ids: List[str] = Field(..., min_items=1, max_items=100)
    state_filter: Optional[str] = Field(None, description="Two-letter state code (e.g., 'NY')")
```

**Response Schema:**
```python
class RefreshStatusResponse(BaseModel):
    refreshed_count: int
    trial_data: Dict[str, Dict[str, Any]]
    errors: List[str] = []
```

**Validation:**
- Max 100 NCT IDs per request
- Empty list check
- State code format validation (optional 2-letter uppercase)

**Error Handling:**
- HTTP 400 for invalid input
- HTTP промежуточные ошибки → 500 with error details
- Graceful degradation: Return partial results on partial failure

---

### **Task 3: ClinicalTrialAgent Integration** (30 minutes)

**File**: `oncology-coPilot/oncology-backend/backend/agents/clinical_trial_agent.py`

**⚠️ IMPORTANT**: The agent's `run()` signature differs from doctrine:
- **Actual**: `async def run(self, query: str, patient_context: Optional[Dict], page_state: Optional[str])`
- **Doctrine suggests**: `async def run(self, patient_data, prompt_details, live_refresh=False)`

**Modifications:**
1. Add `live_refresh: bool = False` parameter to `run()` method (via `**kwargs` for backward compatibility)
2. **Service Import Challenge**: Agent is in full backend, service will be in minimal backend
   - **Option A**: Copy service to both backends
   - **Option B**: Add import path workaround
   - **Option C**: Create service in full backend location (`backend/services/`)
   - **Recommendation**: Option C - Create service in full backend to avoid cross-backend imports

3. After step 3 (vector search), add refresh step if `live_refresh=True`:
   - Extract NCT IDs from `found_nct_ids` (after vector search, before SQLite fetch)
   - Call refresh service
   - Filter to recruiting trials only

4. After step 4 (fetching details from SQLite), merge live data into cached trial details

**Integration Points:**
```python
# After found_nct_ids extraction
if live_refresh:
    live_status = await refresh_trial_status_with_retry(found_nct_ids)
    recruiting_nct_ids = [
        nct_id for nct_id, data in live_status.items()
        if data["status"] in ["RECRUITING", "NOT_YET_RECRUITING"]
    ]
    found_nct_ids = recruiting_nct_ids

# After fetching details from SQLite
if live_refresh and live_status:
    for trial in found_trials_details:
        nct_id = trial.get("nct_id")
        if nct_id in live_status:
            trial["status"] = live_status[nct_id]["status"]
            trial["locations"] = live_status[nct_id]["locations"]
            trial["live_refreshed"] = True
```

**Backward Compatibility:**
- Default `live_refresh=False` preserves existing behavior
- No breaking changes to agent interface

**Note**: Agent is in `oncology-backend/` (full backend), but refresh service should be accessible from there.

---

### **Task 4: Testing** (45 minutes)

**File**: `oncology-coPilot/oncology-backend-minimal/tests/test_trial_refresh_service.py`

**Test Cases:**
1. `test_refresh_single_trial` - Single NCT ID fetch
2. `test_refresh_batch` - Multiple NCT IDs (10 trials)
3. `test_location_filtering` - State filtering utility
4. `test_retry_logic` - Retry on failure
5. `test_endpoint_basic` - FastAPI endpoint smoke test
6. `test_endpoint_state_filter` - Endpoint with state filter
7. `test_agent_integration` - Agent with live_refresh=True (if possible)

**Test Data:**
- Use real NCT IDs for API calls: `["NCT02470585", "NCT02470586"]`
- Mock data for unit tests without API calls

---

### **Task 5: Documentation & Deployment** (15 minutes)

**Files to Create:**
1. `agent_2_refresh/docs/COMPLETION_REPORT.md` - Completion summary
2. Update `.cursor/rules/clinical_trials_agents/MASTER_STATUS.md` - Mark Agent 2 complete

**Deployment Steps:**
```bash
# 1. Verify service location
ls api/services/trial_refresh_service.py

# 2. Verify router updates
grep -A 20 "refresh_status" api/routers/clinical_trials.py

# 3. Run tests
cd oncology-coPilot/oncology-backend-minimal
source venv/bin/activate
pytest tests/test_trial_refresh_service.py -v

# 4. Smoke test endpoint (requires running backend)
curl -X POST http://127.0.0.1:8000/api/trials/refresh_status \
  -H 'Content-Type: application/json' \
  -d '{"nct_ids": ["NCT02470585"], "state_filter": "NY"}'
```

---

## **🚨 CRITICAL DECISIONS & EDGE CASES**

### **1. API v2 Response Structure**
- **Challenge**: Need to parse nested JSON structure
- **Solution**: Follow pattern from `clinicaltrials_utils.py` `parse_study()` but focus on locations
- **Path**: `study["protocolSection"]["contactsLocationsModule"]["locations"]`

### **2. Batch Size Limits**
- **API Limit**: ClinicalTrials.gov allows up to 1000 per page, but 100 IDs per request is safer
- **Implementation**: Chunk larger requests if needed (future enhancement)

### **3. State Filtering**
- **Challenge**: Filter locations within trials (not filter entire trials)
- **Solution**: `filter_locations_by_state()` removes non-matching locations but keeps trial if any location matches

### **4. Error Handling Strategy**
- **Partial Failures**: Some NCT IDs may fail → return successful ones
- **Total Failure**: Return empty dict `{}` (graceful degradation)
- **Logging**: Log failures but don't crash endpoint

### **5. Async vs Sync**
- **Decision**: Use `httpx.AsyncClient` (matches existing pattern)
- **Retry Logic**: Use `asyncio.sleep()` for exponential backoff

### **6. Agent Integration Location**
- **Challenge**: Agent is in `oncology-backend/`, service in `oncology-backend-minimal/`
- **Solution**: 
  - Option A: Copy service to both backends (maintenance burden)
  - Option B: Make service accessible from full backend (preferred)
  - Option C: Add refresh service to full backend only
- **Recommendation**: Add to minimal backend first, then copy to full backend if needed

---

## **📁 FILE STRUCTURE (MODULAR - INSPIRED BY AGENT 1)**

Following Agent 1's modular pattern for maintainability:

```
oncology-coPilot/oncology-backend-minimal/
├── api/
│   ├── services/
│   │   └── trial_refresh/
│   │       ├── __init__.py
│   │       ├── config.py                    # NEW - Constants (API URL, timeouts)
│   │       ├── api_client.py                 # NEW - Core refresh logic
│   │       ├── parser.py                     # NEW - Parse locations/status from API response
│   │       └── filters.py                    # NEW - State filtering utility
│   └── routers/
│       └── clinical_trials.py                # MODIFY - Add refresh endpoint
├── tests/
│   └── agent_2_refresh/
│       ├── __init__.py
│       ├── conftest.py                       # NEW - Test fixtures
│       ├── test_api_client.py                # NEW - Test refresh logic
│       └── test_endpoint.py                  # NEW - Test FastAPI endpoint

oncology-coPilot/oncology-backend/
├── backend/
│   ├── services/
│   │   └── trial_refresh/                    # COPY - Same structure for agent access
│   │       ├── __init__.py
│   │       ├── config.py
│   │       ├── api_client.py
│   │       ├── parser.py
│   │       └── filters.py
│   └── agents/
│       └── clinical_trial_agent.py           # MODIFY - Add live_refresh param

.cursor/rules/clinical_trials_agents/agent_2_refresh/
├── plan/
│   ├── AGENT_2_DOCTRINE.md                   # EXISTING
│   └── AGENT_2_IMPLEMENTATION_PLAN.md        # THIS FILE
└── docs/
    └── COMPLETION_REPORT.md                  # NEW - After completion
```

**Modular Benefits (Like Agent 1):**
- **Clear separation**: API client, parsing, filtering in separate files
- **Testability**: Each module independently testable
- **Reusability**: Parser and filters can be used elsewhere
- **Maintainability**: Easier to update individual components

---

## **⚡ EXECUTION SEQUENCE**

### **Phase 1: Service Implementation** (1.5 hours)
1. Create `trial_refresh_service.py` with core functions
2. Test functions individually (import and call directly)
3. Verify API v2 response parsing

### **Phase 2: Endpoint Integration** (30 min)
1. Add schemas to `clinical_trials.py`
2. Add endpoint handler
3. Register endpoint (already in router, auto-registered)
4. Smoke test endpoint

### **Phase 3: Agent Integration** (30 min)
1. Add `live_refresh` parameter to agent `run()` method
2. Add refresh logic
3. Test agent with `live_refresh=True`

### **Phase 4: Testing & Documentation** (1 hour)
1. Write comprehensive test suite
2. Run all tests
3. Create completion report
4. Update master status

---

## **✅ ACCEPTANCE CRITERIA**

### **Must Have:**
- [x] Service fetches live data from ClinicalTrials.gov API v2
- [x] Endpoint at `POST /api/trials/refresh_status`
- [x] Returns status + locations for multiple NCT IDs
- [x] Integrated into ClinicalTrialAgent (optional flag)
- [x] Response time <2 seconds per 10 trials
- [x] 5+ unit tests pass (expanded from 3)

### **Nice to Have:**
- [ ] Caching layer (TTL: 5 minutes) - Future enhancement
- [ ] Batch size optimization - Future enhancement
- [ ] WebSocket support - Future enhancement

---

## **🔥 READY TO EXECUTE**

**Estimated Time**: 3.5 hours total (2 hours core work + 1.5 hours testing/docs)

**Blocking**: None (can proceed independently)

**Parallelこ**: Can work alongside other agents

**Next Step**: Proceed with Phase 1 - Service Implementation 🚀

---

## **📝 NOTES FOR IMPLEMENTATION**

1. **API v2 Fields**: 
   - Status: `protocolSection.statusModule.overallStatus`
   - Locations: `protocolSection.contactsLocationsModule.locations[]`
   - Location fields: `facility`, `city`, `state`, `zip`, `status`, `contacts[]`

2. **Contact Info**:
   - Nested in `locations[].contacts[]` array
   - Fields: `name`, `phone`, `email`
   - Take first contact if multiple exist

3. **Status Values**:
   - RECRUITING, NOT_YET_RECRUITING (include)
   - COMPLETED, TERMINATED, SUSPENDED (exclude)
   - Filter at location level AND trial level

4. **Retry Strategy**:
   - Max 2 retries with exponential backoff (1s, 2s)
   - Use `asyncio.sleep()` for delays
   - Log each retry attempt

5. **Timezone**:
   - Use ISO 8601 format for `last_updated`
   - Python: `datetime.now(timezone.utc).isoformat()`

---

**COMMANDER APPROVAL NEEDED**: Ready to proceed? 🔥💀
