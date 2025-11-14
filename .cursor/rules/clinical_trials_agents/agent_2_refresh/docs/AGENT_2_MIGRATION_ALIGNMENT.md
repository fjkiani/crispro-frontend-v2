# 🔄 AGENT 2: MIGRATION ALIGNMENT WITH CLINICAL TRIALS SEARCH

## **📋 MISSION**

Align Agent 2 (Trial Refresh Service) with the same migration pattern used for Clinical Trials Search migration, ensuring `oncology-backend-minimal` is the **main production backend**.

---

## **✅ CURRENT STATUS AUDIT**

### **Agent 2 Implementation (Following Migration Pattern):**

**✅ Service in Minimal Backend:**
- Location: `oncology-backend-minimal/api/services/trial_refresh/`
- Modular architecture (config, api_client, parser, filters)
- Self-contained, no dependencies on full backend
- Production-ready with comprehensive tests

**✅ Endpoint in Minimal Backend:**
- Location: `oncology-backend-minimal/api/routers/clinical_trials.py`
- Endpoint: `POST /api/trials/refresh_status`
- Lightweight service, no AgentInterface bloat

**⚠️ Agent Integration Issue:**
- Full backend agent (`oncology-backend/backend/agents/clinical_trial_agent.py`) tries to import `backend.services.trial_refresh`
- This service doesn't exist in full backend (correctly - it's in minimal backend now)

---

## **🎯 ALIGNMENT WITH MIGRATION PATTERN**

### **Pattern from Clinical Trials Search Migration:**

1. ✅ **Self-Contained Service** - No dependencies on main backend
2. ✅ **Lightweight Service** - No AgentInterface bloat (259 lines vs 567 lines)
3. ✅ **Modular Architecture** - Clear separation (config, client, parser, filters)
4. ✅ **Production-Ready** - Comprehensive tests, error handling, retry logic
5. ✅ **Database Pattern** - Uses centralized `database_connections.py` (if needed)

### **Agent 2 Compliance:**

| Criterion | Status | Notes |
|-----------|--------|-------|
| **Service in Minimal Backend** | ✅ Yes | `api/services/trial_refresh/` |
| **No Main Backend Dependencies** | ✅ Yes | Fully self-contained |
| **Lightweight Service** | ✅ Yes | 4 modular files, ~300 lines total |
| **Comprehensive Tests** | ✅ Yes | 18+ tests in `tests/agent_2_refresh/` |
| **Production-Ready** | ✅ Yes | Retry logic, error handling, logging |
| **Follows Same Pattern** | ✅ Yes | Matches ClinicalTrialSearchService structure |

---

## **🔧 REQUIRED ADJUSTMENTS**

### **1. Agent Integration Strategy**

**Problem:** Full backend agent needs refresh service but it's in minimal backend.

**Solution Options:**

#### **Option A: HTTP API Call (Recommended - Matches Migration Pattern)**

Full backend agent calls minimal backend API endpoint:

```python
# In full backend agent
async def run(self, query: str, ..., live_refresh: bool = False, **kwargs):
    # ... existing vector search logic ...
    
    if live_refresh and nct_ids:
        try:
            # Call minimal backend API (HTTP)
            minimal_backend_url = os.getenv("MINIMAL_BACKEND_URL", "http://localhost:8000")
            async with httpx.AsyncClient() as client:
                response = await client.post(
                    f"{minimal_backend_url}/api/trials/refresh_status",
                    json={"nct_ids": nct_ids},
                    timeout=10.0
                )
                response.raise_for_status()
                live_status = response.json()["trial_data"]
                
                # Filter to recruiting trials
                recruiting_nct_ids = [
                    nct_id for nct_id, data in live_status.items()
                    if data.get("status", "").upper() in ["RECRUITING", "NOT_YET_RECRUITING"]
                ]
                
                if recruiting_nct_ids:
                    nct_ids = recruiting_nct_ids
        except Exception as e:
            logger.error(f"Live refresh API call failed: {e}")
            # Graceful degradation - continue with cached data
```

**Benefits:**
- ✅ Maintains clean separation (minimal backend = API, full backend = agent)
- ✅ Matches migration pattern (full backend calls minimal backend API)
- ✅ No code duplication
- ✅ Single source of truth (service only in minimal backend)

#### **Option B: Copy Service (Not Recommended)**

Copy service to full backend for agent access.

**Downsides:**
- ❌ Code duplication
- ❌ Two places to maintain
- ❌ Violates DRY principle
- ❌ Doesn't match migration pattern

---

### **2. Router Consolidation Check**

**Current State:**
- `clinical_trials.py` - Has refresh endpoint ✅
- `trials.py` - Also exists (need to check for duplication)

**Action:** Verify no duplicate endpoints, consolidate if needed.

---

### **3. Update Documentation**

Update Agent 2 documentation to reflect:
- ✅ Service is in minimal backend (production)
- ✅ Full backend agent calls minimal backend API
- ✅ Follows same pattern as Clinical Trials Search migration

---

## **📊 ARCHITECTURE COMPARISON**

### **Clinical Trials Search (Migrated):**

```
Frontend → Minimal Backend → ClinicalTrialSearchService → AstraDB + SQLite
```

**Agent (Full Backend):**
- Deprecated search endpoint (backward compat only)
- AgentOrchestrator may still use it, but migration path clear

### **Trial Refresh (Agent 2 - Aligned):**

```
Frontend → Minimal Backend → TrialRefreshService → ClinicalTrials.gov API v2
                                                      (External API, no local storage needed)

Full Backend Agent (if needed):
  → HTTP call to Minimal Backend /api/trials/refresh_status
```

**Key Difference:** Refresh service doesn't need local storage (SQLite/AstraDB), it's a pure API proxy.

---

## **✅ VALIDATION CHECKLIST**

- [x] Service is in minimal backend
- [x] Endpoint is in minimal backend router
- [x] No dependencies on full backend
- [x] Comprehensive test suite exists
- [x] Follows modular architecture pattern
- [ ] Agent integration updated (Option A: HTTP API call)
- [ ] Documentation updated to reflect migration alignment
- [ ] Router consolidation verified (no duplicates)

---

## **🚀 RECOMMENDED NEXT STEPS**

### **Immediate (Today):**

1. **Update Full Backend Agent** - Change to HTTP API call (Option A)
2. **Verify Router Structure** - Check for duplicate endpoints
3. **Update Documentation** - Reflect migration alignment

### **Before Production:**

1. **Set Environment Variable** - `MINIMAL_BACKEND_URL` in full backend `.env`
2. **Test Integration** - Verify agent → minimal backend API call works
3. **Add Retry Logic** - HTTP calls should have retry (already in service, add in agent)

---

## **📝 CODE CHANGES NEEDED**

### **Full Backend Agent Update:**

Replace direct import with HTTP API call:

```python
# BEFORE (current - won't work):
from backend.services.trial_refresh import refresh_trial_status_with_retry

# AFTER (recommended - HTTP API call):
import httpx
minimal_backend_url = os.getenv("MINIMAL_BACKEND_URL", "http://localhost:8000")
async with httpx.AsyncClient(timeout=10.0) as client:
    response = await client.post(
        f"{minimal_backend_url}/api/trials/refresh_status",
        json={"nct_ids": nct_ids}
    )
    live_status = response.json()["trial_data"]
```

---

## **🎯 SUCCESS CRITERIA**

- ✅ Agent 2 service matches Clinical Trials Search migration pattern
- ✅ Full backend agent uses HTTP API call (not direct import)
- ✅ Minimal backend is single source of truth
- ✅ No code duplication between backends
- ✅ Production deployment ready

---

**STATUS:** ✅ **ALIGNED** - Just needs agent integration update to use HTTP API call instead of direct import.










