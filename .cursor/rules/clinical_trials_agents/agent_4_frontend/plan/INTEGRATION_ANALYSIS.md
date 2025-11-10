# AGENT 4: FRONTEND INTEGRATION ANALYSIS

## **Current State Assessment**

---

## **📄 EXISTING PAGES**

### **1. ResearchPortal.jsx** (Simpler)
**Location:** `oncology-frontend/src/pages/ResearchPortal/ResearchPortal.jsx`

**Current Features:**
- ✅ Basic search interface
- ✅ PubMed + ClinicalTrials search types
- ✅ Uses SearchBar and ResultsDisplay components
- ⚠️ Uses **MOCK endpoint**: `/api/research/clinicaltrials/search` (returns placeholder data)
- ⚠️ Simple structure (108 lines) - easier to enhance

**API Endpoint:**
```javascript
POST /api/research/clinicaltrials/search
Body: { criteria: { "query.cond": "ovarian cancer" } }
Response: Mock data (placeholder)
```

### **2. Research.jsx** (More Complex)
**Location:** `oncology-frontend/src/pages/Research.jsx`

**Current Features:**
- ✅ Advanced search with patient context
- ✅ Uses **REAL endpoint**: `/api/search-trials` (uses ClinicalTrialAgent)
- ✅ Patient data integration
- ✅ Task management (Kanban)
- ✅ Multiple views (Search, Patient Match, Cohort Lab)
- ⚠️ More complex (417 lines) - higher risk of breaking changes

**API Endpoint:**
```javascript
POST /api/search-trials
Body: { query: string, patient_context: {...}, page_state?: string }
Response: { success: true, data: { found_trials: [...], next_page_state: "..." } }
```

---

## **🔍 BACKEND ENDPOINTS ANALYSIS**

### **Existing Endpoints:**
1. **`/api/search-trials`** (REAL) ✅
   - Uses `ClinicalTrialAgent` with ChromaDB/SQLite search
   - Returns trials with full metadata
   - Supports pagination

2. **`/api/research/clinicaltrials/search`** (MOCK) ⚠️
   - Returns placeholder data
   - Should be replaced with real search

3. **`/api/trial-details/{trial_id}`** ✅
   - Returns single trial details

### **Missing Endpoints (Need to Create):**
1. **`/api/trials/refresh_status`** ⏳
   - Service exists: `backend/services/trial_refresh/`
   - Need to expose as FastAPI endpoint
   - Function: `refresh_trial_status(nct_ids: List[str])`

---

## **🎯 INTEGRATION STRATEGY**

### **Option A: Enhance ResearchPortal.jsx (RECOMMENDED)**

**Pros:**
- ✅ Simpler starting point (108 lines)
- ✅ Clean structure for adding features
- ✅ Can switch to real `/api/search-trials` endpoint easily
- ✅ Lower risk of breaking existing complex workflows

**Plan:**
1. Switch `ResearchPortal.jsx` from mock endpoint → real `/api/search-trials`
2. Add Agent 4 modules (Filters, Refresh, Locations, PDF)
3. Keep it focused on clinical trials only (remove PubMed if needed)
4. Use as dedicated "Clinical Trial Finder" page

**Cons:**
- Need to update endpoint (one-line change)
- Two research pages (Research.jsx stays complex, ResearchPortal becomes simple)

---

### **Option B: Enhance Research.jsx**

**Pros:**
- ✅ Already uses real endpoint
- ✅ Already has patient context integration
- ✅ Single page for all research needs

**Cons:**
- ⚠️ More complex (417 lines) - higher breakage risk
- ⚠️ Already has multiple features (tasks, patient match, cohort lab)
- ⚠️ Requires careful integration to not break existing workflows

---

## **✅ RECOMMENDATION: Option A (Enhance ResearchPortal.jsx)**

**Rationale:**
- Clean separation: ResearchPortal = Clinical Trial Finder (Agent 4 focus)
- Research.jsx = Advanced research with patient context (keep separate)
- Easier to test and deploy incrementally
- Lower risk of regression

---

## **🔧 IMPLEMENTATION PLAN FOR RESEARCHPORTAL.JSX**

### **Phase 1: Switch to Real Endpoint (15 min)**

**Current:**
```javascript
url = `${API_BASE_URL}/clinicaltrials/search`;
body = { criteria: parseClinicalTrialsQuery(query) };
```

**New:**
```javascript
url = `${API_BASE_URL}/search-trials`;  // Use real endpoint
body = { query: query };  // Simple query string, backend handles parsing
```

**Changes:**
- Update `API_BASE_URL` to point to main backend (not `/api/research`)
- Simplify query format (backend handles parsing)
- Update response parsing (expects `{ success, data: { found_trials } }`)

---

### **Phase 2: Add Agent 4 Modules (Following Step-by-Step Plan)**

Same as documented in `IMPLEMENTATION/step_by_step.md`:

1. Module 7: Refresh Hook
2. Module 3: Location Card
3. Module 1: Trial Filters
4. Module 2: Refresh Button
5. Module 4: Enhanced ResultsDisplay
6. Module 6: PDF Export
7. Module 5: Integration

---

### **Phase 3: Backend Endpoint (If Missing)**

**Need to Create:** `/api/trials/refresh_status`

**File:** Add to `oncology-backend/main.py` or create new router

**Code:**
```python
from backend.services.trial_refresh import refresh_trial_status_with_retry

@app.post("/api/trials/refresh_status")
async def refresh_status_endpoint(request: Dict[str, Any] = Body(...)):
    """
    Refresh live trial status and locations from ClinicalTrials.gov API v2.
    
    Body: {
        "nct_ids": ["NCT12345", "NCT67890"],
        "state_filter": "NY"  # Optional
    }
    """
    nct_ids = request.get("nct_ids", [])
    state_filter = request.get("state_filter")
    
    if not nct_ids:
        raise HTTPException(status_code=400, detail="nct_ids required")
    
    # Refresh trials
    refreshed_data = await refresh_trial_status_with_retry(nct_ids)
    
    # Filter by state if requested
    if state_filter:
        # Apply state filter logic
        filtered_data = {
            nct_id: data for nct_id, data in refreshed_data.items()
            if any(loc.get("state") == state_filter for loc in data.get("locations", []))
        }
        refreshed_data = filtered_data
    
    return {
        "refreshed_count": len(refreshed_data),
        "trial_data": refreshed_data
    }
```

---

## **📊 DATA FLOW INTEGRATION**

### **Search Flow (Current → Enhanced):**
```
1. User enters query → SearchBar
2. handleSearch() → POST /api/search-trials
3. Backend: ClinicalTrialAgent → ChromaDB/SQLite search
4. Returns: { success: true, data: { found_trials: [...] } }
5. ResultsDisplay shows trials
```

### **Enhanced Flow (With Agent 4):**
```
1. User enters query → SearchBar
2. handleSearch() → POST /api/search-trials
3. Results → Apply filters (Disease/Phase/State)
4. ResultsDisplay shows trials with LocationCard
5. User clicks "Refresh Status" → POST /api/trials/refresh_status
6. Backend: refresh_trial_status_with_retry()
7. Results updated with live status + locations
8. User clicks "Export PDF" → Browser print dialog
```

---

## **🔗 COMPONENT INTEGRATION MAP**

### **ResearchPortal.jsx Structure:**
```jsx
<div>
  <h1>Clinical Trials Research</h1>
  
  {/* NEW: Filters */}
  <TrialFilters 
    filters={filters}
    onFiltersChange={handleFiltersChange}
  />
  
  {/* Existing: Search Bar */}
  <SearchBar onSearch={handleSearch} />
  
  {/* NEW: Action Buttons */}
  <RefreshStatusButton ... />
  <PDFExportButton ... />
  
  {/* Enhanced: Results Display */}
  <ResultsDisplay 
    results={filteredResults}  // Changed from searchResults
    loading={isLoading}
    error={error}
  />
</div>
```

---

## **📋 DATA STRUCTURE EXPECTATIONS**

### **Trial Object from `/api/search-trials`:**
```typescript
{
  nct_id: string;
  title: string;
  status: string;
  phase: string;  // "PHASE2, PHASE3" format
  disease_category?: string;  // From Agent 1
  disease_subcategory?: string;  // From Agent 1
  locations_data?: string;  // JSON string from Agent 1
  // ... other fields
}
```

### **Trial Object After Refresh:**
```typescript
{
  ...originalTrial,
  status: "RECRUITING",  // Updated
  locations_data: JSON.stringify([...newLocations]),  // Updated
  live_refreshed: true  // Flag
}
```

---

## **🚀 QUICK START INTEGRATION**

### **Step 1: Update ResearchPortal.jsx Endpoint (5 min)**
```javascript
// Change from:
const API_BASE_URL = 'http://localhost:8000/api/research';

// To:
const API_BASE_URL = import.meta.env.VITE_API_ROOT || 'http://localhost:8000';

// Update handleSearch for clinicaltrials:
if (type === 'clinicaltrials') {
  url = `${API_BASE_URL}/api/search-trials`;  // Real endpoint
  body = { query: query };  // Simple query
  // Remove parseClinicalTrialsQuery logic
}
```

### **Step 2: Add Filter State (10 min)**
```javascript
const [filters, setFilters] = useState({
  diseaseCategory: '',
  phase: [],
  state: ''
});
```

### **Step 3: Build Modules (Follow Step-by-Step Plan)**
- Start with Module 7 (Refresh Hook) - no dependencies
- Then Module 3 (Location Card) - standalone
- Continue in dependency order

---

## **✅ ACCEPTANCE CRITERIA**

### **Must Have:**
- [ ] ResearchPortal.jsx uses real `/api/search-trials` endpoint
- [ ] Filters apply to search results
- [ ] Refresh button calls `/api/trials/refresh_status`
- [ ] Locations display in ResultsDisplay
- [ ] PDF export works
- [ ] No breaking changes to existing search

### **Backend Requirements:**
- [ ] `/api/trials/refresh_status` endpoint created
- [ ] Endpoint uses existing `refresh_trial_status_with_retry` function
- [ ] Endpoint supports state filtering

---

## **🎯 DECISION: Use ResearchPortal.jsx**

**Why:**
- Cleaner starting point
- Easier to test incrementally
- Can become dedicated "Clinical Trial Finder"
- Lower risk than modifying complex Research.jsx

**Next Steps:**
1. Review this analysis
2. Confirm using ResearchPortal.jsx
3. Start with Phase 1 (endpoint switch)
4. Then follow modular implementation plan

---

**STATUS: READY FOR INTEGRATION** ✅

