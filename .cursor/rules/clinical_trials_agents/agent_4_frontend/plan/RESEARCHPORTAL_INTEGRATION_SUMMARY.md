# RESEARCHPORTAL.JSX INTEGRATION SUMMARY

## **Quick Reference: What We're Building**

---

## **📄 TARGET PAGE: ResearchPortal.jsx**

**Current State:**
- Simple search page (108 lines)
- Uses **mock endpoint** `/api/research/clinicaltrials/search`
- Basic SearchBar + ResultsDisplay

**After Integration:**
- Real search endpoint `/api/search-trials`
- **TrialFilters** (Disease/Phase/State)
- **RefreshStatusButton** (live status updates)
- **LocationCard** (in ResultsDisplay)
- **PDF Export** button
- Full Agent 4 feature set

---

## **🔄 ENDPOINT SWITCH (Required First Step)**

**Current:**
```javascript
// ResearchPortal.jsx line 44
url = `${API_BASE_URL}/clinicaltrials/search`;
body = { criteria: parseClinicalTrialsQuery(query) };
```

**New:**
```javascript
// Switch to real endpoint
const API_ROOT = import.meta.env.VITE_API_ROOT || 'http://localhost:8000';

// In handleSearch for clinicaltrials:
if (type === 'clinicaltrials') {
  url = `${API_ROOT}/api/search-trials`;  // Real endpoint
  body = { query: query };  // Simple query string
  // Remove parseClinicalTrialsQuery logic
}
```

**Response Format Change:**
```javascript
// Old (mock): Direct array
setSearchResults(data);

// New (real): Wrapped response
if (result.success && result.data?.found_trials) {
  setSearchResults(result.data.found_trials);
}
```

---

## **📦 COMPONENTS TO ADD**

### **1. TrialFilters**
- Disease category dropdown
- Phase multi-select with chips
- State dropdown
- Clear Filters button

### **2. RefreshStatusButton**
- Calls `/api/trials/refresh_status`
- Shows loading state
- Updates trial status + locations

### **3. LocationCard**
- Displays facility, city, state, zip
- Contact info (phone/email)
- Status badge (RECRUITING/NOT_YET_RECRUITING)

### **4. Enhanced ResultsDisplay**
- Integrates LocationCard component
- Shows locations below trial title
- Shows live refresh indicator

### **5. PDF Export Utility**
- Browser-based print dialog
- Includes top 10 trials
- Trial summaries + locations

---

## **🔗 BACKEND ENDPOINT (Need to Create)**

**File:** `oncology-backend/main.py`

**Endpoint:** `POST /api/trials/refresh_status`

**Uses:** Existing `backend/services/trial_refresh/refresh_trial_status_with_retry()`

**See:** `COMPONENTS/08_backend_endpoint.md` for full implementation

---

## **📋 DATA EXPECTATIONS**

### **Trial Object from `/api/search-trials`:**
```typescript
{
  nct_id: "NCT12345",
  title: "Trial Title",
  status: "RECRUITING",
  phase: "PHASE2, PHASE3",  // Comma-separated
  disease_category: "gynecologic_oncology",  // From Agent 1
  disease_subcategory: "ovarian_cancer",  // From Agent 1
  locations_data: "[{\"facility\":\"MSK\",\"state\":\"NY\"}]"  // JSON string
}
```

### **After Refresh:**
```typescript
{
  ...originalTrial,
  status: "RECRUITING",  // Updated from API
  locations_data: "[{...new locations...}]",  // Updated
  live_refreshed: true  // Flag for UI
}
```

---

## **🎯 IMPLEMENTATION ORDER**

1. **Switch Endpoint** (5 min) - Update ResearchPortal.jsx to use `/api/search-trials`
2. **Module 8** (20 min) - Create backend `/api/trials/refresh_status` endpoint
3. **Module 7** (15 min) - useTrialRefresh hook
4. **Module 3** (30 min) - LocationCard component
5. **Module 1** (45 min) - TrialFilters component
6. **Module 2** (30 min) - RefreshStatusButton
7. **Module 4** (30 min) - Enhance ResultsDisplay
8. **Module 6** (30 min) - PDF export utility
9. **Module 5** (30 min) - Wire everything into ResearchPortal.jsx

**Total: ~3.5 hours**

---

## **✅ ACCEPTANCE CRITERIA**

- [ ] ResearchPortal uses real `/api/search-trials` endpoint
- [ ] Filters work (Disease/Phase/State)
- [ ] Refresh button updates live status
- [ ] Locations display correctly
- [ ] PDF export works
- [ ] No breaking changes to existing search

---

## **🚨 RISKS & MITIGATIONS**

**Risk:** Breaking existing search functionality  
**Mitigation:** Test endpoint switch separately before adding other features

**Risk:** Missing `locations_data` field  
**Mitigation:** Check Agent 1 seeded data has `locations_data` populated

**Risk:** Backend endpoint not created  
**Mitigation:** Create Module 8 (backend endpoint) first, test with curl

---

**READY TO START** ✅

