# MODULE 5: RESEARCH.JSX INTEGRATION

## **Purpose**

Wire filters, refresh button, and PDF export into main `Research.jsx` page.

---

## **File Location**

`oncology-coPilot/oncology-frontend/src/pages/Research.jsx` (modify existing)

---

## **Changes Required**

### **1. Import New Components:**
```jsx
import { TrialFilters } from '../components/research/TrialFilters';
import { RefreshStatusButton } from '../components/research/RefreshStatusButton';
import { exportTrialsPDF } from '../utils/exportTrialsPDF';
```

### **2. Add Filter State:**
```jsx
const [filters, setFilters] = useState({
  diseaseCategory: '',
  phase: [],
  state: ''
});

const [filteredResults, setFilteredResults] = useState([]);
```

### **3. Add Filter Application Logic:**
```jsx
const applyFilters = useCallback((filterState, allResults) => {
  if (!allResults || allResults.length === 0) return [];

  return allResults.filter(trial => {
    // Disease category
    if (filterState.diseaseCategory && 
        trial.disease_category !== filterState.diseaseCategory) {
      return false;
    }

    // Phase (multi-select)
    if (filterState.phase.length > 0) {
      const trialPhases = trial.phase?.split(', ').map(p => p.trim()) || [];
      const hasMatch = filterState.phase.some(selectedPhase => 
        trialPhases.includes(selectedPhase)
      );
      if (!hasMatch) return false;
    }

    // State filter
    if (filterState.state) {
      try {
        const locations = JSON.parse(trial.locations_data || '[]');
        const hasStateMatch = locations.some(loc => loc.state === filterState.state);
        if (!hasStateMatch) return false;
      } catch (e) {
        return false; // Invalid JSON = no match
      }
    }

    return true;
  });
}, []);
```

### **4. Update Search Handler:**
```jsx
const handleSearch = async (query) => {
  // ... existing search logic ...
  
  // After setSearchResults(result.data.found_trials):
  const filtered = applyFilters(filters, result.data.found_trials);
  setFilteredResults(filtered);
};
```

### **5. Add Filter Change Handler:**
```jsx
const handleFiltersChange = (newFilters) => {
  setFilters(newFilters);
  const filtered = applyFilters(newFilters, searchResults);
  setFilteredResults(filtered);
};
```

### **6. Add Refresh Handler:**
```jsx
const handleRefreshComplete = (refreshedData) => {
  const updated = searchResults.map(trial => {
    const live = refreshedData.trial_data[trial.nct_id];
    if (live) {
      return {
        ...trial,
        status: live.status,
        locations_data: JSON.stringify(live.locations),
        live_refreshed: true
      };
    }
    return trial;
  });
  
  setSearchResults(updated);
  const filtered = applyFilters(filters, updated);
  setFilteredResults(filtered);
};
```

### **7. Add UI Components:**
```jsx
return (
  <div>
    <h1>Clinical Trials Research</h1>

    {/* NEW: Filters */}
    <TrialFilters
      filters={filters}
      onFiltersChange={handleFiltersChange}
      onClearFilters={() => {
        setFilters({ diseaseCategory: '', phase: [], state: '' });
        setFilteredResults(searchResults);
      }}
    />

    {/* Existing search bar */}
    <SearchBar onSearch={handleSearch} />

    {/* NEW: Refresh button */}
    {filteredResults.length > 0 && (
      <RefreshStatusButton
        nctIds={filteredResults.map(t => t.nct_id)}
        stateFilter={filters.state || null}
        onRefreshComplete={handleRefreshComplete}
      />
    )}

    {/* NEW: PDF Export button */}
    {filteredResults.length > 0 && (
      <Button
        variant="outlined"
        onClick={() => exportTrialsPDF(filteredResults)}
        sx={{ ml: 2 }}
      >
        📄 Export PDF
      </Button>
    )}

    {/* Results - use filteredResults instead of searchResults */}
    <ResultsDisplay 
      results={filteredResults}
      loading={isLoading}
      error={error}
    />
  </div>
);
```

---

## **Acceptance Criteria**

- [ ] Filters appear above search bar
- [ ] Filters apply to search results correctly
- [ ] Refresh button appears when results exist
- [ ] Refresh updates trial data correctly
- [ ] PDF export button appears when results exist
- [ ] PDF export works
- [ ] ResultsDisplay shows filtered results
- [ ] No breaking changes to existing search

---

## **Testing Requirements**

- Integration test: Filters + search workflow
- Integration test: Refresh + filter workflow
- E2E test: Complete user flow (search → filter → refresh → export)

---

**ESTIMATED TIME:** 30 minutes

