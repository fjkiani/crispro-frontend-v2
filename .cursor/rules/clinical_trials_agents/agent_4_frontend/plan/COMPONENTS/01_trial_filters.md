# MODULE 1: TRIAL FILTERS COMPONENT

## **Purpose**

Filter component for Disease Category, Phase, and State that integrates with existing search results.

---

## **File Location**

`oncology-coPilot/oncology-frontend/src/components/research/TrialFilters.jsx`

---

## **Component Specification**

### **Props:**
```typescript
interface TrialFiltersProps {
  filters: {
    diseaseCategory: string;
    phase: string[];
    state: string;
  };
  onFiltersChange: (filters: TrialFilters) => void;
  onClearFilters: () => void;
}
```

### **Filter Options:**

**Disease Category:**
- All (empty string)
- Gynecologic Oncology (`gynecologic_oncology`)
- Breast Cancer (`breast_cancer`)
- Lung Cancer (`lung_cancer`)

**Phase (Multi-select):**
- Phase 2 (`PHASE2`)
- Phase 3 (`PHASE3`)
- Phase 4 (`PHASE4`)

**State:**
- All (empty string)
- New York (`NY`)
- New Jersey (`NJ`)
- Connecticut (`CT`)
- California (`CA`)

---

## **Implementation Details**

### **Dependencies:**
```jsx
import React from 'react';
import {
  Card,
  CardContent,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  Button,
  Chip,
  Box
} from '@mui/material';
```

### **Component Structure:**
- Card container with padding
- Three FormControl Selects (Disease, Phase, State)
- Phase uses `multiple` prop for multi-select
- Clear Filters button
- Selected phase chips display

### **Filter Logic:**
- Filters apply to `searchResults` array in parent component
- Filter by `disease_category` field (from database)
- Filter by `phase` field (contains any selected phase)
- Filter by `locations_data` JSON (state field)

---

## **Acceptance Criteria**

- [ ] Component renders with three filter dropdowns
- [ ] Disease category filter works (single select)
- [ ] Phase filter works (multi-select with chips)
- [ ] State filter works (single select)
- [ ] Clear Filters button resets all filters
- [ ] `onFiltersChange` callback fires on filter change
- [ ] Selected phases display as chips
- [ ] Styling matches existing MUI theme

---

## **Example Usage**

```jsx
import { TrialFilters } from '../components/research/TrialFilters';

const Research = () => {
  const [filters, setFilters] = useState({
    diseaseCategory: '',
    phase: [],
    state: ''
  });

  const handleFiltersChange = (newFilters) => {
    setFilters(newFilters);
    // Apply filters to searchResults
    applyFilters(newFilters);
  };

  return (
    <TrialFilters
      filters={filters}
      onFiltersChange={handleFiltersChange}
      onClearFilters={() => setFilters({ diseaseCategory: '', phase: [], state: '' })}
    />
  );
};
```

---

## **Filter Application Logic**

**File:** Integrate into `Research.jsx`

```jsx
const applyFilters = (filters, allResults) => {
  return allResults.filter(trial => {
    // Disease category filter
    if (filters.diseaseCategory && trial.disease_category !== filters.diseaseCategory) {
      return false;
    }

    // Phase filter (multi-select)
    if (filters.phase.length > 0) {
      const trialPhases = trial.phase?.split(', ').map(p => p.trim()) || [];
      const hasMatch = filters.phase.some(selectedPhase => 
        trialPhases.includes(selectedPhase)
      );
      if (!hasMatch) return false;
    }

    // State filter
    if (filters.state) {
      const locations = JSON.parse(trial.locations_data || '[]');
      const hasStateMatch = locations.some(loc => loc.state === filters.state);
      if (!hasStateMatch) return false;
    }

    return true;
  });
};
```

---

## **Testing Requirements**

- Unit test: Filter state changes
- Unit test: Clear filters resets state
- Integration test: Filters apply to search results correctly
- E2E test: User can filter trials by disease/phase/state

---

**ESTIMATED TIME:** 45 minutes

