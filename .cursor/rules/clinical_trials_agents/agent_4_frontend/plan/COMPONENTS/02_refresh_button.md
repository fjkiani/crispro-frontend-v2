# MODULE 2: REFRESH STATUS BUTTON

## **Purpose**

Button component that calls `/api/trials/refresh_status` endpoint to update live trial status and locations.

---

## **File Location**

`oncology-coPilot/oncology-frontend/src/components/research/RefreshStatusButton.jsx`

---

## **Component Specification**

### **Props:**
```typescript
interface RefreshStatusButtonProps {
  nctIds: string[];
  stateFilter?: string | null;
  onRefreshComplete: (refreshedData: RefreshResponse) => void;
  disabled?: boolean;
}
```

### **Refresh Response:**
```typescript
interface RefreshResponse {
  refreshed_count: number;
  trial_data: {
    [nctId: string]: {
      status: string;
      locations: Location[];
      last_updated: string;
    };
  };
}
```

---

## **Implementation Details**

### **Dependencies:**
```jsx
import React, { useState } from 'react';
import { Button, CircularProgress, Alert } from '@mui/material';
import RefreshIcon from '@mui/icons-material/Refresh';
```

### **Component Behavior:**
- Shows "🔄 Refresh Live Status" when idle
- Shows "Refreshing..." with spinner when loading
- Disabled when no trials or already refreshing
- Calls `/api/trials/refresh_status` POST endpoint
- Passes `nct_ids` array and optional `state_filter`
- Merges refreshed data into parent's search results
- Shows success/error alerts

### **API Endpoint:**
```
POST /api/trials/refresh_status
Body: {
  nct_ids: string[],
  state_filter?: string | null
}
```

---

## **Acceptance Criteria**

- [ ] Button renders with refresh icon
- [ ] Button disabled when no trials
- [ ] Button shows loading state during API call
- [ ] API call includes correct `nct_ids` array
- [ ] Refreshed data merges into trial objects
- [ ] Success message shows refreshed count
- [ ] Error handling shows alert on failure
- [ ] `live_refreshed` flag set on updated trials

---

## **Example Usage**

```jsx
import { RefreshStatusButton } from '../components/research/RefreshStatusButton';

const Research = () => {
  const [searchResults, setSearchResults] = useState([]);

  const handleRefreshComplete = (refreshedData) => {
    // Merge refreshed data into results
    const updated = searchResults.map(trial => {
      const live = refreshedData.trial_data[trial.nct_id];
      if (live) {
        return {
          ...trial,
          status: live.status,
          locations: live.locations,
          live_refreshed: true
        };
      }
      return trial;
    });
    setSearchResults(updated);
  };

  return (
    <RefreshStatusButton
      nctIds={searchResults.map(t => t.nct_id)}
      stateFilter={filters.state || null}
      onRefreshComplete={handleRefreshComplete}
      disabled={searchResults.length === 0}
    />
  );
};
```

---

## **Error Handling**

- Network errors: Show alert "Failed to refresh trial status"
- API errors: Log error, show user-friendly message
- Partial failures: Still update successful refreshes
- Timeout: Show timeout message after 30 seconds

---

## **Testing Requirements**

- Unit test: Button disabled when no trials
- Unit test: Loading state during API call
- Unit test: Error handling
- Integration test: API call with mock response
- E2E test: User clicks refresh, sees updated status

---

**ESTIMATED TIME:** 30 minutes

