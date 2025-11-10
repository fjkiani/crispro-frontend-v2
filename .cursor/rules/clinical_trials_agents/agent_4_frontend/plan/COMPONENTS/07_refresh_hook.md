# MODULE 7: REFRESH STATUS HOOK

## **Purpose**

Custom React hook for calling `/api/trials/refresh_status` endpoint.

---

## **File Location**

`oncology-coPilot/oncology-frontend/src/hooks/useTrialRefresh.js`

---

## **Hook Specification**

### **Signature:**
```javascript
const {
  refreshing,
  error,
  refreshStatus
} = useTrialRefresh();
```

### **Method:**
```javascript
refreshStatus(nctIds: string[], stateFilter?: string | null): Promise<RefreshResponse>
```

---

## **Implementation Details**

### **Hook Structure:**
```javascript
import { useState, useCallback } from 'react';

export const useTrialRefresh = () => {
  const [refreshing, setRefreshing] = useState(false);
  const [error, setError] = useState(null);

  const refreshStatus = useCallback(async (nctIds, stateFilter = null) => {
    setRefreshing(true);
    setError(null);

    try {
      const response = await fetch('/api/trials/refresh_status', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          nct_ids: nctIds,
          state_filter: stateFilter
        })
      });

      if (!response.ok) {
        throw new Error(`HTTP ${response.status}`);
      }

      const data = await response.json();
      return data;

    } catch (err) {
      setError(err.message);
      throw err;
    } finally {
      setRefreshing(false);
    }
  }, []);

  return {
    refreshing,
    error,
    refreshStatus
  };
};
```

---

## **Acceptance Criteria**

- [ ] Hook returns `refreshing`, `error`, `refreshStatus`
- [ ] `refreshStatus` calls correct API endpoint
- [ ] Loading state managed correctly
- [ ] Error handling works
- [ ] Returns promise with refresh data

---

## **Example Usage**

```jsx
import { useTrialRefresh } from '../hooks/useTrialRefresh';

const RefreshStatusButton = ({ nctIds, onRefreshComplete }) => {
  const { refreshing, error, refreshStatus } = useTrialRefresh();

  const handleClick = async () => {
    try {
      const data = await refreshStatus(nctIds);
      onRefreshComplete(data);
    } catch (err) {
      alert('Failed to refresh: ' + err.message);
    }
  };

  return (
    <Button onClick={handleClick} disabled={refreshing}>
      {refreshing ? 'Refreshing...' : '🔄 Refresh Status'}
    </Button>
  );
};
```

---

## **Testing Requirements**

- Unit test: Hook returns correct shape
- Unit test: API call with mock
- Unit test: Error handling
- Integration test: Used in RefreshStatusButton

---

**ESTIMATED TIME:** 15 minutes

