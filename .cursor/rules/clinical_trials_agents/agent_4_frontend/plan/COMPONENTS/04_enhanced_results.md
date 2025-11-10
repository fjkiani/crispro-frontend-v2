# MODULE 4: ENHANCED RESULTSDISPLAY

## **Purpose**

Enhance existing `ResultsDisplay.jsx` to integrate LocationCard and show live refresh indicators.

---

## **File Location**

`oncology-coPilot/oncology-frontend/src/components/research/ResultsDisplay.jsx` (modify existing)

---

## **Changes Required**

### **1. Import LocationCard:**
```jsx
import { LocationCard } from './LocationCard';
```

### **2. Add Location Display Section:**
Add after trial title/basic info, before eligibility criteria:

```jsx
{/* Location Display */}
{trial.locations_data && (() => {
  try {
    const locations = JSON.parse(trial.locations_data);
    if (locations && locations.length > 0) {
      return (
        <div style={{ marginTop: '15px', marginBottom: '15px' }}>
          <strong>📍 Locations ({locations.length}):</strong>
          {locations.slice(0, 3).map((loc, idx) => (
            <LocationCard key={idx} location={loc} />
          ))}
          {locations.length > 3 && (
            <div style={{ marginTop: '5px', fontSize: '14px', color: '#666' }}>
              + {locations.length - 3} more locations
            </div>
          )}
        </div>
      );
    }
  } catch (e) {
    console.error('Error parsing locations_data:', e);
    return null;
  }
  return null;
})()}
```

### **3. Add Live Refresh Indicator:**
Add after locations section:

```jsx
{trial.live_refreshed && (
  <Alert severity="success" sx={{ mt: 1, mb: 1 }}>
    ✅ Live status refreshed
  </Alert>
)}
```

---

## **Acceptance Criteria**

- [ ] Locations display below trial title
- [ ] Shows first 3 locations, "+ X more" if more exist
- [ ] Handles invalid/missing `locations_data` gracefully
- [ ] Live refresh indicator appears when `live_refreshed = true`
- [ ] No breaking changes to existing functionality
- [ ] Styling consistent with existing design

---

## **Testing Requirements**

- Integration test: Locations display correctly
- Integration test: Live refresh indicator
- Regression test: Existing ResultsDisplay functionality still works

---

**ESTIMATED TIME:** 30 minutes

