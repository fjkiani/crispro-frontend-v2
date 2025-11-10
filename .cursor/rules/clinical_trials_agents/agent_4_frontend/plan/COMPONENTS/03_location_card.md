# MODULE 3: LOCATION CARD COMPONENT

## **Purpose**

Display component for trial location data with facility, contact info, and status badge.

---

## **File Location**

`oncology-coPilot/oncology-frontend/src/components/research/LocationCard.jsx`

---

## **Component Specification**

### **Props:**
```typescript
interface LocationCardProps {
  location: {
    facility: string;
    city: string;
    state: string;
    zip: string;
    status?: string;
    contact_name?: string;
    contact_phone?: string;
    contact_email?: string;
  };
  index?: number;
  showMax?: number; // Show first N locations, then "+ X more"
}
```

### **Location Data Format:**
```json
{
  "facility": "Memorial Sloan Kettering Cancer Center",
  "city": "New York",
  "state": "NY",
  "zip": "10065",
  "status": "RECRUITING",
  "contact_name": "Dr. Smith",
  "contact_phone": "212-555-1234",
  "contact_email": "smith@mskcc.org"
}
```

---

## **Implementation Details**

### **Dependencies:**
```jsx
import React from 'react';
import { Card, CardContent, Chip, Box, Typography } from '@mui/material';
import LocationOnIcon from '@mui/icons-material/LocationOn';
import PhoneIcon from '@mui/icons-material/Phone';
import EmailIcon from '@mui/icons-material/Email';
```

### **Component Structure:**
- Card with facility name (bold)
- City, State, ZIP on second line
- Status badge (RECRUITING = green, NOT_YET_RECRUITING = orange)
- Contact phone with icon (if available)
- Contact email with icon (if available)
- Gray background (#f5f5f5)
- Rounded corners

### **Status Badge Colors:**
- `RECRUITING` → Green (#4caf50)
- `NOT_YET_RECRUITING` → Orange (#ff9800)
- Other → Gray

---

## **Acceptance Criteria**

- [ ] Displays facility name prominently
- [ ] Shows city, state, zip
- [ ] Status badge with correct color
- [ ] Contact phone displays if available
- [ ] Contact email displays if available
- [ ] Handles missing fields gracefully
- [ ] Styling matches design spec

---

## **Example Usage**

```jsx
import { LocationCard } from '../components/research/LocationCard';

const ResultsDisplay = ({ trial }) => {
  const locations = JSON.parse(trial.locations_data || '[]');

  return (
    <div>
      <h3>{trial.title}</h3>
      
      {locations.length > 0 && (
        <div style={{ marginTop: '15px' }}>
          <strong>📍 Locations ({locations.length}):</strong>
          {locations.slice(0, 3).map((loc, idx) => (
            <LocationCard key={idx} location={loc} />
          ))}
          {locations.length > 3 && (
            <div>+ {locations.length - 3} more locations</div>
          )}
        </div>
      )}
    </div>
  );
};
```

---

## **Edge Cases**

- Missing `locations_data` → Show nothing (don't crash)
- Empty locations array → Show nothing
- Missing contact info → Show location only (no contact section)
- Invalid JSON → Show error message or skip

---

## **Testing Requirements**

- Unit test: Renders location data correctly
- Unit test: Handles missing contact info
- Unit test: Status badge colors
- Unit test: Invalid JSON handling
- Integration test: Displays in ResultsDisplay

---

**ESTIMATED TIME:** 30 minutes

