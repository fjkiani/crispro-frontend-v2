# MODULE 6: PDF EXPORT UTILITY

## **Purpose**

Simple browser-based PDF export for trial summaries (top 10 trials).

---

## **File Location**

`oncology-coPilot/oncology-frontend/src/utils/exportTrialsPDF.js`

---

## **Function Specification**

### **Signature:**
```javascript
export const exportTrialsPDF = (trials: Trial[]) => void;
```

### **Behavior:**
- Creates HTML content with trial summaries
- Opens new window with HTML
- Triggers browser print dialog
- User can save as PDF from print dialog

---

## **Implementation Details**

### **HTML Template:**
- Header: "Clinical Trials Summary"
- Trial count: "Found X matching trials"
- Top 10 trials (limit for PDF size)
- Per trial:
  - Title (h3)
  - NCT ID, Status, Phase
  - Eligibility summary (if available)
  - Top 3 locations with contact info
- Footer: Generated date, source, RUO disclaimer

### **Styling:**
```css
body { font-family: Arial, sans-serif; margin: 40px; }
h1 { color: #1976d2; }
.trial { margin-bottom: 30px; page-break-inside: avoid; }
.trial-header { background: #f5f5f5; padding: 10px; }
.location { margin-left: 20px; margin-top: 10px; }
.footer { margin-top: 50px; font-size: 12px; color: #666; }
```

---

## **Acceptance Criteria**

- [ ] Function accepts trials array
- [ ] Generates HTML with trial summaries
- [ ] Opens print dialog in new window
- [ ] Contains top 10 trials
- [ ] Includes locations and contact info
- [ ] Includes RUO disclaimer
- [ ] Handles empty trials array gracefully

---

## **Example Usage**

```jsx
import { exportTrialsPDF } from '../utils/exportTrialsPDF';

<Button
  variant="outlined"
  onClick={() => exportTrialsPDF(filteredResults)}
  disabled={filteredResults.length === 0}
>
  📄 Export PDF
</Button>
```

---

## **Future Enhancement**

- Use `jsPDF` or `pdfmake` for client-side PDF generation
- Add patient context header
- Add filters summary section

---

## **Testing Requirements**

- Unit test: HTML generation
- Unit test: Handles empty array
- Integration test: Print dialog opens
- E2E test: User exports PDF successfully

---

**ESTIMATED TIME:** 30 minutes

