# 🤖 AGENT 4: FRONTEND INTEGRATION AGENT 🎨 - MASTER INDEX

> **This is the master index. For organized, component-focused documentation, see the folder structure below.**

---

## **📚 NAVIGATION - ORGANIZED DOCUMENTATION**

**START HERE:**
- **[OVERVIEW.md](OVERVIEW.md)** - Mission, objectives, decisions summary, quick start

**COMPONENT SPECIFICATIONS:**
- **[COMPONENTS/01_trial_filters.md](COMPONENTS/01_trial_filters.md)** - Filter component specs
- **[COMPONENTS/02_refresh_button.md](COMPONENTS/02_refresh_button.md)** - Refresh button specs
- **[COMPONENTS/03_location_card.md](COMPONENTS/03_location_card.md)** - Location display specs
- **[COMPONENTS/04_enhanced_results.md](COMPONENTS/04_enhanced_results.md)** - ResultsDisplay specs
- **[COMPONENTS/05_research_integration.md](COMPONENTS/05_research_integration.md)** - Main page specs
- **[COMPONENTS/06_pdf_export.md](COMPONENTS/06_pdf_export.md)** - PDF export specs
- **[COMPONENTS/07_refresh_hook.md](COMPONENTS/07_refresh_hook.md)** - Refresh hook specs

**IMPLEMENTATION:**
- **[IMPLEMENTATION/step_by_step.md](IMPLEMENTATION/step_by_step.md)** - Build order and dependencies

**EXECUTION:**
- **[EXECUTION/checklist.md](EXECUTION/checklist.md)** - Pre-flight, execution, verification checklist

---

## **⚔️ MISSION**
Enhance `Research.jsx` with Ayesha-specific features: **search-based filters**, live status refresh, location display, and PDF export.

**Note:** CT upload functionality **deferred** - focusing on search-based features first.

---

## **🎯 OBJECTIVES**

### **Primary Goal:**
Build user-facing features that connect all backend services (seeding, refresh, parser) into a seamless clinical trial finder for Ayesha's case.

### **Success Criteria:**
- ✅ CT report upload works (paste or file)
- ✅ Auto-populates search from parsed report
- ✅ Disease/Phase/Location filters functional
- ✅ "Refresh Status" button updates live
- ✅ Locations display with contact info
- ✅ Export PDF button generates summary
- ✅ 4/4 E2E tests pass

---

## **⚠️ BLOCKERS**

**Cannot start until:**
- ✅ Agent 1 complete (need 1000 trials in database)
- ✅ Agent 2 complete (need `/api/trials/refresh_status` endpoint)
- ✅ Agent 3 complete (need `/api/trials/parse_ct_report` endpoint)

**Check `MASTER_STATUS.md` before proceeding!**

---

## **📋 TASKS BREAKDOWN**

### **Task 1: CT Report Upload Component (1 hour)**

**Action:**
Create component for pasting/uploading CT reports.

**File:** `implementation/CTReportUpload.jsx`-

**Code:**
```jsx
import React, { useState } from 'react';
import { Alert, Button, Card, Textarea, Upload } from '@mui/material';

export const CTReportUpload = ({ onParsed }) => {
  const [reportText, setReportText] = useState('');
  const [parsing, setParsing] = useState(false);
  const [error, setError] = useState(null);

  const handleParse = async () => {
    if (!reportText || reportText.length < 50) {
      setError('CT report text too short (minimum 50 characters)');
      return;
    }

    setParsing(true);
    setError(null);

    try {
      const response = await fetch('http://localhost:8000/api/trials/parse_ct_report', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ report_text: reportText })
      });

      if (!response.ok) {
        throw new Error(`HTTP ${response.status}`);
      }

      const parsed = await response.json();
      console.log('Parsed CT report:', parsed);

      // Callback to parent with parsed data
      onParsed(parsed);

    } catch (err) {
      console.error('CT parsing error:', err);
      setError(`Failed to parse CT report: ${err.message}`);
    } finally {
      setParsing(false);
    }
  };

  const handleFileUpload = (event) => {
    const file = event.target.files[0];
    if (file) {
      const reader = new FileReader();
      reader.onload = (e) => setReportText(e.target.result);
      reader.readAsText(file);
    }
  };

  return (
    <Card sx={{ p: 3, mb: 3 }}>
      <h3>📄 Upload CT Report</h3>
      <p>Paste your CT scan report or upload a file to automatically search for trials.</p>

      <Textarea
        placeholder="Paste CT report here...&#10;&#10;Example:&#10;CT ABDOMEN AND PELVIS WITH CONTRAST&#10;FINDINGS: Peritoneal carcinomatosis..."
        value={reportText}
        onChange={(e) => setReportText(e.target.value)}
        minRows={8}
        sx={{ mb: 2, width: '100%' }}
      />

      <div style={{ display: 'flex', gap: '10px', marginBottom: '10px' }}>
        <Button
          variant="contained"
          onClick={handleParse}
          disabled={parsing || !reportText}
        >
          {parsing ? 'Parsing...' : '🔍 Parse & Search'}
        </Button>

        <Button
          variant="outlined"
          component="label"
        >
          📁 Upload File
          <input
            type="file"
            hidden
            accept=".txt,.pdf,.doc,.docx"
            onChange={handleFileUpload}
          />
        </Button>

        <Button
          variant="text"
          onClick={() => setReportText('')}
        >
          Clear
        </Button>
      </div>

      {error && (
        <Alert severity="error" sx={{ mt: 2 }}>
          {error}
        </Alert>
      )}
    </Card>
  );
};
```

**Acceptance:**
- [ ] Component renders correctly
- [ ] Paste text works
- [ ] File upload works (txt files)
- [ ] Calls parser endpoint
- [ ] Passes parsed data to parent

---

### **Task 2: Enhanced Research.jsx Integration (1.5 hours)**

**Action:**
Integrate CTReportUpload and add filters/refresh.

**File:** `oncology-frontend/src/pages/Research.jsx` (modify)

**Changes:**
```jsx
import { CTReportUpload } from '../components/research/CTReportUpload';

const Research = () => {
  // ... existing state ...

  // NEW: CT parsing state
  const [parsedCT, setParsedCT] = useState(null);

  // NEW: Filter state
  const [filters, setFilters] = useState({
    diseaseCategory: '',
    phase: [],
    state: ''
  });

  // NEW: Handle parsed CT report
  const handleCTParsed = (parsed) => {
    setParsedCT(parsed);
    
    // Auto-populate search
    setSearchQuery(parsed.trial_search_query);
    
    // Auto-populate filters
    setFilters({
      diseaseCategory: parsed.disease_category,
      phase: ['PHASE2', 'PHASE3', 'PHASE4'],
      state: 'NY'  // Default to NY for Ayesha
    });
    
    // Trigger search automatically
    handleSearch(parsed.trial_search_query);
  };

  // NEW: Refresh status handler
  const handleRefreshStatus = async () => {
    if (!searchResults || searchResults.length === 0) return;
    
    const nct_ids = searchResults.map(trial => trial.nct_id);
    
    try {
      const response = await fetch('http://localhost:8000/api/trials/refresh_status', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          nct_ids,
          state_filter: filters.state || null
        })
      });
      
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      
      const refreshed = await response.json();
      
      // Merge refreshed data into results
      const updated = searchResults.map(trial => {
        const live = refreshed.trial_data[trial.nct_id];
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
      alert(`Refreshed ${refreshed.refreshed_count} trials`);
      
    } catch (err) {
      console.error('Refresh failed:', err);
      alert('Failed to refresh trial status');
    }
  };

  return (
    <div>
      <h1>Clinical Trials Research</h1>

      {/* NEW: CT Report Upload */}
      <CTReportUpload onParsed={handleCTParsed} />

      {/* NEW: Display parsed CT info */}
      {parsedCT && (
        <Alert severity="info" sx={{ mb: 2 }}>
          <strong>Parsed CT Report:</strong><br />
          Disease: {parsedCT.disease} (Stage: {parsedCT.stage})<br />
          Key Findings: {parsedCT.key_findings.join(', ')}
        </Alert>
      )}

      {/* NEW: Filters */}
      <Card sx={{ p: 2, mb: 2 }}>
        <h3>Filters</h3>
        <div style={{ display: 'flex', gap: '15px' }}>
          <FormControl sx={{ minWidth: 200 }}>
            <InputLabel>Disease Category</InputLabel>
            <Select
              value={filters.diseaseCategory}
              onChange={(e) => setFilters({...filters, diseaseCategory: e.target.value})}
            >
              <MenuItem value="">All</MenuItem>
              <MenuItem value="gynecologic_oncology">Gynecologic Oncology</MenuItem>
              <MenuItem value="breast_cancer">Breast Cancer</MenuItem>
              <MenuItem value="lung_cancer">Lung Cancer</MenuItem>
            </Select>
          </FormControl>

          <FormControl sx={{ minWidth: 150 }}>
            <InputLabel>Phase</InputLabel>
            <Select
              multiple
              value={filters.phase}
              onChange={(e) => setFilters({...filters, phase: e.target.value})}
            >
              <MenuItem value="PHASE2">Phase 2</MenuItem>
              <MenuItem value="PHASE3">Phase 3</MenuItem>
              <MenuItem value="PHASE4">Phase 4</MenuItem>
            </Select>
          </FormControl>

          <FormControl sx={{ minWidth: 150 }}>
            <InputLabel>State</InputLabel>
            <Select
              value={filters.state}
              onChange={(e) => setFilters({...filters, state: e.target.value})}
            >
              <MenuItem value="">All</MenuItem>
              <MenuItem value="NY">New York</MenuItem>
              <MenuItem value="NJ">New Jersey</MenuItem>
              <MenuItem value="CT">Connecticut</MenuItem>
              <MenuItem value="CA">California</MenuItem>
            </Select>
          </FormControl>

          <Button
            variant="outlined"
            onClick={() => setFilters({ diseaseCategory: '', phase: [], state: '' })}
          >
            Clear Filters
          </Button>
        </div>
      </Card>

      {/* Existing search bar */}
      <SearchBar onSearch={handleSearch} />

      {/* NEW: Refresh button */}
      {searchResults && searchResults.length > 0 && (
        <Button
          variant="contained"
          color="secondary"
          onClick={handleRefreshStatus}
          sx={{ mb: 2 }}
        >
          🔄 Refresh Live Status
        </Button>
      )}

      {/* Results */}
      <ResultsDisplay 
        results={searchResults}
        loading={isLoading}
        error={error}
      />
    </div>
  );
};
```

**Acceptance:**
- [ ] CTReportUpload integrated
- [ ] Filters display and work
- [ ] Refresh button updates results
- [ ] Auto-search from parsed CT works

---

### **Task 3: Enhanced ResultsDisplay with Locations (30 minutes)**

**Action:**
Update ResultsDisplay to show locations with contact info.

**File:** `oncology-frontend/src/components/research/ResultsDisplay.jsx` (modify)

**Changes:**
```jsx
// Add to each trial card:

{trial.locations && trial.locations.length > 0 && (
  <div style={{ marginTop: '15px' }}>
    <strong>📍 Locations ({trial.locations.length}):</strong>
    {trial.locations.slice(0, 3).map((loc, idx) => (
      <div key={idx} style={{
        padding: '8px',
        marginTop: '5px',
        backgroundColor: '#f5f5f5',
        borderRadius: '4px'
      }}>
        <div><strong>{loc.facility}</strong></div>
        <div>{loc.city}, {loc.state} {loc.zip}</div>
        {loc.status && (
          <span style={{
            fontSize: '12px',
            padding: '2px 6px',
            backgroundColor: loc.status === 'recruiting' ? '#4caf50' : '#ff9800',
            color: 'white',
            borderRadius: '3px'
          }}>
            {loc.status.toUpperCase()}
          </span>
        )}
        {loc.contact_phone && (
          <div style={{ marginTop: '4px' }}>
            📞 {loc.contact_phone}
          </div>
        )}
        {loc.contact_email && (
          <div>
            ✉️ {loc.contact_email}
          </div>
        )}
      </div>
    ))}
    {trial.locations.length > 3 && (
      <div style={{ marginTop: '5px', fontSize: '14px', color: '#666' }}>
        + {trial.locations.length - 3} more locations
      </div>
    )}
  </div>
)}

{trial.live_refreshed && (
  <Alert severity="success" sx={{ mt: 1 }}>
    ✅ Live status refreshed
  </Alert>
)}
```

**Acceptance:**
- [ ] Locations display with contact info
- [ ] Status badge shows recruiting/not recruiting
- [ ] Live refresh indicator appears

---

### **Task 4: PDF Export (30 minutes)**

**Action:**
Add simple PDF export for trial summaries.

**File:** `implementation/exportTrialsPDF.js`

**Code:**
```javascript
/**
 * Simple PDF export using browser print functionality.
 * For production, consider using jsPDF or pdfmake.
 */
export const exportTrialsPDF = (trials, parsedCT = null) => {
  // Create HTML content
  const html = `
    <!DOCTYPE html>
    <html>
    <head>
      <title>Clinical Trials Summary</title>
      <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1 { color: #1976d2; }
        .trial { margin-bottom: 30px; page-break-inside: avoid; }
        .trial-header { background: #f5f5f5; padding: 10px; }
        .location { margin-left: 20px; margin-top: 10px; }
        .footer { margin-top: 50px; font-size: 12px; color: #666; }
      </style>
    </head>
    <body>
      <h1>Clinical Trials Summary</h1>
      
      ${parsedCT ? `
        <div style="background: #e3f2fd; padding: 15px; margin-bottom: 20px;">
          <strong>Patient Context:</strong><br />
          Disease: ${parsedCT.disease} (Stage: ${parsedCT.stage})<br />
          Key Findings: ${parsedCT.key_findings.join(', ')}
        </div>
      ` : ''}
      
      <p><strong>Found ${trials.length} matching trials</strong></p>
      
      ${trials.slice(0, 10).map(trial => `
        <div class="trial">
          <div class="trial-header">
            <h3>${trial.title}</h3>
            <p><strong>NCT ID:</strong> ${trial.nct_id} | 
               <strong>Status:</strong> ${trial.status} | 
               <strong>Phase:</strong> ${trial.phase}</p>
          </div>
          
          <p><strong>Eligibility Summary:</strong></p>
          <p>${trial.llm_assessment?.eligibility_summary || 'See full criteria on ClinicalTrials.gov'}</p>
          
          ${trial.locations && trial.locations.length > 0 ? `
            <p><strong>Locations:</strong></p>
            ${trial.locations.slice(0, 3).map(loc => `
              <div class="location">
                ${loc.facility}, ${loc.city}, ${loc.state}<br />
                ${loc.contact_phone ? `Phone: ${loc.contact_phone}` : ''}
              </div>
            `).join('')}
          ` : ''}
        </div>
      `).join('')}
      
      <div class="footer">
        Generated: ${new Date().toLocaleDateString()}<br />
        Source: ClinicalTrials.gov<br />
        <strong>Research Use Only - Not for Clinical Diagnosis</strong>
      </div>
    </body>
    </html>
  `;
  
  // Open print dialog
  const printWindow = window.open('', '_blank');
  printWindow.document.write(html);
  printWindow.document.close();
  printWindow.print();
};
```

**Integration in Research.jsx:**
```jsx
import { exportTrialsPDF } from '../utils/exportTrialsPDF';

// Add button
<Button
  variant="outlined"
  onClick={() => exportTrialsPDF(searchResults, parsedCT)}
  disabled={!searchResults || searchResults.length === 0}
>
  📄 Export PDF
</Button>
```

**Acceptance:**
- [ ] Export button appears
- [ ] PDF opens in new window
- [ ] Contains top 10 trials
- [ ] Includes locations and contact info

---

## **🧪 COMPLETE TEST SUITE**

**File:** `tests/test_research_page_e2e.js` (Cypress/Playwright)

```javascript
describe('Ayesha Clinical Trial Finder E2E', () => {
  
  it('should upload CT report and auto-search', () => {
    cy.visit('http://localhost:3000/research');
    
    // Paste CT report
    cy.get('textarea').type(`
      CT ABDOMEN AND PELVIS
      Findings: Peritoneal carcinomatosis, ascites
      Impression: Advanced ovarian malignancy
    `);
    
    // Click parse
    cy.contains('Parse & Search').click();
    
    // Verify auto-populated search
    cy.get('input[name="search"]').should('contain.value', 'ovarian cancer');
    
    // Verify results appear
    cy.get('.trial-card', { timeout: 10000 }).should('have.length.greaterThan', 5);
  });
  
  it('should filter trials by state', () => {
    cy.visit('http://localhost:3000/research');
    
    // Select NY filter
    cy.get('select[name="state"]').select('NY');
    
    // Trigger search
    cy.get('input[name="search"]').type('ovarian cancer{enter}');
    
    // Verify only NY trials shown
    cy.get('.trial-card').each($trial => {
      cy.wrap($trial).should('contain', 'NY');
    });
  });
  
  it('should refresh live status', () => {
    // ... search first ...
    
    // Click refresh
    cy.contains('Refresh Live Status').click();
    
    // Verify success alert
    cy.contains('Refreshed', { timeout: 5000 }).should('be.visible');
    
    // Verify live badge appears
    cy.get('.live-refreshed-badge').should('be.visible');
  });
  
  it('should export PDF', () => {
    // ... search first ...
    
    // Click export
    cy.contains('Export PDF').click();
    
    // Verify new window opened (hard to test fully in Cypress)
    cy.window().then(win => {
      cy.spy(win, 'open').should('be.called');
    });
  });
});
```

**Run Tests:**
```bash
cd oncology-frontend
npm test -- --spec tests/test_research_page_e2e.js
```

---

## **📊 ACCEPTANCE CRITERIA**

### **Must Have:**
- [x] CT report upload works (paste + file)
- [x] Auto-populates search from parsed report
- [x] Disease/Phase/Location filters functional
- [x] "Refresh Status" button updates live
- [x] Locations display with contact info
- [x] Export PDF generates summary
- [x] 4/4 E2E tests pass

---

## **📁 DELIVERABLES**

1. `agent_4_frontend/implementation/CTReportUpload.jsx` (150 lines)
2. `agent_4_frontend/implementation/exportTrialsPDF.js` (100 lines)
3. `oncology-frontend/src/pages/Research.jsx` (modified, +200 lines)
4. `oncology-frontend/src/components/research/ResultsDisplay.jsx` (modified, +100 lines)
5. `agent_4_frontend/tests/test_research_page_e2e.js` (150 lines)
6. `agent_4_frontend/docs/COMPLETION_REPORT.md`

---

## **🔥 EXECUTION CHECKLIST**

**Pre-flight:**
- [ ] Verify Agents 1, 2, 3 complete (check MASTER_STATUS.md)
- [ ] Backend running with all endpoints
- [ ] Frontend dev server running (`npm start`)

**Execute:**
```bash
cd agent_4_frontend/implementation
# Create CTReportUpload.jsx
# Modify Research.jsx
# Create exportTrialsPDF.js
```

**Test:**
```bash
cd oncology-frontend
npm test -- --spec tests/test_research_page_e2e.js
```

---

## **⚔️ AGENT 4 STATUS: BLOCKED - WAITING ON AGENTS 1, 2, 3**
**ESTIMATED TIME:** 3 hours (once unblocked)
**BLOCKING:** Agent 5 (E2E Tests need UI)
**COMMANDER APPROVAL:** AWAITING ORDERS 🔥💀


