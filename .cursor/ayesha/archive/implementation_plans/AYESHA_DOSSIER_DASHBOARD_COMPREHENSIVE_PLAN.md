# ⚔️ AYESHA DOSSIER DASHBOARD - COMPREHENSIVE BUILD PLAN ⚔️

**Date**: January 13, 2025  
**Mission**: Transform Ayesha's trial dossiers from text dumps into an intelligent, agentic dashboard with PDF export, Kanban organization, and Co-Pilot integration

---

## 🎯 EXECUTIVE SUMMARY

### **Current State**
- ✅ 60+ Commander-grade dossiers generated (`.cursor/ayesha/zo_fresh_dossiers/`)
- ✅ Backend API operational (`/api/ayesha/dossiers/*`)
- ✅ Basic frontend pages (AyeshaDossierBrowser, AyeshaDossierDetail)
- ❌ **Problem**: DossierDetail is just markdown dump (no formatting)
- ❌ **Problem**: No PDF export capability
- ❌ **Problem**: No dashboard organization (manual one-by-one viewing)
- ❌ **Problem**: No filtering/organization UI (backend has it, frontend doesn't use it)
- ❌ **Problem**: No agentic management (Co-Pilot can't help organize trials)

### **Target State**
- ✅ Beautifully formatted dossier viewer with sections, tables, decision trees
- ✅ PDF export with react-pdf (professional clinical reports)
- ✅ Intelligent dashboard with Kanban board (organize by status: Interested → Contacted → Enrolled → Rejected)
- ✅ Frontend filtering UI (tier, score, location, phase, LLM status)
- ✅ Agentic Co-Pilot integration ("Move trial to Interested", "Export all Top-Tier", "Find trials with HER2 requirement")
- ✅ Real-time organization and workflow management

---

## 📋 PHASE 1: ENHANCE DOSSIER DETAIL VIEW (2-3 hours)

### **Task 1.1: Install react-pdf** (15 min)
```bash
cd oncology-coPilot/oncology-frontend
npm install @react-pdf/renderer
```

**Documentation**: https://react-pdf.org/

### **Task 1.2: Create Formatted Dossier Viewer Component** (1.5 hours)

**File**: `oncology-coPilot/oncology-frontend/src/components/ayesha/DossierViewer.jsx`

**Features**:
- Parse markdown into structured sections
- Extract key sections:
  - **Header**: Trial name, NCT ID, tier, score
  - **Pipeline Assessment**: Stage-by-stage breakdown with icons
  - **Eligibility Table**: Formatted table with color-coded status (✅/⚠️/❌)
  - **LLM Analysis**: Expandable sections (Drug Mechanism, Location Logistics, Risk-Benefit, SOC Comparison)
  - **Decision Trees**: Render ASCII flowcharts (if present in markdown)
  - **Strategic Scenarios**: Best/Likely/Challenge scenarios
  - **WIN-WIN Analysis**: Visual comparison cards
  - **Critical Gates**: Color-coded gate analysis

**Component Structure**:
```jsx
<DossierViewer>
  <DossierHeader />          // Trial metadata, tier badge, score
  <PipelineAssessment />     // Stage 1-6 breakdown
  <EligibilityTable />       // Formatted eligibility checklist
  <LLMAnalysis />            // Expandable LLM sections
  <DecisionTrees />          // ASCII flowchart rendering
  <StrategicScenarios />     // Best/Likely/Challenge cards
  <WinWinAnalysis />         // Comparison visualization
  <CriticalGates />          // Gate analysis
  <ActionButtons />          // Export PDF, Share, Add to Kanban
</DossierViewer>
```

### **Task 1.3: Create PDF Export Component** (1 hour)

**File**: `oncology-coPilot/oncology-frontend/src/components/ayesha/DossierPDFExporter.jsx`

**Using react-pdf**:
```jsx
import { Document, Page, Text, View, StyleSheet, PDFDownloadLink } from '@react-pdf/renderer';

// Convert markdown sections to PDF components
// Professional clinical report format
// Include header, footer, page numbers
// Preserve tables, decision trees, scenarios
```

**Features**:
- Professional clinical report styling
- Header with patient name, date, NCT ID
- Footer with page numbers, "Research Use Only" disclaimer
- Preserve all markdown sections (tables, lists, code blocks)
- Export button in DossierDetail page

### **Task 1.4: Update AyeshaDossierDetail Page** (30 min)

**Replace**: Current markdown dump with `DossierViewer` component

**Add**:
- PDF export button (using `PDFDownloadLink`)
- Share functionality (copy link)
- "Add to Kanban" button (Phase 2 integration)

---

## 📋 PHASE 2: DASHBOARD WITH KANBAN ORGANIZATION (4-5 hours)

### **Task 2.1: Create Trial Kanban Board Component** (2 hours)

**File**: `oncology-coPilot/oncology-frontend/src/components/ayesha/TrialKanbanBoard.jsx`

**Reuse/Improve**: Existing `KanbanBoard.jsx` component

**Columns**:
1. **🔍 Interested** (default for Top-Tier)
2. **📞 Contacted** (site coordinator reached)
3. **📋 Screening** (eligibility review in progress)
4. **✅ Enrolled** (patient enrolled)
5. **❌ Rejected** (not eligible or declined)

**Trial Cards**:
- Display: NCT ID, title (truncated), tier badge, match score, phase
- Quick actions: View dossier, Export PDF, Move to column, Remove
- Color coding: Top-Tier (green), Good-Tier (blue), Acceptable (gray)

**Features**:
- Drag-and-drop between columns (using `@dnd-kit`)
- Persist state to localStorage
- Auto-populate "Interested" with Top-Tier trials on load
- Filter by tier, score, phase (sidebar filters)

### **Task 2.2: Create Dashboard Page** (1.5 hours)

**File**: `oncology-coPilot/oncology-frontend/src/pages/AyeshaTrialDashboard.jsx`

**Layout**:
```
┌─────────────────────────────────────────────────────────┐
│  Header: Ayesha's Trial Intelligence Dashboard         │
│  Stats: 60 dossiers | 10 Top-Tier | 50 Good-Tier      │
├─────────────────────────────────────────────────────────┤
│  [Filters Sidebar]  │  [Kanban Board - Main View]     │
│  - Tier filter      │  ┌─────┐ ┌─────┐ ┌─────┐       │
│  - Score range      │  │Int. │ │Cont.│ │Scrn.│       │
│  - Phase filter     │  │     │ │     │ │     │       │
│  - Location filter  │  │[T1] │ │[T2] │ │[T3] │       │
│  - LLM status       │  │[T4] │ │     │ │     │       │
│  - Search           │  └─────┘ └─────┘ └─────┘       │
│                     │  ┌─────┐ ┌─────┐                │
│                     │  │Enr. │ │Rej. │                │
│                     │  │     │ │     │                │
│                     │  │[T5] │ │[T6] │                │
│                     │  └─────┘ └─────┘                │
└─────────────────────────────────────────────────────────┘
```

**Features**:
- Sidebar filters (tier, score, phase, location, LLM status, search)
- Kanban board (main view)
- Stats cards (total, by tier, by column)
- Export all button (batch PDF export)
- Co-Pilot integration (Phase 3)

### **Task 2.3: Create Trial Card Component** (1 hour)

**File**: `oncology-coPilot/oncology-frontend/src/components/ayesha/TrialKanbanCard.jsx`

**Props**:
- `trial`: Dossier data
- `columnId`: Current column
- `onMove`: Callback for drag-and-drop
- `onView`: Navigate to detail page
- `onExport`: Export PDF

**Display**:
- NCT ID (clickable → detail page)
- Title (truncated, hover for full)
- Tier badge (color-coded)
- Match score (progress bar)
- Phase badge
- LLM analysis indicator (✅ if has LLM)
- Quick actions (View, Export, Remove)

### **Task 2.4: Integrate Backend Filtering** (30 min)

**Connect**: Frontend filters to backend `/api/ayesha/dossiers/list` endpoint

**Query Parameters**:
- `tier`: TOP_TIER | GOOD_TIER | ALL
- `min_score`: 0.0-1.0
- `has_llm`: true | false
- `limit`: number

**Update**: Dashboard to use filtered results from backend

---

## 📋 PHASE 3: AGENTIC CO-PILOT INTEGRATION (2-3 hours)

### **Task 3.1: Add Trial Management Actions to Co-Pilot** (1.5 hours)

**File**: `oncology-coPilot/oncology-frontend/src/components/CoPilot/Actions/trialManagementActions.js`

**New Actions**:
1. **"Move trial to [column]"**
   - Intent: "Move NCT06619236 to Interested"
   - Action: Update Kanban board state
   - Endpoint: None (local state)

2. **"Export all Top-Tier trials"**
   - Intent: "Export all top-tier dossiers as PDF"
   - Action: Batch PDF export
   - Endpoint: None (client-side)

3. **"Find trials with [requirement]"**
   - Intent: "Find trials that require HER2 testing"
   - Action: Filter dossiers by eligibility criteria
   - Endpoint: `/api/ayesha/dossiers/list` with filters

4. **"Show me trials in [location]"**
   - Intent: "Show trials at MSK"
   - Action: Filter by location keywords
   - Endpoint: `/api/ayesha/dossiers/list` with search

5. **"What's the best trial for Ayesha?"**
   - Intent: General trial recommendation
   - Action: Return top-tier trial with highest score
   - Endpoint: `/api/ayesha/dossiers/list?tier=TOP_TIER&limit=1`

### **Task 3.2: Update Q2C Router** (30 min)

**File**: `oncology-coPilot/oncology-frontend/src/components/CoPilot/Q2CRouter/intents.js`

**Add Intent Patterns**:
```javascript
trial_management: {
  patterns: [
    /move.*trial.*to/i,
    /organize.*trial/i,
    /export.*(trial|dossier)/i,
    /find.*trial.*(with|requiring|needing)/i,
    /show.*trial.*(at|in|near)/i,
    /best.*trial.*(for|for ayesha)/i
  ],
  endpoint: null,  // Client-side actions
  description: 'Manage and organize clinical trials',
  confidence: 'high'
}
```

### **Task 3.3: Create Co-Pilot Trial Context** (1 hour)

**File**: `oncology-coPilot/oncology-frontend/src/components/CoPilot/context/TrialContext.jsx`

**Context Provides**:
- Current dashboard state (Kanban columns, filters)
- Available trials (from API)
- Actions: `moveTrial`, `exportTrial`, `filterTrials`, `getBestTrial`

**Integration**: Co-Pilot can read/write trial organization state

---

## 📋 PHASE 4: ADVANCED FILTERING UI (2 hours)

### **Task 4.1: Create Filter Sidebar Component** (1 hour)

**File**: `oncology-coPilot/oncology-frontend/src/components/ayesha/TrialFilterSidebar.jsx`

**Filters**:
1. **Tier**: Radio buttons (All / Top-Tier / Good-Tier)
2. **Score Range**: Slider (0.0 - 1.0)
3. **Phase**: Checkboxes (Phase I / Phase II / Phase III / N/A)
4. **Location**: Multi-select (NYC Metro / MSK / Mount Sinai / etc.)
5. **LLM Status**: Toggle (Has LLM Analysis / All)
6. **Search**: Text input (NCT ID, keywords)

**Features**:
- Real-time filtering (debounced)
- Filter count badges
- Clear all button
- Save filter presets (localStorage)

### **Task 4.2: Create Stats Cards Component** (30 min)

**File**: `oncology-coPilot/oncology-frontend/src/components/ayesha/TrialStatsCards.jsx`

**Display**:
- Total dossiers
- By tier (Top/Good/Acceptable)
- By column (Interested/Contacted/etc.)
- Average score
- LLM-enhanced count

### **Task 4.3: Integrate with Backend Pipeline Metadata** (30 min)

**Extract**: Pipeline stage metadata from dossiers

**Display**:
- Stage-by-stage breakdown (how many passed each stage)
- Rejection reasons (if available)
- Composite score breakdown (weights visualization)

---

## 📋 PHASE 5: PDF EXPORT ENHANCEMENT (1-2 hours)

### **Task 5.1: Create Professional PDF Template** (1 hour)

**File**: `oncology-coPilot/oncology-frontend/src/components/ayesha/DossierPDFTemplate.jsx`

**Using react-pdf**:
- Professional clinical report styling
- Header: Patient name, date, NCT ID, tier
- Footer: Page numbers, "Research Use Only" disclaimer
- Sections: All markdown sections preserved
- Tables: Formatted eligibility tables
- Decision trees: ASCII flowcharts (monospace font)
- Color coding: Tier badges, status indicators

### **Task 5.2: Batch PDF Export** (30 min)

**Feature**: Export multiple dossiers as single PDF or separate files

**Implementation**:
- Select trials from Kanban
- Generate PDF for each
- Zip download or single combined PDF

---

## 📋 PHASE 6: TESTING & POLISH (1 hour)

### **Task 6.1: E2E Testing**
- Test dossier detail rendering
- Test PDF export (single and batch)
- Test Kanban drag-and-drop
- Test Co-Pilot actions
- Test filtering

### **Task 6.2: UI Polish**
- Loading states
- Error handling
- Empty states
- Responsive design (mobile/tablet/desktop)

---

## 🎯 TECHNICAL ARCHITECTURE

### **Component Hierarchy**
```
AyeshaTrialDashboard (Page)
├── TrialFilterSidebar (Filters)
├── TrialStatsCards (Stats)
└── TrialKanbanBoard (Main View)
    └── TrialKanbanCard[] (Trial Cards)
        └── Quick Actions (View, Export, Remove)

AyeshaDossierDetail (Page)
├── DossierViewer (Formatted Display)
│   ├── DossierHeader
│   ├── PipelineAssessment
│   ├── EligibilityTable
│   ├── LLMAnalysis
│   ├── DecisionTrees
│   ├── StrategicScenarios
│   └── CriticalGates
└── DossierPDFExporter (PDF Export)
```

### **State Management**
- **Kanban State**: localStorage (persist across sessions)
- **Filter State**: URL query params (shareable links)
- **Trial Data**: React state (from API)
- **Co-Pilot Context**: React Context (shared with Co-Pilot)

### **API Integration**
- **List Dossiers**: `GET /api/ayesha/dossiers/list?tier=...&min_score=...`
- **Get Detail**: `GET /api/ayesha/dossiers/detail/{nct_id}`
- **Export**: `GET /api/ayesha/dossiers/export/{nct_id}?format=markdown`
- **Stats**: `GET /api/ayesha/dossiers/stats`

---

## 📊 SUCCESS CRITERIA

### **Phase 1 (Dossier Viewer)**
- ✅ Dossier detail page shows formatted sections (not raw markdown)
- ✅ PDF export works (single dossier)
- ✅ All markdown sections render correctly (tables, lists, code blocks)

### **Phase 2 (Dashboard)**
- ✅ Kanban board displays trials in 5 columns
- ✅ Drag-and-drop works (persist to localStorage)
- ✅ Filters work (tier, score, phase, location, LLM)
- ✅ Stats cards show accurate counts

### **Phase 3 (Co-Pilot)**
- ✅ Co-Pilot can move trials between columns
- ✅ Co-Pilot can export trials
- ✅ Co-Pilot can filter/search trials
- ✅ Co-Pilot can recommend best trial

### **Phase 4 (Advanced Features)**
- ✅ Filter sidebar with all options
- ✅ Stats cards with breakdowns
- ✅ Pipeline metadata visualization

### **Phase 5 (PDF)**
- ✅ Professional PDF export (single)
- ✅ Batch PDF export (multiple)
- ✅ All sections preserved in PDF

---

## 🚀 EXECUTION ORDER

1. **Phase 1** (2-3 hours) - Fix dossier viewer + PDF export
2. **Phase 2** (4-5 hours) - Build dashboard + Kanban
3. **Phase 3** (2-3 hours) - Co-Pilot integration
4. **Phase 4** (2 hours) - Advanced filtering
5. **Phase 5** (1-2 hours) - PDF enhancement
6. **Phase 6** (1 hour) - Testing & polish

**Total**: 12-16 hours

---

## 📚 REFERENCES

- **react-pdf**: https://react-pdf.org/
- **Kanban Component**: `oncology-coPilot/oncology-frontend/src/components/KanbanBoard.jsx`
- **Co-Pilot**: `oncology-coPilot/oncology-frontend/src/components/CoPilot/`
- **Backend API**: `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_dossiers.py`
- **Pipeline Config**: `oncology-coPilot/oncology-backend-minimal/api/services/trial_intelligence/config.py`
- **Sample Dossiers**: `.cursor/ayesha/zo_fresh_dossiers/`

---

**READY TO EXECUTE** ⚔️

