# Frontend Development Plan: MBD4+TP53 Clinical Dossier

**Date**: January 28, 2025  
**Purpose**: Comprehensive frontend plan for scaling MBD4+TP53 analysis capabilities to clinicians  
**Status**: Planning Phase

---

## 🔍 Manager Review (January 28, 2025)

### Overall Assessment: ✅ **APPROVED** with corrections below

| Aspect | Status | Notes |
|--------|--------|-------|
| **Scope Alignment** | ✅ Correct | Honest framing: mechanism alignment, NOT outcome predictions |
| **Safety Mechanisms** | ✅ Excellent | Error boundaries, validation, fallbacks all covered |
| **Disclaimers** | ✅ Present | Section 2.5 has the critical disclaimer |
| **Architecture** | ✅ Sound | Clean component hierarchy, proper data flow |
| **Clinician-Centric** | ✅ Good | Quick scan layout, action-oriented design |
| **Timeline** | ✅ Realistic | 5 weeks with clear deliverables |

### ⚠️ Required Corrections (Before Implementation)

1. **Terminology Fix**: Line 633 uses `efficacy_score` - must transform to `alignment_score` on frontend
2. **Terminology Fix**: Line 335 says "Efficacy prediction" - must say "Mechanism alignment"
3. **Add Scope Banner**: Add prominent scope clarification at top of every dossier

### 📋 Recommended Additions

1. **Glossary Tooltips**: Explain "Mechanism Alignment Score" on hover
2. **Guideline Comparison**: Show how recommendations align with NCCN guidelines
3. **Value Proposition Banner**: Explain why mechanism alignment is valuable (rare cases like MBD4)

### ✅ Connection to Benchmark Work

This frontend displays the **same S/P/E framework** that the benchmark agent validates:
- **MBD4 Agent**: Proved mechanism alignment for 1 patient (depth)
- **Benchmark Agent**: Validates mechanism accuracy at scale (breadth)
- **Frontend**: Presents mechanism alignment to clinicians

---

## Executive Summary

This document outlines the frontend development approach for presenting the MBD4+TP53 clinical dossier (`MBD4_TP53_CLINICAL_DOSSIER.md`) in a beautiful, scalable, and safe interface for clinicians. The goal is to transform the existing markdown dossier into an interactive, user-friendly web application that enables clinicians to:

1. **Run comprehensive genomic analyses** on patient variants
2. **View mechanism-based drug recommendations** with clear confidence indicators
3. **Access clinical decision support** with proper disclaimers and limitations
4. **Export and share** analysis results for tumor board discussions

**Key Principles**:
- **Honest Framing**: Mechanism alignment scores, NOT outcome predictions
- **Safety First**: Built-in failure mechanisms, validation, and error handling
- **Clinician-Centric**: Designed for busy oncologists who need quick, actionable insights
- **Scalable**: Architecture supports multiple concurrent analyses

---

## 1. Architecture Overview

### 1.1 Component Hierarchy

```
ClinicalDossierView (Main Container)
├── DossierHeader (Patient info, analysis metadata)
├── ExecutiveSummaryCard (Key findings at a glance)
├── VariantImpactSection
│   ├── VariantCard (MBD4)
│   └── VariantCard (TP53)
├── PathwayDisruptionSection
│   ├── PathwayVisualization (DDR, TP53 pathways)
│   └── DNARepairCapacityGauge
├── TherapeuticRecommendationsSection
│   ├── DrugRankingCard (Top 5 drugs)
│   ├── DrugDetailModal (Expandable drug info)
│   └── EvidenceTierBadges
├── ClinicalTrialMatchingSection
│   └── TrialListCard (With mechanism fit scores)
├── ResistanceSurveillanceSection
│   └── ResistanceRiskCard (DNA repair capacity trends)
├── ImmunotherapyEligibilitySection
│   ├── TMBIndicator
│   ├── MSIStatusBadge
│   └── IORecommendationCard
├── ClinicalActionPlanSection
│   └── ActionPlanTimeline (Priority 1, 2, 3)
├── EvidenceQualitySection
│   └── ConfidenceBreakdownChart
└── ExportActionsBar (PDF, JSON, Share)
```

### 1.2 Data Flow

```
User Input (Variants + Disease)
    ↓
API Call: POST /api/efficacy/predict
    ↓
Backend Processing (S/P/E Framework)
    ↓
Response Transformation (Dossier Format)
    ↓
Frontend State Management
    ↓
Component Rendering (Cards, Charts, Tables)
    ↓
User Interaction (Expand, Export, Share)
```

### 1.3 Technology Stack

- **Framework**: React (existing codebase)
- **UI Library**: Material-UI (MUI) v5+ (existing)
- **State Management**: React Hooks (useState, useEffect, useContext)
- **Data Fetching**: Custom hooks (`useEfficacy`, `useToxicity`, etc.)
- **Charts**: Recharts or Chart.js (for pathway visualizations)
- **PDF Export**: jsPDF + html2canvas (for clinical reports)
- **Error Handling**: React Error Boundaries + Toast notifications

---

## 2. Component Specifications

### 2.1 ClinicalDossierView (Main Container)

**Purpose**: Orchestrates the entire dossier display

**Props**:
```typescript
interface ClinicalDossierViewProps {
  mutations: Mutation[];
  disease: string;
  tumorContext?: TumorContext;
  dossierData?: DossierData; // Pre-computed or fetched
  onExport?: (format: 'pdf' | 'json' | 'markdown') => void;
}
```

**Features**:
- Loading states with progress indicators
- Error boundaries for graceful failure
- Scroll-to-section navigation
- Print-friendly layout toggle
- Responsive design (mobile, tablet, desktop)

**Safety Mechanisms**:
- Validates mutations array before API call
- Timeout handling (30s default, configurable)
- Retry logic with exponential backoff
- Fallback to cached results if API fails

---

### 2.2 ExecutiveSummaryCard

**Purpose**: Quick overview of key findings (Section 0 from dossier)

**Design**:
- **Hero Section**: Large, prominent display
- **Key Metrics Grid**: 
  - DDR Pathway Burden: 1.00/1.00 (with visual gauge)
  - Top Drug: Olaparib (80% alignment)
  - TMB: 25.0 (HIGH badge)
  - Actionability: HIGH (color-coded badge)
- **Quick Actions**: Jump to recommendations, export, share

**Visual Elements**:
- Color-coded severity indicators (red/yellow/green)
- Progress bars for pathway scores
- Icon badges for clinical significance

**Safety**:
- Validates all metrics exist before rendering
- Shows "N/A" for missing data (never crashes)
- Tooltips explain each metric

---

### 2.3 VariantImpactSection

**Purpose**: Display variant-level analysis (Section 1 from dossier)

**Components**:
- `VariantCard`: Reusable card for each variant
  - Variant name and HGVS notation
  - Classification badge (Pathogenic, VUS, etc.)
  - Inheritance indicator (Germline/Somatic)
  - Functional impact scores (Functionality, Essentiality, Regulatory)
  - Biological rationale (expandable)
  - Affects Drug Response count

**Design**:
- Side-by-side cards for multiple variants
- Expandable sections for detailed rationale
- Visual indicators for driver status (HIGH/MODERATE/LOW)
- Hotspot badges (for TP53 R175H)

**Safety**:
- Handles missing variant data gracefully
- Validates HGVS notation format
- Shows "Unknown" for unclassified variants

---

### 2.4 PathwayDisruptionSection

**Purpose**: Visualize pathway disruption scores (Section 3 from dossier)

**Components**:
- `PathwayVisualization`: Interactive pathway diagram
  - DDR Pathway: 1.00/1.00 (MAXIMUM) - Red indicator
  - TP53 Pathway: 0.80/1.00 (HIGH) - Orange indicator
  - Other pathways (if present) - Gray indicators
- `DNARepairCapacityGauge`: Circular gauge showing 0.60/1.00
  - Color zones: Green (0.8-1.0), Yellow (0.5-0.8), Red (<0.5)
  - Formula breakdown (expandable)
  - Clinical interpretation tooltip

**Design**:
- Sankey diagram showing variant → pathway → drug connections
- Interactive tooltips on hover
- Export pathway diagram as PNG

**Safety**:
- Validates pathway scores are 0-1 range
- Handles missing pathway data (shows "Not assessed")
- Prevents division by zero in calculations

---

### 2.5 TherapeuticRecommendationsSection

**Purpose**: Display drug rankings with mechanism alignment scores (Section 4 from dossier)

**Components**:
- `DrugRankingCard`: List of top drugs
  - Drug name, class, mechanism
  - Alignment score (80% = 0.80) - **RENAMED from "efficacy_score"**
  - Confidence level (40% = MODERATE)
  - Evidence tier badge (CONSIDER/SUPPORTED/INSUFFICIENT)
  - Clinical badges (ClinVar-Moderate, PathwayAligned)
  - Expandable rationale breakdown
- `DrugDetailModal`: Full drug information
  - FDA approval status
  - Off-label rationale
  - Clinical action recommendations
  - S/P/E breakdown (Sequence, Pathway, Evidence)

**Design**:
- Tiered display (Tier 1, Tier 2, Tier 3)
- Color-coded by evidence tier
- Sortable/filterable table view option
- Comparison view (side-by-side drugs)

**Safety**:
- **CRITICAL**: Displays "Mechanism Alignment Score" NOT "Efficacy Score"
- Shows disclaimer: "Scores reflect biological plausibility, NOT predicted response rates"
- Validates drug data structure before rendering
- Handles missing rationale gracefully

**Disclaimers** (Required per Production Readiness):
```jsx
<Alert severity="info" sx={{ mt: 2 }}>
  <strong>Mechanism Alignment, Not Outcome Prediction</strong>
  <br />
  These scores reflect how well each drug targets the disrupted pathways in this tumor.
  They do NOT predict response rates or survival outcomes.
  <br />
  <strong>Drug ranking accuracy:</strong> 100% Top-5 (validated)
  <br />
  <strong>Outcome prediction:</strong> NOT VALIDATED (r=0.037 with PFS)
</Alert>
```

---

### 2.6 ClinicalTrialMatchingSection

**Purpose**: Display matched clinical trials (Section 5 from dossier)

**Components**:
- `TrialListCard`: List of trials
  - NCT ID, Phase, Title
  - Eligibility score (0-1)
  - Mechanism fit score (0-1)
  - Combined score (0.7 × eligibility + 0.3 × mechanism)
  - Match reasoning (expandable)
- `TrialDetailModal`: Full trial information
  - Inclusion/exclusion criteria
  - Primary endpoints
  - Contact information

**Design**:
- Sortable by score (highest first)
- Filterable by phase, status, location
- Empty state: "No trials matched" with helpful suggestions

**Safety**:
- Handles empty trial list gracefully
- Validates trial data structure
- Shows "Search in progress" during API call

---

### 2.7 ResistanceSurveillanceSection

**Purpose**: Display resistance risk assessment (Section 6 from dossier)

**Components**:
- `ResistanceRiskCard`:
  - Current resistance signals (NONE DETECTED badge)
  - DNA repair capacity gauge (0.60/1.00)
  - Risk level indicator (MODERATE - color-coded)
  - Risk factors list
  - Monitoring strategy (expandable)
  - Alert thresholds (expandable)

**Design**:
- Visual risk meter (Low/Moderate/High)
- Timeline view for capacity trends (if historical data available)
- Alert badges for threshold breaches

**Safety**:
- Validates capacity score is 0-1 range
- Handles missing resistance data

---

### 2.8 ImmunotherapyEligibilitySection

**Purpose**: Display IO eligibility assessment (Section 7 from dossier)

**Components**:
- `TMBIndicator`: Large, prominent TMB display
  - Score: 25.0 mutations/Mb
  - Status: HIGH (≥20) - Green badge
  - FDA threshold comparison (2.5× above threshold)
- `MSIStatusBadge`: MSS/MSI-H indicator
- `IORecommendationCard`:
  - Overall eligibility: YES ✅
  - Rationale (expandable)
  - Therapeutic options list
  - Clinical recommendation

**Design**:
- Prominent TMB display (large number, color-coded)
- Comparison chart (patient TMB vs. FDA threshold)
- IO eligibility checklist

**Safety**:
- Validates TMB is numeric
- Handles missing TMB/MSI data

---

### 2.9 ClinicalActionPlanSection

**Purpose**: Actionable treatment recommendations (Section 9 from dossier)

**Components**:
- `ActionPlanTimeline`: Priority-based timeline
  - Priority 1 (Immediate Actions): Olaparib, Carboplatin, Combination
  - Priority 2 (Secondary Considerations): Pembrolizumab, Clinical Trials
  - Priority 3 (Monitoring Strategy): Response assessment, Resistance monitoring
- `ActionItemCard`: Individual action item
  - Therapy name
  - Rationale
  - Mechanism alignment assessment (with disclaimer)
  - Consideration notes

**Design**:
- Timeline visualization (vertical or horizontal)
- Checkbox list for tracking actions
- Export action plan as PDF

**Safety**:
- Validates priority levels (1, 2, 3)
- Handles missing action items

---

### 2.10 EvidenceQualitySection

**Purpose**: Display confidence breakdown and limitations (Section 10 from dossier)

**Components**:
- `ConfidenceBreakdownChart`: Bar chart
  - High Confidence (≥70%): Pathway disruption, Variant classification, TMB
  - Moderate Confidence (40-69%): Drug predictions, DNA repair capacity
  - Low Confidence (<40%): Clinical trials, Nutritional therapy
- `EvidenceStrengthsList`: ✅ Strengths
- `EvidenceLimitationsList`: ⚠️ Limitations

**Design**:
- Color-coded confidence levels
- Expandable sections for details
- Link to full methodology

**Safety**:
- Validates confidence scores are 0-1 range
- Handles missing evidence data

---

## 3. Safety & Failure Mechanisms

### 3.1 Input Validation

**Before API Call**:
```typescript
function validateInputs(mutations: Mutation[], disease: string): ValidationResult {
  const errors: string[] = [];
  
  if (!mutations || mutations.length === 0) {
    errors.push("At least one mutation is required");
  }
  
  mutations.forEach((mut, idx) => {
    if (!mut.gene) errors.push(`Mutation ${idx + 1}: Missing gene symbol`);
    if (!mut.hgvs_p && !mut.chrom || !mut.pos) {
      errors.push(`Mutation ${idx + 1}: Missing variant coordinates`);
    }
  });
  
  if (!disease) errors.push("Disease type is required");
  
  return {
    valid: errors.length === 0,
    errors
  };
}
```

**After API Response**:
```typescript
function validateDossierData(data: any): DossierData | null {
  try {
    // Validate required fields
    if (!data.executive_summary) throw new Error("Missing executive summary");
    if (!data.variants || data.variants.length === 0) throw new Error("Missing variants");
    if (!data.drugs || data.drugs.length === 0) throw new Error("Missing drug recommendations");
    
    // Validate numeric ranges
    data.pathway_scores?.forEach((score: number) => {
      if (score < 0 || score > 1) throw new Error("Pathway score out of range");
    });
    
    return data as DossierData;
  } catch (error) {
    console.error("Dossier data validation failed:", error);
    return null;
  }
}
```

---

### 3.2 Error Handling

**API Error Handling**:
```typescript
async function fetchDossier(mutations: Mutation[], disease: string): Promise<DossierData> {
  try {
    const response = await fetch('/api/efficacy/predict', {
      method: 'POST',
      body: JSON.stringify({ mutations, disease }),
      timeout: 30000 // 30s timeout
    });
    
    if (!response.ok) {
      throw new Error(`API error: ${response.status} ${response.statusText}`);
    }
    
    const data = await response.json();
    const validated = validateDossierData(data);
    
    if (!validated) {
      throw new Error("Invalid dossier data structure");
    }
    
    return validated;
  } catch (error) {
    // Log error for debugging
    console.error("Dossier fetch failed:", error);
    
    // Show user-friendly error
    toast.error("Analysis failed. Please check your inputs and try again.");
    
    // Return fallback data or null
    return null;
  }
}
```

**Component Error Boundaries**:
```typescript
class DossierErrorBoundary extends React.Component {
  state = { hasError: false, error: null };
  
  static getDerivedStateFromError(error: Error) {
    return { hasError: true, error };
  }
  
  componentDidCatch(error: Error, errorInfo: React.ErrorInfo) {
    // Log to error tracking service
    console.error("Dossier component error:", error, errorInfo);
  }
  
  render() {
    if (this.state.hasError) {
      return (
        <Alert severity="error">
          <Typography variant="h6">Something went wrong</Typography>
          <Typography variant="body2">
            The analysis could not be displayed. Please try refreshing the page.
          </Typography>
          <Button onClick={() => window.location.reload()}>
            Refresh Page
          </Button>
        </Alert>
      );
    }
    
    return this.props.children;
  }
}
```

---

### 3.3 Loading States

**Progressive Loading**:
- Show skeleton loaders for each section
- Display sections as they load (not all-or-nothing)
- Show progress indicator: "Loading pathway analysis... (2/11 sections)"

**Timeout Handling**:
- 30s default timeout (configurable)
- Show "Analysis taking longer than expected" message after 15s
- Offer "Cancel" option
- Retry with exponential backoff

---

### 3.4 Data Fallbacks

**Missing Data Handling**:
- Show "Not assessed" for missing pathway scores
- Show "N/A" for missing drug recommendations
- Show "No trials matched" for empty trial list
- Never show undefined/null values (always have fallback)

**Cached Results**:
- Cache dossier results in localStorage (keyed by mutation hash)
- Show "Using cached results" badge if API fails
- Allow "Force refresh" to bypass cache

---

## 4. User Experience Design

### 4.1 Clinician-Centric Design

**Quick Scan Layout**:
- Executive summary at top (most important info first)
- Color-coded severity indicators (red/yellow/green)
- Expandable sections (don't overwhelm with details)
- "Jump to" navigation for long dossiers

**Action-Oriented**:
- Clear "Primary Recommendation" section
- Action items with checkboxes
- Export buttons prominently placed
- Share functionality for tumor board

**Information Hierarchy**:
1. **Critical**: Executive summary, top drug recommendations
2. **Important**: Pathway disruption, resistance risk
3. **Supporting**: Evidence quality, technical details

---

### 4.2 Visual Design Principles

**Color Coding**:
- **Red**: Critical/high risk (DDR 1.00, resistance risk)
- **Orange**: Moderate risk (TP53 0.80, moderate confidence)
- **Green**: Low risk/safe (TMB-High eligibility, low resistance)
- **Gray**: Neutral/informational (supporting data)

**Typography**:
- **Headings**: Bold, larger font (h4, h5)
- **Body**: Readable size (14px minimum)
- **Metrics**: Large, prominent numbers (24px+)
- **Labels**: Smaller, muted color (12px)

**Spacing**:
- Generous padding (16px minimum)
- Clear section separation (24px between sections)
- Card padding (16px internal)

---

### 4.3 Responsive Design

**Breakpoints**:
- **Mobile** (< 600px): Single column, stacked cards
- **Tablet** (600-960px): Two columns, side-by-side variants
- **Desktop** (> 960px): Full layout, all sections visible

**Mobile Optimizations**:
- Collapsible sections (accordion style)
- Bottom navigation for quick access
- Swipe gestures for navigation
- Touch-friendly button sizes (44px minimum)

---

## 5. API Integration

### 5.1 Backend Endpoint

**Endpoint**: `POST /api/efficacy/predict` (existing)

**Request Format**:
```json
{
  "mutations": [
    {
      "gene": "MBD4",
      "hgvs_p": "p.Ile413Serfs*2",
      "chrom": "3",
      "pos": 129430456,
      "ref": "A",
      "alt": "",
      "build": "GRCh37"
    },
    {
      "gene": "TP53",
      "hgvs_p": "p.Arg175His",
      "chrom": "17",
      "pos": 7577120,
      "ref": "G",
      "alt": "A",
      "build": "GRCh37"
    }
  ],
  "disease": "ovarian_cancer",
  "germline_status": "positive",
  "tumor_context": {
    "disease": "ovarian_cancer",
    "tmb": 25.0,
    "msi_status": "MSS"
  }
}
```

**Response Transformation**:
```typescript
function transformApiResponse(apiResponse: EfficacyResponse): DossierData {
  return {
    executive_summary: {
      ddr_pathway_burden: apiResponse.provenance?.confidence_breakdown?.pathway_disruption?.ddr || 0,
      top_drug: apiResponse.efficacy?.drugs?.[0]?.name || "N/A",
      top_drug_alignment: apiResponse.efficacy?.drugs?.[0]?.alignment_score || 0,
      tmb: apiResponse.tumor_context?.tmb || 0,
      actionability: calculateActionability(apiResponse)
    },
    variants: transformVariants(apiResponse),
    pathway_disruption: transformPathways(apiResponse),
    drugs: transformDrugs(apiResponse),
    // ... etc
  };
}
```

---

### 5.2 Custom Hooks

**useClinicalDossier Hook**:
```typescript
function useClinicalDossier(mutations: Mutation[], disease: string) {
  const [dossier, setDossier] = useState<DossierData | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  
  const fetchDossier = async () => {
    setLoading(true);
    setError(null);
    
    try {
      const validated = validateInputs(mutations, disease);
      if (!validated.valid) {
        throw new Error(validated.errors.join(", "));
      }
      
      const response = await fetch('/api/efficacy/predict', {
        method: 'POST',
        body: JSON.stringify({ mutations, disease }),
        timeout: 30000
      });
      
      if (!response.ok) throw new Error(`API error: ${response.status}`);
      
      const data = await response.json();
      const transformed = transformApiResponse(data);
      const validatedData = validateDossierData(transformed);
      
      if (!validatedData) throw new Error("Invalid dossier data");
      
      setDossier(validatedData);
    } catch (err) {
      setError(err.message);
      console.error("Dossier fetch failed:", err);
    } finally {
      setLoading(false);
    }
  };
  
  useEffect(() => {
    if (mutations.length > 0 && disease) {
      fetchDossier();
    }
  }, [mutations, disease]);
  
  return { dossier, loading, error, refetch: fetchDossier };
}
```

---

## 6. Export & Sharing

### 6.1 PDF Export

**Library**: jsPDF + html2canvas

**Features**:
- Full dossier as PDF
- Executive summary only (quick share)
- Customizable sections (user selects)
- Watermark: "Research Use Only"
- Footer: Report ID, analysis date, disclaimers

**Implementation**:
```typescript
async function exportToPDF(dossier: DossierData, sections: string[] = ['all']) {
  const pdf = new jsPDF();
  
  // Add header
  pdf.setFontSize(18);
  pdf.text("Clinical Genomic Analysis Dossier", 20, 20);
  
  // Add disclaimer
  pdf.setFontSize(10);
  pdf.text("Research Use Only - Mechanism Alignment, Not Outcome Prediction", 20, 30);
  
  // Add sections
  if (sections.includes('all') || sections.includes('executive')) {
    addExecutiveSummary(pdf, dossier.executive_summary);
  }
  
  // ... add other sections
  
  // Add footer
  pdf.setFontSize(8);
  pdf.text(`Report ID: ${dossier.report_id}`, 20, pdf.internal.pageSize.height - 10);
  pdf.text(`Analysis Date: ${dossier.analysis_date}`, 20, pdf.internal.pageSize.height - 5);
  
  pdf.save(`clinical_dossier_${dossier.report_id}.pdf`);
}
```

---

### 6.2 JSON Export

**Format**: Full dossier data as JSON

**Use Cases**:
- Integration with other systems
- Programmatic access
- Data backup

---

### 6.3 Share Functionality

**Options**:
- Copy link (with dossier ID)
- Email (pre-filled template)
- Tumor board integration (if available)

---

## 7. Sprint Planning

### Sprint Overview

**Sprint Duration**: 1 week per sprint  
**Total Sprints**: 5 sprints  
**Sprint Goal Pattern**: Each sprint delivers a working, testable increment

---

### Sprint 1: Foundation & Core Infrastructure (Week 1)

**Sprint Goal**: Establish core architecture, API integration, and basic dossier display

**User Story**: "As a clinician, I want to see a basic dossier view with executive summary so that I can quickly understand key findings for a patient."

**Backlog Items**:
- [ ] **SP1-1**: Set up project structure and component scaffolding
  - Create `ClinicalDossierView` container component
  - Set up routing and navigation
  - Initialize TypeScript interfaces for dossier data structures
- [ ] **SP1-2**: Implement `useClinicalDossier` custom hook
  - API integration with `/api/efficacy/predict`
  - Input validation (`validateInputs`)
  - Response transformation (`transformApiResponse`)
  - Error handling with retry logic
- [ ] **SP1-3**: Create `DossierHeader` component
  - Patient info display
  - Analysis metadata (date, report ID)
  - Scope banner: "Mechanism Alignment, Not Outcome Prediction"
- [ ] **SP1-4**: Build `ExecutiveSummaryCard` component
  - Key metrics grid (DDR burden, top drug, TMB, actionability)
  - Visual gauges and progress bars
  - Quick action buttons
- [ ] **SP1-5**: Implement error boundaries and loading states
  - `DossierErrorBoundary` component
  - Skeleton loaders for each section
  - Progressive loading indicators
  - Timeout handling (30s default)

**Definition of Done**:
- ✅ All components render without errors
- ✅ API integration works with MBD4+TP53 test case
- ✅ Error boundaries catch and display errors gracefully
- ✅ Loading states show during API calls
- ✅ Scope banner visible at top of dossier
- ✅ Unit tests for validation functions (80%+ coverage)

**Acceptance Criteria**:
- Can display basic dossier for MBD4+TP53 case
- All safety mechanisms in place (validation, error boundaries)
- Error handling works (network errors, invalid data, timeouts)
- Terminology correct: "alignment_score" not "efficacy_score"

---

### Sprint 2: Variant Analysis & Drug Recommendations (Week 2)

**Sprint Goal**: Display variant-level analysis and therapeutic recommendations with proper disclaimers

**User Story**: "As a clinician, I want to see variant impact and drug recommendations with alignment scores so that I can understand the biological rationale for treatment options."

**Backlog Items**:
- [ ] **SP2-1**: Build `VariantImpactSection` component
  - `VariantCard` reusable component
  - Display HGVS notation, classification, inheritance
  - Functional impact scores (Functionality, Essentiality, Regulatory)
  - Expandable biological rationale
  - Hotspot badges (for TP53 R175H)
- [ ] **SP2-2**: Create `TherapeuticRecommendationsSection` component
  - `DrugRankingCard` component with tiered display
  - **CRITICAL**: Display "Mechanism Alignment Score" (not "Efficacy Score")
  - Evidence tier badges (CONSIDER/SUPPORTED/INSUFFICIENT)
  - Clinical badges (ClinVar-Moderate, PathwayAligned)
  - Expandable rationale breakdown
- [x] **SP2-3**: Implement `DrugDetailModal` component ✅ **COMPLETE**
  - Full drug information modal
  - FDA approval status
  - Off-label rationale
  - S/P/E breakdown (Sequence, Pathway, Evidence)
- [x] **SP2-4**: Add required disclaimers ✅ **COMPLETE**
  - Mechanism alignment disclaimer in TherapeuticRecommendationsSection ✅
  - Tooltips explaining "Mechanism Alignment Score" ✅
  - Glossary tooltips for key terms ✅
  - Value Proposition Banner added to ClinicalDossierView ✅

**Definition of Done**:
- ✅ Variant cards display all variant information correctly
- ✅ Drug recommendations show alignment scores (not efficacy scores)
- ✅ Disclaimers are prominent and clear
- ✅ Expandable sections work (expand/collapse)
- ✅ Modal interactions work smoothly
- ✅ All terminology fixes applied (efficacy_score → alignment_score)

**Acceptance Criteria**:
- Variant impact section displays MBD4 and TP53 variants correctly
- Drug recommendations show top 5 drugs with alignment scores
- Disclaimers are visible and explain mechanism alignment vs. outcome prediction
- Tooltips explain "Mechanism Alignment Score" on hover
- No terminology errors (checked via grep for "efficacy_score")

---

### Sprint 3: Pathway Visualization & Responsive Design (Week 3)

**Sprint Goal**: Add pathway visualizations, charts, and ensure mobile responsiveness

**User Story**: "As a clinician, I want to see pathway disruption visualizations and use the interface on mobile devices so that I can understand pathway context and access information anywhere."

**Backlog Items**:
- [x] **SP3-1**: Build `PathwayDisruptionSection` component ✅ **COMPLETE**
  - `PathwayVisualization` with interactive pathway diagram ✅
  - Pathway scores display (DDR: 1.00/1.00, TP53: 0.80/1.00) ✅
  - Color-coded indicators (Red/Orange/Gray) ✅
  - Interactive tooltips on hover ✅
- [x] **SP3-2**: Create `DNARepairCapacityGauge` component ✅ **COMPLETE**
  - Circular gauge showing 0.60/1.00 ✅
  - Color zones (Green/Yellow/Red) ✅
  - Formula breakdown (expandable) ✅
  - Clinical interpretation tooltip ✅
- [x] **SP3-3**: Implement `ClinicalTrialMatchingSection` component ✅ **COMPLETE**
  - `TrialListCard` with sortable/filterable trials ✅
  - Eligibility and mechanism fit scores ✅
  - Combined score calculation (0.7 × eligibility + 0.3 × mechanism) ✅
  - Empty state handling ✅
- [x] **SP3-4**: Responsive design implementation ✅ **COMPLETE**
  - Mobile breakpoints (< 600px): Single column, stacked cards ✅
  - Tablet breakpoints (600-960px): Two columns ✅
  - Desktop (> 960px): Full layout ✅
  - Collapsible sections (accordion style) for mobile ✅
  - Touch-friendly button sizes (44px minimum) ✅

**Definition of Done**:
- ✅ Pathway visualizations render correctly
- ✅ Charts display accurate data (Recharts or Chart.js)
- ✅ Responsive design works on all breakpoints
- ✅ Mobile optimizations implemented (accordion, touch-friendly)
- ✅ All sections accessible on mobile devices

**Acceptance Criteria**:
- Pathway disruption section displays DDR and TP53 pathways correctly
- DNA repair capacity gauge shows correct value with color zones
- Clinical trial matching section works (even if empty)
- Interface is fully functional on mobile devices
- All interactive elements work on touch devices

---

### Sprint 4: Advanced Features & Export (Week 4)

**Sprint Goal**: Complete all remaining dossier sections and implement export functionality

**User Story**: "As a clinician, I want to see complete dossier information including resistance risk, IO eligibility, action plans, and export reports so that I have full clinical context and can share findings."

**Backlog Items**:
- [ ] **SP4-1**: Build `ResistanceSurveillanceSection` component
  - `ResistanceRiskCard` with visual risk meter
  - DNA repair capacity gauge integration
  - Risk level indicators (Low/Moderate/High)
  - Monitoring strategy (expandable)
- [ ] **SP4-2**: Create `ImmunotherapyEligibilitySection` component
  - `TMBIndicator` with prominent TMB display
  - `MSIStatusBadge` (MSS/MSI-H)
  - `IORecommendationCard` with eligibility assessment
  - Comparison chart (patient TMB vs. FDA threshold)
- [ ] **SP4-3**: Implement `ClinicalActionPlanSection` component
  - `ActionPlanTimeline` with priority-based timeline
  - `ActionItemCard` for individual actions
  - Checkbox list for tracking actions
  - Priority visualization (Priority 1, 2, 3)
- [ ] **SP4-4**: Build `EvidenceQualitySection` component
  - `ConfidenceBreakdownChart` (bar chart)
  - `EvidenceStrengthsList` and `EvidenceLimitationsList`
  - Color-coded confidence levels
- [ ] **SP4-5**: Implement export functionality
  - PDF export (jsPDF + html2canvas)
  - JSON export
  - Share functionality (copy link, email template)
  - Watermarks and disclaimers in PDF

**Definition of Done**:
- ✅ All dossier sections implemented and displayed
- ✅ Export functionality works (PDF and JSON)
- ✅ All disclaimers present in exported PDFs
- ✅ Share functionality works
- ✅ Complete dossier displays for MBD4+TP53 case

**Acceptance Criteria**:
- All 10 dossier sections render correctly
- PDF export includes all sections with watermarks
- JSON export contains full dossier data
- Share functionality works (copy link, email)
- All disclaimers visible in exports

---

### Sprint 5: Polish, Testing & Accessibility (Week 5)

**Sprint Goal**: Polish user experience, complete testing, and ensure accessibility compliance

**User Story**: "As a clinician, I want a polished, accessible interface that works reliably so that I can use it confidently in clinical practice."

**Backlog Items**:
- [ ] **SP5-1**: User testing with clinicians
  - Recruit 5-10 oncologists for testing
  - Collect feedback on usability, clarity, missing features
  - Iterate based on feedback
- [ ] **SP5-2**: Performance optimization
  - Bundle size optimization
  - Code splitting for lazy loading
  - Memoization of expensive computations
  - Image optimization
  - Target: < 3s load time for full dossier
- [ ] **SP5-3**: Accessibility improvements (WCAG 2.1 AA)
  - ARIA labels for all interactive elements
  - Keyboard navigation support
  - Screen reader compatibility
  - Color contrast compliance
  - Focus indicators
- [ ] **SP5-4**: Comprehensive testing
  - Unit tests: 80%+ coverage target
  - Integration tests for API calls
  - Error scenario testing
  - Cross-browser testing
  - Mobile device testing
- [ ] **SP5-5**: Documentation and bug fixes
  - User guide documentation
  - API documentation
  - Component documentation
  - Bug fixes from testing
  - Final code review

**Definition of Done**:
- ✅ User testing completed with positive feedback
- ✅ Performance targets met (< 3s load time)
- ✅ WCAG 2.1 AA compliance verified
- ✅ Test coverage ≥ 80%
- ✅ All critical bugs fixed
- ✅ Documentation complete

**Acceptance Criteria**:
- Clinicians can use interface effectively (positive feedback)
- Load time < 3s for full dossier
- WCAG 2.1 AA compliant (verified with accessibility tools)
- Test coverage ≥ 80%
- No critical bugs (P0/P1 issues)
- Documentation available for users and developers

---

### Sprint Retrospective & Planning

**After Each Sprint**:
- Sprint review (demo working features)
- Sprint retrospective (what went well, what to improve)
- Backlog refinement for next sprint
- Update sprint goals based on learnings

**Sprint Ceremonies**:
- **Daily Standups**: 15 minutes (progress, blockers, plans)
- **Sprint Planning**: 2 hours (backlog refinement, sprint goal)
- **Sprint Review**: 1 hour (demo, stakeholder feedback)
- **Sprint Retrospective**: 1 hour (improvements, action items)

---

## 8. Testing Strategy

### 8.1 Unit Tests

**Components**:
- Test each card component with mock data
- Test error boundaries
- Test validation functions
- Test data transformation functions

**Coverage Target**: 80%+

---

### 8.2 Integration Tests

**API Integration**:
- Test successful API calls
- Test API error handling
- Test timeout handling
- Test data validation

---

### 8.3 User Acceptance Testing

**Clinician Testing**:
- 5-10 oncologists test the interface
- Feedback on:
  - Usability
  - Information clarity
  - Missing features
  - Confusing elements

---

## 9. Success Metrics

### 9.1 Technical Metrics

- **Load Time**: < 3s for full dossier
- **Error Rate**: < 1% of analyses fail
- **API Success Rate**: > 99%
- **Mobile Usability**: 100% of features accessible on mobile

---

### 9.2 User Metrics

- **Time to Insight**: < 30s to find top recommendation
- **Export Usage**: > 50% of users export dossier
- **Error Recovery**: 100% of errors show helpful message
- **Clinician Satisfaction**: > 4/5 rating

---

## 10. Future Enhancements

### 10.1 Short-Term (Next Quarter)

- Historical comparison (track DNA repair capacity over time)
- Multi-patient comparison
- Customizable dossier sections
- Integration with EMR systems

---

### 10.2 Long-Term (Next Year)

- Real-time updates (if new evidence emerges)
- AI-powered insights (explain pathway connections)
- Collaborative features (annotate, discuss)
- Outcome tracking (if validation studies complete)

---

## 11. Critical Requirements Checklist

### Before Production Launch

- [ ] **Disclaimers**: All outputs show "Mechanism Alignment, Not Outcome Prediction"
- [ ] **Terminology**: "efficacy_score" renamed to "alignment_score" everywhere
- [ ] **Error Handling**: All API failures handled gracefully
- [ ] **Validation**: All inputs validated before API calls
- [ ] **Loading States**: All sections show loading indicators
- [ ] **Responsive**: Works on mobile, tablet, desktop
- [ ] **Accessibility**: WCAG 2.1 AA compliant
- [ ] **Performance**: < 3s load time
- [ ] **Testing**: 80%+ test coverage
- [ ] **Documentation**: User guide and API docs

---

## 12. References

- **Dossier Content**: `.cursor/ayesha/MBD4_TP53_CLINICAL_DOSSIER.md`
- **Production Readiness**: `.cursor/ayesha/PRODUCTION_READINESS_ASSESSMENT.md`
- **Existing Components**: `oncology-frontend/src/components/ClinicalGenomicsCommandCenter/`
- **API Endpoints**: `oncology-backend-minimal/api/routers/efficacy.py`

---

**End of Frontend Development Plan**

*This plan is a living document and will be updated as development progresses and requirements evolve.*


