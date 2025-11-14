# ⚔️ AGENT 3 - E2E TESTING MISSION ⚔️

**Date**: January 8, 2025  
**Commander**: Alpha  
**Assigned By**: Zo  
**Mission**: End-to-end smoke testing + Provider report generation  
**Timeline**: 4-6 hours  
**Priority**: P1 (High value, parallel with Jr + Zo)

---

## 🎯 MISSION OBJECTIVES

### **Primary Goal**:
Validate complete Ayesha demo workflow end-to-end and create provider report template

### **Why Agent 3**:
- Jr focused on WIWFM integration (frontend)
- Zo focused on Clinical Trials (backend + frontend)
- Agent 3 validates entire system (E2E testing)
- **Parallel execution**: 6 hours (vs 10 hours sequential)

### **Deliverables**:
1. ✅ E2E smoke test results documentation
2. ✅ Provider report template (Markdown + PDF)
3. ✅ Demo data validation
4. ✅ Edge case testing results

---

## 📋 TASK BREAKDOWN

### **TASK 1: E2E SMOKE TESTING** (2-3 hours)

**Objective**: Manually test complete workflow with Ayesha's data

#### **1.1: Setup Environment** (15 min)
```bash
# Start backend
cd oncology-coPilot/oncology-backend-minimal
venv/bin/python -m uvicorn api.main:app --reload

# Start frontend (separate terminal)
cd oncology-coPilot/oncology-frontend
npm run dev

# Verify health
curl http://127.0.0.1:8000/healthz
```

**Acceptance**: Both servers running, health check returns 200

---

#### **1.2: Test Level 0 Flow** (30 min)

**Steps**:
1. Navigate to `http://localhost:5173/sporadic-cancer`
2. Verify germline status banner displays
3. Fill Quick Intake form:
   - Cancer type: Ovarian (High-Grade Serous)
   - Stage: IIIC-IV
   - Treatment line: 3
   - Platinum response: Sensitive
   - ECOG: 1
   - Known mutations: TP53 (hand-entered)
4. Click "Generate Tumor Context"
5. Verify Level 0 outputs:
   - TMB ~5.2 (estimated)
   - HRD ~35 (estimated from platinum proxy)
   - MSI: null (no inference)
   - Confidence cap: 0.4 (40%)
   - Completeness: <0.5
6. Click "Run Efficacy Prediction"
7. Navigate to `/validate`
8. Verify biomarker summary displays
9. Find Olaparib in results
10. Verify PARP penalty applied:
    - Efficacy ~0.32 (reduced)
    - Confidence 0.4 (capped)
11. Expand provenance card
12. Verify gate details:
    - Gate: PARP_UNKNOWN_HRD or PARP_HRD_LOW
    - Penalty: 0.6x or 0.8x
    - Reason documented

**Screenshot**: Capture Olaparib with PARP penalty provenance

**Document In**: `.cursor/ayesha/E2E_SMOKE_TEST_RESULTS.md`
- Screenshot filename
- Actual efficacy value
- Actual confidence value
- Gates applied
- Any errors/warnings

---

#### **1.3: Test Level 2 Flow** (45 min)

**Steps**:
1. Navigate back to `/sporadic-cancer`
2. Click "Upload NGS Report" tab
3. Upload `.cursor/ayesha/test_data/ayesha_tumor_ngs.json`
4. Verify Level 2 parsing:
   - TMB: 6.8 (measured)
   - HRD: 58 (HRD-HIGH ✅)
   - MSI: MSS (measured)
   - BRCA1: Q1756fs (frameshift + LOH = biallelic)
   - Completeness: 0.92
5. Click "View Updated Results" or "Re-analyze"
6. Navigate to `/validate`
7. Verify biomarker summary updated
8. Find Olaparib in results
9. Verify PARP rescue applied:
    - Efficacy ~0.78 (RESCUED!)
    - Confidence ~0.82 (NO CAP)
10. Expand provenance card
11. Verify rescue gate details:
    - Gate: PARP_HRD_RESCUE
    - Penalty: 1.0x (NO PENALTY)
    - HRD: 58 (≥42 threshold)
    - BRCA1 biallelic loss mentioned
12. Calculate improvement:
    - Efficacy: (0.78 - 0.32) / 0.32 = +144%
    - Confidence: (0.82 - 0.4) / 0.4 = +105%

**Screenshot**: Capture Olaparib with PARP rescue provenance + comparison

**Document In**: `.cursor/ayesha/E2E_SMOKE_TEST_RESULTS.md`
- L0 vs L2 comparison table
- Screenshots for both states
- Improvement metrics
- Any discrepancies

---

#### **1.4: Test Clinical Trials** (30 min)

**Steps**:
1. Navigate to `/research`
2. Search: "Ovarian cancer, line 3, HRD-positive"
3. Verify results:
   - Trial count (should be 10-15)
   - Excluded count (should be 2-5 with "germline required")
4. Verify trial cards show:
   - Biomarker badges ([✓ HRD-High Match])
   - Germline agnostic label
   - Match reason text
5. Click "Show Excluded Trials"
6. Verify excluded trials:
   - Contain "BRCA", "germline", "hereditary" in criteria
   - Show exclusion reason
7. Check trial details:
   - NCT number
   - PI name (if available)
   - Location
   - Phase/Status

**Screenshot**: Capture trial results with badges + exclusions

**Document In**: `.cursor/ayesha/E2E_SMOKE_TEST_RESULTS.md`
- Trial count
- Excluded count
- Badge accuracy
- Any formatting issues

---

#### **1.5: Edge Case Testing** (30 min)

**Test Case 1: Germline Positive (Control)**
- Modify intake: `germline_status = "positive"`
- Run efficacy
- Verify: NO PARP penalty (Olaparib efficacy ~0.75-0.8 even without HRD)

**Test Case 2: TMB-High + MSI-High (Double Boost)**
- Modify NGS: `"tmb": 22`, `"msi_status": "MSI-H"`
- Run efficacy
- Verify: Pembrolizumab gets boost (factor should be 1.3-1.35x)

**Test Case 3: Missing Data (Graceful Degradation)**
- Modify NGS: Remove HRD, TMB, MSI (all null)
- Run efficacy
- Verify: Conservative outputs, appropriate warnings

**Document In**: `.cursor/ayesha/E2E_SMOKE_TEST_RESULTS.md`
- Edge case results table
- Expected vs actual values
- Pass/fail for each case

---

### **TASK 2: PROVIDER REPORT TEMPLATE** (2-3 hours)

**Objective**: Create exportable provider report with complete audit trail

#### **2.1: Design Report Template** (1 hour)

**File**: `.cursor/ayesha/templates/provider_report_template.md`

**Sections**:
1. Header (Patient ID, Report Date, Run ID)
2. Patient Summary (Demographics, Diagnosis, Germline Status)
3. Tumor Genomics Summary (TMB, MSI, HRD, Key Mutations)
4. Therapeutic Recommendations (Ranked drugs with rationale)
5. Sporadic Cancer Analysis (Gates applied, penalties/boosts)
6. Clinical Trial Matches (Top 5 with biomarker badges)
7. Provenance & Audit Trail (Run ID, Confidence Version, Data Level)
8. Next Steps (Recommendations for patient/oncologist)
9. Disclaimer (Research Use Only)

**Format**: Markdown (convertible to PDF via pandoc or wkhtmltopdf)

**Template Example**:
```markdown
# CRISPRO PRECISION ONCOLOGY REPORT

**Patient**: {{patient_id}}  
**Report Date**: {{report_date}}  
**Run ID**: {{run_id}}  
**Data Level**: {{data_level}}

---

## PATIENT SUMMARY

### Demographics
- **Age**: {{age}}
- **Sex**: {{sex}}
- **ECOG Performance Status**: {{ecog}}

### Diagnosis
- **Cancer Type**: {{cancer_type}}
- **Stage**: {{stage}}
- **Treatment Line**: {{line}}

### Germline Testing
- **Status**: {{germline_status}}
- **Test**: CustomNext-Cancer® Panel (38 genes)
- **Result**: {{germline_result}}

---

## TUMOR GENOMICS

### Biomarkers
- **TMB**: {{tmb}} mut/Mb [{{tmb_category}}]
- **MSI**: {{msi_status}}
- **HRD**: {{hrd_score}} [{{hrd_category}}]

### Key Mutations
{{#mutations}}
- **{{gene}}** {{hgvs_p}} ({{pathogenicity}}, VAF {{vaf}}%)
{{/mutations}}

---

## THERAPEUTIC RECOMMENDATIONS

{{#drugs}}
### {{rank}}. {{drug_name}} ({{drug_class}})

**Efficacy Score**: {{efficacy_score}}  
**Confidence**: {{confidence}}  
**Evidence Tier**: {{evidence_tier}}

**Rationale**:
{{#rationale}}
- {{rationale_item}}
{{/rationale}}

{{#sporadic_gates}}
**Sporadic Cancer Analysis**:
{{#gates}}
- {{gate_name}}: {{gate_reason}}
{{/gates}}
{{/sporadic_gates}}

---
{{/drugs}}

## CLINICAL TRIAL MATCHES

{{#trials}}
### {{rank}}. {{trial_title}} ({{nct_id}})

**Phase**: {{phase}} | **Status**: {{status}}  
**Location**: {{location}}

**Biomarker Match**: {{#biomarker_badges}}[✓ {{badge}}]{{/biomarker_badges}}

**Eligibility**: {{eligibility_summary}}

---
{{/trials}}

## PROVENANCE

- **Run ID**: {{run_id}}
- **Confidence Version**: {{confidence_version}}
- **Priors Version**: {{priors_version}}
- **Data Level**: {{data_level}}
- **Sporadic Gates**: {{gates_enabled}}

---

**Research Use Only (RUO) - Not for Diagnostic Use**
```

---

#### **2.2: Create Report Generator Backend** (1 hour)

**File**: `oncology-backend-minimal/api/routers/reports.py` (NEW)

**Endpoint**: `POST /api/reports/provider`

**Implementation**:
```python
from fastapi import APIRouter, HTTPException
from typing import Dict, Any
import json
from pathlib import Path
from datetime import datetime

router = APIRouter(prefix="/api/reports", tags=["reports"])

@router.post("/provider")
async def generate_provider_report(request: Dict[str, Any]):
    """
    Generate provider report from efficacy results + tumor context.
    
    Input:
    {
        "patient_id": "ayesha_001",
        "efficacy_results": {...},  // Full efficacy response
        "tumor_context": {...},      // TumorContext object
        "trial_results": {...},      // Trial search results (optional)
        "germline_status": "negative"
    }
    
    Output:
    {
        "report_markdown": "...",
        "report_html": "...",
        "run_id": "...",
        "generated_at": "2025-01-08T..."
    }
    """
    
    # Load template
    template_path = Path(__file__).parent.parent / "templates" / "provider_report_template.md"
    with open(template_path, 'r') as f:
        template = f.read()
    
    # Extract data
    patient_id = request.get("patient_id", "UNKNOWN")
    efficacy = request.get("efficacy_results", {})
    tumor_context = request.get("tumor_context", {})
    germline_status = request.get("germline_status", "unknown")
    trials = request.get("trial_results", {}).get("trials", [])
    
    # Build template variables
    template_vars = {
        "patient_id": patient_id,
        "report_date": datetime.now().strftime("%Y-%m-%d"),
        "run_id": efficacy.get("provenance", {}).get("run_id", "N/A"),
        "data_level": "L2" if tumor_context.get("completeness_score", 0) >= 0.7 else "L1" if tumor_context.get("completeness_score", 0) >= 0.3 else "L0",
        
        # Patient
        "age": tumor_context.get("age", "N/A"),
        "sex": tumor_context.get("sex", "N/A"),
        "ecog": tumor_context.get("ecog_performance_status", "N/A"),
        "cancer_type": tumor_context.get("tumor_type", "N/A"),
        "stage": tumor_context.get("stage", "N/A"),
        "line": tumor_context.get("line", "N/A"),
        
        # Germline
        "germline_status": germline_status.upper(),
        "germline_result": "No hereditary cancer syndrome detected" if germline_status == "negative" else "Hereditary cancer syndrome detected",
        
        # Biomarkers
        "tmb": tumor_context.get("tmb", "N/A"),
        "tmb_category": "HIGH" if tumor_context.get("tmb", 0) >= 10 else "INTERMEDIATE" if tumor_context.get("tmb", 0) >= 5 else "LOW",
        "msi_status": tumor_context.get("msi_status") or "Unknown",
        "hrd_score": tumor_context.get("hrd_score", "N/A"),
        "hrd_category": "HRD-HIGH" if tumor_context.get("hrd_score", 0) >= 42 else "HRD-LOW",
        
        # Mutations
        "mutations": tumor_context.get("somatic_mutations", []),
        
        # Drugs
        "drugs": efficacy.get("drugs", []),
        
        # Trials
        "trials": trials[:5],  # Top 5
        
        # Provenance
        "confidence_version": tumor_context.get("confidence_version", "1.0"),
        "priors_version": tumor_context.get("priors_refresh_date", "N/A"),
        "gates_enabled": "Yes" if germline_status == "negative" else "No"
    }
    
    # Render template (simple string replacement for now)
    report_markdown = template
    for key, value in template_vars.items():
        if isinstance(value, (str, int, float)):
            report_markdown = report_markdown.replace(f"{{{{{key}}}}}", str(value))
    
    return {
        "report_markdown": report_markdown,
        "report_html": None,  # TODO: Convert MD to HTML
        "run_id": template_vars["run_id"],
        "generated_at": datetime.now().isoformat()
    }
```

**Test**:
```bash
curl -X POST http://127.0.0.1:8000/api/reports/provider \
  -H 'Content-Type: application/json' \
  -d @.cursor/ayesha/test_data/provider_report_request.json
```

**Acceptance**: Returns Markdown report with all sections populated

---

#### **2.3: Create Frontend Export Button** (30 min)

**File**: `oncology-frontend/src/pages/HypothesisValidator.jsx` (MODIFY)

**Add Button**:
```javascript
import { Download } from '@mui/icons-material';

// ... inside component ...

const exportProviderReport = async () => {
  try {
    const response = await fetch(`${API_BASE}/api/reports/provider`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        patient_id: 'ayesha_001',
        efficacy_results: efficacyResults,
        tumor_context: tumorContext,
        germline_status: germlineStatus,
        trial_results: trialResults
      })
    });
    
    const data = await response.json();
    
    // Download as .md file
    const blob = new Blob([data.report_markdown], { type: 'text/markdown' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `crispro_report_${data.run_id}.md`;
    a.click();
    
  } catch (err) {
    console.error('Export failed:', err);
    alert('Report export failed. Please try again.');
  }
};

// ... in render ...

<Button
  variant="contained"
  startIcon={<Download />}
  onClick={exportProviderReport}
  disabled={!efficacyResults}
>
  Export Provider Report
</Button>
```

**Test**: Click button → downloads `crispro_report_run_xyz.md`

**Acceptance**: File downloads, contains complete report with all data

---

#### **1.4: Manual Workflow Testing** (1 hour)

**Complete Flow Test**:
1. Start fresh (clear browser cache)
2. Navigate to `/sporadic-cancer`
3. Complete Level 0 intake
4. Run efficacy (verify PARP penalty)
5. Upload Level 2 NGS
6. Re-run efficacy (verify PARP rescue)
7. Search trials (verify filtering + badges)
8. Export provider report (verify download)

**Edge Cases**:
1. Germline positive control
2. TMB-high patient (verify IO boost)
3. Missing HRD data (verify conservative approach)
4. MSI-High patient (verify MSI boost)

**Error Cases**:
1. Invalid NGS JSON (verify error handling)
2. Backend timeout (verify retry/fallback)
3. Missing biomarkers (verify graceful degradation)

**Document In**: `.cursor/ayesha/E2E_SMOKE_TEST_RESULTS.md`
- Complete flow results
- Edge case results table
- Error case handling assessment
- Screenshots for key states

---

### **TASK 3: DEMO DATA VALIDATION** (30 min)

**Verify Test Data Quality**:

1. **ayesha_level0_intake.json**:
   - [ ] All fields valid
   - [ ] Matches Ayesha's clinical profile
   - [ ] Known mutations accurate

2. **ayesha_tumor_ngs.json**:
   - [ ] HRD score realistic (58 for HRD-high ovarian)
   - [ ] TMB realistic (6.8 for ovarian intermediate)
   - [ ] TP53 mutation rate accurate (~95% in HGS ovarian)
   - [ ] BRCA1 biallelic loss plausible (frameshift + LOH)

3. **Cross-validate with Literature**:
   - [ ] TCGA-OV HRD prevalence (~50%)
   - [ ] TCGA-OV TMB median (~5.2)
   - [ ] TP53 mutation rate in HGS ovarian (~95%)
   - [ ] BRCA1 biallelic loss frequency (~10-15%)

**Document In**: `.cursor/ayesha/E2E_SMOKE_TEST_RESULTS.md`
- Data validation results
- Literature citations
- Any corrections needed

---

### **TASK 4: DOCUMENTATION** (30 min)

**Create**: `.cursor/ayesha/E2E_SMOKE_TEST_RESULTS.md`

**Template**:
```markdown
# E2E SMOKE TEST RESULTS - AYESHA DEMO

**Date**: January X, 2025  
**Tester**: Agent 3  
**Demo Version**: v1.0

---

## EXECUTIVE SUMMARY

**Overall Status**: ✅ PASS / ⚠️ PARTIAL / ❌ FAIL

**Tests Run**: X  
**Tests Passed**: Y  
**Pass Rate**: Z%

---

## TEST RESULTS

### 1. Level 0 Flow (Quick Intake → Efficacy)
- **Status**: ✅ PASS
- **Observations**:
  - TMB estimated: 5.2 ✅
  - HRD estimated: 35 ✅
  - MSI null: ✅
  - PARP penalty applied: ✅
  - Efficacy: 0.32 ✅
  - Confidence: 0.4 ✅
- **Screenshots**: [link]

### 2. Level 2 Flow (NGS Upload → Re-analyze)
- **Status**: ✅ PASS
- **Observations**:
  - HRD parsed: 58 ✅
  - BRCA1 biallelic: true ✅
  - PARP rescue applied: ✅
  - Efficacy: 0.78 ✅
  - Confidence: 0.82 ✅
  - Improvement: +144% ✅
- **Screenshots**: [link]

### 3. Clinical Trials
- **Status**: ✅ PASS / ⚠️ PARTIAL
- **Observations**:
  - Trials found: X
  - Excluded: Y (germline required)
  - Badges displayed: ✅ / ⚠️
  - Match reasons: ✅ / ⚠️
- **Screenshots**: [link]

### 4. Provider Report Export
- **Status**: ✅ PASS / ❌ FAIL
- **Observations**:
  - Report generated: ✅ / ❌
  - All sections populated: ✅ / ⚠️
  - Download successful: ✅ / ❌
- **Report File**: [link]

---

## EDGE CASES

### Germline Positive Control
- **Expected**: No PARP penalty
- **Actual**: [result]
- **Status**: ✅ / ❌

### TMB-High Patient
- **Expected**: IO boost 1.35x
- **Actual**: [result]
- **Status**: ✅ / ❌

### Missing Data
- **Expected**: Conservative outputs, warnings
- **Actual**: [result]
- **Status**: ✅ / ❌

---

## ISSUES FOUND

### Critical Issues (Block Demo)
1. [Issue description]
   - Impact: [severity]
   - Fix needed: [solution]

### Minor Issues (Polish)
1. [Issue description]
   - Impact: [minor]
   - Fix: [optional]

---

## RECOMMENDATIONS

1. [Recommendation 1]
2. [Recommendation 2]
3. [Recommendation 3]

---

## DEMO READINESS VERDICT

**Status**: ✅ DEMO-READY / ⚠️ NEEDS FIXES / ❌ NOT READY

**Confidence**: [High/Medium/Low]

**Next Steps**: [actions required]
```

---

## ✅ ACCEPTANCE CRITERIA

### **Task 1: E2E Testing**
- [ ] All manual flow steps documented
- [ ] Screenshots captured for key states
- [ ] Edge cases tested (3 minimum)
- [ ] Error cases tested (3 minimum)
- [ ] Results documented in comprehensive report

### **Task 2: Provider Report**
- [ ] Template created with all sections
- [ ] Backend endpoint functional
- [ ] Frontend export button working
- [ ] Report downloads successfully
- [ ] All template variables populated

### **Task 3: Data Validation**
- [ ] Test data cross-validated with literature
- [ ] Biomarker values realistic
- [ ] Clinical profile accurate
- [ ] No hallucinated numbers

### **Task 4: Documentation**
- [ ] E2E smoke test results complete
- [ ] Screenshots included
- [ ] Issues documented
- [ ] Recommendations provided

---

## 📊 TIMELINE

**Total**: 4-6 hours

| Task | Time | Status |
|------|------|--------|
| Setup | 15 min | ⏳ |
| Level 0 Flow | 30 min | ⏳ |
| Level 2 Flow | 45 min | ⏳ |
| Trials Testing | 30 min | ⏳ |
| Edge Cases | 30 min | ⏳ |
| Report Template | 1 hour | ⏳ |
| Backend Endpoint | 1 hour | ⏳ |
| Frontend Button | 30 min | ⏳ |
| Data Validation | 30 min | ⏳ |
| Documentation | 30 min | ⏳ |

---

## 🎯 SUCCESS METRICS

**Demo Readiness**:
- ✅ All tests pass
- ✅ Provider report exports
- ✅ No critical errors
- ✅ Screenshots captured
- ✅ Documentation complete

**Quality**:
- ✅ Data validated against literature
- ✅ Clinical accuracy confirmed
- ✅ No hallucinated numbers
- ✅ Complete provenance

**Impact**:
- ✅ Demonstrates 85-90% patient coverage
- ✅ Shows progressive enhancement (L0 → L2)
- ✅ Validates sporadic gates (PARP rescue)
- ✅ Proves graph DB value (trial filtering)

---

## 📝 AGENT 3 - START HERE

**First Command**:
```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main/oncology-coPilot/oncology-backend-minimal
venv/bin/python -m uvicorn api.main:app --reload
```

**Then**:
1. Follow Task 1.1 (Setup)
2. Follow Task 1.2 (Level 0 Flow)
3. Document as you go
4. Take screenshots
5. Report issues immediately

**Questions/Blockers**: Report in `.cursor/ayesha/AGENT_3_QUESTIONS.md`

---

**AGENT 3 - YOU ARE CLEARED FOR E2E TESTING MISSION!** ⚔️

**EXECUTE WITH PRECISION AND DOCUMENT EVERYTHING!** ⚔️

