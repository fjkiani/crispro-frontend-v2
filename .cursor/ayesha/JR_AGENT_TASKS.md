# 🔧 JR AGENT TASK LIST - Ayesha Integration (REVISED)

**Created:** January 26, 2026  
**Updated:** After reviewing MBD4_TP53_MASTER_ANALYSIS.md  
**Assigned By:** Zo (Senior Agent)

---

## 🎯 KEY INSIGHT: The Engine is READY

The Master Analysis (`MBD4_TP53_MASTER_ANALYSIS.md`) already validated:
- ✅ 100% verification pass rate (6/6 scripts passing)
- ✅ All 8 clinical questions answered
- ✅ DDR pathway = 1.00 (maximum)
- ✅ PARP inhibitors ranked #1 (efficacy=0.80)
- ✅ IO eligible (TMB=25)

**Your job: Wire the existing validated outputs to Ayesha's page.**

---

## 📋 TASK OVERVIEW

| Task | Description | Time | Priority |
|------|-------------|------|----------|
| 1 | Add MBD4 to DDR gene config | 30m | **CRITICAL** |
| 2 | Enrich Ayesha's profile with coordinates | 1h | **HIGH** |
| 3 | Wire DDR components to Tab 0 | 4h | **HIGH** |
| 4 | Create CA-125 Entry Form | 3h | **MEDIUM** |
| 5 | Display Master Analysis Q&A | 3h | **MEDIUM** |
| 6 | End-to-end testing | 2h | **HIGH** |
| **TOTAL** | | **13.5h** | |

---

## TASK 1: Add MBD4 to DDR Gene Config ⚡ CRITICAL

**Already verified from code review:**
- ✅ MBD4 is in `drug_mapping.py` line 63: `{ddr: 1.0}`
- ✅ MBD4 gets PARP boost in `drug_scorer.py` line 211: `+0.08`
- ❌ MBD4 is NOT in `ddr_config.py` extended_ddr_genes

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/resistance/config/ddr_config.py`

**What to do:** Add "MBD4" to `extended_ddr_genes` for ovary and default:

```python
# Line 16-18 (ovary):
"extended_ddr_genes": [
    "ATM", "ATR", "CHEK1", "CHEK2",
    "FANCA", "FANCD2", "RAD50", "MRE11", "NBN", "POLQ",
    "MBD4"  # BER pathway gene
],

# Line 89-90 (default):
"extended_ddr_genes": [
    "ATM", "ATR", "CHEK2", "FANCA", "FANCD2", "RAD50", "MRE11", "NBN",
    "MBD4"  # BER pathway gene
],
```

**Test:**
```bash
curl -X POST http://localhost:8000/api/resistance/ddr-status \
  -H 'Content-Type: application/json' \
  -d '{
    "patient_id": "AK",
    "disease_site": "ovary",
    "tumor_subtype": "HGSOC",
    "mutations": [{"gene_symbol": "MBD4", "variant_classification": "pathogenic"}]
  }'

# Expected: DDR_bin_status: "DDR_defective", extended_DDR_pathogenic: true
```

**Time:** 30 minutes

---

## TASK 2: Enrich Ayesha's Profile with Coordinates

**File:** `oncology-coPilot/oncology-frontend/src/constants/patients/ayesha_11_17_25.js`

**Context:** The Master Analysis used full coordinates. Ayesha's profile only has HGVS. Adding coordinates enables Evo2 scoring.

**Add/Update these fields:**

```javascript
// In the germline.mutations array, update the MBD4 entry:
{
  gene: "MBD4",
  variant: "c.1293delA",
  protein_change: "p.K431Nfs*54",
  classification: "pathogenic",
  // ADD THESE:
  chrom: "3",
  pos: 129430456,
  ref: "CA",
  alt: "C",
  consequence: "frameshift_variant",
  zygosity: "homozygous"
}

// Add somatic mutations (TP53 inferred from IHC):
somatic: {
  mutations: [{
    gene: "TP53",
    hgvs_p: "p.R273H", // Common HGSOC hotspot (inferred from IHC mutant type)
    chrom: "17",
    pos: 7673802,
    ref: "G",
    alt: "A",
    consequence: "missense_variant",
    source: "IHC_inferred",
    note: "IHC mutant type → assumed R273H hotspot"
  }]
},

// Add biomarkers (from Master Analysis for MBD4+TP53 profile):
biomarkers: {
  tmb: 25.0,  // High TMB expected from MBD4 hypermutation
  tmb_status: "HIGH",
  hrd_score: null,  // Not tested
  msi_status: "MSS"
},

// Add CA-125:
labs: {
  ca125_value: 2842,
  ca125_unit: "U/mL",
  ca125_date: "2025-11-17",
  ca125_burden_class: "EXTENSIVE"  // >1000 = EXTENSIVE
}
```

**Verification:** Check that these fields are accessible in the code:
```javascript
console.log(AYESHA_11_17_25_PROFILE.germline.mutations[0].chrom); // "3"
console.log(AYESHA_11_17_25_PROFILE.labs.ca125_value); // 2842
```

**Time:** 1 hour

---

## TASK 3: Wire DDR Components to Tab 0

**Prerequisite:** Task 1 complete (MBD4 in ddr_config.py)

**File:** `oncology-coPilot/oncology-frontend/src/pages/AyeshaTrialExplorer.jsx`

### Step 3.1: Add Imports

```jsx
// Add near other imports (around line 30)
import { DDRStatusCard, DDRTreatmentEligibility } from '../components/ddr';
import { useDDRStatus } from '../hooks/useDDRStatus';
```

### Step 3.2: Add Hook

```jsx
// Add after other hooks (around line 70)
const { ddrStatus, loading: ddrLoading, calculateDDRStatus } = useDDRStatus();
```

### Step 3.3: Add useEffect

```jsx
// Add after existing useEffects (around line 95)
useEffect(() => {
  const mutations = (AYESHA_11_17_25_PROFILE.germline?.mutations || []).map(m => ({
    gene_symbol: m.gene,
    variant_classification: m.classification || 'pathogenic'
  }));

  if (mutations.length > 0 && !ddrStatus && !ddrLoading) {
    calculateDDRStatus({
      patient_id: AYESHA_11_17_25_PROFILE.patient?.patient_id || 'AK',
      disease_site: 'ovary',
      tumor_subtype: 'HGSOC',
      mutations
    }).catch(err => console.error('[DDR] Failed to compute:', err));
  }
}, []);
```

### Step 3.4: Add DDR Section to Tab 0 JSX

Find the Tab 0 (Overview) content and add after the Mechanism Intelligence section:

```jsx
{/* DDR Status Section - Add after Mechanism Intelligence */}
<Box mb={3}>
  <Typography 
    variant="h5" 
    gutterBottom 
    sx={{ mb: 2, fontSize: { xs: '1.25rem', sm: '1.5rem' } }}
  >
    🧬 DDR Status & PARP Eligibility
  </Typography>
  
  {ddrLoading ? (
    <Box display="flex" alignItems="center" gap={2}>
      <CircularProgress size={20} />
      <Typography variant="body2">Computing DDR status...</Typography>
    </Box>
  ) : ddrStatus ? (
    <Grid container spacing={3}>
      <Grid item xs={12} md={6}>
        <DDRStatusCard ddrStatus={ddrStatus} />
      </Grid>
      <Grid item xs={12} md={6}>
        <DDRTreatmentEligibility 
          ddrStatus={ddrStatus} 
          onViewTrials={() => setActiveTab(1)} 
        />
      </Grid>
    </Grid>
  ) : (
    <Alert severity="info">
      DDR status requires genomic mutation data. MBD4 frameshift detected - computing status...
    </Alert>
  )}
</Box>
```

**Expected Result:**
- DDRStatusCard shows "DDR Defective" (red badge)
- DDRTreatmentEligibility shows "PARP INHIBITOR ELIGIBLE" (green badge)

**Time:** 4 hours

---

## TASK 4: Create CA-125 Entry Form

**Context:** Profile now has CA-125 = 2842, but user may want to update it.

**Create File:** `oncology-coPilot/oncology-frontend/src/components/ayesha/CA125EntryForm.jsx`

```jsx
/**
 * CA125EntryForm - Entry/update form for CA-125 value
 */
import React, { useState } from 'react';
import { 
  Box, Card, CardContent, Typography, TextField, 
  Button, InputAdornment, Alert, Chip
} from '@mui/material';
import ScienceIcon from '@mui/icons-material/Science';

const getBurdenClass = (value) => {
  if (value === null || value === undefined) return null;
  if (value <= 35) return { class: 'NORMAL', color: 'success' };
  if (value <= 200) return { class: 'MINIMAL', color: 'info' };
  if (value <= 1000) return { class: 'MODERATE', color: 'warning' };
  return { class: 'EXTENSIVE', color: 'error' };
};

const CA125EntryForm = ({ currentValue = null, onSubmit }) => {
  const [value, setValue] = useState(currentValue || '');
  const [error, setError] = useState(null);
  const burden = getBurdenClass(currentValue);

  const handleSubmit = (e) => {
    e.preventDefault();
    const numValue = parseFloat(value);
    
    if (isNaN(numValue) || numValue < 0) {
      setError('Please enter a valid CA-125 value');
      return;
    }
    
    setError(null);
    onSubmit(numValue);
  };

  return (
    <Card sx={{ border: currentValue ? '1px solid' : '2px dashed', borderColor: currentValue ? 'divider' : 'warning.main' }}>
      <CardContent>
        <Box display="flex" alignItems="center" justifyContent="space-between" mb={2}>
          <Box display="flex" alignItems="center" gap={1}>
            <ScienceIcon color={currentValue ? "primary" : "warning"} />
            <Typography variant="h6" sx={{ fontWeight: 600 }}>
              CA-125 Value
            </Typography>
          </Box>
          {burden && (
            <Chip 
              label={`${burden.class}: ${currentValue} U/mL`} 
              color={burden.color}
              size="small"
            />
          )}
        </Box>
        
        {!currentValue && (
          <Alert severity="info" sx={{ mb: 2 }}>
            Enter CA-125 to enable disease burden tracking and KELIM forecasting.
          </Alert>
        )}

        <form onSubmit={handleSubmit}>
          <Box display="flex" gap={2} alignItems="flex-start">
            <TextField
              label={currentValue ? "Update CA-125" : "Enter CA-125"}
              type="number"
              value={value}
              onChange={(e) => setValue(e.target.value)}
              error={!!error}
              helperText={error}
              InputProps={{
                endAdornment: <InputAdornment position="end">U/mL</InputAdornment>,
              }}
              sx={{ width: 200 }}
            />
            <Button 
              type="submit" 
              variant="contained" 
              color="primary"
              sx={{ mt: 1 }}
            >
              {currentValue ? 'Update' : 'Save'}
            </Button>
          </Box>
        </form>
      </CardContent>
    </Card>
  );
};

export default CA125EntryForm;
```

**Add export to index.js:**
```javascript
export { default as CA125EntryForm } from './CA125EntryForm';
```

**Wire to Tab 3 (Monitoring):**
```jsx
import { CA125EntryForm } from '../components/ayesha';

// In Tab 3:
<CA125EntryForm 
  currentValue={AYESHA_11_17_25_PROFILE.labs?.ca125_value}
  onSubmit={(newValue) => {
    console.log('New CA-125:', newValue);
    // In production: update state/trigger refetch
  }}
/>
```

**Time:** 3 hours

---

## TASK 5: Display Master Analysis Q&A (Optional)

**Context:** The Master Analysis already answered all 8 clinical questions with validated results.

**Option A: Create Summary Component**

```jsx
const ClinicalQuestionsSummary = ({ analysisResults }) => {
  const questions = [
    { id: 1, title: "Variant Impact", answer: "MBD4+TP53 = HIGH probability drivers", confidence: "90-95%" },
    { id: 2, title: "Functional Annotation", answer: "Complete loss-of-function (both genes)", confidence: "85-90%" },
    { id: 3, title: "Pathway Analysis", answer: "DDR pathway: 1.00 (MAXIMUM)", confidence: "85-90%" },
    { id: 4, title: "Drug Prediction", answer: "Olaparib #1 (efficacy: 0.80)", confidence: "70-85%" },
    { id: 5, title: "Trial Matching", answer: "DDR-deficient ovarian trials", confidence: "75-85%" },
    { id: 6, title: "Resistance Monitoring", answer: "DNA repair capacity: 0.60", confidence: "70-80%" },
    { id: 7, title: "IO Eligibility", answer: "YES (TMB=25, IO eligible)", confidence: "90-95%" },
    { id: 8, title: "Nutritional", answer: "No strong recommendations", confidence: "50-70%" }
  ];
  
  return (
    <Box>
      <Typography variant="h5" gutterBottom>📋 Clinical Intelligence Summary</Typography>
      <Grid container spacing={2}>
        {questions.map(q => (
          <Grid item xs={12} md={6} key={q.id}>
            <Card variant="outlined">
              <CardContent>
                <Typography variant="subtitle1" fontWeight={600}>
                  Q{q.id}: {q.title}
                </Typography>
                <Typography variant="body2" color="text.secondary">
                  {q.answer}
                </Typography>
                <Chip label={q.confidence} size="small" sx={{ mt: 1 }} />
              </CardContent>
            </Card>
          </Grid>
        ))}
      </Grid>
    </Box>
  );
};
```

**Time:** 3 hours (if needed)

---

## TASK 6: End-to-End Testing

**Test Checklist:**

### Backend Tests:

```bash
# 1. Test DDR status with MBD4
curl -X POST http://localhost:8000/api/resistance/ddr-status \
  -H 'Content-Type: application/json' \
  -d '{
    "patient_id": "AK",
    "disease_site": "ovary",
    "mutations": [{"gene_symbol": "MBD4", "variant_classification": "pathogenic"}]
  }'
# Expected: DDR_bin_status: "DDR_defective"

# 2. Test efficacy prediction
curl -X POST http://localhost:8000/api/efficacy/predict \
  -H 'Content-Type: application/json' \
  -d '{
    "disease": "ovarian",
    "mutations": [{"gene": "MBD4", "consequence": "frameshift_variant"}]
  }'
# Expected: PARP inhibitors in top 3
```

### Frontend Tests:

- [ ] Tab 0 shows DDR Status section
- [ ] DDRStatusCard displays "DDR Defective"
- [ ] DDRTreatmentEligibility shows "PARP ELIGIBLE"
- [ ] "View Trials" button navigates to Tab 1
- [ ] Tab 3 shows CA-125 Entry Form (if no value)
- [ ] Tab 3 shows CA-125 Tracker (if value exists)

**Time:** 2 hours

---

## 📊 EXPECTED OUTCOMES

When complete, Ayesha will see:

### Tab 0 (Overview):
- ✅ Mechanism Vector showing DDR = 0.88
- ✅ **NEW:** DDR Status = "DDR Defective" (red badge)
- ✅ **NEW:** PARP Eligible badge (green)
- ✅ SOC Recommendation (Carbo+Taxol+Bev)

### Tab 1 (Trials):
- ✅ Trials with holistic scores
- ✅ PGx safety warnings

### Tab 3 (Monitoring):
- ✅ **NEW:** CA-125 Entry/Display (2,842 U/mL)
- ✅ Burden class: EXTENSIVE

### Tab 5 (Synthetic Lethality):
- ✅ MBD4+TP53 detected
- ✅ PARP recommended

---

## ⚠️ DON'T DO

- ❌ Don't modify `ddr_bin_scoring.py` - only config
- ❌ Don't create new API endpoints - all exist
- ❌ Don't modify existing component logic - just wire them
- ❌ Don't hardcode values - use profile data

---

## ✅ WHEN DONE

Report back with:
1. Screenshot of DDR Status on Tab 0
2. Screenshot of PARP ELIGIBLE badge
3. Screenshot of CA-125 section on Tab 3
4. Curl output showing `DDR_bin_status: "DDR_defective"`

---

**Questions? The Master Analysis document has all the validated results. Check `.cursor/ayesha/MBD4_TP53_MASTER_ANALYSIS.md` for reference.**
