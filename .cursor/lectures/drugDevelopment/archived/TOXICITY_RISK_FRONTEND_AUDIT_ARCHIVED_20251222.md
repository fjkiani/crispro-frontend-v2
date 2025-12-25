# ⚠️ Toxicity Risk Assessment - Frontend Audit & Integration Plan

**Purpose:** Audit current toxicity risk capabilities and plan standalone page + Advanced Care Plan integration  
**Date:** January 28, 2025  
**Source of Truth:** `.cursor/MOAT/ADVANCED_CARE_PLAN_TOXCITY.md`  
**Reference:** `.cursor/lectures/drugDevelopment/toxicity_risk_contribution.mdc`

---

## 📊 CURRENT STATE AUDIT

### ✅ **Backend Implementation (100% Complete)**

| Component | Status | Location | Notes |
|-----------|--------|----------|-------|
| **API Endpoint** | ✅ Complete | `api/routers/safety.py` | `/api/safety/toxicity_risk` - Fully operational |
| **Safety Service** | ✅ Complete | `api/services/safety_service.py` | `compute_toxicity_risk()` - Three-factor model implemented |
| **Pathway Mappings** | ✅ Complete | `api/services/toxicity_pathway_mappings.py` | 30+ pharmacogenes, 11 MoA patterns, 3 pathways |
| **Schemas** | ✅ Complete | `api/schemas/safety.py` | `ToxicityRiskRequest`, `ToxicityRiskResponse`, `ToxicityFactor` |
| **Mitigating Foods** | ✅ Complete | `toxicity_pathway_mappings.py` | `get_mitigating_foods()` - THE MOAT implemented |

**Backend Capabilities:**
- ✅ Risk score calculation (0-1)
- ✅ Risk level classification (HIGH/MODERATE/LOW)
- ✅ Contributing factors (pharmacogene, pathway, tissue)
- ✅ Confidence adjustment (conservative for high-risk)
- ✅ Complete provenance (run_id, profile, methods, timestamp)
- ✅ Mitigating foods mapping (DNA repair, inflammation, cardiometabolic)

---

### ⚠️ **Frontend Implementation (60% Complete)**

| Component | Status | Location | Notes |
|-----------|--------|----------|-------|
| **ToxicityChip** | ⚠️ Placeholder | `components/vus/ToxicityChip.jsx` | NOT wired to API - shows placeholder only |
| **ToxicityRiskCard** | ✅ Complete | `components/ClinicalGenomicsCommandCenter/cards/ToxicityRiskCard.jsx` | Fully wired, displays risk score, factors, confidence |
| **useToxicity Hook** | ✅ Complete | `components/ClinicalGenomicsCommandCenter/hooks/useToxicity.js` | Calls `/api/safety/toxicity_risk` correctly |
| **CoPilot Integration** | ✅ Complete | `integrations/ClinicalGenomicsCoPilotIntegration.jsx` | 3 quick actions, suggested questions |
| **Standalone Page** | ❌ Missing | N/A | No dedicated `/toxicity-risk` route or page |

**Current Usage:**
- ✅ Used in `MechanisticEvidenceTab.jsx` (ClinicalGenomicsCommandCenter)
- ✅ Used in `AnalysisResults.jsx` (VUS analysis - placeholder chip)
- ❌ NOT used in Complete Care Plan pages
- ❌ NOT a standalone page

---

## 🎯 WHAT'S MISSING FOR STANDALONE PAGE

### **1. Standalone Page Component** ❌

**Required:**
- New page: `pages/ToxicityRiskAssessment.jsx`
- Route: `/toxicity-risk` or `/safety/toxicity-risk`
- Full-page layout with:
  - Patient input form (germline variants, drug selection)
  - Real-time assessment
  - Detailed results display
  - Export functionality

**Current Gap:**
- No standalone page exists
- Only embedded in ClinicalGenomicsCommandCenter

---

### **2. Patient Input Form** ❌

**Required:**
- Germline variant input (VCF upload or manual entry)
- Drug selection (dropdown with MoA mapping)
- Disease context selection
- Treatment line selection (optional)

**Current Gap:**
- No dedicated input form
- Relies on parent component passing data

---

### **3. Multi-Drug Assessment** ❌

**Required:**
- Assess multiple drugs simultaneously
- Compare toxicity risks across drug options
- Show risk ranking (lowest to highest)

**Current Gap:**
- Only single-drug assessment supported
- No comparison view

---

### **4. Mitigating Foods Display** ⚠️ Partial

**Current State:**
- ✅ Backend returns `mitigating_foods` in response
- ❌ Frontend `ToxicityRiskCard` does NOT display mitigating foods
- ✅ Food validator shows toxicity mitigation badge (THE MOAT)

**Required:**
- Display mitigating foods in ToxicityRiskCard
- Show timing guidance ("post-chemo, not during")
- Link to food validation for detailed recommendations

---

### **5. Advanced Care Plan Integration** ⚠️ Partial

**Current State:**
- ✅ Complete Care Plan Universal endpoint exists (`/api/complete_care/universal`)
- ❌ Does NOT call toxicity risk assessment
- ✅ Food validator integrates toxicity mitigation (THE MOAT)

**Required:**
- Add toxicity risk assessment to Complete Care Plan flow
- Display toxicity risk in care plan summary
- Show mitigating foods recommendations
- Link toxicity risk to drug recommendations

---

## 🔗 HOW TOXICITY RISK SUPPORTS ADVANCED CARE PLAN

### **Integration Points (Per ADVANCED_CARE_PLAN_TOXCITY.md)**

#### **1. Pre-Enrollment Toxicity Screening** ✅

**From ADVANCED_CARE_PLAN_TOXCITY.md:**
> "Pre-enrollment toxicity screening flags patients with germline toxicity factors before enrollment. HIGH risk (≥0.5) → enhanced monitoring or alternative therapies."

**Current Implementation:**
- ✅ Backend supports this (`compute_toxicity_risk()`)
- ⚠️ Frontend needs standalone page for pre-enrollment screening
- ❌ Not integrated into Complete Care Plan flow

**Integration Required:**
```javascript
// In Complete Care Plan flow:
const toxicityRisk = await assessToxicityRisk(
  patient.germlineVariants,
  recommendedDrugs.map(d => d.moa),
  patient.disease
);

// Flag HIGH risk drugs
const highRiskDrugs = recommendedDrugs.filter(d => 
  toxicityRisk[d.moa]?.risk_score >= 0.5
);
```

---

#### **2. Toxicity-Aware Nutrition (THE MOAT)** ✅

**From ADVANCED_CARE_PLAN_TOXCITY.md:**
> "When you validate a food, the system now checks: What drugs are you on? What's your germline profile? What toxicity pathways are stressed? Does this food help those pathways?"

**Current Implementation:**
- ✅ Backend: `get_mitigating_foods()` implemented
- ✅ Food validator: Shows toxicity mitigation badge
- ❌ ToxicityRiskCard: Does NOT display mitigating foods

**Integration Required:**
- Display mitigating foods in ToxicityRiskCard
- Show timing guidance ("post-chemo, not during")
- Link to food validation for detailed recommendations

---

#### **3. Pharmacogene Risk Flagging** ✅

**From ADVANCED_CARE_PLAN_TOXCITY.md:**
> "DPYD variant + 5-FU: Can't break down 5-FU → Toxic levels → Severe diarrhea, death (5-10% mortality)"

**Current Implementation:**
- ✅ Backend: Pharmacogene detection (DPYD, TPMT, UGT1A1, etc.)
- ✅ Risk weights: High-impact (0.4), CYP (0.3), Others (0.2)
- ⚠️ Frontend: Shows in ToxicityRiskCard but not prominently flagged

**Integration Required:**
- Prominent warning for high-impact pharmacogenes (DPYD, TPMT)
- Dose adjustment recommendations
- Alternative drug suggestions

---

#### **4. Complete Care Plan Integration** ❌

**From ADVANCED_CARE_PLAN_UNIVERSAL.md:**
> "Complete care plan includes: Drug recommendations, Trials, Food validation, Monitoring, **Toxicity screening**"

**Current State:**
- ❌ Complete Care Plan does NOT call toxicity risk assessment
- ❌ No toxicity risk section in care plan output
- ✅ Food validator has toxicity mitigation (separate flow)

**Integration Required:**
```python
# In complete_care_universal.py:
async def _assess_toxicity_risks(
    patient_profile: Dict[str, Any],
    recommended_drugs: List[Dict[str, Any]]
) -> Dict[str, Any]:
    """Assess toxicity risks for all recommended drugs."""
    toxicity_results = {}
    
    for drug in recommended_drugs:
        moa = drug.get("moa")
        if not moa:
            continue
        
        # Call toxicity risk API
        request = ToxicityRiskRequest(
            patient=PatientContext(
                germlineVariants=patient_profile.get("germlineVariants", [])
            ),
            candidate=TherapeuticCandidate(type="drug", moa=moa),
            context=ClinicalContext(disease=patient_profile.get("disease"))
        )
        
        result = await safety_service.compute_toxicity_risk(request)
        toxicity_results[drug["name"]] = result
    
    return toxicity_results
```

---

## 📋 IMPLEMENTATION ROADMAP

### **Phase 1: Standalone Toxicity Risk Page** (Priority: HIGH)

**Deliverables:**
1. **New Page Component** (`pages/ToxicityRiskAssessment.jsx`)
   - Patient input form (germline variants, drug selection)
   - Real-time assessment
   - Results display (ToxicityRiskCard)
   - Export functionality

2. **Route Addition** (`App.jsx`)
   ```javascript
   <Route path="/toxicity-risk" element={<ToxicityRiskAssessment />} />
   <Route path="/toxicity-risk/:patientId" element={<ToxicityRiskAssessment />} />
   ```

3. **Multi-Drug Support**
   - Assess multiple drugs simultaneously
   - Comparison table (risk scores, factors, mitigating foods)
   - Risk ranking (lowest to highest)

**Estimated Time:** 8-12 hours

---

### **Phase 2: ToxicityRiskCard Enhancement** (Priority: HIGH)

**Deliverables:**
1. **Display Mitigating Foods**
   - Show mitigating foods from `result.mitigating_foods`
   - Display timing guidance ("post-chemo, not during")
   - Link to food validation for detailed recommendations

2. **Prominent Pharmacogene Warnings**
   - High-impact pharmacogenes (DPYD, TPMT) → Red alert
   - Dose adjustment recommendations
   - Alternative drug suggestions

3. **Export Functionality**
   - PDF export
   - JSON export
   - Shareable link

**Estimated Time:** 4-6 hours

---

### **Phase 3: Complete Care Plan Integration** (Priority: MEDIUM)

**Deliverables:**
1. **Backend Integration** (`complete_care_universal.py`)
   - Add `_assess_toxicity_risks()` function
   - Call for all recommended drugs
   - Include in response

2. **Frontend Display** (`AyeshaCompleteCare.jsx` or new UniversalCarePlan page)
   - Toxicity risk section in care plan
   - Risk chips for each drug
   - Mitigating foods recommendations
   - Link to detailed toxicity assessment

**Estimated Time:** 6-8 hours

---

### **Phase 4: ToxicityChip Wiring** (Priority: LOW)

**Deliverables:**
1. **Wire ToxicityChip to API**
   - Replace placeholder with actual API call
   - Show risk level chip (HIGH/MODERATE/LOW)
   - Tooltip with details

**Estimated Time:** 2-3 hours

---

## 🎯 STANDALONE PAGE SPECIFICATION

### **Page Structure**

```
/toxicity-risk
├── Header
│   ├── Title: "Toxicity Risk Assessment (RUO)"
│   ├── Subtitle: "Germline-based toxicity prediction for precision safety"
│   └── RUO Disclaimer
│
├── Input Section
│   ├── Patient Selection (if patientId in URL)
│   ├── Germline Variants Input
│   │   ├── VCF Upload
│   │   ├── Manual Entry
│   │   └── Load from Patient Profile
│   ├── Drug Selection
│   │   ├── Single Drug (dropdown)
│   │   ├── Multiple Drugs (multi-select)
│   │   └── MoA Auto-Detection
│   └── Clinical Context
│       ├── Disease Selection
│       ├── Treatment Line (optional)
│       └── Tissue Context (optional)
│
├── Assessment Results
│   ├── Single Drug View (if one drug selected)
│   │   └── ToxicityRiskCard (enhanced)
│   │       ├── Risk Score Visualization
│   │       ├── Risk Level Chip (HIGH/MODERATE/LOW)
│   │       ├── Confidence Chip
│   │       ├── Helper Text
│   │       ├── Contributing Factors
│   │       ├── Mitigating Foods (NEW)
│   │       └── Provenance
│   │
│   └── Multi-Drug Comparison (if multiple drugs selected)
│       ├── Comparison Table
│       │   ├── Drug Name
│       │   ├── Risk Score
│       │   ├── Risk Level
│       │   ├── Key Factors
│       │   └── Mitigating Foods
│       └── Risk Ranking (lowest to highest)
│
└── Actions
    ├── Export PDF
    ├── Export JSON
    ├── Share Link
    └── Add to Care Plan
```

---

### **Enhanced ToxicityRiskCard Specification**

**New Fields to Display:**

1. **Mitigating Foods Section** (NEW)
   ```jsx
   {result.mitigating_foods && result.mitigating_foods.length > 0 && (
     <Box sx={{ mt: 2 }}>
       <Typography variant="subtitle2" gutterBottom>
         Mitigating Foods (THE MOAT):
       </Typography>
       {result.mitigating_foods.map((food, idx) => (
         <Card key={idx} sx={{ mb: 1 }}>
           <CardContent>
             <Typography variant="body1" fontWeight="bold">
               {food.compound}
             </Typography>
             <Typography variant="body2" color="text.secondary">
               Dose: {food.dose}
             </Typography>
             <Typography variant="body2" color="text.secondary">
               Timing: {food.timing}
             </Typography>
             <Typography variant="body2">
               {food.mechanism}
             </Typography>
             <Chip 
               label={food.evidence_tier}
               size="small"
               color={food.evidence_tier === "SUPPORTED" ? "success" : "default"}
             />
           </CardContent>
         </Card>
       ))}
     </Box>
   )}
   ```

2. **Prominent Pharmacogene Warnings** (ENHANCED)
   ```jsx
   {factors.some(f => f.type === "germline" && f.weight >= 0.4) && (
     <Alert severity="error" sx={{ mt: 2 }}>
       <AlertTitle>High-Impact Pharmacogene Detected</AlertTitle>
       {factors
         .filter(f => f.type === "germline" && f.weight >= 0.4)
         .map(f => (
           <Typography key={f.detail}>
             {f.detail} - Consider dose reduction or alternative therapy
           </Typography>
         ))}
     </Alert>
   )}
   ```

---

## 🔗 ADVANCED CARE PLAN INTEGRATION SPECIFICATION

### **Backend Integration** (`complete_care_universal.py`)

**Add to `get_complete_care_v2()`:**

```python
async def _assess_toxicity_risks(
    patient_profile: Dict[str, Any],
    recommended_drugs: List[Dict[str, Any]],
    safety_service: SafetyService
) -> Dict[str, Any]:
    """
    Assess toxicity risks for all recommended drugs.
    
    Returns:
        {
            "drug_name": {
                "risk_score": float,
                "risk_level": str,  # HIGH/MODERATE/LOW
                "factors": List[Dict],
                "mitigating_foods": List[Dict],
                "confidence": float
            }
        }
    """
    toxicity_results = {}
    
    # Extract germline variants
    germline_variants = []
    if "germlineVariants" in patient_profile:
        germline_variants = patient_profile["germlineVariants"]
    elif "mutations" in patient_profile:
        # Extract germline mutations
        germline_variants = [
            m for m in patient_profile["mutations"]
            if m.get("type") == "germline"
        ]
    
    # Assess each drug
    for drug in recommended_drugs:
        drug_name = drug.get("name", "Unknown")
        moa = drug.get("moa")
        
        if not moa:
            continue
        
        try:
            request = ToxicityRiskRequest(
                patient=PatientContext(
                    germlineVariants=germline_variants
                ),
                candidate=TherapeuticCandidate(type="drug", moa=moa),
                context=ClinicalContext(
                    disease=patient_profile.get("disease", "cancer")
                ),
                options={"profile": "baseline"}
            )
            
            result = await safety_service.compute_toxicity_risk(request)
            
            # Derive risk level
            risk_level = "HIGH" if result.risk_score >= 0.5 else \
                        "MODERATE" if result.risk_score >= 0.3 else "LOW"
            
            toxicity_results[drug_name] = {
                "risk_score": result.risk_score,
                "risk_level": risk_level,
                "confidence": result.confidence,
                "reason": result.reason,
                "factors": [f.dict() for f in result.factors],
                "mitigating_foods": result.mitigating_foods,
                "provenance": result.provenance
            }
        except Exception as e:
            logger.warning(f"Toxicity assessment failed for {drug_name}: {e}")
            continue
    
    return toxicity_results
```

**Add to response assembly:**

```python
# In get_complete_care_v2():
toxicity_risks = await _assess_toxicity_risks(
    patient_profile,
    results.get("wiwfm", {}).get("drugs", []),
    get_safety_service()
)

results["toxicity_risks"] = toxicity_risks
```

---

### **Frontend Integration** (`AyeshaCompleteCare.jsx` or new UniversalCarePlan page)

**Add Toxicity Risk Section:**

```jsx
{carePlan.toxicity_risks && Object.keys(carePlan.toxicity_risks).length > 0 && (
  <Box sx={{ mt: 3 }}>
    <Typography variant="h5" gutterBottom>
      Toxicity Risk Assessment
    </Typography>
    
    {Object.entries(carePlan.toxicity_risks).map(([drugName, risk]) => (
      <ToxicityRiskCard
        key={drugName}
        result={{
          risk_score: risk.risk_score,
          confidence: risk.confidence,
          reason: risk.reason,
          factors: risk.factors,
          mitigating_foods: risk.mitigating_foods  // NEW
        }}
        drugName={drugName}
      />
    ))}
  </Box>
)}
```

---

## 📊 CAPABILITY MATRIX

| Capability | Backend | Frontend | Standalone Page | Care Plan Integration |
|------------|---------|----------|-----------------|----------------------|
| **Risk Score Calculation** | ✅ | ✅ | ❌ | ❌ |
| **Risk Level Classification** | ✅ | ✅ | ❌ | ❌ |
| **Contributing Factors** | ✅ | ✅ | ❌ | ❌ |
| **Mitigating Foods** | ✅ | ⚠️ Partial | ❌ | ❌ |
| **Multi-Drug Assessment** | ✅ | ❌ | ❌ | ❌ |
| **Patient Input Form** | N/A | ❌ | ❌ | N/A |
| **Export Functionality** | N/A | ❌ | ❌ | N/A |
| **Pharmacogene Warnings** | ✅ | ⚠️ Partial | ❌ | ❌ |
| **Complete Care Plan Integration** | ❌ | ❌ | N/A | ❌ |

**Legend:**
- ✅ Complete
- ⚠️ Partial (needs enhancement)
- ❌ Missing
- N/A Not applicable

---

## 🎯 SUCCESS CRITERIA

### **Standalone Page:**
- [ ] User can input germline variants (VCF or manual)
- [ ] User can select single or multiple drugs
- [ ] Real-time toxicity assessment
- [ ] Risk level chips (HIGH/MODERATE/LOW) with color coding
- [ ] Contributing factors displayed
- [ ] Mitigating foods displayed with timing guidance
- [ ] Export functionality (PDF, JSON)
- [ ] Shareable link generation

### **Care Plan Integration:**
- [ ] Complete Care Plan calls toxicity risk assessment
- [ ] Toxicity risks displayed for all recommended drugs
- [ ] Mitigating foods shown in care plan summary
- [ ] High-risk drugs flagged prominently
- [ ] Link to detailed toxicity assessment page

### **Enhanced ToxicityRiskCard:**
- [ ] Displays mitigating foods section
- [ ] Shows timing guidance ("post-chemo, not during")
- [ ] Prominent warnings for high-impact pharmacogenes
- [ ] Export functionality
- [ ] Link to food validation for detailed recommendations

---

## 📝 IMPLEMENTATION PRIORITY

### **P0 (Critical - Blocks Product Launch):**
1. ✅ Backend implementation (DONE)
2. ⚠️ ToxicityRiskCard enhancement (display mitigating foods)
3. ❌ Standalone page creation
4. ❌ Complete Care Plan integration

### **P1 (Important - Product Enhancement):**
5. Multi-drug comparison view
6. Export functionality
7. Prominent pharmacogene warnings
8. ToxicityChip wiring

### **P2 (Nice to Have):**
9. Advanced filtering (by risk level, pharmacogene type)
10. Historical tracking (risk scores over time)
11. Patient-specific recommendations based on toxicity risk

---

## 🔗 REFERENCES

- **Source of Truth:** `.cursor/MOAT/ADVANCED_CARE_PLAN_TOXCITY.md`
- **Contribution Document:** `.cursor/lectures/drugDevelopment/toxicity_risk_contribution.mdc`
- **Concept Document:** `.cursor/rules/research/toxicity_risk_concept.mdc`
- **Backend API:** `api/routers/safety.py` - `/api/safety/toxicity_risk`
- **Frontend Components:** 
  - `components/ClinicalGenomicsCommandCenter/cards/ToxicityRiskCard.jsx`
  - `components/vus/ToxicityChip.jsx`
  - `components/ClinicalGenomicsCommandCenter/hooks/useToxicity.js`

---

**Last Updated:** January 28, 2025  
**Status:** Audit Complete - Implementation Plan Ready



