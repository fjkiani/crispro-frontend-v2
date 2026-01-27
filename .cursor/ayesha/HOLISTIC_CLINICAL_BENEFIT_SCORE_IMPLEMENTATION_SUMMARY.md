# 🚀 HOLISTIC CLINICAL BENEFIT SCORE - IMPLEMENTATION SUMMARY

**Date**: January 29, 2025  
**Status**: ✅ **ALL QUESTIONS RESOLVED - READY FOR IMPLEMENTATION**  
**Priority**: **P0** - Integrator task after engines stable  
**Timeline**: 5-6 days (same sprint after DDR_bin, Timing, CSI v0, KELIM MVP complete)

---

## 📋 MANAGER DECISIONS SUMMARY

### **✅ Q5: Weight Normalization (Trial Enrollment)**
**Decision**: Option B + C
- For `use_case = "trial_enrollment"`: Set `w_T = 0` and renormalize other weights (D, P, M, S only)
- **Still compute and return T** in breakdown for transparency
- Do not let T affect overall enrollment score

**Rationale**: Pre-treatment trial matching should not depend on on-treatment kinetics.

### **✅ Q7: CSI Integration**
**Decision**: Option C (CSI-first, with graceful fallback)
- **When CSI available**: Use CSI's mechanism-fit for M, CSI's baseline-risk for P_contrib
- **When CSI not available**: Fall back to 7D mechanism vectors and Timing/BRI model
- **Always expose M and P explicitly** in output
- **Document in provenance** whether M/P came from CSI or separate calculations

### **✅ Q8: Integration Point**
**Decision**: Option C (Complement now, replace after validation)

**Phase 1 (Current Sprint):**
- Complement existing ranking logic (keep `mechanism_fit_ranker` and `holistic_score`)
- Add Holistic Clinical Benefit Score as **parallel metric**
- New endpoint: `POST /api/resistance/holistic-clinical-benefit`
- Return alongside existing scores (don't replace yet)

**Phase 2 (After Validation):**
- Promote to primary ranking metric for specific flows
- Deprecate older rankers where redundant

---

## 🎯 IMPLEMENTATION PLAN

### **Phase 1: Component Functions** (2-3 days)

#### **1.1 Create Service Structure**
```
api/services/holistic_clinical_benefit/
├── __init__.py
├── service.py           # Main orchestration
├── diagnostic_fit.py    # Component D
├── prognostic_risk.py   # Component P
├── therapeutic_dynamics.py  # Component T
└── config.py            # Weights & thresholds
```

#### **1.2 Component Functions**

**D: Diagnostic Fit** (`diagnostic_fit.py`)
- Inputs: DDR features, disease_site, clinical_features
- Logic: Start from 1.0, subtract penalties for missing info
- Returns: float [0, 1]

**P: Prognostic Risk** (`prognostic_risk.py`)
- Inputs: timing_features, disease_site, regimen_class, clinical_features
- Logic: Simple disease-specific model predicting PFS probability
- CSI integration: Use CSI's baseline-risk if available
- Returns: float [0, 1]

**M: Predictive Mechanism Fit** (`service.py` - reuse existing)
- Import `compute_mechanism_fit()` from `holistic_score.mechanism_fit`
- CSI integration: Use CSI's mechanism-fit if available, else 7D vectors
- Returns: float [0, 1]

**T: Therapeutic Dynamics** (`therapeutic_dynamics.py`)
- Inputs: kinetic_features, regimen_start_date, disease_site, regimen_class
- Logic: Map KELIM/CA-125 metrics to T via disease-specific rules
- Neutral prior (0.5) if insufficient data
- Returns: float [0, 1]

**S: Safety/Tolerability** (`service.py` - reuse existing)
- Import `compute_pgx_safety()` from `holistic_score.pgx_safety`
- Returns: float [0, 1]

#### **1.3 Regimen MoA Vector Computation**
**Reuse Trial MoA Tagging Pattern:**
```python
def compute_regimen_moa_vector(regimen_drugs: List[str], regimen_type: Optional[str] = None) -> List[float]:
    """Compute 7D MoA vector using same pattern as trials."""
    from api.services.client_dossier.dossier_generator import get_drug_mechanism
    from api.services.pathway_to_mechanism_vector import convert_moa_dict_to_vector
    
    moa_dict = {"ddr": 0.0, "mapk": 0.0, "pi3k": 0.0, "vegf": 0.0, "her2": 0.0, "io": 0.0, "efflux": 0.0}
    
    # Map drugs to mechanisms using keyword matching (same as trials)
    for drug in regimen_drugs:
        mechanism = get_drug_mechanism(drug)
        if mechanism and mechanism.get("mechanism"):
            mech_text = mechanism["mechanism"].lower()
            # PARP/ATR → DDR, MEK/BRAF → MAPK, etc.
            if any(kw in mech_text for kw in ["parp", "atr", "chk", "wee1", "dna repair"]):
                moa_dict["ddr"] = max(moa_dict["ddr"], 0.9)
            # ... (repeat for other pathways)
    
    return convert_moa_dict_to_vector(moa_dict, use_7d=True)
```

---

### **Phase 2: Configuration System** (1 day)

#### **2.1 Create Config File**
**File**: `api/services/resistance/config/holistic_clinical_benefit_config.py`

```python
# Use-case default weights (disease-agnostic)
DEFAULT_WEIGHTS = {
    "trial_enrollment": {
        "D": 0.20,
        "P": 0.10,
        "M": 0.45,
        "T": 0.00,  # Set to 0 for trial enrollment
        "S": 0.25
    },
    "next_line": {
        "D": 0.10,
        "P": 0.25,
        "M": 0.35,
        "T": 0.20,
        "S": 0.10
    },
    "monitoring": {
        "D": 0.05,
        "P": 0.20,
        "M": 0.15,
        "T": 0.45,
        "S": 0.15
    }
}

# Disease-specific overrides (optional)
DISEASE_WEIGHTS = {
    "ovary": {
        "trial_enrollment": {...},
        ...
    }
}

def get_weights(use_case: str, disease_site: str, regimen_class: str) -> Dict[str, float]:
    """Get weights per disease_site + regimen_class + use_case."""
    # Check for disease/regimen-specific weights first
    # Fall back to use-case defaults
    # Normalize if components missing
    pass
```

---

### **Phase 3: Main Service** (1-2 days)

#### **3.1 Main Function**
**File**: `api/services/holistic_clinical_benefit/service.py`

```python
async def compute_holistic_clinical_benefit_score(
    patient_id: str,
    regimen_id: str,
    use_case: str,  # "trial_enrollment" | "next_line" | "monitoring"
    ddr_features: Optional[Dict] = None,        # Auto-fetched if None
    timing_features: Optional[Dict] = None,     # Auto-fetched if None
    kinetic_features: Optional[Dict] = None,    # Auto-fetched if None
    safety_features: Optional[Dict] = None,     # Auto-fetched if None
    clinical_features: Optional[Dict] = None,   # Auto-fetched if None
    csi_outputs: Optional[Dict] = None,         # CSI outputs if available
    config: Optional[Dict] = None
) -> Dict[str, Any]:
    """
    Compute Holistic Clinical Benefit Score for patient-regimen pair.
    
    Returns:
    {
        "D": float,      # Diagnostic fit
        "P": float,      # Prognostic risk
        "M": float,      # Predictive mechanism fit
        "T": float,      # Therapeutic response dynamics
        "S": float,      # Safety/tolerability
        "overall": float,  # Overall score
        "weights": Dict[str, float],  # Normalized weights
        "component_available": Dict[str, bool],
        "breakdown": Dict[str, Any],
        "interpretation": str,  # "HIGH" | "MEDIUM" | "LOW"
        "recommendation": str,
        "provenance": Dict[str, Any]
    }
    """
    # 1. Auto-fetch features from engines if not provided
    # 2. Get weights from config (per disease_site + regimen_class + use_case)
    # 3. Normalize weights if components missing (or T=0 for trial_enrollment)
    # 4. Compute each component (with CSI integration for M/P)
    # 5. Compute overall score using normalized weights
    # 6. Generate interpretation
    # 7. Return with provenance
```

#### **3.2 Feature Fetching**
**Import engine functions directly** (matches existing patterns):
```python
from api.services.resistance.biomarkers.diagnostic.ddr_bin_scoring import assign_ddr_status
from api.services.resistance.biomarkers.therapeutic.timing_chemo_features import build_timing_chemo_features
from api.services.holistic_score.pgx_safety import compute_pgx_safety
```

#### **3.3 Weight Normalization**
**Special handling for trial enrollment:**
```python
def _normalize_weights_for_missing_components(weights: Dict, component_available: Dict, use_case: str) -> Dict:
    """Normalize weights, with special handling for trial_enrollment."""
    # If trial_enrollment: set T weight to 0, renormalize others
    if use_case == "trial_enrollment":
        weights["T"] = 0.0
    
    # Set weight to 0 for missing components
    adjusted = {k: weights[k] if component_available[k] else 0.0 for k in ["D", "P", "M", "T", "S"]}
    
    # Renormalize to sum to 1.0
    total = sum(adjusted.values())
    return {k: v / total for k, v in adjusted.items()} if total > 0 else {k: 0.2 for k in ["D", "P", "M", "T", "S"]}
```

#### **3.4 CSI Integration**
**CSI-first with graceful fallback:**
```python
# M: Mechanism Fit
if csi_outputs and csi_outputs.get("mechanism_fit"):
    M = csi_outputs["mechanism_fit"]
    m_source = "csi"
else:
    # Fall back to 7D mechanism vectors
    patient_vector = await _fetch_patient_mechanism_vector(patient_id)
    regimen_vector = await _fetch_regimen_moa_vector(regimen_id)
    M = compute_mechanism_fit(patient_vector, regimen_vector)
    m_source = "mechanism_vectors"

# P: Prognostic Risk
if csi_outputs and csi_outputs.get("baseline_risk"):
    P = csi_outputs["baseline_risk"]
    p_source = "csi"
else:
    # Fall back to Timing/BRI model
    P = compute_prognostic_risk(timing_features, ...)
    p_source = "timing_model"
```

---

### **Phase 4: API Integration** (1 day)

#### **4.1 Create Endpoint**
**File**: `api/routers/resistance.py`

```python
@router.post("/api/resistance/holistic-clinical-benefit")
async def get_holistic_clinical_benefit(request: HolisticClinicalBenefitRequest):
    """
    Compute Holistic Clinical Benefit Score for patient-regimen pair.
    
    Phase 1: Complement existing ranking (parallel metric)
    Phase 2: Promote to primary ranking after validation
    """
    result = await compute_holistic_clinical_benefit_score(
        patient_id=request.patient_id,
        regimen_id=request.regimen_id,
        use_case=request.use_case,
        ...
    )
    return HolisticClinicalBenefitResponse(**result)
```

#### **4.2 Pydantic Models**
```python
class HolisticClinicalBenefitRequest(BaseModel):
    patient_id: str
    regimen_id: str
    use_case: Literal["trial_enrollment", "next_line", "monitoring"]
    # Optional: override features (for testing)
    ddr_features: Optional[Dict] = None
    timing_features: Optional[Dict] = None
    ...

class HolisticClinicalBenefitResponse(BaseModel):
    D: float
    P: float
    M: float
    T: float
    S: float
    overall: float
    weights: Dict[str, float]
    component_available: Dict[str, bool]
    breakdown: Dict[str, Any]
    interpretation: str
    recommendation: str
    provenance: Dict[str, Any]
```

---

### **Phase 5: Testing** (1 day)

#### **5.1 Unit Tests**
- Each component function (D, P, M, T, S)
- Weight normalization logic
- CSI integration (with/without CSI)
- Missing component handling

#### **5.2 Integration Tests**
- End-to-end with real engine outputs
- Synthetic patient-regimen pairs
- All three use cases (trial_enrollment, next_line, monitoring)

#### **5.3 Validation Tests**
- Known patient-regimen pairs
- Compare to expected scores
- Side-by-side with existing holistic_score

---

## 🎯 IMPLEMENTATION CHECKLIST

### **Day 1-2: Component Functions**
- [ ] Create service structure
- [ ] Implement D (Diagnostic Fit)
- [ ] Implement P (Prognostic Risk) with CSI integration
- [ ] Implement T (Therapeutic Dynamics)
- [ ] Reuse M (Mechanism Fit) with CSI integration
- [ ] Reuse S (Safety/Tolerability)
- [ ] Implement regimen MoA vector computation

### **Day 3: Configuration & Main Service**
- [ ] Create config file with weights
- [ ] Implement main service function
- [ ] Implement feature auto-fetching
- [ ] Implement weight normalization (with trial_enrollment special case)
- [ ] Implement CSI integration (CSI-first, graceful fallback)
- [ ] Implement interpretation logic

### **Day 4: API Integration**
- [ ] Create Pydantic models
- [ ] Create API endpoint
- [ ] Add to router
- [ ] Test endpoint

### **Day 5-6: Testing & Documentation**
- [ ] Unit tests for all components
- [ ] Integration tests
- [ ] Validation tests
- [ ] Update architecture doc with implementation details
- [ ] Document CSI integration approach

---

## 🔗 KEY INTEGRATION POINTS

### **Reuse Existing Components:**
1. **Mechanism Fit**: `holistic_score.mechanism_fit.compute_mechanism_fit()`
2. **PGx Safety**: `holistic_score.pgx_safety.compute_pgx_safety()`
3. **Regimen MoA Vector**: Reuse trial MoA tagging pattern
4. **Interpretation Thresholds**: Same as holistic_score (HIGH ≥0.8, MEDIUM ≥0.6, LOW ≥0.4)

### **Engine Integration:**
1. **DDR Features**: `assign_ddr_status()` from DDR_bin engine
2. **Timing Features**: `build_timing_chemo_features()` from Timing engine
3. **Kinetic Features**: KELIM/CA-125 from Kinetic engine
4. **CSI Outputs**: CSI v0 outputs (when available)

### **Default Values:**
- **D/P/M/T**: 0.5 (neutral) if missing
- **S**: 1.0 (assume safe) if missing
- **Weights**: Normalize if components missing

---

## 📊 SUCCESS CRITERIA

### **Functionality:**
- ✅ All 5 components (D, P, M, T, S) computed correctly
- ✅ Overall score computed with normalized weights
- ✅ Trial enrollment sets T weight to 0, still computes T
- ✅ CSI integration works (CSI-first, graceful fallback)
- ✅ Missing components handled gracefully

### **Integration:**
- ✅ API endpoint works
- ✅ Returns alongside existing scores (complement, not replace)
- ✅ Clear provenance (M/P source: CSI vs separate)

### **Testing:**
- ✅ Unit tests pass
- ✅ Integration tests pass
- ✅ Validation tests match expected scores

---

**Last Updated**: January 29, 2025  
**Status**: ✅ **READY FOR IMPLEMENTATION**  
**All Questions Resolved - Manager Decisions Received**
