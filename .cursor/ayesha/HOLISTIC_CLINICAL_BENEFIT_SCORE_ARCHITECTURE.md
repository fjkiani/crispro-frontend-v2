# 🎯 HOLISTIC CLINICAL BENEFIT SCORE - ARCHITECTURE & INTEGRATION PLAN

**Date**: January 29, 2025  
**Status**: 📋 **DESIGN DOCUMENT**  
**Purpose**: Map existing CrisPRO engines to D-P-M-T-S framework, following validated Holistic Score pattern

---

## 🎯 EXECUTIVE SUMMARY

This document shows how to apply the **validated Holistic Score framework** (from trial matching, AUROC=0.714) to create a **Holistic Clinical Benefit Score** for **patient-regimen pairs** (not just trials). The framework combines five role-specific sub-scores:

1. **D: Diagnostic Fit** → Is this the right diagnostic context?
2. **P: Prognostic Risk** → What is my expected outlook on this regimen?
3. **M: Predictive Mechanism Fit** → How well does tumor biology align with drug mechanism?
4. **T: Therapeutic Response Dynamics** → Is the current regimen working now?
5. **S: Safety/Tolerability** → Can this patient safely receive this regimen?

**Key Insight**: We already have **80% of the pieces** - we just need to orchestrate them into this unified score.

---

## 📊 MAPPING: EXISTING ENGINES → D-P-M-T-S COMPONENTS

### **Component D: Diagnostic Fit**

**Goal**: "Is this the right diagnostic context?"

**Existing Engine**: ✅ **DDR_bin Engine**

**What We Have**:
- `assign_ddr_status()` → Returns `DDR_bin_status`, `HRD_status_inferred`, `DDR_score`
- Disease-specific configurations (ovary, breast, pancreas, prostate)
- Diagnostic biomarkers detection (BRCA, HRD, core HRR, extended DDR)

**How to Map**:
```python
def compute_diagnostic_fit(
    ddr_features: Dict[str, Any],
    disease_site: str,
    tumor_subtype: Optional[str],
    clinical_features: Dict[str, Any]
) -> float:
    """
    Compute diagnostic fit score D in [0, 1].
    
    Start from 1.0, subtract penalties:
    - Missing essential diagnostic info (unknown disease_site): -0.5
    - Discordant markers (CUP if flagged): -0.3
    - Unknown DDR_bin_status: -0.2
    - DDR_bin_status = "unknown" and required biomarkers missing: -0.3
    """
    score = 1.0
    
    # Penalty: Missing disease_site
    if not disease_site or disease_site == "unknown":
        score -= 0.5
    
    # Penalty: Discordant markers (CUP)
    if clinical_features.get("cup_discordant"):
        score -= 0.3
    
    # Penalty: Unknown DDR status when it's required for disease
    ddr_status = ddr_features.get("DDR_bin_status", "unknown")
    if ddr_status == "unknown":
        # Check if DDR status is required for this disease
        if disease_site in ["ovary", "breast", "pancreas", "prostate"]:
            # DDR status is important for these diseases
            score -= 0.2
    
    # Bonus: High-confidence diagnostic markers
    if ddr_status in ["DDR_defective", "DDR_proficient"]:
        score += 0.1  # Boost for known status
    
    return max(0.0, min(1.0, score))
```

**Input Data Flow**:
```
DDR_bin Engine → DDRStatusResponse → {
    "DDR_bin_status": "DDR_defective" | "DDR_proficient" | "unknown",
    "HRD_status_inferred": "HRD_positive" | "HRD_negative" | "unknown",
    "DDR_score": 0.0-8.0,
    "BRCA_pathogenic": bool,
    "core_HRR_pathogenic": bool,
    "extended_DDR_pathogenic": bool
}
```

**Status**: ✅ **READY** - DDR_bin engine provides all needed inputs

---

### **Component P: Prognostic Risk**

**Goal**: "What is my expected outlook on this regimen, regardless of exact mechanism?"

**Existing Engine**: ✅ **Timing & Chemosensitivity Engine**

**What We Have**:
- `build_timing_chemo_features()` → Returns `PFI_days`, `PTPI_days`, `TFI_days`, `PFS_days`, `OS_days`, `line_of_therapy`
- Disease-specific PFI cutpoints (ovary: 6 months, others configurable)
- Per-regimen PFS/OS computation

**How to Map**:
```python
def compute_prognostic_risk(
    timing_features: Dict[str, Any],
    disease_site: str,
    regimen_class: str,
    baseline_prognostic: Optional[Dict[str, Any]] = None,
    config: Dict[str, Any]
) -> float:
    """
    Compute prognostic risk score P in [0, 1].
    
    Higher P = better prognosis (higher probability of 6-12 month PFS).
    
    Inputs:
    - PFI_days, PTPI_days, TFI_days
    - line_of_therapy
    - PFS_days (if available from prior regimen)
    - baseline_prognostic: stage, performance_status, etc.
    
    Implementation:
    - Use simple disease-specific logistic/Cox model
    - Predict 6-month PFS probability
    - Normalize to [0, 1] (P = predicted_probability)
    """
    # Extract timing features
    pfi_days = timing_features.get("PFI_days")
    pfi_category = timing_features.get("PFI_category")  # "resistant", "sensitive", "partially_sensitive"
    line_of_therapy = timing_features.get("line_of_therapy", 1)
    pfs_days = timing_features.get("PFS_days")  # From prior regimen
    
    # Get disease-specific prognostic model
    prognostic_model = config.get("prognostic_models", {}).get(disease_site, {}).get(regimen_class, {})
    
    # Simple logistic model coefficients (configurable)
    intercept = prognostic_model.get("intercept", 0.0)
    pfi_coef = prognostic_model.get("pfi_coef", 0.5)
    line_coef = prognostic_model.get("line_coef", -0.3)
    
    # Compute linear predictor
    logit = intercept
    
    # PFI contribution (higher PFI = better prognosis)
    if pfi_days is not None:
        # Normalize PFI to months, then scale
        pfi_months = pfi_days / 30.0
        logit += pfi_coef * min(pfi_months / 12.0, 1.0)  # Cap at 12 months
    elif pfi_category:
        # Categorical mapping
        if pfi_category == "sensitive":
            logit += pfi_coef * 1.0
        elif pfi_category == "partially_sensitive":
            logit += pfi_coef * 0.5
        else:  # resistant
            logit += pfi_coef * 0.0
    
    # Treatment line penalty (later lines = worse prognosis)
    if line_of_therapy > 1:
        logit += line_coef * min((line_of_therapy - 1) / 3.0, 1.0)
    
    # Baseline prognostic factors (stage, PS)
    if baseline_prognostic:
        stage = baseline_prognostic.get("stage")
        if stage and "IV" in str(stage):
            logit -= 0.2  # Stage IV penalty
    
    # Convert logit to probability (sigmoid)
    import numpy as np
    pfs_probability = 1.0 / (1.0 + np.exp(-logit))
    
    # Normalize to [0, 1] (P score)
    return float(pfs_probability)
```

**Input Data Flow**:
```
Timing Engine → build_timing_chemo_features() → {
    "PFI_days": int,
    "PFI_category": "resistant" | "sensitive" | "partially_sensitive",
    "PTPI_days": int,
    "TFI_days": int,
    "PFS_days": int,  # Per-regimen
    "OS_days": int,   # Per-regimen
    "line_of_therapy": int,
    "regimen_id": str
}
```

**Status**: ✅ **READY** - Timing engine provides PFI/PTPI/TFI/PFS/OS per regimen

---

### **Component M: Predictive Mechanism Fit**

**Goal**: "How well does the tumor's biology align with this drug's mechanism?"

**Existing Engine**: ✅ **Holistic Score Service (Mechanism Fit Component)**

**What We Have**:
- `compute_mechanism_fit()` → Returns `mechanism_fit_score` (0.0-1.0), `mechanism_alignment` (per-pathway breakdown)
- 7D mechanism vectors: [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
- Cosine similarity between patient vector and drug/regimen vector
- **Validated**: TOPACIO study shows mechanism fit is the primary driver (50% weight)

**How to Map**:
```python
def compute_predictive_mechanism_fit(
    patient_mechanism_vector: List[float],
    regimen_moa_vector: List[float],
    ddr_features: Dict[str, Any],
    disease_site: str,
    regimen_class: str,
    config: Dict[str, Any]
) -> float:
    """
    Compute predictive mechanism fit score M in [0, 1].
    
    Uses cosine similarity of 7D mechanism vectors (validated approach).
    Supports disease-specific calibration curves.
    """
    from api.services.holistic_score.mechanism_fit import compute_mechanism_fit
    
    # Build patient profile and regimen dict
    patient_profile = {"mechanism_vector": patient_mechanism_vector}
    regimen = {"moa_vector": regimen_moa_vector}
    
    # Compute mechanism fit (same as holistic score service)
    mechanism_fit_score, mechanism_alignment = compute_mechanism_fit(
        patient_profile, regimen
    )
    
    if mechanism_fit_score is None:
        return 0.5  # Default neutral
    
    # Disease-specific calibration (optional)
    calibration = config.get("mechanism_calibration", {}).get(disease_site, {}).get(regimen_class, {})
    
    if calibration:
        # Apply logistic transform if specified
        # Example: M = sigmoid(alpha * (raw_score - beta))
        alpha = calibration.get("alpha", 1.0)
        beta = calibration.get("beta", 0.5)
        calibrated_score = 1.0 / (1.0 + np.exp(-alpha * (mechanism_fit_score - beta)))
        return float(calibrated_score)
    
    return float(mechanism_fit_score)
```

**Input Data Flow**:
```
Holistic Score Service → compute_mechanism_fit() → {
    "mechanism_fit_score": 0.0-1.0,
    "mechanism_alignment": {
        "DDR": 0.792,
        "MAPK": 0.012,
        "PI3K": 0.015,
        ...
    }
}

Patient Mechanism Vector: From tumor_context OR SAE features OR defaults
Regimen MoA Vector: Pre-tagged per regimen (e.g., carboplatin → [0.85, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
```

**Status**: ✅ **READY** - Mechanism fit component already exists and validated

---

### **Component T: Therapeutic Response Dynamics**

**Goal**: "Is the current regimen working now?"

**Existing Engine**: ✅ **Kinetic Biomarker Engine (KELIM/CA-125/PSA)**

**What We Have**:
- `CA125KELIMOvarian.compute_k_value()` → Returns `k_value`, `category` (favorable/unfavorable)
- `build_timing_chemo_features()` → Can compute KELIM on-the-fly from raw CA-125 measurements
- Early response indicators (cycle 3 normalization, % decline)

**How to Map**:
```python
def compute_therapeutic_dynamics(
    kinetic_features: Optional[Dict[str, Any]],
    regimen_start_date: datetime,
    current_date: datetime,
    disease_site: str,
    regimen_class: str,
    config: Dict[str, Any]
) -> float:
    """
    Compute therapeutic dynamics score T in [0, 1].
    
    Higher T = regimen is working (favorable kinetics).
    
    Inputs:
    - kelim_k_value: KELIM elimination rate constant
    - kelim_category: "favorable" | "unfavorable" | "intermediate"
    - ca125_early_decline_pct: % decline by cycle 3
    - ca125_normalized_by_cycle3: bool
    
    Implementation:
    - Early in regimen (few data points): T = 0.5 (neutral prior)
    - After enough data: Map kinetic metrics to T via disease-specific rules
    """
    # Check if we have enough data
    days_on_regimen = (current_date - regimen_start_date).days
    min_days = config.get("therapeutic_dynamics", {}).get(disease_site, {}).get("min_days_for_assessment", 90)
    
    if not kinetic_features or days_on_regimen < min_days:
        # Too early - use neutral prior
        return 0.5
    
    # Extract kinetic metrics
    kelim_k_value = kinetic_features.get("kelim_k_value")
    kelim_category = kinetic_features.get("kelim_category")
    ca125_decline_pct = kinetic_features.get("ca125_early_decline_pct")
    ca125_normalized = kinetic_features.get("ca125_normalized_by_cycle3", False)
    
    # Disease-specific rules
    rules = config.get("therapeutic_dynamics", {}).get(disease_site, {}).get(regimen_class, {})
    
    score = 0.5  # Start neutral
    
    # KELIM contribution (ovarian cancer)
    if disease_site == "ovary" and kelim_k_value is not None:
        if kelim_category == "favorable":
            score = 0.9  # High confidence favorable
        elif kelim_category == "unfavorable":
            score = 0.1  # High confidence unfavorable
        elif kelim_category == "intermediate":
            score = 0.5  # Neutral
    
    # CA-125 early decline (ovarian cancer)
    if disease_site == "ovary" and ca125_decline_pct is not None:
        if ca125_decline_pct >= 70:
            score = max(score, 0.9)  # Excellent response
        elif ca125_decline_pct >= 50:
            score = max(score, 0.7)  # Good response
        elif ca125_decline_pct >= 30:
            score = max(score, 0.5)  # Moderate response
        else:
            score = min(score, 0.3)  # Poor response
    
    # Cycle 3 normalization bonus
    if ca125_normalized:
        score = min(1.0, score + 0.1)  # Boost for early normalization
    
    return float(max(0.0, min(1.0, score)))
```

**Input Data Flow**:
```
Kinetic Engine → CA125KELIMOvarian.compute_k_value() → {
    "kelim_k_value": float,  # e.g., 1.2 (favorable), 0.3 (unfavorable)
    "kelim_category": "favorable" | "unfavorable" | "intermediate",
    "measurements_used": int,
    "warnings": List[str],
    "error": Optional[str]
}

Timing Engine → Can compute on-the-fly from raw CA-125 measurements
```

**Status**: ✅ **READY** - Kinetic engine provides KELIM scores and early response indicators

---

### **Component S: Safety/Tolerability**

**Goal**: "Can this patient safely receive this regimen?"

**Existing Engine**: ✅ **Holistic Score Service (PGx Safety Component)**

**What We Have**:
- `compute_pgx_safety()` → Returns `pgx_safety_score` (0.0-1.0), `pgx_details` (contraindications, dose adjustments)
- PGx screening service integration (DPYD, TPMT, UGT1A1, CYP2D6, etc.)
- Contraindication detection (threshold ≤0.1)
- Dose adjustment recommendations

**How to Map**:
```python
async def compute_safety_tolerability(
    safety_features: Dict[str, Any],
    regimen_drugs: List[str],
    baseline_organ_function: Optional[Dict[str, Any]] = None,
    prior_toxicity: Optional[Dict[str, Any]] = None,
    config: Dict[str, Any]
) -> float:
    """
    Compute safety/tolerability score S in [0, 1].
    
    Higher S = safer (no contraindications, no dose adjustments needed).
    
    Start at 1.0, subtract risk contributions:
    - Serious PGx contraindication: -1.0 (S = 0.0)
    - Moderate PGx concern (dose adjustment): -0.3
    - Organ function issues: -0.2
    - Prior severe toxicity: -0.2
    """
    from api.services.holistic_score.pgx_safety import compute_pgx_safety
    from api.services.pgx_screening_service import get_pgx_screening_service
    
    score = 1.0
    pgx_service = get_pgx_screening_service()
    
    # PGx safety (per drug in regimen)
    pgx_scores = []
    for drug in regimen_drugs:
        pharmacogenes = safety_features.get("germline_variants", [])
        
        # Use existing PGx safety computation
        pgx_score, pgx_details = await compute_pgx_safety(
            pharmacogenes, drug, {"interventions": [{"drug_names": [drug]}]}, pgx_service
        )
        
        pgx_scores.append(pgx_score)
        
        # Hard contraindication → S = 0.0
        if pgx_details.get("contraindicated"):
            return 0.0
    
    # Aggregate PGx scores (worst drug determines overall)
    min_pgx_score = min(pgx_scores) if pgx_scores else 1.0
    score = min(score, min_pgx_score)
    
    # Organ function penalties
    if baseline_organ_function:
        if baseline_organ_function.get("creatinine_clearance", 100) < 30:
            score -= 0.2  # Renal impairment
        if baseline_organ_function.get("bilirubin", 0) > 2.0:
            score -= 0.2  # Hepatic impairment
    
    # Prior toxicity penalties
    if prior_toxicity:
        regimen_class = safety_features.get("regimen_class", "")
        if prior_toxicity.get("severe_toxicity_with_same_class"):
            score -= 0.2
    
    return float(max(0.0, min(1.0, score)))
```

**Input Data Flow**:
```
PGx Safety Component → compute_pgx_safety() → {
    "pgx_safety_score": 0.0-1.0,
    "pgx_details": {
        "drug": str,
        "contraindicated": bool,
        "dose_adjustments": List[str],
        "variants_screened": List[Dict],
        "toxicity_tier": "LOW" | "MODERATE" | "HIGH"
    }
}
```

**Status**: ✅ **READY** - PGx safety component already exists

---

## 🏗️ OVERALL HOLISTIC CLINICAL BENEFIT SCORE

### **Formula**

```python
Overall = w_D * D + w_P * P + w_M * M + w_T * T + w_S * S
```

### **Use-Case-Specific Weights** ✅ **MANAGER-APPROVED**

**Weights depend on use-case and disease** (configurable, with defaults):

**Manager Guidance (January 29, 2025):**
- Keep them simple and override-able per disease/regimen later
- All weights must live in config per `disease_site + regimen_class + use_case`
- If a component is missing (e.g., no KELIM), renormalize remaining weights

```python
# Trial enrollment (patient-trial pair)
weights_trial_enrollment = {
    "D": 0.20,  # Diagnostic fit important
    "P": 0.10,  # Prognosis less critical (trial is exploratory)
    "M": 0.45,  # Mechanism fit is KEY (validated in TOPACIO) - increased from 0.40
    "T": 0.00,  # Therapeutic dynamics not relevant (pre-treatment)
    "S": 0.25   # Safety is critical - increased from 0.20
}

# Next-line regimen selection (patient-regimen pair)
weights_next_line = {
    "D": 0.10,  # Diagnostic fit assumed
    "P": 0.25,  # Prognosis important (PFI, line of therapy)
    "M": 0.35,  # Mechanism fit is KEY
    "T": 0.20,  # Therapeutic dynamics from prior regimens (neutral 0.5 if no kinetics yet)
    "S": 0.10   # Safety important but less critical
}

# On-treatment monitoring / continuation decision
weights_monitoring = {
    "D": 0.05,  # Diagnostic fit not changing
    "P": 0.20,  # Prognosis based on current response - reduced from 0.25
    "M": 0.15,  # Mechanism fit already known
    "T": 0.45,  # Therapeutic dynamics is PRIMARY (DRI/KELIM/early response)
    "S": 0.15   # Safety monitoring - increased from 0.10
}
```

**Normalization Rules:**
- If T is missing (no kinetics): Set T=0.5 (neutral), renormalize other weights
- If any component is unavailable: Set to neutral value, renormalize remaining weights
- Always ensure weights sum to 1.0

---

## 🔌 INTEGRATION INTERFACE

### **Main Function** ✅ **MANAGER-APPROVED INTERFACE**

**Manager Guidance:** Expose one endpoint/service: `compute_holistic_score(patient_id, regimen_id, use_case)`

```python
async def compute_holistic_clinical_benefit_score(
    patient_id: str,
    regimen_id: str,
    use_case: str,  # "trial_enrollment" | "next_line" | "monitoring"
    ddr_features: Optional[Dict[str, Any]] = None,        # from DDR_bin engine (pulled automatically)
    timing_features: Optional[Dict[str, Any]] = None,     # from timing engine (pulled automatically)
    kinetic_features: Optional[Dict[str, Any]] = None,    # from KELIM/CA-125/PSA engine (pulled automatically)
    safety_features: Optional[Dict[str, Any]] = None,     # from PGx engine (pulled automatically)
    clinical_features: Optional[Dict[str, Any]] = None,   # baseline clinicopathologic features (pulled automatically)
    patient_mechanism_vector: Optional[List[float]] = None,  # 7D mechanism vector (pulled automatically)
    regimen_moa_vector: Optional[List[float]] = None,        # 7D MoA vector (pulled automatically)
    config: Optional[Dict[str, Any]] = None               # scoring weights & disease-specific params (optional override)
) -> Dict[str, Any]:
    """
    Compute Holistic Clinical Benefit Score for patient-regimen pair.
    
    Manager Guidance: This is CSI-plus orchestration layer.
    - CSI is the Predictive core (mostly M + part of P)
    - Holistic Score is the clinically-facing layer that exposes D/P/M/T/S explicitly
    
    All scores must be decomposable (always return D,P,M,T,S separately) so they map
    cleanly to clinical-benefit categories (Diagnostic, Prognostic, Predictive, 
    Therapeutic, Safety, Monitoring) in docs and UI.
    
    Returns:
    {
        "D": float,      # Diagnostic fit
        "P": float,      # Prognostic risk
        "M": float,      # Predictive mechanism fit
        "T": float,      # Therapeutic response dynamics
        "S": float,      # Safety/tolerability
        "overall": float,  # Overall score
        "weights": Dict[str, float],  # Weights used (normalized)
        "breakdown": Dict[str, Any],  # Detailed breakdown per component
        "interpretation": str,        # "HIGH" | "MEDIUM" | "LOW"
        "recommendation": str,        # Human-readable recommendation
        "component_available": Dict[str, bool],  # Which components were available
        "provenance": Dict[str, Any]
    }
    """
    # Pull features from engines if not provided (automatic fetch)
    if not ddr_features:
        ddr_features = await _fetch_ddr_features(patient_id)
    if not timing_features:
        timing_features = await _fetch_timing_features(patient_id, regimen_id)
    if not kinetic_features:
        kinetic_features = await _fetch_kinetic_features(patient_id, regimen_id)
    if not safety_features:
        safety_features = await _fetch_safety_features(patient_id)
    if not clinical_features:
        clinical_features = await _fetch_clinical_features(patient_id)
    if not patient_mechanism_vector:
        patient_mechanism_vector = await _fetch_mechanism_vector(patient_id)
    if not regimen_moa_vector:
        regimen_moa_vector = await _fetch_regimen_moa_vector(regimen_id)
    
    disease_site = clinical_features.get("disease_site")
    regimen_class = clinical_features.get("regimen_class", "unknown")
    
    # Get use-case-specific weights from config (per disease_site + regimen_class + use_case)
    weights_config = config or {}
    weights = weights_config.get("weights", {}).get(use_case, {}).get(disease_site, {}).get(regimen_class)
    
    # Fallback to use-case defaults if no disease/regimen-specific weights
    if not weights:
        if use_case == "trial_enrollment":
            weights = {"D": 0.20, "P": 0.10, "M": 0.45, "T": 0.00, "S": 0.25}
        elif use_case == "next_line":
            weights = {"D": 0.10, "P": 0.25, "M": 0.35, "T": 0.20, "S": 0.10}
        elif use_case == "monitoring":
            weights = {"D": 0.05, "P": 0.20, "M": 0.15, "T": 0.45, "S": 0.15}
        else:
            weights = {"D": 0.15, "P": 0.25, "M": 0.30, "T": 0.20, "S": 0.10}  # Default
    
    # Check component availability and normalize weights if needed
    component_available = {
        "D": ddr_features is not None,
        "P": timing_features is not None,
        "M": patient_mechanism_vector is not None and regimen_moa_vector is not None,
        "T": kinetic_features is not None,
        "S": safety_features is not None
    }
    
    # Normalize weights if components are missing (set missing components to neutral, renormalize)
    normalized_weights = _normalize_weights_for_missing_components(weights, component_available)
    
    # Compute each component (use neutral value if missing)
    D = compute_diagnostic_fit(ddr_features, disease_site, clinical_features, config) if component_available["D"] else 0.5
    P = compute_prognostic_risk(timing_features, disease_site, regimen_class, clinical_features, config) if component_available["P"] else 0.5
    M = compute_predictive_mechanism_fit(patient_mechanism_vector, regimen_moa_vector, ddr_features, disease_site, regimen_class, config) if component_available["M"] else 0.5
    T = compute_therapeutic_dynamics(kinetic_features, clinical_features.get("regimen_start_date"), datetime.now(), disease_site, regimen_class, config) if component_available["T"] else 0.5
    S = await compute_safety_tolerability(safety_features, clinical_features.get("regimen_drugs", []), clinical_features.get("baseline_organ_function"), clinical_features.get("prior_toxicity"), config) if component_available["S"] else 0.5
    
    # Compute overall score using normalized weights
    overall = (
        normalized_weights["D"] * D +
        normalized_weights["P"] * P +
        normalized_weights["M"] * M +
        normalized_weights["T"] * T +
        normalized_weights["S"] * S
    )
    
    # Generate interpretation
    interpretation, recommendation = interpret_holistic_clinical_benefit_score(
        overall, D, P, M, T, S, use_case
    )
    
    return {
        "D": round(D, 3),
        "P": round(P, 3),
        "M": round(M, 3),
        "T": round(T, 3),
        "S": round(S, 3),
        "overall": round(overall, 3),
        "weights": normalized_weights,  # Return normalized weights (sum to 1.0)
        "component_available": component_available,  # Which components were available
        "breakdown": {
            "diagnostic_fit": {"score": D, "details": ddr_features, "available": component_available["D"]},
            "prognostic_risk": {"score": P, "details": timing_features, "available": component_available["P"]},
            "mechanism_fit": {"score": M, "details": {"mechanism_alignment": M}, "available": component_available["M"]},
            "therapeutic_dynamics": {"score": T, "details": kinetic_features, "available": component_available["T"]},
            "safety": {"score": S, "details": safety_features, "available": component_available["S"]}
        },
        "interpretation": interpretation,
        "recommendation": recommendation,
        "provenance": {
            "service": "HolisticClinicalBenefitScore",
            "version": "1.0",
            "use_case": use_case,
            "disease_site": disease_site,
            "regimen_class": regimen_class,
            "formula": f"{normalized_weights['D']}×D + {normalized_weights['P']}×P + {normalized_weights['M']}×M + {normalized_weights['T']}×T + {normalized_weights['S']}×S",
            "ruo": "Research Use Only",
            "note": "CSI-plus orchestration: CSI is Predictive core (M+part of P), Holistic Score exposes D/P/M/T/S explicitly"
        }
    }


def _normalize_weights_for_missing_components(
    weights: Dict[str, float],
    component_available: Dict[str, bool]
) -> Dict[str, float]:
    """
    Normalize weights if components are missing.
    
    If a component is unavailable, set its weight to 0 and renormalize remaining weights.
    """
    # Set weight to 0 for missing components
    adjusted_weights = {
        k: weights[k] if component_available[k] else 0.0
        for k in ["D", "P", "M", "T", "S"]
    }
    
    # Renormalize to sum to 1.0
    total = sum(adjusted_weights.values())
    if total > 0:
        normalized = {k: v / total for k, v in adjusted_weights.items()}
    else:
        # Fallback: equal weights if all missing (shouldn't happen)
        normalized = {k: 0.2 for k in ["D", "P", "M", "T", "S"]}
    
    return normalized
```

### **Batch Version**

```python
async def compute_holistic_clinical_benefit_scores_batch(
    patient_id: str,
    regimens: List[Dict[str, Any]],  # List of regimen dicts with regimen_id
    # ... same inputs as above, but regimens is a list
) -> List[Dict[str, Any]]:
    """Compute holistic scores for multiple regimens for same patient."""
    results = []
    for regimen in regimens:
        score = await compute_holistic_clinical_benefit_score(
            patient_id=patient_id,
            regimen_id=regimen["regimen_id"],
            ddr_features=ddr_features,  # Same for all regimens
            timing_features=get_timing_features_for_regimen(regimen["regimen_id"]),
            kinetic_features=get_kinetic_features_for_regimen(regimen["regimen_id"]),
            # ... etc
        )
        results.append(score)
    
    # Sort by overall score (descending)
    results.sort(key=lambda x: x["overall"], reverse=True)
    return results
```

---

## 📋 IMPLEMENTATION ROADMAP ✅ **MANAGER-GUIDED**

### **Priority Guidance (Manager - January 29, 2025):**

**P0 Integrator Task** - Execute after engines are stable:
1. ✅ Finish: DDR_bin, Timing (PFI/PTPI/TFI), CSI v0, KELIM/CA-125 MVP
2. ✅ Then execute this Holistic Score wrapper in the same sprint (5-6 days acceptable)

**Guardrails:**
- Do not block engine hardening
- Should be the first orchestration layer once CSI v0 is in place
- Do not build new ML inside this wrapper; just apply calibrated linear/logistic mappings from engines
- All scores must be decomposable (always return D,P,M,T,S separately)

### **Phase 1: Component Functions** (2-3 days)

1. ✅ **D: Diagnostic Fit** → Implement `compute_diagnostic_fit()`
   - Use DDR_bin engine output
   - Add penalties for missing/unknown data
   - Test with synthetic examples

2. ✅ **P: Prognostic Risk** → Implement `compute_prognostic_risk()`
   - Use Timing Engine output (PFI, PFS, line_of_therapy)
   - Add simple logistic model
   - Test with known PFI categories

3. ✅ **M: Predictive Mechanism Fit** → **ALREADY EXISTS**
   - Use `compute_mechanism_fit()` from holistic_score service
   - Add disease-specific calibration (optional)

4. ✅ **T: Therapeutic Dynamics** → Implement `compute_therapeutic_dynamics()`
   - Use Kinetic Engine output (KELIM, CA-125 decline)
   - Add neutral prior for early regimens
   - Test with favorable/unfavorable KELIM

5. ✅ **S: Safety/Tolerability** → **ALREADY EXISTS**
   - Use `compute_pgx_safety()` from holistic_score service
   - Add organ function and prior toxicity logic

### **Phase 2: Orchestration** (1-2 days)

1. ✅ Implement `compute_holistic_clinical_benefit_score()`
2. ✅ Implement batch version
3. ✅ Add interpretation function
4. ✅ Create configuration system (weights, thresholds, models)

### **Phase 3: Testing** (1-2 days)

1. ✅ Unit tests for each component
2. ✅ Integration tests with real engine outputs
3. ✅ Test with missing data (kinetic_features empty, etc.)
4. ✅ Test different use-cases (trial enrollment, next line, monitoring)

### **Phase 4: Integration** (1-2 days) ✅ **MANAGER-GUIDED**

**Integration Point (Manager Guidance):**

1. ✅ **Standalone Service First**: 
   - Expose one endpoint/service: `compute_holistic_score(patient_id, regimen_id, use_case)`
   - Pulls D/P/M/T/S from existing engines automatically
   - Returns: overall score + the five sub-scores + feature contributions

2. ✅ **Then Integrate into Care Plan**:
   - Have `complete_care` call it as a sub-step when:
     - (a) proposing trials (use_case="trial_enrollment")
     - (b) ranking next systemic regimens (use_case="next_line")
     - (c) evaluating whether to continue/stop a regimen (use_case="monitoring")

3. ✅ Add API endpoint: `POST /api/resistance/holistic-clinical-benefit`

4. ✅ Add frontend display components

**Total Estimated Effort**: 5-6 days (same sprint after engines are stable)

---

## 🔍 KEY INSIGHTS FROM HOLISTIC SCORE VALIDATION

### **What Worked in Trial Matching (TOPACIO)**

1. ✅ **Mechanism Fit is Primary Driver** (50% weight)
   - AUROC=0.714, p=0.023
   - Mechanism fit tracked with genomic features (BRCA-mutant 0.849 vs HRD-negative 0.579)

2. ✅ **Simple Transparent Models Work**
   - Cosine similarity (not complex ML)
   - Linear weighted combination (not deep learning)
   - Configurable weights (not hard-coded)

3. ✅ **Component Breakdown is Critical**
   - Shows which dimension drives the score
   - Enables explainability
   - Helps clinicians understand recommendations

### **How to Apply to Clinical Benefit Score**

1. **M (Mechanism Fit) Should Remain High Weight** (30-40%)
   - Proven to be predictive
   - Aligns with biomarker slide emphasis

2. **Use-Case-Specific Weights Are Essential**
   - Trial enrollment: M+D+S emphasis
   - Next line: M+P+T emphasis
   - Monitoring: T+P emphasis

3. **Transparency Over Complexity**
   - Simple logistic models for P
   - Rule-based for D and T
   - Reuse validated components (M, S)

---

## ✅ READINESS ASSESSMENT

| Component | Engine Status | Integration Status | Effort |
|-----------|---------------|-------------------|--------|
| **D: Diagnostic Fit** | ✅ DDR_bin ready | ⚠️ Needs wrapper function | 4 hours |
| **P: Prognostic Risk** | ✅ Timing engine ready | ⚠️ Needs prognostic model | 8 hours |
| **M: Mechanism Fit** | ✅ Complete | ✅ Ready to use | 0 hours |
| **T: Therapeutic Dynamics** | ✅ Kinetic engine ready | ⚠️ Needs mapping function | 6 hours |
| **S: Safety/Tolerability** | ✅ Complete | ✅ Ready to use | 0 hours |
| **Orchestration** | ⚠️ Not started | ⚠️ Needs implementation | 16 hours |
| **Testing** | ⚠️ Not started | ⚠️ Needs unit/integration tests | 8 hours |

**Total**: ~42 hours (~5-6 days)

---

## 🎯 NEXT STEPS ✅ **MANAGER-APPROVED**

### **Prerequisites (Must Complete First):**
1. ✅ DDR_bin engine stable
2. ✅ Timing engine stable (PFI/PTPI/TFI)
3. ✅ CSI v0 in place
4. ✅ KELIM/CA-125 MVP complete

### **Implementation (5-6 days in same sprint):**

1. **Create new service**: `api/services/holistic_clinical_benefit/`
   - Structure similar to `holistic_score` service
   - Component functions: `diagnostic_fit.py`, `prognostic_risk.py`, `therapeutic_dynamics.py`
   - Main service: `service.py` with `compute_holistic_score(patient_id, regimen_id, use_case)`

2. **Reuse existing components**:
   - Import `compute_mechanism_fit()` from `holistic_score.mechanism_fit`
   - Import `compute_pgx_safety()` from `holistic_score.pgx_safety`

3. **Create configuration system**:
   - `config/holistic_clinical_benefit_config.py`
   - Weights per `disease_site + regimen_class + use_case` (with defaults)
   - Disease-specific thresholds, models

4. **Add API endpoint**:
   - `POST /api/resistance/holistic-clinical-benefit`
   - Request: `patient_id`, `regimen_id`, `use_case` (optionally override features/config)
   - Response: `D, P, M, T, S, overall, weights, component_available, breakdown, interpretation`

5. **Integrate into care plan**:
   - Make it a standalone service first
   - Then have `complete_care` call it as sub-step for:
     - (a) proposing trials → `use_case="trial_enrollment"`
     - (b) ranking next regimens → `use_case="next_line"`
     - (c) monitoring decisions → `use_case="monitoring"`

6. **Display in UI**:
   - Map D/P/M/T/S to clinical-benefit categories in docs and UI
   - Show component availability flags
   - Display normalized weights used

---

## 📝 MANAGER GUIDANCE SUMMARY

**Key Points:**
1. ✅ **Priority**: P0 integrator task after engines stable (don't block engine hardening)
2. ✅ **Weights**: Simple defaults, override-able per disease/regimen, renormalize if components missing
3. ✅ **Integration**: Standalone service first, then integrate into `complete_care` as sub-step
4. ✅ **Guardrails**: No new ML, just orchestration layer, always return D/P/M/T/S separately
5. ✅ **Relationship to CSI**: CSI is Predictive core (M+part of P), Holistic Score is clinically-facing layer

---

**Last Updated**: January 29, 2025 (Manager Guidance Incorporated)  
**Status**: ✅ **APPROVED FOR IMPLEMENTATION**  
**Priority**: **P0** - Integrator task after engines stable  
**Timeline**: 5-6 days in same sprint after DDR_bin, Timing, CSI v0, KELIM MVP complete
