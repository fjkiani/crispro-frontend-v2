# ⚔️ SAE-WIWFM INTEGRATION REVIEW - COMPLETE ANALYSIS

**Date**: January 13, 2025  
**Status**: 🔍 **COMPREHENSIVE REVIEW COMPLETE**  
**Purpose**: Review current SAE integration with WIWFM (Will It Work For Me) drug efficacy prediction

---

## 📋 EXECUTIVE SUMMARY

### **Current State: SAE is "Display Only"**

**Finding**: SAE features are computed but **NOT used to modulate drug efficacy scores or confidence** in WIWFM.

**Architecture Gap**:
- ✅ SAE features computed in `sae_feature_service.py` (Phase 2)
- ✅ SAE features extracted in `efficacy_orchestrator.py` (optional, line 330-380)
- ❌ **SAE features NOT used in `drug_scorer.py`** (no confidence modulation)
- ❌ **No SAE lifts/penalties applied to drug scores**
- ❌ **SAE lives in Ayesha orchestrator, not integrated into S/P/E pipeline**

**Manager's Vision** (from `.cursorrules`):
> "SAE must live inside S/P/E and modulate confidence"  
> "S/P/E base ± SAE lifts/penalties"

**Current Reality**:
> SAE is isolated in Ayesha orchestrator (display only, no drug influence)

---

## 🔍 DETAILED ANALYSIS

### **1. CURRENT IMPLEMENTATION STATE**

#### **A. Efficacy Orchestrator (`orchestrator.py`)**

**Location**: `oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/orchestrator.py`

**Current Behavior** (Lines 330-380):
```python
# 8) Extract SAE features if requested (P2 feature)
if (request.options or {}).get("include_sae_features"):
    try:
        from api.services.sae_service import extract_sae_features_from_real_data
        
        # Collect all real data sources for SAE
        sae_bundle = extract_sae_features_from_real_data(...)
        
        # Add SAE features to response
        response.sae_features = sae_features_to_dict(sae_bundle)
        
        # Add SAE attribution to confidence breakdown
        if "confidence_breakdown" in response.provenance:
            response.provenance["confidence_breakdown"]["sae_attribution"] = {
                "boosting_features": sae_bundle.boosting_features,
                "limiting_features": sae_bundle.limiting_features,
                "overall_impact": sae_bundle.overall_impact
            }
        
        response.provenance["sae_enabled"] = True
    except Exception as e:
        response.provenance["sae_error"] = str(e)
        response.provenance["sae_enabled"] = False
```

**Analysis**:
- ✅ SAE features are extracted and added to response
- ✅ SAE attribution is tracked in provenance
- ❌ **SAE features are NOT passed to `drug_scorer.py`**
- ❌ **No confidence modulation based on SAE**
- ❌ **No drug score adjustments based on SAE**

**Gap**: SAE is computed but never used in drug scoring logic.

---

#### **B. Drug Scorer (`drug_scorer.py`)**

**Location**: `oncology-coPilot/oncology-backend-minimal/api/services/efficacy_orchestrator/drug_scorer.py`

**Current Behavior** (Lines 24-208):
```python
async def score_drug(
    self,
    drug: Dict[str, Any],
    seq_scores: List[SeqScore],
    pathway_scores: Dict[str, float],
    evidence_result: Any,
    clinvar_result: Any,
    insights: InsightsBundle,
    confidence_config,
    disease: str = "",
    include_fda_badges: bool = False,
) -> DrugScoreResult:
    # ... S/P/E scoring ...
    
    # Base confidence
    confidence = compute_confidence(tier, seq_pct, path_pct, insights_dict, confidence_config)
    
    # ClinVar-based confidence bump
    if clinvar_prior > 0 and path_pct >= 0.2:
        confidence += min(0.1, clinvar_prior)
    
    # Deterministic gene-drug MoA tie-breaker
    if primary_gene == "BRAF" and drug_name == "BRAF inhibitor":
        confidence += 0.01
    
    # Efficacy score (likelihood of benefit)
    raw_lob = 0.3 * seq_pct + 0.4 * path_pct + 0.3 * s_evd + clinvar_prior
```

**Analysis**:
- ✅ S/P/E framework implemented (30/40/30 weighting)
- ✅ Confidence computed from tier, sequence, pathway, evidence
- ✅ ClinVar prior boosts confidence
- ✅ Deterministic gene-drug MoA tie-breakers
- ❌ **NO SAE parameter in `score_drug()` method signature**
- ❌ **NO SAE-based confidence modulation**
- ❌ **NO SAE lifts/penalties applied**

**Gap**: Drug scorer has no knowledge of SAE features.

---

#### **C. SAE Feature Service (`sae_feature_service.py`)**

**Location**: `oncology-coPilot/oncology-backend-minimal/api/services/sae_feature_service.py`

**Current Behavior** (Lines 123-448):
```python
def compute_sae_features(
    self,
    insights_bundle: Dict[str, Any],
    pathway_scores: Dict[str, float],
    tumor_context: Dict[str, Any],
    treatment_history: Optional[List[Dict]] = None,
    ca125_intelligence: Optional[Dict] = None,
    previous_hrd_score: Optional[float] = None,
    previous_dna_repair_capacity: Optional[float] = None
) -> SAEFeatures:
    """
    Compute SAE Features (Post-NGS Only)
    
    Returns:
        SAEFeatures dataclass with:
        - dna_repair_capacity: float  # 0-1, Manager's formula (C5)
        - pathway_burden_ddr: float   # 0-1
        - pathway_burden_mapk: float  # 0-1
        - pathway_burden_pi3k: float  # 0-1
        - pathway_burden_vegf: float  # 0-1
        - pathway_burden_her2: float  # 0-1
        - io_eligible: bool           # TMB >= 20 OR MSI-High
        - mechanism_vector: List[float]  # [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
        - resistance_signals: Dict[str, Any]
        - hotspot_mutation: bool
        - ...
    """
```

**Analysis**:
- ✅ SAE features computed correctly (Manager's C1-C10 policy)
- ✅ DNA repair capacity formula implemented (0.6/0.2/0.2 weights)
- ✅ Mechanism vector computed (7D for cosine similarity)
- ✅ Hotspot detection implemented (P0 Fix #3)
- ✅ Resistance signals detected (2-of-3 trigger logic)
- ❌ **SAE features NOT passed to drug scoring**
- ❌ **No integration with `drug_scorer.py`**

**Gap**: SAE service computes features but they're not consumed by drug scoring.

---

#### **D. Ayesha Orchestrator (`ayesha_orchestrator_v2.py`)**

**Location**: `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_orchestrator_v2.py`

**Current Behavior** (Lines 440-513):
```python
# Get SAE features from WIWFM if available
sae_features = None
if results["wiwfm"] and isinstance(results["wiwfm"], dict):
    sae_features = results["wiwfm"].get("sae_features")

# Phase 2 SAE Services
if request.tumor_context:
    results["sae_features"] = compute_sae_features(...)
    results["resistance_alert"] = detect_resistance(...)
```

**Analysis**:
- ✅ SAE features extracted from WIWFM response (if available)
- ✅ SAE features computed independently in Ayesha orchestrator
- ✅ SAE features displayed in UI (Mechanism Map, Hint Tiles)
- ❌ **SAE features NOT fed back into WIWFM drug scoring**
- ❌ **SAE is "display only", not affecting drug recommendations**

**Gap**: SAE lives in Ayesha orchestrator but doesn't influence drug efficacy.

---

### **2. MANAGER'S POLICY REVIEW**

#### **A. Previous Guidance (from `.cursorrules` lines 1426-1456)**

**Manager's Previous Guidance (Q2c, lines 332-354)**:
> **"Short‑Term Role (Until Validation + Refactor):**  
> Keep SAE **display + Resistance Playbook only** (no direct changes to drug scores in `/api/efficacy/predict`).  
> No SAE‑based lifts/caps in WIWFM until:
> 1. HRD/platinum validation is running, and  
> 2. We have a written SAE policy for lifts/gates"

**Status**: ✅ **CURRENT STATE ALIGNS WITH THIS GUIDANCE**
- SAE is display only
- No SAE hooks in `/api/efficacy/predict` drug scoring
- SAE powers Resistance Playbook

---

#### **B. Manager's Vision (from Gap #30 audit)**

**Manager's Plan (19.1, 19.2)**:
> SAE is PART OF S/P/E pipeline: "S/P/E base ± SAE lifts/penalties"  
> Data flow: Evo2 → Insights → Pathway → **SAE extractor** → Confidence computation  
> SAE should modulate S/P/E scores, not replace them

**Current Architecture**:
```
Ayesha Orchestrator:
  - Trials (separate)
  - SOC (separate)
  - CA-125 (separate)
  - SAE Features (separate)  ← ISOLATED, NOT INTEGRATED

Should Be:
Drug Efficacy Router (WIWFM):
  S → P → E → SAE → Confidence → Final Score
```

**Gap**: Manager's vision not implemented yet (blocked by validation requirement).

---

#### **C. SAE Lift/Gate Policy (from `SAE_LIFT_GATE_POLICY_V1.md`)**

**Documented Policy** (Lines 253-313):
- ✅ Policy document exists (`.cursor/ayesha/SAE_LIFT_GATE_POLICY_V1.md`)
- ✅ Lift/penalty rules defined (PARP +0.10, MEK +0.15, Taxane -0.20)
- ✅ Safety guardrails documented (max lift +0.15, max penalty -0.20)
- ⚠️ **Status: AWAITING VALIDATION + MANAGER APPROVAL**

**Implementation Checklist** (Lines 265-280):
- [ ] Create feature flag `EFFICACY_ENABLE_SAE=true`
- [ ] Add SAE computation to `/api/efficacy/predict` (after pathway, before confidence)
- [ ] Import `sae_modulation_rules.py`
- [ ] Apply lifts/penalties using `modulate_drug_confidence()`
- [ ] Update `DrugEfficacyResponse` schema to include `sae_features`
- [ ] Add `confidence_breakdown` with SAE contribution
- [ ] Surface SAE reasoning in drug responses
- [ ] Update frontend to display SAE delta and reasoning
- [ ] Log all lifts/penalties to provenance
- [ ] Run E2E tests with Ayesha's profile
- [ ] Deploy behind feature flag

**Status**: ⚔️ **POLICY DOCUMENTED - AWAITING VALIDATION + MANAGER APPROVAL**

**Do NOT implement until**:
1. Validation running (≥200 TCGA patients, AUROC/AUPRC computed)
2. Manager explicitly approves implementation

---

### **3. INTEGRATION GAPS IDENTIFIED**

#### **Gap #1: SAE Not in Drug Scoring Pipeline** 🔥

**Current**: SAE features computed but never passed to `drug_scorer.py`

**Required**: 
- Add `sae_features: Optional[SAEFeatures]` parameter to `score_drug()`
- Pass SAE features from orchestrator to drug scorer
- Apply SAE-based confidence modulation

**Impact**: CRITICAL - SAE doesn't influence drug recommendations

---

#### **Gap #2: No SAE Confidence Modulation** 🔥

**Current**: Confidence computed from S/P/E only, no SAE contribution

**Required**:
- Implement `modulate_drug_confidence()` function
- Apply SAE lifts/penalties based on:
  - DNA repair capacity <0.40 → PARP +0.10
  - Hotspot mutation (KRAS/BRAF) → MEK/RAF +0.15
  - Cross-resistance risk → Taxane -0.20
- Track SAE contribution in `confidence_breakdown`

**Impact**: CRITICAL - Manager's vision not realized

---

#### **Gap #3: SAE Isolated in Ayesha Orchestrator** ⚠️

**Current**: SAE computed in Ayesha orchestrator, not in efficacy pipeline

**Required**:
- Move SAE computation into efficacy orchestrator
- Compute SAE alongside S/P/E (not separately)
- Feed SAE into confidence computation

**Impact**: MEDIUM - Architectural mismatch with Manager's vision

---

#### **Gap #4: No SAE Attribution in Drug Response** ⚠️

**Current**: SAE attribution tracked in provenance but not in drug response

**Required**:
- Add `sae_context` to each drug in `DrugEfficacyResponse`
- Show SAE lift/penalty in drug rationale
- Display SAE contribution in frontend

**Impact**: MEDIUM - Transparency issue

---

### **4. WHAT WORKS (CURRENT STATE)**

#### **✅ SAE Feature Computation**
- DNA repair capacity formula correct (0.6/0.2/0.2 weights)
- Mechanism vector computed (7D for cosine similarity)
- Hotspot detection working (P0 Fix #3)
- Resistance signals detected (2-of-3 trigger logic)

#### **✅ SAE Display in UI**
- Mechanism Map shows SAE pathway burdens (post-NGS)
- Hint Tiles use SAE features for recommendations
- Resistance Playbook uses SAE for resistance detection

#### **✅ SAE Attribution Tracking**
- SAE attribution tracked in provenance
- `sae_enabled` flag in response
- Error handling for SAE failures

---

### **5. RECOMMENDATIONS**

#### **Option A: Wait for Validation (RECOMMENDED - Aligns with Manager's Policy)**

**Status**: ⏸️ **BLOCKED - AWAITING VALIDATION**

**Requirements**:
1. ✅ HRD/platinum validation running (≥200 TCGA patients)
2. ✅ AUROC/AUPRC computed and reviewed
3. ✅ Manager explicitly approves implementation
4. ✅ Written SAE policy for lifts/gates finalized

**Action**: Do NOT implement SAE→WIWFM integration until validation complete.

**Timeline**: TBD (depends on Jr2 HRD extraction + validation results)

---

#### **Option B: Hybrid Integration (Phased Approach)**

**Phase 1 (Now)**: Keep current state (SAE display only)
- ✅ SAE in Ayesha orchestrator
- ✅ SAE powers Mechanism Map, Hint Tiles, Resistance Playbook
- ❌ No SAE in drug scoring

**Phase 2 (After Validation)**: Add SAE module to efficacy orchestrator
- Add SAE computation inside `/api/efficacy/predict` (gated by feature flag)
- Compute SAE alongside S/P/E
- Apply SAE lifts/penalties (behind feature flag)

**Phase 3 (Final)**: Make SAE-enhanced efficacy default
- Remove feature flag
- Keep "Baseline (no SAE)" profile for comparison

**Timeline**: Phase 2 starts after validation approval

---

#### **Option C: Full Integration Now (NOT RECOMMENDED - Violates Manager's Policy)**

**Action**: Implement SAE→WIWFM integration immediately

**Risk**: 
- ❌ Violates Manager's explicit guidance (Q2c)
- ❌ No validation data to support lifts/penalties
- ❌ Could break existing WIWFM functionality

**Status**: ❌ **NOT APPROVED - Manager explicitly said "wait for validation"**

---

### **6. NEXT STEPS**

#### **Immediate (Do Now)**:
1. ✅ **Document current state** (this review)
2. ✅ **Confirm alignment with Manager's policy** (current state is correct)
3. ⏸️ **Wait for validation** (Jr2 HRD extraction + validation script)

#### **After Validation**:
1. [ ] Run validation script on ≥200 TCGA patients
2. [ ] Compute AUROC/AUPRC metrics
3. [ ] Present results to Manager
4. [ ] Get explicit approval for SAE→WIWFM integration
5. [ ] Implement SAE module in efficacy orchestrator (behind feature flag)
6. [ ] Apply SAE lifts/penalties using documented policy
7. [ ] Update drug response schema to include SAE context
8. [ ] E2E test with Ayesha's profile
9. [ ] Deploy behind feature flag
10. [ ] Monitor for 1 week before making default

---

## 📊 SUMMARY

### **Current State**: ✅ **ALIGNED WITH MANAGER'S POLICY**

- SAE is "display only" (not affecting drug scores)
- No SAE hooks in `/api/efficacy/predict` drug scoring
- SAE powers Mechanism Map, Hint Tiles, Resistance Playbook
- Manager's guidance: "Wait for validation + written policy"

### **Gap**: ⚠️ **ARCHITECTURAL MISMATCH WITH MANAGER'S VISION**

- Manager's vision: SAE should modulate S/P/E scores
- Current reality: SAE is isolated, not integrated
- **This is INTENTIONAL** (per Manager's policy to wait for validation)

### **Next Action**: ⏸️ **WAIT FOR VALIDATION**

- Jr2 must extract HRD scores from cBioPortal
- Run validation script on ≥200 TCGA patients
- Present results to Manager
- Get explicit approval before implementing SAE→WIWFM integration

---

**Status**: ✅ **REVIEW COMPLETE - CURRENT STATE DOCUMENTED**  
**Recommendation**: ⏸️ **MAINTAIN CURRENT STATE - WAIT FOR VALIDATION APPROVAL**

---

**Last Updated**: January 13, 2025  
**Next Review**: After validation results available



