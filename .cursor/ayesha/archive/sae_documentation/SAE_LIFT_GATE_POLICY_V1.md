# ⚔️ SAE MASTER DOCUMENT - COMPLETE REFERENCE

**Date:** January 13, 2025  
**Owner:** Zo  
**Status:** 📝 **POLICY DOCUMENTED - AWAITING VALIDATION**  
**Single Source of Truth:** This document contains all SAE context, policy, and status

**⚠️ CRITICAL:** SAE lifts/penalties are DOCUMENTATION ONLY. Do NOT implement until:
1. HRD/platinum validation is running (≥200 TCGA patients)
2. Manager explicitly approves implementation

---

## 📋 TABLE OF CONTENTS

1. [What We Wanted to Do](#what-we-wanted-to-do)
2. [Current State](#current-state)
3. [Manager's Policy](#managers-policy)
4. [Lift/Gate Rules](#liftgate-rules)
5. [Implementation Plan](#implementation-plan)
6. [Status & Next Steps](#status--next-steps)

---

## 🎯 WHAT WE WANTED TO DO

### **Original Intent: SAE Integration with WIWFM**

**Goal:** Integrate SAE (Sparse Autoencoder) features into the WIWFM (Will It Work For Me) drug efficacy prediction pipeline to improve recommendations.

**Vision:**
- **S/P/E → SAE → Confidence**: Make SAE part of the S/P/E pipeline: "S → P → E → SAE → Confidence → Final Score"
- **Confidence Modulation**: Apply SAE-based lifts/penalties to drug confidence scores
  - PARP inhibitors: +0.10 if DNA repair capacity <0.40 (HRD rescue)
  - MEK/RAF inhibitors: +0.15 if hotspot mutation detected (KRAS/BRAF/NRAS)
  - Taxanes: -0.20 if cross-resistance risk ≥0.70
- **Transparency**: Show SAE contribution in confidence breakdown and drug rationale
- **Provenance**: Track all SAE lifts/penalties with full audit trail

**Expected Impact:**
- More accurate drug recommendations (confidence 0.48-0.51 → 0.58-0.65 with SAE)
- Better PARP rescue for sporadic cancer patients (germline-negative + HRD≥42)
- Hotspot-driven MEK/RAF recommendations
- Early resistance detection (3-6 weeks faster than imaging)

---

## 📊 CURRENT STATE

### **Architecture: SAE is "Display Only"**

**What Works:**
- ✅ SAE features computed in `sae_feature_service.py` (DNA repair capacity, mechanism vector, hotspot detection)
- ✅ SAE features extracted in `efficacy_orchestrator.py` (optional, line 330-380)
- ✅ SAE displayed in UI (Mechanism Map, Hint Tiles, Resistance Playbook)
- ✅ SAE attribution tracked in provenance

**What's Missing:**
- ❌ **SAE features NOT used in `drug_scorer.py`** (no confidence modulation)
- ❌ **No SAE lifts/penalties applied to drug scores**
- ❌ **SAE lives in Ayesha orchestrator, not integrated into S/P/E pipeline**

**Current Data Flow:**
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

**Why Not Integrated Yet:**
- Manager's explicit guidance: "Wait for validation + written policy"
- Current state is **INTENTIONAL** (per Manager's policy)
- Validation blocked on Jr2 HRD extraction from cBioPortal

---

## 📜 MANAGER'S POLICY

### **Approved Strategy (Q2a, Q2b, Q2c)**

**Phase 1 (NOW - Current State):**
- ✅ Keep SAE in `ayesha_orchestrator_v2.py` (display + Resistance Playbook)
- ✅ Do NOT touch `/api/efficacy/predict` during P0 triage
- ✅ SAE powers: Next-Test, Hint Tiles, Mechanism Map, Resistance Playbook
- ❌ NO SAE-based lifts/caps in WIWFM ranking

**Phase 2 (AFTER Validation + Written Policy):**
- Add SAE module inside `/api/efficacy/predict` behind feature flag
- Compute SAE features alongside S/P/E
- Feed SAE into confidence/lift logic ONLY after:
  1. HRD/platinum validation is running (blocked on Jr2)
  2. Written SAE policy for lifts/gates exists ✅ (this document)

**Phase 3 (Final - Once Stable):**
- Make SAE-enhanced efficacy the default
- Keep "Baseline (no SAE)" profile for comparison

### **Validation Gates**

**"Validation is running" means:**
1. HRD scores successfully extracted for TCGA-OV cohort (via cBioPortal)
2. Validation script executes end-to-end on ≥200 patients
3. Initial AUROC/AUPRC computed for platinum response
4. DNA repair capacity ↔ HRD correlation calculated

**Performance Thresholds:**
- Platinum response AUROC ≥0.60 (baseline)
- DNA repair ↔ HRD correlation r ≥0.50
- No worse than baseline for HRD-negative patients

---

## ⚔️ LIFT/GATE RULES

## 1. PARP LIFT/PENALTY RULES

### **Lift Conditions (+0.10 confidence)**

**When to Apply:**
- DNA repair capacity <0.40 (LOW - indicates HRD)
- AND HRD score ≥42 (PARP eligibility gate)
- AND No resistance signals detected

**Rationale:** Low DNA repair capacity suggests homologous recombination deficiency, indicating strong PARP sensitivity.

**Implementation:**
```python
if sae_features.dna_repair_capacity < 0.40 and tumor_context.hrd_score >= 42:
    if not sae_features.resistance_detected:
        drug.confidence += 0.10
        drug.sae_reasoning.append("Low DNA repair capacity (<0.40) supports PARP efficacy")
```

### **Penalty Conditions (-0.15 confidence)**

**When to Apply:**
- HR restoration pattern detected (2-of-3 triggers):
  - HRD drop ≥10 points vs baseline
  - DNA repair capacity drop ≥0.15 vs baseline
  - CA-125 <50% drop by cycle 3 OR on-therapy rise

**Rationale:** HR restoration indicates reversion mutations or alternative DNA repair pathway activation, leading to PARP resistance.

**Implementation:**
```python
if sae_features.resistance_detected and sae_features.triggers_count >= 2:
    if drug.mechanism == "PARP":
        drug.confidence -= 0.15
        drug.sae_reasoning.append(f"HR restoration detected ({sae_features.triggers_count} of 3 triggers)")
```

---

## 2. MEK/RAF HOTSPOT GATES

### **Lift Conditions (+0.15 confidence)**

**When to Apply:**
- KRAS/BRAF/NRAS hotspot mutation detected (COSMIC verified)
- AND MAPK pathway burden ≥0.40

**Rationale:** Hotspot mutations are activating mutations with high pathway burden, indicating strong MEK/RAF inhibitor efficacy.

**Implementation:**
```python
if sae_features.hotspot_mutation and sae_features.pathway_burden_mapk >= 0.40:
    if drug.mechanism in ["MEK", "RAF"]:
        drug.confidence += 0.15
        hotspot = sae_features.hotspot_details["mutation"]
        drug.sae_reasoning.append(f"MAPK hotspot ({hotspot}) with burden {sae_features.pathway_burden_mapk:.2f}")
```

### **No Boost Conditions (Show trials but no monotherapy boost)**

**When to Apply:**
- Hotspot present BUT MAPK burden <0.40

**Rationale:** Hotspot mutation present but low pathway burden suggests other dominant pathways; MEK monotherapy unlikely effective.

**Implementation:**
```python
if sae_features.hotspot_mutation and sae_features.pathway_burden_mapk < 0.40:
    # Show trials but no confidence boost
    drug.sae_reasoning.append(f"MAPK hotspot detected but burden low ({sae_features.pathway_burden_mapk:.2f}) - combo preferred")
```

### **Deprioritize Monotherapy (-0.15 confidence)**

**When to Apply:**
- MEK monotherapy candidate
- AND MAPK burden <0.40
- AND No hotspot mutation

**Rationale:** Without hotspot or high pathway burden, MEK monotherapy ineffective in ovarian cancer.

**Implementation:**
```python
if drug.mechanism == "MEK" and drug.type == "monotherapy":
    if sae_features.pathway_burden_mapk < 0.40 and not sae_features.hotspot_mutation:
        drug.confidence -= 0.15
        drug.sae_reasoning.append("Low MAPK burden without hotspot - MEK monotherapy not recommended")
```

---

## 3. HER2 PATHWAY THRESHOLDS

### **Lift Conditions (+0.12 confidence)**

**When to Apply:**
- HER2 pathway burden ≥0.70 (HIGH)
- HER2 amplification detected (CNV >4) OR HER2 IHC 3+

**Rationale:** High HER2 pathway burden indicates HER2 addiction, supporting HER2-targeted therapy efficacy.

**Implementation:**
```python
if sae_features.pathway_burden_her2 >= 0.70:
    if tumor_context.her2_status in ["amplified", "IHC_3+"]:
        if drug.mechanism == "HER2":
            drug.confidence += 0.12
            drug.sae_reasoning.append(f"High HER2 pathway burden ({sae_features.pathway_burden_her2:.2f})")
```

---

## 4. CROSS-RESISTANCE PENALTIES

### **Penalty Conditions (-0.20 confidence)**

**When to Apply:**
- Cross-resistance risk ≥0.70 (HIGH)
- Prior taxane exposure with progression ≤6 months
- OR ABCB1 CNV >4 (efflux pump overexpression)

**Rationale:** High ABCB1 or prior taxane failure indicates efflux-mediated resistance, affecting all taxane substrates.

**Implementation:**
```python
if sae_features.cross_resistance_risk >= 0.70:
    if drug.substrate_class == "taxane":
        drug.confidence -= 0.20
        drug.sae_reasoning.append(f"Cross-resistance risk ({sae_features.cross_resistance_risk:.2f}) - taxane inefficacy likely")
        drug.sae_reasoning.append("Consider non-substrates: platinum, PARP, ATR/CHK1")
```

### **Non-Substrate Alternatives (No Penalty)**

**Drugs NOT Affected:**
- Platinum agents (carboplatin, cisplatin)
- PARP inhibitors (olaparib, niraparib, rucaparib)
- ATR/CHK1/WEE1 inhibitors
- Immunotherapy (PD-1/PD-L1)

---

## 5. CONFIDENCE CAPS

### **Cap at 0.60 (60%) Confidence**

**When to Apply:**
- Mechanism vector all gray (<0.40 across ALL pathways)
- AND No hotspot mutations detected
- AND DNA repair capacity in gray zone (0.40-0.70)

**Rationale:** Weak signals across all pathways suggest uncertain efficacy; cap confidence to avoid over-recommendation.

**Implementation:**
```python
if all(burden < 0.40 for burden in sae_features.mechanism_vector.values()):
    if not sae_features.hotspot_mutation and 0.40 <= sae_features.dna_repair_capacity <= 0.70:
        if drug.confidence > 0.60:
            drug.confidence = 0.60
            drug.sae_reasoning.append("Confidence capped at 60% - weak mechanism signals across all pathways")
```

---

## 6. PROVENANCE REQUIREMENTS

### **Required Logging for All Lifts/Penalties:**

```python
provenance = {
    "sae_version": "v1.0",
    "policy_version": "v1.0",
    "policy_source": "SAE_LIFT_GATE_POLICY_V1.md",
    "manager": "SR",
    "date": "2025-01-13",
    "feature_flags": {
        "sae_enabled": True,
        "lifts_enabled": True,
        "gates_enabled": True
    },
    "thresholds_used": {
        "dna_repair_high": 0.70,
        "dna_repair_low": 0.40,
        "mapk_burden_threshold": 0.40,
        "her2_high_threshold": 0.70,
        "cross_resistance_high": 0.70
    },
    "lifts_applied": [
        {
            "drug": "Olaparib",
            "mechanism": "PARP",
            "lift": +0.10,
            "reason": "Low DNA repair capacity (0.32)",
            "confidence_before": 0.75,
            "confidence_after": 0.85
        }
    ],
    "penalties_applied": [
        {
            "drug": "Paclitaxel",
            "substrate_class": "taxane",
            "penalty": -0.20,
            "reason": "Cross-resistance risk (0.85)",
            "confidence_before": 0.65,
            "confidence_after": 0.45
        }
    ]
}
```

### **Transparency Requirements:**

1. **All lifts/penalties must be visible in UI**
   - Show SAE delta in drug card
   - Display reasoning in expandable section
   - Include RUO label

2. **Confidence breakdown must show SAE contribution**
   ```python
   confidence_breakdown = {
       "evo2_contribution": 0.40,
       "pathway_contribution": 0.30,
       "sae_contribution": 0.30
   }
   ```

3. **Provenance logs must be versioned and timestamped**
   - Track SAE version, policy version, thresholds used
   - Store in database for audit trail
   - Surface in "View Details" modal

---

## 7. VALIDATION GATES

### **Before Implementation:**

**Must Complete:**
1. ✅ Extract HRD scores for TCGA-OV cohort (≥200 patients)
2. ✅ Run validation script end-to-end
3. ✅ Achieve baseline AUROC for platinum response
4. ✅ Correlate DNA repair capacity with HRD scores
5. ✅ Manager review and explicit approval

**Performance Thresholds:**
- Platinum response AUROC ≥0.60 (baseline)
- DNA repair ↔ HRD correlation r ≥0.50
- No worse than baseline for HRD-negative patients

---

## 8. IMPLEMENTATION CHECKLIST

**When Manager Approves Implementation:**

- [ ] Create feature flag `EFFICACY_ENABLE_SAE=true`
- [ ] Add SAE computation to `/api/efficacy/predict` (after pathway, before confidence)
- [ ] Import `sae_modulation_rules.py`
- [ ] Apply lifts/penalties using `modulate_drug_confidence()`
- [ ] Update `DrugEfficacyResponse` schema to include `sae_features`
- [ ] Add `confidence_breakdown` with SAE contribution
- [ ] Surface SAE reasoning in drug responses
- [ ] Update frontend to display SAE delta and reasoning
- [ ] Log all lifts/penalties to provenance
- [ ] Run E2E tests with Ayesha's profile (BRCA1, KRAS G12D)
- [ ] Deploy behind feature flag
- [ ] Monitor for 1 week before making default

---

## 9. SAFETY GUARDRAILS

**Hard Limits:**
- **Max Lift:** +0.15 (MEK/RAF hotspot)
- **Max Penalty:** -0.20 (cross-resistance)
- **Min Confidence:** 0.10 (never drop below 10%)
- **Max Confidence:** 0.95 (never exceed 95%, even with lifts)

**Override Conditions:**
- If SAE service fails → graceful fallback (no lifts/penalties)
- If tumor_context missing → SAE disabled (pre-NGS mode)
- If any threshold ambiguous → default to no lift/penalty

**RUO Labeling:**
- All SAE-derived lifts/penalties must include "RUO: Research Use Only"
- UI must display disclaimer: "SAE features are investigational and not validated for clinical use"

---

---

## 📊 STATUS & NEXT STEPS

### **Current Status: ⚔️ POLICY DOCUMENTED - AWAITING VALIDATION**

**What's Complete:**
- ✅ SAE feature computation service (`sae_feature_service.py`)
- ✅ DNA repair capacity formula (0.6/0.2/0.2 weights - Manager's C1)
- ✅ Hotspot mutation detection (KRAS/BRAF/NRAS - P0 Fix #3)
- ✅ Mechanism vector computation (7D for cosine similarity)
- ✅ Resistance detection (2-of-3 trigger logic)
- ✅ SAE lift/gate policy documented (this document)
- ✅ SAE display in UI (Mechanism Map, Hint Tiles, Resistance Playbook)

**What's Blocked:**
- ⏸️ SAE→WIWFM integration (waiting for validation)
- ⏸️ HRD score extraction (Jr2 mission - cBioPortal)
- ⏸️ Validation script execution (needs HRD ground truth)

**What's NOT Implemented (By Design):**
- ❌ SAE confidence modulation in `/api/efficacy/predict`
- ❌ SAE lifts/penalties applied to drug scores
- ❌ SAE integration into S/P/E pipeline

**Reason:** Manager's explicit guidance - "Wait for validation + written policy"

---

### **Implementation Checklist (When Approved)**

**Before Implementation:**
- [ ] Jr2 extracts HRD scores from cBioPortal (≥200 TCGA patients)
- [ ] Validation script runs end-to-end
- [ ] AUROC/AUPRC computed and reviewed
- [ ] Manager explicitly approves implementation

**During Implementation:**
- [ ] Create feature flag `EFFICACY_ENABLE_SAE=true`
- [ ] Add SAE computation to `/api/efficacy/predict` (after pathway, before confidence)
- [ ] Import `sae_modulation_rules.py`
- [ ] Apply lifts/penalties using `modulate_drug_confidence()`
- [ ] Update `DrugEfficacyResponse` schema to include `sae_features`
- [ ] Add `confidence_breakdown` with SAE contribution
- [ ] Surface SAE reasoning in drug responses
- [ ] Update frontend to display SAE delta and reasoning
- [ ] Log all lifts/penalties to provenance
- [ ] Run E2E tests with Ayesha's profile (BRCA1, KRAS G12D)
- [ ] Deploy behind feature flag
- [ ] Monitor for 1 week before making default

---

### **Next Steps (Immediate)**

1. **Wait for Jr2** - HRD extraction from cBioPortal
2. **Run Validation** - Execute validation script on ≥200 patients
3. **Present Results** - Show AUROC/AUPRC to Manager
4. **Get Approval** - Manager explicitly approves SAE→WIWFM integration
5. **Implement** - Follow implementation checklist above

---

**Status:** ⚔️ **POLICY DOCUMENTED - AWAITING VALIDATION + MANAGER APPROVAL** ⚔️

**Last Updated:** January 13, 2025  
**Owner:** Zo  
**Reference:** This is the single source of truth for all SAE context







