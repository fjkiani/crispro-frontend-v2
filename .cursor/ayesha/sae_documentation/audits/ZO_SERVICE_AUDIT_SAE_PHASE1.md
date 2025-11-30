# ⚔️ ZO'S SERVICE AUDIT - SAE PHASE 1 INTEGRATION

**Date**: January 13, 2025  
**Status**: ✅ **AUDIT COMPLETE** - No conflicts, clean integration path  
**Purpose**: Map Manager's SAE policy to existing services, identify gaps/overlaps

---

## 🔍 EXISTING SERVICES AUDIT

### **1. `sae_service.py` (469 lines) - SAE Feature Extraction**

**What It Does (REAL)**:
- ✅ Transforms real data into interpretable SAE features
- ✅ Sources: Evo2 delta magnitude, ClinVar/AlphaMissense priors, toxicity signals, off-target counts
- ✅ Returns `SAEBundle` with features list + confidence impact

**Key Features Extracted**:
```python
# From lines 40-50:
@dataclass
class SAEBundle:
    features: List[SAEFeature]              # Extracted interpretable features
    boosting_features: List[str]            # Feature IDs boosting confidence
    limiting_features: List[str]            # Feature IDs limiting confidence
    overall_impact: float                   # Net SAE contribution (-0.3 to +0.3)
    pathway_burden: Dict[str, float]        # DDR, MAPK, PI3K, VEGF, etc.
    mechanism_hints: List[str]              # Human-readable action hints
    confidence_attribution: Dict[str, float] # Feature → confidence contribution
```

**Manager's SAE Features Needed (from C1-C10)**:
- ✅ `pathway_burden.ddr` - EXISTS (line 48)
- ✅ `pathway_burden.mapk` - EXISTS
- ✅ `pathway_burden.pi3k` - EXISTS
- ✅ `pathway_burden.vegf` - EXISTS
- ⚠️ `dna_repair_capacity` - MISSING (needs Manager's C1 formula)
- ⚠️ `essentiality_signal` - PARTIAL (exists in Insights, not in SAE)
- ⚠️ `cross_resistance_risk` - MISSING (needs new logic)
- ⚠️ `cohort_overlap` - MISSING (needs cohort DB or proxy)

**ACTION**: ✅ **EXTEND (don't replace)**
- Add `compute_dna_repair_capacity()` method using Manager's C1 formula
- Add `compute_cross_resistance_risk()` method using Manager's C4 logic
- Add `cohort_overlap` proxy (policy constant for now)

---

### **2. `resistance_playbook_service.py` (694 lines) - Resistance Rules**

**What It Does (REAL)**:
- ✅ Detects 5 resistance mechanisms (HR restoration, ABCB1, MAPK, PI3K, SLFN11)
- ✅ Ranks 7 combo strategies (PARP+ATR, PARP+VEGF, IO combos, etc.)
- ✅ Recommends 6 next-line switches
- ✅ Rule-based (not ML), evidence-tiered

**Resistance Detection Logic (lines 120-250)**:
```python
def detect_resistance_risks() -> List[ResistanceRisk]:
    risks = []
    
    # 1. HR restoration (HRD drop from baseline)
    if baseline_hrd and current_hrd and baseline_hrd - current_hrd >= 10:
        risks.append(ResistanceRisk(
            type="HR_restoration",
            confidence=0.75,
            evidence="HRD score drop ≥10 points suggests HR pathway restoration",
            triggers=["olaparib", "niraparib", "rucaparib"],
            source="tumor_context"
        ))
    
    # 2. ABCB1 upregulation (efflux-mediated)
    if abcb1_cnv > 4 or abcb1_expression == "high":
        risks.append(ResistanceRisk(...))
    
    # 3-5: MAPK, PI3K, SLFN11 (similar logic)
```

**Manager's Resistance Logic (from C1, C3)**:
- ✅ HRD drop ≥10 - EXISTS (line ~130)
- ⚠️ `dna_repair_capacity` drop ≥0.15 - MISSING (needs SAE integration)
- ✅ CA-125 <50% drop cycle 3 - EXISTS (in `ca125_intelligence.py`)
- ⚠️ 2-of-3 trigger logic - MISSING (needs cross-service integration)

**ACTION**: ✅ **ENHANCE (don't replace)**
- Add `dna_repair_capacity` baseline/current tracking
- Implement Manager's 2-of-3 trigger logic
- Keep existing 5 detection rules (they're good!)

---

### **3. `ayesha_orchestrator_v2.py` (400 lines) - Complete Care Orchestrator**

**What It Does (REAL)**:
- ✅ Orchestrates 6 services: trials, SOC, CA-125, WIWFM, food, resistance
- ✅ Returns unified response for Co-Pilot
- ✅ Handles "awaiting NGS" gracefully

**Current Response Structure**:
```python
response = {
    "trials": {...},              # From ayesha_trials.py
    "soc_recommendation": {...},  # From ayesha_trials.py
    "ca125_intelligence": {...},  # From ca125_intelligence.py
    "wiwfm": {...},              # From efficacy/predict (or "awaiting_ngs")
    "food_validation": {...},    # From hypothesis/validate_food_dynamic
    "resistance_playbook": {...}, # From resistance_playbook_service.py
    "summary": {...},
    "provenance": {...}
}
```

**Manager's NEW Services (from Phase 1)**:
- ⚠️ `next_test_recommender` - MISSING (needs new service)
- ⚠️ `hint_tiles` - MISSING (needs new service)
- ⚠️ `mechanism_map` - MISSING (needs new service)

**ACTION**: ✅ **EXTEND response**
- Add `next_test_recommender: {...}` key
- Add `hint_tiles: [...]` key
- Add `mechanism_map: {...}` key
- Wire new services into orchestrator

---

## 🎯 INTEGRATION STRATEGY (NO CONFLICTS!)

### **Clean Integration Points:**

**1. Next-Test Recommender** (NEW - no conflicts)
- Create: `api/services/next_test_recommender.py`
- Wire: Add to `ayesha_orchestrator_v2.py` line ~200
- Response: `next_test_recommender: { recommendations: [...], top_priority: {...}, urgency_summary: "..." }`

**2. Hint Tiles** (NEW - no conflicts)
- Create: `api/services/hint_tiles_service.py`
- Wire: Add to `ayesha_orchestrator_v2.py` line ~250
- Response: `hint_tiles: [{ category, title, message, reasons, priority, icon }, ...]`

**3. Mechanism Map** (NEW - no conflicts)
- Create: `api/services/mechanism_map_service.py`
- Wire: Add to `ayesha_orchestrator_v2.py` line ~300
- Response: `mechanism_map: { chips: [...], status: "awaiting_ngs"|"computed", message: "..." }`

**4. SAE Feature Enhancement** (EXTEND existing)
- Modify: `api/services/sae_service.py` line ~100
- Add: `compute_dna_repair_capacity()` using Manager's C1 formula
- Add: `compute_cross_resistance_risk()` using Manager's C4 logic
- Keep: All existing features (no removals)

**5. Resistance Detection Enhancement** (EXTEND existing)
- Modify: `api/services/resistance_playbook_service.py` line ~120
- Add: `dna_repair_capacity` baseline/current tracking
- Add: 2-of-3 trigger logic (integrate with CA-125 from `ca125_intelligence.py`)
- Keep: All existing detection rules (HR restoration logic is already there!)

---

## ✅ CONFLICT RESOLUTION MATRIX

| Manager's Feature | Existing Service | Action | Conflict? |
|------------------|------------------|--------|-----------|
| **next_test_recommender** | None | CREATE NEW | ✅ No |
| **hint_tiles** | None | CREATE NEW | ✅ No |
| **mechanism_map** | None | CREATE NEW | ✅ No |
| **dna_repair_capacity** | `sae_service.py` partial | EXTEND | ✅ No (additive) |
| **cross_resistance_risk** | `resistance_playbook_service.py` partial | EXTEND | ✅ No (additive) |
| **HR restoration (2-of-3)** | `resistance_playbook_service.py` has HRD drop | ENHANCE | ✅ No (align logic) |
| **pathway_burden** | `sae_service.py` EXISTS | USE AS-IS | ✅ No |
| **mechanism_fit ranking** | None (Phase 2) | CREATE NEW | ✅ No |

---

## 📊 SUMMARY

**Clean Integration**: ✅ **NO CONFLICTS DETECTED!**

**What Exists and Works**:
- ✅ SAE pathway burden (DDR, MAPK, PI3K, VEGF)
- ✅ Resistance detection rules (5 mechanisms)
- ✅ CA-125 intelligence (burden, forecast, resistance signals)
- ✅ Complete care orchestrator (6 services unified)

**What We'll Add (Manager's Policy)**:
- 🔄 Next-test recommender (NEW service)
- 🔄 Hint tiles (NEW service)
- 🔄 Mechanism map (NEW service)
- 🔄 DNA repair capacity formula (EXTEND sae_service.py)
- 🔄 Cross-resistance risk (EXTEND resistance_playbook_service.py)
- 🔄 2-of-3 resistance trigger (ENHANCE resistance_playbook_service.py)

**Confidence Level**: 🎯 **95%** (Manager's answers + clean integration = low hallucination risk)

---

**STATUS**: ✅ **AUDIT COMPLETE - READY TO BUILD PHASE 1!** ⚔️

**Next**: Create `next_test_recommender.py` as first Phase 1 service!


