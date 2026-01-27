# ❓ HOLISTIC CLINICAL BENEFIT SCORE - IMPLEMENTATION QUESTIONS

**Date**: January 29, 2025  
**Purpose**: Questions for Manager to guide implementation  
**Status**: 🔴 **BLOCKERS & CLARIFICATIONS NEEDED**

---

## 🔴 CRITICAL QUESTIONS (Blockers)

### **Q1: Feature Fetching - How Do We Pull D/P/M/T/S from Engines?** ✅ **ANSWERED**

**Issue**: Manager said "pulls D/P/M/T/S from existing engines" but unclear on implementation.

**✅ ANSWER (From Codebase):**
- **Pattern Found**: `trial_service.py` (line 164-243) and `holistic_score/service.py` import engine functions **directly**
- **Implementation**: Import functions directly (no HTTP calls):
  ```python
  from api.services.resistance.biomarkers.diagnostic.ddr_bin_scoring import assign_ddr_status
  from api.services.resistance.biomarkers.therapeutic.timing_chemo_features import build_timing_chemo_features
  from api.services.holistic_score.pgx_safety import compute_pgx_safety
  ```
- **Caching**: Not required initially (engines are fast), can add later if needed
- **Fallbacks**: Use defaults (neutral values 0.5) if engines unavailable

**✅ RESOLVED**: Import functions directly (matches existing patterns)

---

### **Q2: Regimen MoA Vector - Where Does It Come From?** ✅ **RESOLVED**

**Issue**: Component M requires `regimen_moa_vector` (7D vector), but source is unclear.

**✅ ANSWER (From Trial MoA Tagging Pattern):**

**Found Existing Capabilities:**
1. **Trial MoA Tagging Pattern** (`trial_data_enricher.py`, lines 111-147):
   - Runtime keyword matching maps drug mechanisms → 7D MoA dict → 7D vector
   - Uses `get_drug_mechanism()` from `DRUG_MECHANISM_DB` to get mechanism descriptions
   - Then uses keyword matching to map mechanism text → 7D MoA dict (ddr, mapk, pi3k, vegf, her2, io, efflux)
   - Uses `convert_moa_dict_to_vector()` to convert dict → 7D list

2. **Keyword Mapping Pattern** (lines 129-142):
   ```python
   # PARP/ATR/CHK/WEE1 → DDR = 0.7
   # MEK/BRAF/KRAS/MAPK → MAPK = 0.7
   # PI3K/AKT/MTOR → PI3K = 0.7
   # VEGF/angiogenesis → VEGF = 0.7
   # HER2/trastuzumab/pertuzumab → HER2 = 0.7
   # PD-1/PD-L1/CTLA-4/immuno → IO = 0.7
   # ABCB1/P-gp/MDR → Efflux = 0.7
   ```

3. **Functions Available:**
   - `get_drug_mechanism(drug_name)` from `dossier_generator.py` → returns mechanism description
   - `convert_moa_dict_to_vector(moa_dict, use_7d=True)` from `pathway_to_mechanism_vector.py` → converts dict to 7D list
   - `DRUG_TO_MOA` from `toxicity_pathway_mappings.py` → maps drug → MoA category (e.g., "carboplatin" → "platinum_agent")
   - `REGIMEN_TYPE_CLASSIFICATIONS` from `timing_config.py` → maps regimen_type → regimen_class

**✅ SOLUTION: Reuse Trial MoA Tagging Pattern**

**Implementation Approach:**
1. **Extract regimen drugs** from `regimen_drugs` list (or `regimen_type` string)
2. **For each drug**: Use `get_drug_mechanism(drug)` to get mechanism description
3. **Map mechanism text → 7D MoA dict** using same keyword matching as trials:
   - If mechanism contains "parp/atr/chk/wee1" → `moa_dict["ddr"] = 0.9` (use 0.9 for regimens, vs 0.7 for trials - higher confidence)
   - If mechanism contains "mek/braf/kras/mapk" → `moa_dict["mapk"] = 0.9`
   - If mechanism contains "pi3k/akt/mtor" → `moa_dict["pi3k"] = 0.9`
   - If mechanism contains "vegf/angiogenesis" → `moa_dict["vegf"] = 0.9`
   - If mechanism contains "her2/trastuzumab" → `moa_dict["her2"] = 0.9`
   - If mechanism contains "pd-1/pd-l1/ctla-4/immuno" → `moa_dict["io"] = 0.9`
   - If mechanism contains "abcb1/p-gp/mdr/efflux" → `moa_dict["efflux"] = 0.9`
4. **Convert to 7D vector**: Use `convert_moa_dict_to_vector(moa_dict, use_7d=True)`
5. **Fallback**: If no drugs found, use `REGIMEN_TYPE_CLASSIFICATIONS` + simple inference:
   - "platinum" → DDR=0.9
   - "PARPi" → DDR=0.95
   - "IO" → IO=0.9
   - etc.

**Function to Create:**
```python
def compute_regimen_moa_vector(
    regimen_drugs: List[str],
    regimen_type: Optional[str] = None
) -> List[float]:
    """
    Compute 7D MoA vector for a regimen from drugs or regimen_type.
    
    Reuses trial MoA tagging pattern for consistency.
    """
    from api.services.client_dossier.dossier_generator import get_drug_mechanism
    from api.services.pathway_to_mechanism_vector import convert_moa_dict_to_vector
    
    moa_dict = {"ddr": 0.0, "mapk": 0.0, "pi3k": 0.0, "vegf": 0.0, "her2": 0.0, "io": 0.0, "efflux": 0.0}
    
    # Try drugs first
    for drug in regimen_drugs:
        mechanism = get_drug_mechanism(drug)
        if mechanism and mechanism.get("mechanism"):
            mech_text = mechanism["mechanism"].lower()
            # Use same keyword matching as trials (but with 0.9 instead of 0.7)
            if any(kw in mech_text for kw in ["parp", "atr", "chk", "wee1", "dna repair"]):
                moa_dict["ddr"] = max(moa_dict["ddr"], 0.9)
            # ... (repeat for other pathways)
    
    # Fallback to regimen_type if no matches
    if all(v == 0.0 for v in moa_dict.values()) and regimen_type:
        # Use REGIMEN_TYPE_CLASSIFICATIONS + simple inference
        ...
    
    return convert_moa_dict_to_vector(moa_dict, use_7d=True)
```

**✅ RESOLVED**: Reuse existing trial MoA tagging pattern - no new mapping table needed, just reuse keyword matching logic from `trial_data_enricher.py`.

---

### **Q3: Patient Mechanism Vector - Source Priority?** ✅ **ANSWERED**

**Issue**: Component M requires `patient_mechanism_vector` (7D vector), but source could vary.

**✅ ANSWER (From Codebase):**
- **Pattern Found**: `trial_service.py` (line 164-243) shows the exact priority order:
  ```python
  1. From parameter (if provided)
  2. Compute from tumor_context using convert_pathway_scores_to_mechanism_vector()
  3. Default DDR-high vector [0.88, 0.12, 0.15, 0.10, 0.05, 0.2, 0.0]
  ```
- **Function**: `convert_pathway_scores_to_mechanism_vector()` in `pathway_to_mechanism_vector.py` converts pathway scores → 7D vector
- **CSI Integration**: If CSI outputs available, extract pathway_scores from CSI response provenance

**✅ RESOLVED**: 
1. CSI v0 outputs (if available) → extract pathway_disruption
2. SAE features → mechanism_vector (if available)
3. tumor_context → convert_pathway_scores_to_mechanism_vector()
4. Default DDR-high vector [0.88, 0.12, ...]

**Implementation**: Reuse `trial_service._compute_mechanism_vector_from_tumor_context()` pattern

---

### **Q4: Regimen Classification - How Do We Determine `regimen_class`?** ✅ **ANSWERED**

**Issue**: Weights depend on `regimen_class` (e.g., "platinum", "PARPi", "ATRi"), but source is unclear.

**✅ ANSWER (From Codebase):**
- **Function Found**: `get_regimen_biomarker_class()` in `timing_config.py` (line 144-172) ✅
- **Mapping Table**: `REGIMEN_TYPE_CLASSIFICATIONS` maps regimen types to classes:
  ```python
  "platinum": ["platinum", "carboplatin", "cisplatin", "oxaliplatin"],
  "PARPi": ["PARPi", "olaparib", "niraparib", "rucaparib", "talazoparib"],
  "ATR_inhibitor": ["ATRi", "ATR_inhibitor", "berzosertib", "ceralasertib"],
  ```
- **Usage**: `get_regimen_biomarker_class(regimen_type)` returns class or None
- **Default**: If classification is unknown, use "unknown" as regimen_class

**✅ RESOLVED**: Infer from `regimen_type` using `get_regimen_biomarker_class()` function. Default to "unknown" if not found.

---

## ⚠️ DESIGN QUESTIONS (Important)

### **Q5: Weight Normalization - When T Weight is 0.00 in Trial Enrollment?** ✅ **MANAGER DECISION**

**Issue**: Trial enrollment has `T: 0.00` (therapeutic dynamics not relevant pre-treatment).

**✅ MANAGER DECISION: Option B + C**

**For `use_case = "trial_enrollment"`:**
1. **Set `w_T = 0`** and renormalize the other weights to sum to 1 (D, P, M, S only).
2. **Still compute and return T** in the breakdown for transparency, but do not let it affect the overall enrollment score.

**Rationale**: Pre-treatment trial matching should not depend on on-treatment kinetics, even if historical KELIM exists; that belongs in prognostic/predictive models, not in the enrollment feasibility score.

**Implementation**:
- If `use_case == "trial_enrollment"`: set T weight to 0, renormalize D/P/M/S weights
- Compute T component anyway (for breakdown transparency)
- Include T in response but don't multiply by weight in overall score calculation

**✅ RESOLVED**: Option B + C - Renormalize weights, but still compute T for breakdown.

---

### **Q6: Interpretation Thresholds - What Are the Cutoffs?** ✅ **ANSWERED**

**Issue**: Need thresholds for HIGH/MEDIUM/LOW interpretation.

**✅ ANSWER (From Codebase):**
- **Found**: `interpreter.py` (line 51-88) shows exact thresholds:
  ```python
  HIGH: ≥0.8
  MEDIUM: ≥0.6
  LOW: ≥0.4
  VERY_LOW: <0.4
  ```
- **Pattern**: Same thresholds used in `holistic_score` service for consistency

**✅ RESOLVED**: Use same thresholds as holistic_score (HIGH ≥0.8, MEDIUM ≥0.6, LOW ≥0.4, VERY_LOW <0.4) for consistency. No use-case-specific thresholds needed initially.

---

### **Q7: CSI Integration - Should We Pull M and Part of P from CSI?** ✅ **MANAGER DECISION**

**Issue**: Manager said "CSI is the Predictive core (mostly M+part of P)".

**✅ MANAGER DECISION: Option C (CSI-first, with graceful fallback)**

**Treat CSI as the preferred source for the Predictive core:**

**When CSI is available for a given patient-regimen:**
- **M**: Use CSI's mechanism-fit / pathway alignment component (or the calibrated "mechanism fit" probability if CSI exposes it).
- **P_contrib_from_CSI**: Use any explicit baseline-risk/line-adjusted term CSI outputs (e.g., its intercept or non-mechanism features).

**When CSI is not available:**
- **M**: Compute from the 7D mechanism vectors (patient vs regimen) using existing mechanism-fit code.
- **P**: Compute from the Timing/BRI model.

**Implementation Guidance:**
- Keep CSI as a source of inputs, not a black-box override.
- **M and P should still be explicit fields** in the holistic score output.
- **Document clearly in provenance** whether M/P came from CSI or from separate calculations.

**So: use CSI when present, fall back to the separate engines when not.**

**✅ RESOLVED**: Option C - CSI-first with graceful fallback, always expose M/P explicitly.

---

### **Q8: Integration Point - When Should complete_care Call This?** ✅ **MANAGER DECISION**

**Issue**: Manager said "call it as sub-step when proposing trials, ranking regimens, monitoring decisions".

**✅ MANAGER DECISION: Option C – Complement now, plan to replace after validation**

**Phase 1 (Current Sprint):**
- **Complement existing ranking logic:**
  - Keep `mechanism_fit_ranker` and existing `holistic_score` for trials.
  - Add Holistic Clinical Benefit Score as an **additional metric** and expose it via:
    - A new endpoint (e.g., `compute_holistic_clinical_benefit`).
    - Additional fields in responses where appropriate.
  - **Do not remove or change current ranking behavior yet.**

**Phase 2 (After Internal Validation):**
- Once we've run side-by-side comparisons on a set of cases and like the behavior, incrementally:
  - Use Holistic Score as the **primary ranking metric** for specific flows (e.g., trial enrollment suggestions, next-line regimen suggestions).
  - Deprecate older `mechanism_fit_ranker` where redundant.

**Implementation:**
- Implement as a **parallel score now**, with clear hooks so it can be promoted to primary ranking later.
- Add new endpoint: `POST /api/resistance/holistic-clinical-benefit`
- Return Holistic Clinical Benefit Score alongside existing scores (don't replace yet).

**✅ RESOLVED**: Option C - Complement now (parallel metric), plan to make primary after validation.

---

### **Q9: Missing Component Defaults - What Are Neutral Values?** ✅ **ANSWERED**

**Issue**: When components are missing, we need neutral defaults.

**✅ ANSWER (From Codebase):**
- **PGx Safety**: `pgx_safety.py` (line 48-53) returns **1.0** if no pharmacogenes (assumes safe) ✅
- **Mechanism Fit**: Default to **0.5** (neutral alignment) if missing ✅
- **Architecture Doc**: Already specifies **0.5** for missing components ✅

**✅ RESOLVED**: 
- **D**: 0.5 (unknown = neutral)
- **P**: 0.5 (unknown = neutral)
- **M**: 0.5 (neutral alignment)
- **T**: 0.5 (neutral response) - matches code
- **S**: 1.0 (assume safe unless contraindication) - matches PGx safety behavior

---

## 📋 IMPLEMENTATION CLARIFICATIONS (Nice to Have)

### **Q10: Configuration System - File Structure?**

**Question**: Where should config live?
- `api/services/resistance/config/holistic_clinical_benefit_config.py`?
- Or reuse existing config files (`timing_config.py`, `ddr_config.py`)?

**Recommendation**: New file `holistic_clinical_benefit_config.py` for weights, reuse disease-specific configs for thresholds.

**Manager Decision Needed**: ✅ New file OR ✅ Reuse existing?

---

### **Q11: Error Handling - What Should We Return on Failure?**

**Question**: If engine call fails (e.g., DDR_bin unavailable), should we:
- Fail entire score computation?
- Return partial score with available components?
- Return default score with warning flags?

**Recommendation**: Return partial score with available components + `component_available` flags + warnings in provenance.

**Manager Decision Needed**: ✅ Fail OR ✅ Partial OR ✅ Default with warnings?

---

### **Q12: Logging and Monitoring - What Should We Track?**

**Question**: What metrics should we log?
- Component availability rates?
- Score distributions by use_case?
- Missing component patterns?
- Performance (computation time per component)?

**Recommendation**: Log all of the above for monitoring and debugging.

**Manager Decision Needed**: ✅ Confirm tracking OR ✅ Adjust list?

---

### **Q13: Testing Strategy - What Are Success Criteria?**

**Question**: How do we validate the orchestration layer works correctly?
- Unit tests: Each component function?
- Integration tests: End-to-end with real engine outputs?
- Validation tests: Compare to expected scores for known cases?
- Performance tests: Response time < X seconds?

**Recommendation**: All of the above, with synthetic test cases and known patient-regimen pairs.

**Manager Decision Needed**: ✅ Confirm approach OR ✅ Adjust priorities?

---

## 🎯 SUMMARY

### **Critical Blockers (Need Before Implementation):**
1. ✅ **Q1: RESOLVED** - Import functions directly (matches existing patterns)
2. ⚠️ **Q2: PARTIALLY RESOLVED** - Need regimen MoA vector mapping (drug/type → 7D vector)
3. ✅ **Q3: RESOLVED** - Priority: CSI → SAE → tumor_context → defaults (matches trial_service.py)
4. ✅ **Q4: RESOLVED** - Use `get_regimen_biomarker_class()` from timing_config.py

### **Design Decisions (Important):**
5. ❓ **Q5: NEED MANAGER** - Weight normalization when T weight is 0.00 in trial enrollment
6. ✅ **Q6: RESOLVED** - Same thresholds as holistic_score (HIGH ≥0.8, MEDIUM ≥0.6, LOW ≥0.4)
7. ❓ **Q7: NEED MANAGER** - CSI integration approach (use CSI outputs vs separate computation)
8. ❓ **Q8: NEED MANAGER** - Integration point (replace vs complement existing logic)
9. ✅ **Q9: RESOLVED** - Defaults: D/P/M/T=0.5, S=1.0 (matches PGx safety behavior)

### **Implementation Details (Nice to Have):**
10. ✅ **Q10: INFERRED** - New config file `holistic_clinical_benefit_config.py` (reuse timing_config pattern)
11. ✅ **Q11: INFERRED** - Return partial scores with component_available flags + warnings
12. ✅ **Q12: INFERRED** - Log component availability, score distributions, performance metrics
13. ✅ **Q13: INFERRED** - Unit tests (each component), integration tests (end-to-end), synthetic test cases

---

## 📋 **REMAINING QUESTIONS FOR MANAGER**

### **✅ Q2: Regimen MoA Vector Mapping - RESOLVED**
**Question**: How should we map regimen drugs/types to 7D MoA vectors?

**✅ ANSWER**: Reuse existing trial MoA tagging pattern from `trial_data_enricher.py`:
- Use `get_drug_mechanism()` to get mechanism descriptions from `DRUG_MECHANISM_DB`
- Use same keyword matching logic (PARP/ATR → DDR, MEK/BRAF → MAPK, etc.)
- Use `convert_moa_dict_to_vector()` to convert to 7D vector
- **No new mapping table needed** - reuse proven pattern from trials

### **❓ Q5: Weight Normalization When T Weight is 0.00**
**Question**: If T weight is 0.00 for trial enrollment but T data exists, should we:
- **Option A**: Ignore T data (use T weight 0.00, don't include T in overall)
- **Option B**: Set T weight to 0 and renormalize other weights
- **Option C**: Use T data but with 0 weight (include in breakdown but not overall)

### **❓ Q7: CSI Integration**
**Question**: If CSI outputs are available, should we:
- **Option A**: Use CSI for M and part of P, compute D/T/S separately
- **Option B**: Always compute M/P separately (CSI is parallel, not source)
- **Option C**: Use CSI if available, fallback to separate computation

### **❓ Q8: Integration Point**
**Question**: How should this integrate with existing logic?
- **Option A**: Replace existing ranking logic (mechanism_fit_ranker, holistic_score)
- **Option B**: Complement existing logic (add as additional metric)
- **Option C**: Both (complement now, replace later after validation)

---

**Priority**: **P0** - All questions resolved ✅  
**Timeline**: Ready for implementation (5-6 day sprint)  
**Status**: 13/13 questions resolved ✅ - **READY FOR IMPLEMENTATION**

**Key Decisions Summary:**
- ✅ **Q5**: Option B+C - T weight = 0 for trial enrollment, renormalize, but still compute T for breakdown
- ✅ **Q7**: Option C - CSI-first with graceful fallback, always expose M/P explicitly
- ✅ **Q8**: Option C - Complement existing logic now, plan to replace after validation
- ✅ **Q2**: Reuse trial MoA tagging pattern - no new infrastructure needed

**All Manager Decisions Received - Ready to Implement! 🚀**

---

**Last Updated**: January 29, 2025  
**Status**: 🔴 **AWAITING MANAGER GUIDANCE**
