# 🤔 AGENT JR: CLARIFYING QUESTIONS - ✅ ANSWERED

**Date:** January 28, 2025  
**Mission:** MM Resistance Prediction - Next Phase  
**Status:** ✅ **ANSWERS PROVIDED - READY TO BUILD**

---

## 🎯 ARCHITECTURE DECISION: OPTION C (DRY/SHARED)

> **Manager Decision:** Build a shared `ResistancePlaybookService` that handles BOTH Ovarian and Myeloma.
> Same service, disease-specific mappings. DRY principle.

---

## 📋 TASK 1: ResistancePlaybookService

### Q1.1: Service Architecture ✅ ANSWERED
**Question:** Should `ResistancePlaybookService` be standalone or integrated?

**ANSWER:** ✅ **Standalone service** (`api/services/resistance_playbook_service.py`)
- Disease-agnostic design
- Accepts `disease: str` parameter ("ovarian" | "myeloma")
- Routes to disease-specific mappings internally

---

### Q1.2: Multiple Resistance Mechanisms ✅ ANSWERED
**Question:** When multiple resistances detected, how to prioritize?

**ANSWER:** ✅ **Option A: Union of alternatives, weighted by RR**
- Combine all alternatives from all detected genes
- Rank by: (1) priority field, (2) RR of source gene
- Deduplicate (if same drug appears for multiple genes, use highest priority)

---

### Q1.3: Evidence Level Assignment ✅ ANSWERED
**Question:** Evidence level for low-power genes (PSMB5/CRBN)?

**ANSWER:** ✅ **`"LOW_EVIDENCE"`** with note: "n=2-3 in MMRF, mechanism-based"

---

## 📋 TASK 2: Cytogenetics Support

### Q2.1: Data Availability ✅ ANSWERED
**Question:** MMRF has no cytogenetics. Strategy?

**ANSWER:** ✅ **Option A: Implement with graceful fallback**
- Accept cytogenetics as optional input
- Use literature values (del(17p) HR=2.5, t(4;14) HR=1.8)
- Mark as "LITERATURE_BASED" evidence
- If cytogenetics not provided, skip (no error)

---

### Q2.2: Cytogenetics + Gene-Level Risk ✅ ANSWERED
**Question:** How to combine DIS3 + del(17p)?

**ANSWER:** ✅ **Option D: Separate signals, both flagged**
- Both contribute to risk independently
- Combined risk_level = "VERY_HIGH" if 2+ signals
- Show each signal in output with individual RR/HR

---

## 📋 TASK 3: Treatment Line Context

### Q3.1: RR Adjustment Multipliers ✅ ANSWERED
**Question:** Are multipliers validated?

**ANSWER:** ✅ **Mark as "EXPERT_OPINION"**
- 1.2x for 2nd line, 1.4x for 3rd+ are reasonable estimates
- Document in provenance as "EXPERT_OPINION - clone evolution"
- Can validate later if needed

---

### Q3.2: Cross-Resistance Logic ✅ ANSWERED
**Question:** Prior PI → current PI cross-resistance?

**ANSWER:** ✅ **Option A: Additive (base RR × 1.3) + separate flag**
- Multiply base RR by 1.3 for same-class prior exposure
- Also flag "prior_PI_exposure" in output for transparency

---

## 📋 TASK 4: Downstream Agent Handoffs

### Q4.1: Agent Availability ✅ ANSWERED
**Question:** Monitoring agent doesn't exist. What to do?

**ANSWER:** ✅ **Option A: Create handoff contract**
- Define the structure for all 3 agents (Drug Efficacy, Care Plan, Monitoring)
- Return structured requests, orchestrator handles routing
- Monitoring contract is ready for future agent implementation

---

### Q4.2: Handoff Format ✅ ANSWERED
**Question:** Sync vs async handoffs?

**ANSWER:** ✅ **Option B: Async structured request**
- Return `downstream_handoffs: Dict` with structured requests
- Orchestrator decides when/how to call agents
- Follows existing `build_complete_care_plan()` pattern

---

## 📋 TASK 5: API Endpoint

### Q5.1: Router Location ✅ ANSWERED
**Question:** Where to put `/api/resistance/predict`?

**ANSWER:** ✅ **New router: `api/routers/resistance.py`**
- NOT `/api/mm/resistance` - make it disease-agnostic
- Single endpoint: `POST /api/resistance/predict`
- Accept `disease` parameter in request body
- Works for both OV and MM

---

### Q5.2: Authentication ✅ ANSWERED
**Question:** Auth requirements?

**ANSWER:** ✅ **Follow existing pattern** (no auth for internal APIs)
- Check other routers in same directory
- Apply same auth pattern

---

## 📋 TASK 6: Frontend Component (OPTIONAL)

### Q6.1: Component Location ✅ ANSWERED
**Question:** Where to put resistance panel?

**ANSWER:** ✅ **Option C: Component + Integration**
- Create reusable `ResistancePanel.jsx` component
- Integrate into `MyelomaDigitalTwin.jsx`
- Same component can be used for OV later

---

## 📋 TASK 7: NFE2L2/XBP1/IRE1 Genes

### Q7.1: Validation Before Adding ✅ ANSWERED
**Question:** Check MMRF first or add directly?

**ANSWER:** ✅ **Option B: Add with "LITERATURE_BASED" flag**
- Add genes to playbook with literature references
- Check MMRF in parallel (quick search)
- If found in MMRF, upgrade evidence level

---

### Q7.2: Evidence Level ✅ ANSWERED
**Question:** Evidence level for literature genes?

**ANSWER:** ✅ **`"LITERATURE_BASED"`** with PubMed reference

---

## 📋 TASK 8: Validate del(17p) from MMRF

### Q8.1: Data Not Available ✅ ANSWERED
**Question:** MMRF has no cytogenetics. Strategy?

**ANSWER:** ✅ **Option B: Use literature value (HR=2.5)**
- Mark as "LITERATURE_BASED" not "VALIDATED"
- Document source: IMWG consensus guidelines
- Skip MMRF validation (no data)

---

## 📋 TASK 9: Evo2 Delta → Response (OPTIONAL)

### Q9.1: Priority ✅ ANSWERED
**Question:** Skip or implement?

**ANSWER:** ✅ **Option A: Skip entirely**
- Not needed for current scope
- Can add later if outcome prediction needed

---

## 📋 TASK 10: MyelomaDigitalTwin Integration

### Q10.1: Integration Scope ✅ ANSWERED
**Question:** How much to integrate?

**ANSWER:** ✅ **Option A: Minimal integration**
- Add API call to fetch resistance prediction
- Display resistance section in existing layout
- No new tabs, no major UI changes

---

## ✅ ALL QUESTIONS ANSWERED - EXECUTION PLAN

### Build Order (DRY Architecture):

| # | Task | Time | Notes |
|---|------|------|-------|
| 1 | `ResistancePlaybookService` (shared) | 3-4h | Disease-agnostic, OV+MM mappings |
| 2 | Add cytogenetics to `resistance_prophet_service.py` | 2h | Literature-based, graceful fallback |
| 3 | Add treatment line context | 1-2h | Expert opinion multipliers |
| 4 | Define downstream handoff contracts | 2h | Drug, Care Plan, Monitoring |
| 5 | Create `/api/resistance/predict` endpoint | 1-2h | Disease-agnostic |
| 6 | Add NFE2L2/XBP1/IRE1 genes | 1h | Literature-based |
| 7 | Frontend `ResistancePanel.jsx` | 3-4h | Reusable component |
| 8 | MyelomaDigitalTwin integration | 1h | Minimal integration |

**Total: ~14-18 hours**

### Skip (Per Answers):
- ❌ Task 8 (del(17p) validation) - Use literature value instead
- ❌ Task 9 (Evo2 delta) - Not needed

---

**Status:** ✅ **ALL TASKS COMPLETE - IMPLEMENTATION FINISHED**

---

## ✅ IMPLEMENTATION STATUS (January 28, 2025)

### All 8 Tasks Completed:

| # | Task | Status | File Created/Modified |
|---|------|--------|----------------------|
| 1 | `ResistancePlaybookService` (shared) | ✅ DONE | `api/services/resistance_playbook_service.py` |
| 2 | Add cytogenetics to `resistance_prophet_service.py` | ✅ DONE | `api/services/resistance_prophet_service.py` |
| 3 | Add treatment line context | ✅ DONE | `api/services/resistance_prophet_service.py` |
| 4 | Define downstream handoff contracts | ✅ DONE | `api/services/resistance_playbook_service.py` |
| 5 | Create `/api/resistance/predict` endpoint | ✅ DONE | `api/routers/resistance.py` |
| 6 | Add NFE2L2/XBP1/IRE1 genes | ✅ DONE | `api/services/resistance_prophet_service.py` |
| 7 | Frontend `ResistancePanel.jsx` | ✅ DONE | `src/components/myeloma/ResistancePanel.jsx` |
| 8 | MyelomaDigitalTwin integration | ✅ DONE | `src/pages/MyelomaDigitalTwin.jsx` |

### Architecture Delivered: ✅ DRY (Option C)

- ✅ **Single shared service** for OV + MM
- ✅ **Disease-agnostic API** endpoint
- ✅ **Reusable frontend component**
- ✅ **All answers implemented as specified**

### Key Features Delivered:

1. **ResistancePlaybookService**: 
   - MM playbook: DIS3, TP53, CRBN, PSMB5, NFE2L2, XBP1, IRE1
   - OV playbook: NF1, KRAS, PIK3CA, PTEN, BRCA1/2
   - Cytogenetics: del(17p), t(4;14), 1q gain, t(11;14)

2. **Treatment Line Context**:
   - 1st line: 1.0x multiplier
   - 2nd line: 1.2x multiplier
   - 3rd+ line: 1.4x multiplier
   - Cross-resistance: 1.3x for same-class prior exposure

3. **Downstream Handoffs**:
   - Drug Efficacy Agent: Re-rank drugs
   - Care Plan Agent: Update regimen
   - Monitoring Agent: Intensify monitoring

4. **API Endpoint**: `POST /api/resistance/predict`
   - Accepts `disease`, `mutations`, `cytogenetics`, `treatment_line`
   - Returns alternatives, regimen changes, monitoring, handoffs

5. **Frontend Component**: `ResistancePanel.jsx`
   - Auto-fetches on mutation change
   - Displays risk level, alternatives, monitoring
   - Integrated into MyelomaDigitalTwin

---

**All questions answered. All tasks implemented. Ready for production.**

