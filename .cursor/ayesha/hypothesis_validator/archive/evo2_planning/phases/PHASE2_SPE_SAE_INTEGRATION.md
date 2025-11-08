# ⚔️ PHASE 2: S/P/E + SAE INTEGRATION

**Time Estimate:** 2 hours  
**Dependencies:** Phase 1 (Evo2 service)  
**Status:** ⚠️ Pending decisions from [`questions/EXECUTION_DECISIONS.md`](../questions/EXECUTION_DECISIONS.md)

---

## **🎯 GOAL**

Integrate food validator with full S/P/E framework + SAE features for treatment line intelligence.

**Core Innovation:** Multi-modal scoring (Sequence + Pathway + Evidence) with treatment line context.

---

## **📁 FILES TO CREATE**

- `oncology-coPilot/oncology-backend-minimal/api/services/food_spe_integration.py` (250+ lines)

---

## **🔧 IMPLEMENTATION**

See [`../services/food_spe_integration_service.md`](../services/food_spe_integration_service.md) for complete implementation code.

### **Key Components:**

1. **`FoodSPEIntegrationService`** class
   - `compute_spe_score()` - Main method combining S/P/E + SAE
   - `_compute_pathway_score()` - Pathway alignment scoring
   - `_convert_evidence_grade()` - Literature grade → 0-1 score
   - `_aggregate_spe()` - S/P/E weighted aggregation
   - `_compute_confidence()` - Multi-stage confidence modulation
   - `_classify_verdict()` - SUPPORTED/WEAK_SUPPORT/NOT_SUPPORTED

2. **Workflow:**
   ```
   [1] SEQUENCE (S): Evo2 biological plausibility
       → Call evo2_service.compute_biological_plausibility()
       → Extract overall_plausibility as sequence_score
   
   [2] PATHWAY (P): Compound pathways vs. disease pathways
       → Check alignment from evo2_result['pathway_alignment']
       → Score: aligned pathways = 1.0, misaligned = 0.2
   
   [3] EVIDENCE (E): Literature grade conversion
       → Convert STRONG/MODERATE/WEAK to 0.9/0.6/0.3
   
   [4] SAE: Treatment line features (if treatment_history provided)
       → Call compute_treatment_line_features() or supplement-specific logic
       → Extract line_appropriateness, cross_resistance, sequencing_fitness
   
   [5] Aggregate S/P/E:
       → overall_score = 0.4×S + 0.3×P + 0.3×E
   
   [6] Compute confidence with modulation:
       → base_confidence = (S+P+E)/3
       → evo2_boost = +0.05 per HIGH plausibility target
       → sae_boost = (line_app + seq_fit) × 0.05
       → biomarker_boost = +0.05 if HRD+ matches DNA repair pathway
       → final_confidence = min(base + boosts, 0.95)
   
   [7] Classify verdict:
       → SUPPORTED if score ≥0.65 AND confidence ≥0.70
       → WEAK_SUPPORT if score ≥0.45 AND confidence ≥0.50
       → NOT_SUPPORTED otherwise
   ```

3. **Output Schema:**
   ```python
   {
       "overall_score": 0.72,
       "confidence": 0.68,
       "verdict": "SUPPORTED",
       "spe_breakdown": {
           "sequence_score": 0.45,
           "pathway_score": 0.85,
           "evidence_score": 0.60
       },
       "sae_features": {
           "line_appropriateness": 1.0,
           "cross_resistance": 0.0,
           "sequencing_fitness": 0.9
       },
       "confidence_modulation": {
           "base_confidence": 0.60,
           "evo2_boost": 0.10,
           "sae_boost": 0.05,
           "biomarker_boost": 0.05,
           "final_confidence": 0.80
       },
       "evo2_analysis": {...},
       "provenance": {...}
   }
   ```

---

## **❓ DECISIONS NEEDED**

**Before starting, review:**
- [`../questions/EXECUTION_DECISIONS.md`](../questions/EXECUTION_DECISIONS.md) - Q3 (SAE for supplements), Q4 (S/P/E approach), Q6 (Treatment line logic)

### **Critical Questions:**
1. **Q3:** How do we compute treatment line features for supplements? (generic vs. supplement-specific)
2. **Q4:** Reuse efficacy_orchestrator or build new service?
3. **Q6:** What supplement-specific treatment line rules should we use?

### **Recommended Approach:**
- **Q3/Q6:** Build supplement-specific logic (NAC post-platinum, Vitamin D HRD+, etc.)
- **Q4:** Build new `FoodSPEIntegrationService` (simpler, purpose-built)

---

## **🔗 INTEGRATIONS**

### **Dependencies:**
- `api/services/evo2_food_plausibility.py` (Phase 1)
- `api/services/treatment_line_integration.py` (existing) - OR new supplement logic
- `api/services/llm_literature_service.py` (existing) - for evidence grade

### **Import Pattern:**
```python
from api.services.evo2_food_plausibility import get_evo2_food_plausibility_service
from api.services.treatment_line_integration import compute_treatment_line_features  # OR supplement version
from api.services.llm_literature_service import get_llm_service
```

---

## **🧪 TESTING**

### **Test Case: Vitamin D + S/P/E + SAE**
```bash
curl -X POST http://127.0.0.1:8000/api/test/food_spe \
  -H "Content-Type: application/json" \
  -d '{
    "compound": "Vitamin D",
    "targets": ["VDR", "TP53"],
    "pathways": ["VDR signaling", "p53 pathway"],
    "disease_context": {
      "disease": "ovarian_cancer_hgs",
      "mutations": [{"gene": "TP53", "hgvs_p": "R248Q"}],
      "biomarkers": {"HRD": "POSITIVE"},
      "pathways_disrupted": ["DNA repair", "Cell cycle"]
    },
    "evidence_grade": "MODERATE",
    "treatment_history": {
      "current_line": 3,
      "prior_therapies": ["carboplatin", "paclitaxel"]
    }
  }'
```

### **Expected Response:**
```json
{
  "overall_score": 0.61,
  "confidence": 0.75,
  "verdict": "SUPPORTED",
  "spe_breakdown": {
    "sequence_score": 0.45,
    "pathway_score": 0.85,
    "evidence_score": 0.60
  },
  "sae_features": {
    "line_appropriateness": 0.9,
    "cross_resistance": 0.0,
    "sequencing_fitness": 0.9
  }
}
```

---

## **✅ ACCEPTANCE CRITERIA**

- [ ] Evo2 plausibility integrated as Sequence (S) component
- [ ] Pathway alignment scoring works (compound vs. disease)
- [ ] Evidence grade conversion correct (STRONG/MODERATE/WEAK → 0.9/0.6/0.3)
- [ ] S/P/E aggregation formula correct (40% S, 30% P, 30% E)
- [ ] SAE features computed (line appropriateness, cross-resistance, sequencing fitness)
- [ ] Confidence modulation applies all boosts correctly
- [ ] Verdict classification logic correct
- [ ] Provenance tracking complete

---

**Next Phase:** [`PHASE4_ENHANCED_ENDPOINT.md`](./PHASE4_ENHANCED_ENDPOINT.md) (after Phase 3)

