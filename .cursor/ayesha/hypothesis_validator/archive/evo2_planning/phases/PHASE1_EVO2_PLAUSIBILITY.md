# ⚔️ PHASE 1: EVO2 BIOLOGICAL PLAUSIBILITY SERVICE

**Time Estimate:** 3 hours  
**Dependencies:** None (can start immediately)  
**Status:** ⚠️ Pending decisions from [`questions/EXECUTION_DECISIONS.md`](../questions/EXECUTION_DECISIONS.md)

---

## **🎯 GOAL**

Build Evo2 service to score if a food compound can biologically modulate target genes in the patient's disease context.

**Core Innovation:** Use Evo2's sequence-level understanding to predict biological plausibility (baseline vs. intervention).

---

## **📁 FILES TO CREATE**

- `oncology-coPilot/oncology-backend-minimal/api/services/evo2_food_plausibility.py` (300+ lines)

---

## **🔧 IMPLEMENTATION**

See [`../services/evo2_food_plausibility_service.md`](../services/evo2_food_plausibility_service.md) for complete implementation code.

### **Key Components:**

1. **`Evo2FoodPlausibilityService`** class
   - `compute_biological_plausibility()` - Main method
   - `_score_target_gene()` - Score single target (baseline vs. intervention)
   - `_fetch_gene_sequence()` - Get gene sequence from Ensembl
   - `_call_evo2_score()` - Call Evo2 API with prompt
   - `_compute_overall_score()` - Aggregate target scores
   - `_classify_verdict()` - HIGH/MODERATE/LOW classification

2. **Workflow:**
   ```
   For each target gene:
   ├─ Fetch gene sequence (Ensembl REST)
   ├─ Build disease-active context prompt
   ├─ Call Evo2: baseline_score = await evo2.score(disease_prompt)
   ├─ Build intervention context prompt
   ├─ Call Evo2: intervention_score = await evo2.score(intervention_prompt)
   ├─ Compute delta = |baseline - intervention|
   └─ Classify plausibility (HIGH if delta > 0.5, MODERATE if > 0.2, LOW otherwise)
   ```

3. **Output Schema:**
   ```python
   {
       "overall_plausibility": 0.65,  # 0-1 score
       "verdict": "MODERATE",  # HIGH/MODERATE/LOW
       "target_analysis": [
           {
               "gene": "NFKB1",
               "baseline_activity": 0.85,
               "intervention_activity": 0.72,
               "delta": 0.13,
               "plausibility": "LOW",
               "mechanism": "NF-κB inhibition",
               "rationale": "...",
               "confidence": 0.55
           }
       ],
       "mechanisms_validated": ["NF-κB inhibition"],
       "pathway_alignment": {"aligned": [...], "misaligned": [...]},
       "provenance": {...}
   }
   ```

---

## **❓ DECISIONS NEEDED**

**Before starting, review:**
- [`../questions/EXECUTION_DECISIONS.md`](../questions/EXECUTION_DECISIONS.md) - Q1 (Evo2 API approach), Q2 (Gene sequence fetching)

### **Critical Questions:**
1. **Q1:** Can `/api/evo/score` accept text prompts, or do we need `/api/evo/generate`?
2. **Q2:** Use Ensembl REST directly or call `/api/genomic_intel/gene_info` (if it exists)?

### **Recommended Approach:**
- **Q1:** Use `/api/evo/score` with sequence prompts (need verification)
- **Q2:** Use Ensembl REST directly (proven pattern from `evo.py`)

---

## **🧪 TESTING**

### **Test Case 1: Vitamin D + TP53**
```bash
curl -X POST http://127.0.0.1:8000/api/test/evo2_plausibility \
  -H "Content-Type: application/json" \
  -d '{
    "compound": "Vitamin D",
    "targets": ["VDR", "TP53"],
    "disease_context": {
      "disease": "ovarian_cancer_hgs",
      "mutations": [{"gene": "TP53", "hgvs_p": "R248Q"}],
      "biomarkers": {"HRD": "POSITIVE"},
      "pathways_disrupted": ["DNA repair", "Cell cycle"]
    },
    "pathways": ["VDR signaling", "p53 pathway"]
  }'
```

### **Expected Response:**
```json
{
  "overall_plausibility": 0.45,
  "verdict": "MODERATE",
  "target_analysis": [
    {
      "gene": "VDR",
      "delta": 0.35,
      "plausibility": "MODERATE",
      "mechanism": "VDR modulation by Vitamin D"
    },
    {
      "gene": "TP53",
      "delta": 0.25,
      "plausibility": "MODERATE",
      "rationale": "May partially restore p53 function via VDR"
    }
  ],
  "mechanisms_validated": ["VDR signaling", "p53 pathway rescue"]
}
```

---

## **✅ ACCEPTANCE CRITERIA**

- [ ] Service successfully calls Evo2 API with prompts
- [ ] Gene sequences fetched from Ensembl REST
- [ ] Baseline and intervention scores computed correctly
- [ ] Delta calculation produces meaningful plausibility scores
- [ ] Overall score aggregation weights HIGH/MODERATE/LOW appropriately
- [ ] Pathway alignment check works (compound vs. disease pathways)
- [ ] Error handling: graceful degradation when Evo2/Ensembl unavailable
- [ ] Caching: gene sequences cached to avoid repeated Ensembl calls
- [ ] Provenance tracking: full audit trail included

---

## **📝 NOTES**

- **Targets:** For Phase 1, use hardcoded targets from `food_targets.json` (see Q5 in questions doc)
- **Evo2 Model:** Default to `evo2_1b` (fastest, most stable)
- **Timeout:** 60s for Evo2 calls, 30s for Ensembl calls
- **Fallback:** Return neutral score (0.5) if Evo2 unavailable

---

**Next Phase:** [`PHASE3_DYNAMIC_DISCOVERY.md`](./PHASE3_DYNAMIC_DISCOVERY.md) (needed for Phase 2)  
**OR:** [`PHASE2_SPE_SAE_INTEGRATION.md`](./PHASE2_SPE_SAE_INTEGRATION.md) (if using hardcoded targets)

