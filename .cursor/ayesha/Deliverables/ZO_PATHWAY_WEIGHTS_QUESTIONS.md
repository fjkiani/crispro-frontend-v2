# ⚔️ PATHWAY WEIGHTS - QUESTIONS FOR MANAGER

**Author**: Zo (Lead AI Agent)  
**Date**: January 14, 2025  
**Purpose**: Questions about biological rationale and validation of pathway weights

---

## 🎯 CONTEXT

**Current State**: Pathway weights are hardcoded in:
- `api/services/pathway/drug_mapping.py` (gene→pathway weights)
- `api/services/pathway/panel_config.py` (drug→pathway weights)

**Gap**: I understand the structure and implementation, but not the biological rationale or validation methodology.

---

## 📋 CURRENT IMPLEMENTATION

### **Gene→Pathway Weights** (`drug_mapping.py:16-32`)

**Current Mapping**:
```python
# MAPK drivers → ras_mapk pathway
BRAF, KRAS, NRAS, MAP2K1, MAPK1 → {"ras_mapk": 1.0}

# DNA repair/HR genes → tp53 pathway bucket
TP53, MDM2, ATM, ATR, CHEK2, BRCA1, BRCA2, PTEN, RAD51 → {"tp53": 1.0}
```

**Notes**:
- Binary mapping (1.0 for match, 0.0 for no match)
- Comment says: "Minimal MM mapping to unblock pathway signal"
- Comment says: "Extend as needed per disease"

---

### **Drug→Pathway Weights** (`panel_config.py:8-14`)

**Current Weights** (DEFAULT_MM_PANEL):
```python
{
    "BRAF inhibitor": {"ras_mapk": 0.8, "tp53": 0.2},
    "MEK inhibitor": {"ras_mapk": 0.9, "tp53": 0.1},
    "IMiD": {"ras_mapk": 0.2, "tp53": 0.3},
    "Proteasome inhibitor": {"ras_mapk": 0.3, "tp53": 0.4},
    "Anti-CD38": {"ras_mapk": 0.1, "tp53": 0.1}
}
```

**Notes**:
- Weights sum to 1.0 (or close to it)
- Split between primary pathway (0.8-0.9) and secondary pathway (0.1-0.2)
- Disease-specific panel (MM = Multiple Myeloma)

---

## ❓ QUESTIONS FOR MANAGER

### **Question 1: Biological Rationale for Weight Splits**

**What I See**: 
- BRAF inhibitor: 0.8 ras_mapk + 0.2 tp53
- MEK inhibitor: 0.9 ras_mapk + 0.1 tp53

**My Question**:
1. **What is the biological rationale for the 0.8/0.2 split?** (not 0.9/0.1 or 0.7/0.3)
2. **Why does MEK inhibitor have higher ras_mapk weight (0.9) than BRAF inhibitor (0.8)?**
3. **Are these weights based on:**
   - Literature evidence (mechanism of action studies)?
   - Expert opinion (oncologist consensus)?
   - Empirical data (clinical trial outcomes)?
   - First principles (pathway biology)?

**Why This Matters**: 
- Need to justify weights in clinical documentation
- Need to know if weights are evidence-based or heuristic
- Need to understand if weights can be validated/refined

---

### **Question 2: Validation Methodology**

**What I See**: 
- Weights are hardcoded (no validation mentioned)

**My Question**:
1. **Have these weights been validated against clinical outcomes?**
   - Do patients with high pathway alignment (0.8+) have better responses?
   - Do patients with low pathway alignment (<0.5) have worse responses?
2. **If validated, what was the methodology?**
   - Retrospective cohort analysis?
   - Clinical trial data?
   - Literature meta-analysis?
3. **If not validated, should we validate them?**
   - What would validation look like?
   - What metrics would we use? (response rate, PFS, OS?)

**Why This Matters**:
- Clinical credibility requires validation
- Need to know if weights are evidence-based or placeholder
- Need to plan validation studies if missing

---

### **Question 3: Disease-Specificity**

**What I See**: 
- DEFAULT_MM_PANEL (Multiple Myeloma specific)
- Comment says: "Extend as needed per disease"

**My Question**:
1. **Are these weights disease-specific?**
   - Do BRAF inhibitor weights differ for MM vs ovarian vs melanoma?
   - Should we have separate panels per disease?
2. **If disease-specific, what's the rationale?**
   - Different pathway importance per disease?
   - Different drug mechanisms per disease?
3. **If not disease-specific, why?**
   - Are pathways universal?
   - Are weights generalizable?

**Why This Matters**:
- Ayesha has ovarian cancer (not MM)
- Need to know if MM weights apply to ovarian
- Need to know if we should create ovarian-specific panel

---

### **Question 4: Weight Determination Process**

**What I See**: 
- Weights appear to be manually set (hardcoded)

**My Question**:
1. **How were these weights determined?**
   - Expert consensus?
   - Literature review?
   - Data-driven (regression, ML)?
   - First principles (pathway biology)?
2. **Who determined them?**
   - Oncologist input?
   - Computational biologist?
   - Literature synthesis?
3. **Are they documented anywhere?**
   - Internal documentation?
   - Publication?
   - Clinical guidelines?

**Why This Matters**:
- Need to understand decision-making process
- Need to know if weights are reproducible
- Need to know if weights can be updated/refined

---

### **Question 5: Future Enhancement Path**

**What I See**: 
- Weights are hardcoded (not configurable)

**My Question**:
1. **Should we make weights configurable per disease?**
   - Create `ovarian_cancer_panel.json`?
   - Create `breast_cancer_panel.json`?
2. **Should we validate weights against outcomes?**
   - Retrospective validation study?
   - Prospective validation in trials?
3. **Should we refine weights based on feedback?**
   - Update based on clinical outcomes?
   - Update based on new literature?

**Why This Matters**:
- Need to plan future enhancements
- Need to know if weights are static or evolving
- Need to know if we should invest in validation

---

## ✅ PROPOSED ANSWERS FOR AGENT

### Biological rationale (0.8/0.2 vs 0.9/0.1)
- Primary pathway receives the majority weight to reflect direct MoA; a secondary weight captures cross‑talk via hub pathways (e.g., TP53/DDR).
- MEK inhibitor 0.9 ras_mapk vs BRAF inhibitor 0.8: MEK is downstream and more specific; BRAF inhibitors can show paradoxical MAPK activation and broader context effects, warranting a slightly lower specificity weight.

### Validation methodology
- Retrospective outcome correlation:
  - Compute per‑drug pathway alignment (s_path) and correlate with response (ORR/PFS/OS) using Spearman ρ and AUROC (responder vs non‑responder).
  - Within‑patient ranking check: did the top‑aligned therapy perform best?
- Acceptance thresholds: AUROC ≥0.65; top‑1 correct ≥55%; top‑3 contains best ≥80%.

### Disease specificity
- Yes. Maintain per‑disease panels (e.g., ovarian: heavier PARP/DDR weights; MM: MAPK/IMiD/PI).
- Default to MM panel only when a disease‑specific panel is missing.

### Determination process (pragmatic policy)
- Initialize from first principles (MoA→pathways), literature consensus, and expert priors.
- Refine via retrospective data; bounded updates (±0.05 per iteration), keeping weights normalized and audited.

### Future enhancement
- Config‑driven JSON per disease; audit trail of changes.
- Automated calibration job proposes small deltas with guardrails; human approval required before applying.

---

## 📊 CURRENT STATE DOCUMENTATION

**File**: `PATHWAY_WEIGHTS_RATIONALE.md` (to be created)

**Will Document**:
- Current weights for all drugs in DEFAULT_MM_PANEL
- Current gene→pathway mappings
- Usage patterns (where weights are applied)
- Integration points (pathway aggregation, drug scoring)

---

## 🎯 EXPECTED OUTCOMES

**After Manager Answers**:
1. ✅ Understand biological rationale for weight splits
2. ✅ Know if weights are validated or need validation
3. ✅ Know if weights are disease-specific or universal
4. ✅ Understand weight determination process
5. ✅ Plan future enhancements (validation, disease-specific panels)

---

**Status**: ✅ **ANSWERS PROVIDED**  
**Last Updated**: January 14, 2025  
**By**: Zo (Lead AI Agent)  
**Answers Added**: January 14, 2025

