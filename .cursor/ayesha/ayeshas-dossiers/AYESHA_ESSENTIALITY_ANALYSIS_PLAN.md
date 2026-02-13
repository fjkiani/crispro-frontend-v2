# 🧬 AYESHA GENE ESSENTIALITY ANALYSIS PLAN

**Patient:** AK  
**Profile:** MBD4 c.1239delA (homozygous) + TP53 somatic mutation  
**Date:** January 28, 2025  
**Purpose:** Understand pathway dependencies and synthetic lethal opportunities

---

## ✅ MANAGER REVIEW COMPLETE - ALL QUESTIONS ANSWERED

**Status:** ✅ READY FOR EXECUTION

### Summary of Decisions:

| Decision | Answer |
|----------|--------|
| **Scope** | Full analysis: Pathway + Synthetic Lethality + Drug Sensitivity |
| **Approach** | Hybrid: Standalone first → then WIWFM S/P/E |
| **Context** | Disease (ovarian HGSOC) + Mutation (MBD4, TP53) |
| **Output** | Detailed breakdown → Drug-level recommendations |
| **Model** | `evo2_1b` (default, sufficient for clear LoF variants) |
| **Integration** | Part of Phase 2 (Pathway Analysis) |

### Expected Key Findings:

| Gene | Expected Essentiality | Implication |
|------|----------------------|-------------|
| **MBD4** | ≥ 0.8 (frameshift) | BER pathway non-functional |
| **TP53** | ≥ 0.7 (hotspot) | Checkpoint pathway bypassed |
| **Combined** | Synthetic lethality | PARP + ATR high sensitivity |

---

## 📋 CLARIFYING QUESTIONS FOR MANAGER (ANSWERED)

### Question 1: Analysis Scope

**What specific insights do we want from essentiality analysis?**

**Options:**
- **A)** Pathway dependency analysis (which pathways become essential after MBD4/TP53 loss)
- **B)** Synthetic lethal target identification (which genes become essential in this context)
- **C)** Drug sensitivity prediction (which drugs exploit essentiality)
- **D)** All of the above

✅ **MANAGER ANSWER: D) All of the above**

**Rationale:** Ayesha's MBD4+TP53 is a rare combination. We need the full picture:
1. **Pathway dependencies** → Understand what's broken (BER, checkpoint)
2. **Synthetic lethal targets** → Identify druggable vulnerabilities (PARP, ATR, CHK1, WEE1)
3. **Drug sensitivity** → Validate PARP inhibitor recommendation with evidence

This comprehensive approach creates the strongest clinical case for aggressive PARP maintenance.

---

### Question 2: Integration Approach

**How should we run this analysis?**

**Options:**
- **A)** Standalone essentiality endpoint calls (`/api/insights/predict_gene_essentiality`)
  - Separate calls for MBD4 and TP53
  - Manual interpretation of results
  
- **B)** Integrated into full WIWFM S/P/E analysis (`/api/efficacy/predict`)
  - Essentiality becomes part of insights bundle
  - Automatic confidence lift (+0.07 if essentiality ≥ 0.7)
  - Pathway aggregation includes essentiality signals
  
- **C)** Hybrid approach
  - Run standalone for deep dive
  - Then integrate into WIWFM for drug recommendations

✅ **MANAGER ANSWER: C) Hybrid approach**

**Rationale:** 
1. **First: Standalone calls** for MBD4 and TP53 individually
   - Get detailed breakdown (scores, rationale, calibration, flags)
   - Document essentiality scores for clinical record
   - Understand each gene's contribution separately
   
2. **Then: Full WIWFM S/P/E** with both mutations
   - Let system aggregate pathway scores
   - Get confidence-lifted drug recommendations
   - Automatic synthetic lethality detection

This gives us both the scientific deep-dive AND the actionable drug recommendations.

---

### Question 3: Context Specification

**What context should we provide for essentiality analysis?**

From the codebase, essentiality can be context-aware. For Ayesha:

**Options:**
- **A)** Disease context: `{"disease": "ovarian_cancer", "stage": "IVB"}`
- **B)** Mutation context: `{"known_mutations": ["MBD4", "TP53"], "pathway_deficiencies": ["BER", "checkpoint"]}`
- **C)** Minimal: Just gene + variants (no context)
- **D)** Full profile: Include all available clinical data

✅ **MANAGER ANSWER: B) Mutation context + A) Disease context combined**

**Provide this context:**
```json
{
  "disease": "ovarian_cancer",
  "subtype": "high_grade_serous",
  "stage": "IVB",
  "known_mutations": ["MBD4", "TP53"],
  "pathway_deficiencies": ["BER", "checkpoint"],
  "germline_status": "MBD4_positive"
}
```

**Rationale:** 
- Disease context calibrates essentiality to ovarian cancer cell line dependencies (DepMap data)
- Mutation context identifies the "double hit" for synthetic lethality detection
- Pathway deficiencies explicitly flag BER and checkpoint loss for backup pathway analysis

---

### Question 4: Expected Output Format

**What should the output look like?**

**Options:**
- **A)** Simple scores: `{"MBD4": 0.9, "TP53": 0.85}` (essentiality scores)
- **B)** Detailed breakdown: Essentiality score + rationale + calibration + flags
- **C)** Pathway-level: Which pathways become essential (BER, DDR, checkpoint)
- **D)** Drug-level: Which drugs exploit essentiality (PARP, ATR, WEE1)

✅ **MANAGER ANSWER: B) Detailed breakdown → then D) Drug-level**

**Two-phase output:**

**Phase 1 (Essentiality Deep Dive):**
```json
{
  "MBD4": {
    "essentiality_score": 0.9,
    "rationale": "Homozygous frameshift → complete loss-of-function",
    "calibration": {"gene_percentile": 0.95, "z_score": 2.3},
    "flags": {"truncation": true, "frameshift": true, "lof": true},
    "pathway_impact": "BER pathway non-functional"
  },
  "TP53": {
    "essentiality_score": 0.85,
    "rationale": "Hotspot mutation → tumor suppressor inactivation",
    "calibration": {"gene_percentile": 0.88, "z_score": 1.9},
    "flags": {"hotspot": true, "dominant_negative": true},
    "pathway_impact": "Checkpoint pathway bypassed"
  }
}
```

**Phase 2 (Drug-Level Output):**
```json
{
  "synthetic_lethal_drugs": [
    {"drug": "Olaparib", "target": "PARP", "sensitivity": "HIGH"},
    {"drug": "Ceralasertib", "target": "ATR", "sensitivity": "HIGH"},
    {"drug": "Adavosertib", "target": "WEE1", "sensitivity": "MODERATE"}
  ],
  "backup_pathways_exploited": ["HR", "ATR/CHK1", "G2/M checkpoint"]
}
```

---

### Question 5: Synthetic Lethality Focus

**Should we focus on identifying synthetic lethal relationships?**

For Ayesha's MBD4 + TP53 combination:

**Questions:**
- Which genes become essential when BER is deficient?
- Which genes become essential when TP53 checkpoint is lost?
- What are the combined effects (double-hit essentiality)?

✅ **MANAGER ANSWER: YES - All three analyses**

**Execute in this order:**

**Step 1: MBD4 loss → Essential genes**
- BER pathway broken → HR pathway becomes essential → PARP1/PARP2 become synthetic lethal targets
- Expected: BRCA1/2, RAD51, PALB2 become more essential (backup repair)

**Step 2: TP53 loss → Essential genes**
- Checkpoint broken → ATR/CHK1/WEE1 become essential (only remaining checkpoint)
- Expected: ATR, CHK1, WEE1 become synthetic lethal targets

**Step 3: Combined (MBD4 + TP53)**
- Double-hit creates AMPLIFIED synthetic lethality
- Expected findings:
  - PARP sensitivity: VERY HIGH (BER + checkpoint both gone)
  - ATR sensitivity: HIGH (only backup checkpoint remaining)
  - WEE1 sensitivity: HIGH (G2/M checkpoint sole protector)

**This is the key insight:** MBD4+TP53 creates a "perfect storm" for PARP + ATR combination therapy.

---

### Question 6: Integration with Existing Analysis

**How does this relate to the MBD4 analysis plan?**

From `src/services/evo_service/MBD4.mdc`, we have:
- Phase 1: Variant functional annotation (Evo2 scoring)
- Phase 2: Pathway analysis (BER, HRD, DDR)
- Phase 3: Drug predictions (S/P/E framework)
- Phase 4: Clinical trial matching

**Question:** Where does essentiality fit?

**Options:**
- **A)** Part of Phase 1 (variant annotation) - understand functional impact
- **B)** Part of Phase 2 (pathway analysis) - understand pathway dependencies
- **C)** Part of Phase 3 (drug predictions) - inform drug sensitivity
- **D)** New Phase 5 - dedicated essentiality analysis

✅ **MANAGER ANSWER: B) Part of Phase 2 (Pathway Analysis)**

**Rationale:**
Essentiality analysis is fundamentally about **pathway dependencies** - which pathways become essential when others are broken. This fits naturally into Phase 2.

**Updated Phase 2 (Pathway Analysis + Essentiality):**
```
Phase 2: Pathway & Essentiality Analysis
├── 2.1 BER Pathway Assessment (MBD4 impact)
├── 2.2 Checkpoint Pathway Assessment (TP53 impact)
├── 2.3 Essentiality Scoring (MBD4, TP53 individually)
├── 2.4 Synthetic Lethal Target Identification
└── 2.5 Backup Pathway Dependency Mapping
```

**Why NOT Phase 1?** Essentiality goes beyond variant annotation - it's about cellular context and pathway dependencies.

**Why NOT Phase 3?** Drug predictions consume essentiality scores as input - essentiality must be computed first.

**Why NOT Phase 5?** Essentiality is a prerequisite for drug predictions - it shouldn't be after Phase 4.

---

## 🎯 PROPOSED ANALYSIS PLAN (Pending Manager Answers)

### Step 1: Essentiality Scoring

```python
# For MBD4
mbd4_essentiality = await client.post("/api/insights/predict_gene_essentiality", {
    "gene": "MBD4",
    "variants": [{
        "gene": "MBD4",
        "hgvs_p": "p.Ile413Serfs*2",
        "chrom": "3",
        "pos": 129430456,
        "ref": "A",
        "alt": ""
    }],
    "model_id": "evo2_1b"  # or evo2_7b?
})

# For TP53
tp53_essentiality = await client.post("/api/insights/predict_gene_essentiality", {
    "gene": "TP53",
    "variants": [{
        "gene": "TP53",
        "hgvs_p": "p.R175H",  # or actual variant
        "chrom": "17",
        "pos": 7577120,
        "ref": "G",
        "alt": "A"
    }],
    "model_id": "evo2_1b"
})
```

### Step 2: Pathway Essentiality Analysis

**Questions to answer:**
1. After MBD4 loss → Is BER pathway essential? (likely NO - pathway is broken)
2. After TP53 loss → Are checkpoint pathways essential? (likely NO - checkpoint is broken)
3. Combined → Which backup pathways become essential?

**Expected insights:**
- PARP pathway may become essential (synthetic lethality)
- ATR/CHK1 pathways may become essential (checkpoint backup)
- HR pathway may become essential (alternative DNA repair)

### Step 3: Drug Sensitivity Prediction

**Based on essentiality scores:**
- High essentiality in backup pathways → High drug sensitivity
- Low essentiality → Drug resistance risk

**Expected drugs:**
- PARP inhibitors (exploit HR pathway essentiality)
- ATR inhibitors (exploit checkpoint backup essentiality)
- WEE1 inhibitors (exploit G2/M checkpoint essentiality)

---

## 📊 WHAT WE EXPECT TO LEARN

### For MBD4 (BER Deficiency):

| Question | Expected Answer | Clinical Implication |
|----------|----------------|----------------------|
| Is BER pathway essential? | NO (pathway is broken) | Cannot repair base damage |
| Are backup repair pathways essential? | YES (HR, NER) | PARP inhibitors exploit this |
| Is MBD4 itself essential? | NO (loss-of-function) | Gene is non-functional |

### For TP53 (Checkpoint Loss):

| Question | Expected Answer | Clinical Implication |
|----------|----------------|----------------------|
| Are checkpoint pathways essential? | NO (checkpoint is broken) | Cannot stop damaged cells |
| Are backup checkpoints essential? | YES (ATR, CHK1, WEE1) | ATR/CHK1/WEE1 inhibitors exploit this |
| Is TP53 itself essential? | NO (loss-of-function) | Gene is non-functional |

### Combined (MBD4 + TP53):

| Question | Expected Answer | Clinical Implication |
|----------|----------------|----------------------|
| Which pathways become essential? | HR, ATR, CHK1, WEE1 | Multiple drug targets |
| Synthetic lethal targets? | PARP, ATR, CHK1, WEE1 | High drug sensitivity expected |
| Resistance mechanisms? | Pathway restoration | Monitor for backup pathway activation |

---

## 🔬 TECHNICAL DETAILS

### Essentiality Score Computation (from codebase):

```python
# Base score calculation
base_score = 0.2 + 0.15 * len(variants) + 0.5 * evo_magnitude

# Truncation/frameshift boost
if truncation:
    base_score = max(base_score, 0.9)
if frameshift:
    base_score = max(base_score, 0.8)

# Final score (clamped to [0,1])
essentiality_score = max(0.0, min(1.0, round(base_score, 3)))
```

**For Ayesha:**
- MBD4: Frameshift → likely score ≥ 0.8
- TP53: Hotspot mutation → likely score ≥ 0.7

### Confidence Lift (from WIWFM master):

If essentiality ≥ 0.7:
- Legacy confidence: +0.07 lift
- V2 confidence: +0.02 lift (capped at +0.08 total)

---

## ✅ MANAGER ANSWERS COMPLETE

| Question | Answer |
|----------|--------|
| 1. Scope | **D) All of the above** - Full analysis |
| 2. Approach | **C) Hybrid** - Standalone first, then WIWFM |
| 3. Context | **B+A) Mutation + Disease context** |
| 4. Output | **B→D) Detailed breakdown → Drug-level** |
| 5. Synthetic Lethality | **YES** - All three analyses |
| 6. Integration | **B) Part of Phase 2** - Pathway Analysis |
| 7. Model | **evo2_1b** (see below) |
| 8. Follow-up | **Full WIWFM S/P/E → Trial Matching** |

### Question 7: Model Selection

✅ **MANAGER ANSWER: Use `evo2_1b` (default)**

**Rationale:**
- `evo2_1b` is the default for reproducibility (matches submission package)
- MBD4 frameshift and TP53 hotspots are well-characterized - 1B sufficient
- 7B is overkill for clear loss-of-function variants
- Faster execution, lower latency

**Exception:** If essentiality scores seem unexpectedly low, retry with `evo2_7b` for deeper analysis.

### Question 8: Follow-up Steps

✅ **MANAGER ANSWER: Full pipeline after essentiality**

**Execution order:**
```
⚠️ HALLUCINATION AUDIT (Feb 10, 2026)
> **STATUS**: DEBUNKED / LEGACY ARTIFACT
> **CAUSE**: This plan executed a script (`ayesha_mbd4_tp53_hgsoc_analysis.py`) that relied on a hardcoded lookup table.
> **REALITY**: The "Hybrid Approach" described below was actually a deterministic check: `If MBD4 -> Then PARP`.
> **ACTION**: Treat as legacy.

# AYESHA: Essentiality Analysis Plan (MBD4 + TP53) ← WE ARE HERE
        ↓
Step 2: Full WIWFM S/P/E Drug Predictions (Phase 3)
        - Input: MBD4 + TP53 with essentiality scores
        - Output: Ranked drugs with confidence
        ↓
Step 3: Clinical Trial Matching (Phase 4)
        - Keywords: MBD4, TP53, PARP, ATR, synthetic lethality
        - Basket trials for DNA repair deficiency
        ↓
Step 4: Integrate into Ayesha's Care Plan
        - Update AYESHA_CLINICAL_SUMMARY_AND_MONITORING.md
        - Add essentiality scores to profile
        - Document drug recommendations with evidence
```

---

**Status:** ✅ **ALL QUESTIONS ANSWERED - READY FOR EXECUTION**

**Next Step:** Execute essentiality analysis using the parameters specified above.

