# 🎯 OVARIAN CANCER PATIENT - 1 WEEK BEFORE L1 CHEMO

**Scenario:** Ovarian cancer (HGS) patient, 1 week before starting first-line chemotherapy  
**Treatment Line:** L1 (first-line, no prior treatment)  
**Timing:** Prophylactic (before treatment starts)

---

## 📊 SYSTEM RECOMMENDATIONS

### Top Recommendations (High Priority):

| Supplement | Final Score | Line App | Rationale |
|------------|-------------|----------|-----------|
| **Vitamin D** | **0.80** | 0.90 | DNA repair support, immune modulation, appropriate for L1, ovarian cancer boost |
| **Folate** | **0.80** | 0.80 | DNA repair support, critical for HRD+ context, appropriate for L1 |
| **Omega-3** | **0.80** | 0.70 | Anti-inflammatory, immune support, appropriate for L1, ovarian cancer boost |

### Moderate Recommendations:

| Supplement | Final Score | Line App | Rationale |
|------------|-------------|----------|-----------|
| **NAC** | **0.75** | 0.70 | Recovery supplement, LOWER in L1 (no prior treatment), but still gets cancer boost |
| **Green Tea (EGCG)** | **0.70** | 0.65 | Anti-angiogenic, immune support, appropriate for L1 |

### Lower Priority (Not Recommended for L1):

| Supplement | Final Score | Line App | Rationale |
|------------|-------------|----------|-----------|
| **L-glutamine** | **0.65** | 0.65 | Recovery supplement, NOT appropriate for L1 (no prior GI toxicity) |
| **Alpha-lipoic acid** | **0.65** | 0.65 | Neuropathy prevention, NOT needed before treatment starts |

---

## 🔍 HOW THE SYSTEM WORKS

### 1. **Treatment Line Logic (Phase 1)**
- **L1 (Before Treatment):** Prophylactic supplements prioritized
  - Recovery supplements (NAC, L-glutamine) get **reduced scores** (0.65-0.70)
  - General supplements (Vitamin D, Folate) get **appropriate scores** (0.80-0.90)

### 2. **Cancer Type Boost (Phase 2)**
- Ovarian cancer foods get **+0.1 base boost**
- If treatment line matches (L1 in recommended lines), get **+0.05 extra boost**
- **Total boost: +0.15** for matched foods

### 3. **Biomarker Boost (Phase 3)**
- If HRD+ is known: DNA repair foods (Vitamin D, Folate, NAC) get **+0.1 boost**
- If TMB-H: Immune foods (Omega-3, Vitamin D) get **+0.1 boost**
- **Note:** Biomarker status may not be known 1 week before chemo

### 4. **Line-Specific Scores (Phase 4)**
- For supplements with `line_specific_scores` in JSON:
  - **NAC L1:** 0.70 (prophylactic support, but no prior treatment to recover from)
  - **Vitamin D L1:** 0.90 (appropriate across all lines)
  - **L-glutamine L1:** 0.65 (not needed without prior GI toxicity)

---

## 💡 KEY INSIGHTS

### ✅ **What Works Well:**
1. **Differentiates L1 vs L2/L3:** Recovery supplements correctly get lower scores in L1
2. **Cancer type awareness:** Ovarian cancer-specific foods get boosted
3. **Prophylactic focus:** General supplements (Vitamin D, Folate) prioritized for L1

### ⚠️ **Potential Gaps:**
1. **Timing specificity:** "1 week before" isn't explicitly modeled - system treats L1 as "before treatment"
2. **Biomarker uncertainty:** If HRD status unknown, misses biomarker boost
3. **Drug-specific recommendations:** Doesn't know which chemo drugs yet (carboplatin vs. other)

---

## 🎯 ACTUAL RECOMMENDATIONS FOR PATIENT

### **Start Now (1 Week Before):**
1. **Vitamin D:** 2000-4000 IU daily
   - Rationale: DNA repair support, immune modulation, appropriate for L1
   - Score: 0.80 (high confidence)

2. **Folate (5-MTHF):** 400-800mcg daily
   - Rationale: DNA repair support, especially if HRD+ (check biomarker status)
   - Score: 0.80 (high confidence)

3. **Omega-3 (EPA+DHA):** 2-3g daily
   - Rationale: Anti-inflammatory, immune support, reduces chronic inflammation
   - Score: 0.80 (high confidence)

### **Consider (Moderate Priority):**
4. **NAC:** 600mg twice daily
   - Rationale: Prophylactic antioxidant support, but recovery focus means lower L1 score
   - Score: 0.75 (moderate confidence)
   - **Note:** May be more important if starting platinum-based chemo

### **Wait Until After Treatment:**
5. **L-glutamine:** NOT recommended for L1
   - Rationale: GI mucosal protection not needed before treatment starts
   - Score: 0.65 (low appropriateness for L1)
   - **When to start:** After first cycle if GI toxicity develops

6. **Alpha-lipoic acid:** NOT recommended for L1
   - Rationale: Neuropathy prevention not needed before treatment
   - Score: 0.65 (low appropriateness for L1)
   - **When to start:** After first cycle if neuropathy develops

---

## 📋 SYSTEM OUTPUT EXAMPLE

```json
{
  "compound": "Vitamin D",
  "alignment_score": 0.80,
  "overall_score": 0.80,
  "confidence": 0.85,
  "verdict": "SUPPORTED",
  "spe_breakdown": {
    "sequence": 0.5,
    "pathway": 0.75,
    "evidence": 0.70
  },
  "sae_features": {
    "line_fitness": {
      "score": 0.90,
      "status": "appropriate",
      "reason": "Immune modulation, DNA repair support - appropriate across all lines"
    },
    "cross_resistance": {
      "risk": "LOW",
      "score": 0.0,
      "reason": "No prior therapies to assess cross-resistance"
    },
    "sequencing_fitness": {
      "score": 0.85,
      "optimal": true,
      "reason": "Appropriate for L1 prophylactic support"
    }
  },
  "boosts": {
    "base_score": 0.65,
    "cancer_type_boost": 0.15,
    "biomarker_boost": 0.0,
    "total_boost": 0.15,
    "final_score": 0.80,
    "boost_reasons": [
      "Cancer type match (ovarian_cancer_hgs)",
      "Treatment line match (L1)"
    ]
  },
  "dietician_recommendations": {
    "dosage": "2000-4000 IU daily",
    "timing": "continuous, with fatty meal",
    "safety": "Generally safe, monitor calcium levels"
  }
}
```

---

## ✅ BOTTOM LINE

**For an ovarian cancer patient 1 week before L1 chemo, the system would recommend:**

1. ✅ **Vitamin D** (0.80) - Start now
2. ✅ **Folate** (0.80) - Start now  
3. ✅ **Omega-3** (0.80) - Start now
4. ⚠️ **NAC** (0.75) - Consider, but lower priority for L1
5. ❌ **L-glutamine** (0.65) - Wait until after treatment if GI toxicity develops

**The system correctly:**
- Prioritizes prophylactic supplements for L1
- Deprioritizes recovery supplements (appropriate for L2/L3)
- Applies cancer type boosts
- Differentiates treatment lines

**This is the correct behavior!** 🎯

