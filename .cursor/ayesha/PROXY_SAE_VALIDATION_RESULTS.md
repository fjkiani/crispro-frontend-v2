# 💀 PROXY SAE VALIDATION RESULTS

**Date:** January 25, 2026  
**Status:** 🎯 **ANALYSIS COMPLETE**

---

## 📋 SUMMARY: WHAT PROXY SAE CAN AND CANNOT DO

### ❌ WHAT DOESN'T WORK (Validated Today):

| Experiment | Data | Result | Status |
|------------|------|--------|--------|
| **DDR → Platinum Resistance** | TCGA-OV n=469 | AUC **0.537** | ❌ Failed |
| **Composite → Platinum Resistance** | TCGA-OV n=469 | AUC **0.587** | ⚠️ Weak |
| **TP53 → Platinum Resistance** | TCGA-OV n=469 | AUC **0.527** | ❌ Failed |
| **BRCA1/2 → Platinum Resistance** | TCGA-OV n=469 | AUC **0.503** | ❌ Random |
| **DDR → Overall Survival** | TCGA-OV n=203 | p=0.92 | ❌ Failed |

**Why it fails:** 
- Baseline mutation profiles don't predict ACQUIRED resistance
- TCGA samples are taken BEFORE treatment
- Resistance develops DURING treatment (HR restoration, pathway escape)

---

### ✅ WHAT DOES WORK:

| Experiment | Data | Result | Status |
|------------|------|--------|--------|
| **Post-treatment DDR → PFI** | GSE165897 n=11 | ρ=-0.711, p=0.014 | ✅ Significant |
| **Post-treatment PI3K → PFI** | GSE165897 n=11 | ρ=-0.683, p=0.020 | ✅ Significant |
| **Post-treatment composite → Resistance** | GSE165897 n=11 | AUC 0.714 | ✅ Fair |
| **MFAP4 → Resistance** | GSE63885 n=101 | AUC **0.763** | ✅ **Best result** |
| **TIL+Exhaustion → IO Response** | GSE168204 n=11 | AUC **0.714** | ✅ Validated |

---

## 🎯 EXPERIMENT 1: COMBINATION THERAPY VALIDATION

### Status: ⏸️ **BLOCKED - Need Treatment Data**

**The Problem:** TCGA-OV doesn't have detailed treatment annotations (which specific drugs each patient received). We can't test:

> "Do high-DDR + high-VEGF patients who received PARP+Bevacizumab have better outcomes than PARP alone?"

**What Would Be Needed:**
1. **Clinical trial data** with treatment arms (e.g., PAOLA-1)
2. **TCGA clinical supplement** with drug-level treatment data
3. **Real-world data** with treatment regimens

### Alternative Approach: Literature-Based Validation ✅

The combinations ARE validated - just not with PROXY SAE pathway scores:

| Combination | Trial | Result | Status |
|-------------|-------|--------|--------|
| **PARP + Bevacizumab** | PAOLA-1 | HR 0.59, PFS 22.1 vs 16.6 mo | **FDA Approved** |
| **PARP + IO** | TOPACIO | ORR 18% (higher in biomarker+) | Published |
| **PARP + ATR** | CAPRI | Synergy in BRCA-mutant | Phase I/II |

**Conclusion:** We can recommend these combinations based on literature, not our own validation.

---

## 🔬 EXPERIMENT 2: RESISTANCE LEAD TIME

### Status: ⚠️ **POC Complete (n=11), External Validation Pending**

**Current Results (GSE165897, n=11):**

| Feature | ρ vs PFI | p-value | AUC |
|---------|----------|---------|-----|
| **post_ddr** | **-0.711** | **0.014** | 0.714 |
| **post_pi3k** | **-0.683** | **0.020** | 0.750 |
| delta_ddr | <0.3 | NS | — |

**Key Finding:**
> "Post-treatment pathway state (what survives) predicts resistance, not baseline or delta (how it changed)."

**What's Needed for Full Validation:**
1. **BriTROC-1** (n=276 paired) - EGA access
2. **MSK-SPECTRUM** (n=57 paired) - dbGaP access
3. **Compute pathway scores at both timepoints**
4. **Validate AUC > 0.70 with confidence intervals**

---

## 📊 COMPREHENSIVE VALIDATION SUMMARY

### Mutation-Based Biomarkers (PROXY SAE):

| Biomarker | Target | Data | Result | Status |
|-----------|--------|------|--------|--------|
| DDR pathway burden | Resistance | TCGA-OV n=469 | AUC 0.537 | ❌ Failed |
| Composite pathway | Resistance | TCGA-OV n=469 | AUC 0.587 | ⚠️ Weak |
| TP53 mutation | Resistance | TCGA-OV n=469 | AUC 0.527 | ❌ Failed |
| BRCA mutation | Resistance | TCGA-OV n=469 | AUC 0.503 | ❌ Random |
| DDR pathway burden | Survival | TCGA-OV n=203 | p=0.92 | ❌ Failed |

### Expression-Based Biomarkers:

| Biomarker | Target | Data | Result | Status |
|-----------|--------|------|--------|--------|
| **MFAP4** | Resistance | GSE63885 n=101 | **AUC 0.763** | ✅ Best |
| EMT score | Resistance | GSE63885 n=101 | AUC 0.715 (CV) | ✅ Good |
| TIL+Exhaustion | IO Response | GSE168204 n=11 | AUC 0.714 | ✅ Good |

### Serial/Temporal Biomarkers:

| Biomarker | Target | Data | Result | Status |
|-----------|--------|------|--------|--------|
| Post-treatment DDR | Resistance | GSE165897 n=11 | ρ=-0.711, p=0.014 | ✅ POC |
| Post-treatment PI3K | Resistance | GSE165897 n=11 | ρ=-0.683, AUC 0.75 | ✅ POC |

---

## 🎯 WHAT PROXY SAE IS ACTUALLY GOOD FOR

### ✅ Use Cases That Work:

1. **Mechanism-Aware Trial Matching**
   - Match BRCA+ patients to PARP inhibitor trials
   - Match KRAS+ patients to MEK/RAF trials
   - Match HRD-high patients to platinum trials
   - **This is valid because it's biology-based, not prediction-based**

2. **Drug Combination Recommendations**
   - PARP + Bevacizumab for HRD-high + ascites
   - PARP + IO for DDR-high + MSI-H
   - **Based on literature (PAOLA-1, TOPACIO), not our data**

3. **Pathway Profiling for Clinical Decision Support**
   - Display pathway burdens for clinician awareness
   - Highlight mechanism of action alignment
   - **Informational, not predictive**

### ❌ Use Cases That Don't Work:

1. **Predicting platinum resistance from baseline mutations**
2. **Predicting survival from baseline mutations**
3. **Predicting who will benefit from specific drugs**

---

## 📋 RECOMMENDATIONS

### For Publications:

| Paper | Use Case | Biomarker | Status |
|-------|----------|-----------|--------|
| **SAE Resistance (03)** | Resistance prediction | DDR_bin | ❌ Pivot to survival or shelve |
| **Serial SAE** | Early resistance detection | Post-treatment DDR | ⚠️ POC, needs validation |
| **MFAP4 (07)** | Resistance prediction | MFAP4 expression | ✅ Ready to publish |
| **IO Response (06)** | IO prediction | TIL+Exhaustion | ✅ Ready to publish |

### For Production:

| Feature | Status | Recommendation |
|---------|--------|----------------|
| Pathway burden display | ✅ Keep | Informational |
| Trial matching by mechanism | ✅ Keep | Biology-based |
| Resistance prediction | ❌ Remove | No evidence |
| Combination recommendations | ✅ Keep | Literature-backed |

### For Validation:

1. **Submit EGA request for BriTROC-1** (priority for Serial SAE)
2. **Submit dbGaP request for MSK-SPECTRUM** (priority for Serial SAE)
3. **Focus on expression-based biomarkers** (MFAP4, TIL) instead of mutation-based

---

## 💀 BOTTOM LINE

**PROXY SAE (mutation-based) does NOT predict platinum resistance from baseline samples.**

**What works:**
- Expression-based markers (MFAP4: AUC 0.763)
- Post-treatment pathway state (DDR: ρ=-0.711)
- Literature-backed combination recommendations

**What to do:**
1. ✅ Publish MFAP4 paper (expression-based, validated)
2. ⚠️ Validate Serial SAE with BriTROC-1/MSK-SPECTRUM (promising, n=11 too small)
3. ✅ Keep PROXY SAE for trial matching (mechanism-based, not prediction)
4. ❌ Stop claiming mutation-based resistance prediction (no evidence)
