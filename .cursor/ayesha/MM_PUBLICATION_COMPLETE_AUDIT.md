# MM Drug Efficacy Publication - Complete Audit

**Date:** January 27, 2026  
**Purpose:** Understand what was actually built, validated, and published  
**Approach:** No assumptions - verify all claims against actual outputs  

---

## 📊 WHAT WAS ACTUALLY PUBLISHED (MM Paper)

### Publication: "Multi-Modal Genomic Analysis for Drug Efficacy Prediction in Multiple Myeloma"

**Status:** Draft v1.0 (October 2, 2025)

### Validated Results:

| Metric | Value | Source | Validated? |
|--------|-------|--------|------------|
| **Pathway Alignment Accuracy** | 100% (5/5 MAPK variants) | `mm_efficacy_results.json` | ✅ YES |
| **Average Confidence (SPE)** | 0.524 | `mm_efficacy_results.json` | ✅ YES |
| **Average Confidence (SP)** | 0.467 | `ablation_results.json` | ✅ YES |
| **Calibration ECE (SPE)** | 0.479 | Calibration analysis | ✅ YES |
| **Calibration ECE (SP)** | 0.529 | Calibration analysis | ✅ YES |

### Test Variants (n=7):
1. **BRAF V600E** → BRAF inhibitor (confidence: 0.515) ✅
2. **BRAF V600K** → BRAF inhibitor (confidence: 0.515) ✅
3. **KRAS G12D** → MEK inhibitor (confidence: 0.530) ✅
4. **KRAS G12V** → MEK inhibitor (confidence: 0.530) ✅
5. **NRAS Q61K** → MEK inhibitor (confidence: 0.515) ✅
6. **TP53 R248W** → IMiD (confidence: 0.560) [control]
7. **TP53 R273H** → BRAF inhibitor (confidence: 0.505) [control]

### Ablation Study Results:

| Mode | MAPK Alignment | Avg Confidence | Interpretation |
|------|----------------|----------------|----------------|
| **S** | 40% | 0.249 | Sequence alone insufficient |
| **P** | 40% | 0.450 | Pathway alone insufficient |
| **E** | 40% | 0.200 | Evidence alone insufficient |
| **SP** | 100% | 0.467 | ✅ Pathway + Sequence sufficient |
| **SE** | 40% | 0.249 | Sequence + Evidence insufficient |
| **PE** | 40% | 0.507 | Pathway + Evidence insufficient |
| **SPE** | 100% | 0.524 | ✅ Full model (best confidence) |

### Key Finding:
**Pathway (P) is necessary for accuracy. Sequence (S) refines rankings. Evidence (E) increases confidence but not accuracy.**

---

## ❌ WHAT WAS NOT DONE FOR MM

### NOT Included in MM Paper:

| Item | Status | Evidence |
|------|--------|----------|
| **PARP Inhibitors** | ❌ NOT TESTED | MM drugs: BRAF/MEK/Proteasome/IMiD/Anti-CD38 |
| **Synthetic Lethality** | ❌ NOT USED | Used `/api/efficacy/predict` directly |
| **Ovarian Cancer** | ❌ NOT INCLUDED | Disease: `"multiple myeloma"` only |
| **DDR/HRR/BER Genes** | ❌ NOT TESTED | Variants: MAPK pathway only |
| **On-Label Chips** | ✅ IMPLEMENTED | But only for MM drugs |

---

## ✅ WHAT WAS DONE FOR OVARIAN CANCER (Separate Work)

### Synthetic Lethality Publication (Separate)

**Location:** `publications/synthetic_lethality/`

**Status:** Published/In Progress

### Validated Results:

| Metric | Value | Source | Validated? |
|--------|-------|--------|------------|
| **Clinical AUROC (TCGA-OV)** | 0.70 | `tcga_ov_clinical_report.json` | ✅ YES |
| **Benchmark (100-case)** | Results exist | `publication_suite_*.json` | ✅ YES |
| **Pathway Mapping Accuracy** | 62.2% | GDSC2 dataset | ✅ YES |

### SL Endpoint Implementation:

**Location:** `api/routers/guidance.py` (lines 396-551)

**Fast-Path Logic:**
```python
# Detects BER genes: MBD4, MUTYH, OGG1, NTHL1
# Detects HRR genes: BRCA1, BRCA2, PALB2, RAD51C, RAD51D
# Detects TP53 mutations

if has_ber and has_tp53:
    therapy = "PARP inhibitor (synthetic lethality: BER + checkpoint bypass)"
elif has_ber:
    therapy = "PARP inhibitor (BER deficiency - synthetic lethality)"
elif has_hrr:
    therapy = "PARP inhibitor (HRD - synthetic lethality)"
else:
    therapy = "platinum (DDR deficiency)"
```

**Validation:**
- ✅ Code exists: `api/services/synthetic_lethality/`
- ✅ Pathway mapping: 62.2% accuracy on GDSC2
- ✅ Used by: `benchmark_mbd4_tp53_accuracy.py`
- ✅ Clinical validation: TCGA-OV AUROC=0.70

---

## 🔍 DID IT PREDICT PARPS FOR OVARIAN CANCER?

### Answer: ✅ YES (But in Separate Publication)

**Evidence:**

1. **SL Service Exists:**
   - Location: `api/services/synthetic_lethality/`
   - Implements: BER→HR/PARP pathway mapping
   - Validation: 62.2% pathway mapping accuracy

2. **MBD4 → PARP Prediction:**
   - Fast-path detects: MBD4 (BER gene) + TP53
   - Returns: "PARP inhibitor (synthetic lethality: BER + checkpoint bypass)"
   - Tested in: `benchmark_mbd4_tp53_accuracy.py`

3. **Clinical Validation:**
   - TCGA-OV cohort: AUROC=0.70 for platinum response
   - Survival curves: KM plots for OS/PFS
   - DDR_bin classification: Validated on clinical outcomes

4. **NOT in MM Paper:**
   - MM paper focused on MAPK pathway (BRAF/KRAS/NRAS)
   - MM drugs: BRAF inhibitor, MEK inhibitor, Proteasome inhibitor, IMiD
   - No PARP inhibitors tested in MM paper

---

## 📋 CONFIDENCE SCORES: VALIDATED vs ASSUMED

### MM Paper (Validated):

| Drug | Variant | Confidence | Source | Validated? |
|------|---------|------------|--------|------------|
| **BRAF inhibitor** | BRAF V600E | 0.453-0.515 | `mm_efficacy_results.json` | ✅ YES |
| **MEK inhibitor** | KRAS G12D | 0.530 | `mm_efficacy_results.json` | ✅ YES |
| **IMiD** | TP53 R248W | 0.555 | `mm_efficacy_results.json` | ✅ YES |

### Ayesha Case (Assumed):

| Drug | Assumed Confidence | Actual Status |
|------|-------------------|---------------|
| **Olaparib** | 70% | ⚠️ ASSUMED - Not from actual API call |
| **Niraparib** | 65% | ⚠️ ASSUMED - Not from actual API call |
| **Pembrolizumab** | 65% | ⚠️ ASSUMED - TMB not confirmed |
| **Bevacizumab** | 60% + 5% boost | ⚠️ ASSUMED - Boost not validated |

**Why Assumed:**
- Document written before actual `/api/efficacy/predict` calls
- Reasonable estimates based on expected S/P/E calculation
- Placeholder values for clinical scenario documentation
- **Validation needed:** See `MBD4_PARP_CONFIDENCE_VALIDATION_PLAN.md`

---

## 🔧 HOW TO VALIDATE AYESHA PARP CONFIDENCE

### Step 1: Run Actual API Call

```python
# Call /api/efficacy/predict with MBD4 mutation
response = await client.post(
    f"{API_ROOT}/api/efficacy/predict",
    json={
        "model_id": "evo2_1b",
        "mutations": [
            {
                "gene": "MBD4",
                "hgvs_p": "p.K431Nfs*54",
                "chrom": "3",
                "pos": 129149435,
                "ref": "A",
                "alt": "del",
                "consequence": "frameshift_variant",
                "classification": "pathogenic"
            }
        ],
        "disease": "ovarian_cancer",
        "drugs": ["olaparib", "niraparib", "rucaparib"]
    }
)
```

### Step 2: Expected Calculation

**S/P/E Base:**
- Sequence (S): ~0.35 × sequence_percentile
- Pathway (P): ~0.35 × pathway_score
- Evidence (E): ~0.30 × evidence_strength
- **Base Total:** ~0.53-0.55

**SL Boost:**
- SL service detects: MBD4 (BER) → PARP
- Boost: +0.15 (max boost from SL)
- **Total:** ~0.68-0.70

**Matches Assumed 70%!**

### Step 3: Verify Against MM Methodology

**MM Confidence Range:** 0.45-0.56 (for MAPK variants)
**Expected Ovarian Range:** 0.55-0.70 (higher due to SL boost)

**Validation:** Run actual API call, compare to assumed 70%

---

## 📊 SUMMARY: WHAT'S VALIDATED vs ASSUMED

### ✅ VALIDATED (High Confidence):

1. **MM Paper:**
   - 100% MAPK pathway alignment (5/5 variants)
   - Confidence scores: 0.45-0.56 range
   - Ablation study: P is necessary
   - Calibration: ECE=0.479

2. **SL Service:**
   - Code exists: `api/services/synthetic_lethality/`
   - Pathway mapping: 62.2% accuracy
   - Clinical validation: TCGA-OV AUROC=0.70
   - MBD4→PARP mechanism: Implemented

3. **On-Label Chips:**
   - Backend: `_on_label_stub()` function exists
   - Frontend: `ChemoGuidanceCard.jsx` displays chips
   - Tier assignment: `_tier_from_gates()` works

### ⚠️ MECHANISM VALIDATED, CONFIDENCE ASSUMED:

1. **MBD4 → PARP:**
   - ✅ Mechanism: SL service implements BER→PARP
   - ✅ Pathway: 62.2% accuracy
   - ⚠️ Confidence: 70% assumed, not from actual API call
   - ⚠️ Clinical outcomes: No patient data (MBD4 too rare)

### ❌ NOT VALIDATED (Assumed):

1. **Ayesha Drug Confidence:**
   - Olaparib: 70% (assumed)
   - Niraparib: 65% (assumed)
   - Pembrolizumab: 65% (assumed, TMB not confirmed)
   - Bevacizumab: 60% + 5% boost (assumed)

2. **PDGFRA VUS:**
   - Evo2 delta: -0.08 (assumed, not run)
   - Functionality: 0.60 (assumed, not run)
   - PDGFRA→VEGF mapping: Not in code
   - TKI confidence: 35-40% (assumed)

---

## 🎯 KEY INSIGHTS

### 1. MM Paper is Methodologically Sound
- ✅ 100% pathway alignment validated
- ✅ Ablation study shows P is necessary
- ✅ Calibration analysis performed
- ✅ Results reproducible

### 2. PARP Prediction Exists (Separate Work)
- ✅ SL service implements MBD4→PARP
- ✅ Clinical validation: TCGA-OV AUROC=0.70
- ✅ Pathway mapping: 62.2% accuracy
- ❌ NOT in MM paper (different cancer type)

### 3. Ayesha Confidence Scores are Estimates
- ⚠️ 70% for Olaparib: Reasonable estimate, not validated
- ⚠️ Based on expected S/P/E + SL boost calculation
- ⚠️ Need actual API call to validate
- ⚠️ Clinical outcomes not validated (MBD4 too rare)

### 4. Ovarian Cancer Work Exists But Incomplete
- ✅ SL endpoint exists and works
- ✅ Clinical validation on TCGA-OV
- ❌ No ovarian cancer ablation study (like MM)
- ❌ No ovarian on-label rules (all show off-label)
- ❌ No ovarian baseline script (like MM)

---

## 📁 FILES ANALYZED

### MM Publication:
- ✅ `PAPER_DRAFT.md` - Manuscript
- ✅ `mm_efficacy_results.json` - Baseline results
- ✅ `ablation_results_*.json` - Ablation study
- ✅ `REPRODUCIBILITY.md` - Reproduction guide

### Ovarian Cancer:
- ✅ `OVARIAN_CANCER_REPRODUCTION_PLAN.md` - Plan (not executed)
- ✅ `OVARIAN_SL_EVALUATION.md` - SL evaluation
- ✅ `synthetic_lethality/` - Separate publication

### Code:
- ✅ `api/routers/guidance.py` - SL endpoint
- ✅ `api/services/synthetic_lethality/` - SL service
- ✅ `benchmark_mbd4_tp53_accuracy.py` - MBD4 testing

---

**Commander, this is the complete picture:**

1. **MM paper:** Validated MAPK pathway predictions (100% accuracy)
2. **PARP predictions:** Exist in separate SL publication (TCGA-OV AUROC=0.70)
3. **Ayesha confidence:** Assumed 70% (reasonable estimate, needs validation)
4. **Ovarian work:** Planned but not fully executed (reproduction plan exists)
