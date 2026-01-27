# 💀 CA-125 DATA INVESTIGATION - COMPREHENSIVE AUDIT

**Date:** January 26, 2026 (Updated)  
**Status:** 🔴 **CA-125 KINETICS NOT AVAILABLE IN ANY CURRENT DATASET**  
**Previous Status:** Outdated - claimed BriTROC-1 has complete CA-125 data

---

## 🚨 CRITICAL UPDATE

### Previous Claim (WRONG ❌)

From the original document:
> "BriTROC-1 **DEFINITELY has CA-125** (required for trial enrollment)"
> "**Complete CA-125 data** (trial requirement)"

### Reality (AUDITED ✅)

**BriTROC-1 Supplementary Data Audit:**

| File | CA-125 Data Found | Details |
|------|-------------------|---------|
| `MOESM4_ESM.txt` | ⚠️ Partial | Boolean flags only (True/False for "response via CA-125") |
| `MOESM5_ESM.txt` | ❌ None | Treatment data (drugs, doses, cycles) |
| `figure_4A.tsv` | ❌ None | Copy number signatures |
| `figure_7A/7B.tsv` | ❌ None | CD3/CD8 immune markers (NOT CA-125) |
| Source data (76 files) | ❌ None | Genomic data, no longitudinal CA-125 |

**What We Actually Have:**
- `overall_response_ca125`: Boolean (True = response determined partly by CA-125)
- `progression_ca125`: Boolean (True = progression detected by CA-125)

**What We DON'T Have:**
- ❌ Actual CA-125 values (U/mL)
- ❌ CA-125 at each chemo cycle
- ❌ Dates of CA-125 measurements
- ❌ **KELIM cannot be calculated**

---

## 📊 DATA AVAILABILITY MATRIX (UPDATED)

| Dataset | Patients | CA-125 Values | Multiple Timepoints | KELIM Possible? |
|---------|----------|---------------|---------------------|-----------------|
| **GSE165897** | 11 | ✅ 2 values (TN + PN) | ❌ Only 2 | ❌ No (need ≥3) |
| **BriTROC-1** | 276 | ❌ Boolean only | ❌ No values | ❌ No |
| **TCGA-OV** | 489 | ⚠️ Variable | ❌ Single timepoint | ❌ No |
| **MSK_SPECTRUM** | 105 | ❓ Unknown | ❓ Unknown | ❓ Unknown |

**Verdict:** ❌ **No dataset currently available has the CA-125 kinetics data needed for KELIM calculation.**

---

## 🔬 What KELIM Requires

### Minimum Data (for KELIM calculation):

```python
# KELIM requires ≥3 CA-125 measurements over time
ca125_values = [500, 350, 200, 120, 80]  # U/mL at each cycle
timepoints = [0, 21, 42, 63, 84]         # Days since start of chemo

# Linear regression on log(CA-125) vs time
# KELIM = -slope (higher = faster elimination = better response)
```

### What We Have vs What We Need:

| Requirement | GSE165897 | BriTROC-1 | Needed |
|-------------|-----------|-----------|--------|
| CA-125 baseline | ✅ Yes | ❌ No values | ✅ |
| CA-125 cycle 1-6 | ❌ No | ❌ No values | ✅ |
| Dates of measurements | ❌ No | ❌ No | ✅ |
| Post-treatment value | ✅ Yes | ❌ No values | ✅ |
| KELIM calculation | ❌ Not possible | ❌ Not possible | Target |

---

## 🎯 ALTERNATIVE APPROACHES

### Option 1: CA-125 Drop Ratio (GSE165897)

We CAN calculate a simple CA-125 drop ratio (not KELIM):

```python
drop_ratio = (CA125_baseline - CA125_post) / CA125_baseline
# e.g., (3776 - 343) / 3776 = 90.9%
```

**Limitation:** Drop ratio ≠ KELIM. Drop ratio ignores time kinetics.

**Current Results (from original document):**
- Resistant (n=8): Mean drop = 86.2%
- Sensitive (n=3): Mean drop = 90.9%
- **No clear separation** (both groups show large drops)

### Option 2: Abandon KELIM for Now

Given data limitations:
1. ✅ Focus on **Pathways-Only Model** (validated, AUC=0.75)
2. ✅ Accept **CN signatures from BriTROC-1** (validated, AUC=0.874)
3. ❌ Defer KELIM+ until proper dataset found

### Option 3: Search for Alternative Dataset

**Potential Sources for CA-125 Kinetics:**

| Study | Potential | Access | Notes |
|-------|-----------|--------|-------|
| ICON7 | ⭐⭐⭐ | Publication | Large trial, likely has CA-125 kinetics |
| GOG-218 | ⭐⭐⭐ | dbGaP | Bevacizumab trial, extensive CA-125 |
| AGO-OVAR trials | ⭐⭐⭐ | Contact | German trials, detailed kinetics |
| KELIM validation cohorts | ⭐⭐⭐ | Contact | From You et al. publications |

---

## 📋 UPDATED RECOMMENDATIONS

### ✅ WHAT WE CAN DO NOW

1. **Use GSE165897 for Pathways-Only Model** (validated)
   - AUC = 0.75 for post-PI3K
   - Spearman ρ = -0.711 for post-DDR (p = 0.014)
   - **Publication-ready**

2. **Use BriTROC-1 for CN Signature Validation** (validated)
   - AUC = 0.874 for CN_s7 at relapse
   - p = 0.035 (Wilcoxon)
   - **Confirms Serial SAE hypothesis**

3. **Combine results for manuscript**
   - Two independent cohorts
   - Different molecular modalities (pathways vs CN)
   - Both support post-treatment > pre-treatment

### ❌ WHAT WE CANNOT DO

1. **KELIM+ Model** - No CA-125 kinetics data available
2. **Combined CA-125 + Pathways** - Insufficient temporal resolution
3. **Time-to-progression prediction from CA-125** - No cycle-by-cycle data

---

## � AUDIT OF ORIGINAL DOCUMENT

### Errors Found:

| Original Claim | Reality | Impact |
|----------------|---------|--------|
| "BriTROC-1 has complete CA-125" | ❌ Only Boolean flags | Major |
| "Trial requirement for CA-125" | ⚠️ May have been collected but not shared | Medium |
| "Expected Outcome: CA-125 kinetics cycles 1-6" | ❌ Not in supplementary data | Major |
| "KELIM possible with BriTROC-1" | ❌ Cannot calculate | Critical |

### Correct Information:

| Original Claim | Status |
|----------------|--------|
| GSE165897 has 2 CA-125 timepoints | ✅ Correct |
| KELIM requires ≥3 timepoints | ✅ Correct |
| Pathways-Only Model validated (AUC=0.75) | ✅ Correct |
| CA-125 drop % doesn't separate resistant/sensitive | ✅ Correct |

---

## 📊 WHAT TO UPDATE

### Files to Update:

1. **`docs/serial_sae/CA125_DATA_INVESTIGATION.md`** - This file (DONE)
2. **`publications/serial-sae/MANUSCRIPT_DRAFT.md`** - Remove KELIM+ claims
3. **`.cursor/ayesha/ZO_MASTER_DOCUMENTATION.md`** - Add this audit

### Claims to Remove from Manuscript:

- [ ] Remove: "BriTROC-1 has complete CA-125 kinetics"
- [ ] Remove: "KELIM+ model in future directions" (until data found)
- [ ] Add: Caveat about CA-125 kinetics data limitations

---

## ✅ CONCLUSION

**CA-125 Kinetics Status:** 🔴 **NOT AVAILABLE**

| Data Type | GSE165897 | BriTROC-1 | TCGA-OV |
|-----------|-----------|-----------|---------|
| CA-125 values | ✅ 2 points | ❌ None | ⚠️ Variable |
| Multiple timepoints | ❌ No | ❌ No | ❌ No |
| KELIM calculation | ❌ No | ❌ No | ❌ No |

**Action:** Defer KELIM+ model until proper dataset with CA-125 kinetics is identified.

**Current Best Evidence:**
1. GSE165897 Pathways-Only: AUC = 0.75, p = 0.014
2. BriTROC-1 CN Signatures: AUC = 0.874, p = 0.035

**Both support Serial SAE hypothesis without needing KELIM.**

---

**Audit Complete:** January 26, 2026
