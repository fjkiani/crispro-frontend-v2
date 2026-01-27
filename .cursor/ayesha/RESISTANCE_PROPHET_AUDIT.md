# RESISTANCE PROPHET - COMPLETE MECHANISM AUDIT

**Date:** January 8, 2025  
**Patient:** Ayesha (AK) - Stage IVB Ovarian Cancer  
**Issue:** HIGH risk (97.1%) showing due to baseline comparison against population average  
**Goal:** Understand mechanism to decrease false positives

---

## 🔍 **ROOT CAUSE ANALYSIS**

### **The Problem:**
```
Current SAE Features:
  - DNA repair capacity: 0.0
  - Mechanism vector: [0.5, 0.2, 0.2, 0.3, 0.0, 0.0, 0.0]

Population Baseline (used when no real baseline):
  - DNA repair capacity: 0.50
  - Mechanism vector: [0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50]

Result:
  - DNA Repair Restoration: -0.50 change → DETECTED (threshold: -0.15)
  - Pathway Escape: DDR pathway change → DETECTED
  - Probability: 97.1% (weighted average of signal probabilities)
  - Risk Level: HIGH (≥0.70 probability + ≥2 signals)
```

### **Why This Happens:**
1. **Treatment-naive patient** → No real baseline measurements
2. **Population baseline (0.50)** → Arbitrary average, not Ayesha's actual baseline
3. **Current SAE features (0.0)** → Computed from proxy pathway scores (not real NGS)
4. **Large delta (-0.50)** → Triggers restoration signal even though it's a false comparison

---

## ⚙️ **EXISTING MECHANISMS TO DECREASE FALSE POSITIVES**

### **1. Confidence Penalty (Manager Q16)**
**Location:** `_compute_confidence()` (line 1626-1628)

```python
if not baseline_available:
    avg_confidence *= 0.80  # 20% penalty
```

**What it does:**
- Reduces confidence from 0.60 → 0.48 (for DNA repair signal)
- Reduces confidence from 0.85 → 0.68 (for pathway escape signal)

**Limitation:**
- Only affects **confidence**, not **probability**
- Risk level is determined by **probability** (97.1%), not confidence
- HIGH risk still triggers because probability ≥0.70

### **2. Confidence Cap (Manager Q15)**
**Location:** `_compute_confidence()` (line 1631-1635)

```python
if not ca125_available and signal_count < 2:
    if avg_confidence > 0.60:
        avg_confidence = 0.60
        confidence_cap = "MEDIUM"
```

**What it does:**
- Caps confidence at 0.60 (MEDIUM) when CA-125 missing AND <2 signals

**Limitation:**
- Ayesha has **2 signals** (DNA repair + pathway escape)
- So this cap doesn't apply
- Even if it did, it only affects confidence, not probability

### **3. Risk Stratification Cap (Manager Q15)**
**Location:** `_stratify_risk()` (line 1591-1593)

```python
if not ca125_available and signal_count < 2 and probability >= 0.70:
    return ResistanceRiskLevel.MEDIUM
```

**What it does:**
- Caps risk at MEDIUM when CA-125 missing AND <2 signals

**Limitation:**
- Ayesha has **2 signals**, so this doesn't apply
- This is the closest mechanism, but requires <2 signals

---

## 🎯 **MECHANISMS TO DECREASE RESISTANCE METRIC**

### **Option 1: Apply Probability Penalty When Baseline Missing**
**Location:** `_compute_resistance_probability()` (line 1552-1574)

**Current Logic:**
```python
# Weighted average by confidence
total_probability = sum(sig.probability * sig.confidence for sig in active_signals)
total_weight = sum(sig.confidence for sig in active_signals)
overall_probability = total_probability / total_weight
```

**Proposed Fix:**
```python
overall_probability = total_probability / total_weight

# Manager Q16 Extension: Apply probability penalty when baseline missing
if baseline_source == "population_average":
    overall_probability *= 0.70  # 30% penalty (more aggressive than confidence penalty)
    logger.info("Probability penalty applied: baseline SAE missing (30% reduction)")
```

**Effect:**
- 97.1% → 68.0% (drops to MEDIUM risk)
- Still requires ≥2 signals for HIGH, but probability now below 0.70 threshold

### **Option 2: Require More Signals for HIGH Risk When Baseline Missing**
**Location:** `_stratify_risk()` (line 1596)

**Current Logic:**
```python
if probability >= 0.70 and signal_count >= 2:
    return ResistanceRiskLevel.HIGH
```

**Proposed Fix:**
```python
# Manager Q16 Extension: Require ≥3 signals for HIGH when baseline missing
min_signals_for_high = 3 if baseline_source == "population_average" else 2

if probability >= 0.70 and signal_count >= min_signals_for_high:
    return ResistanceRiskLevel.HIGH
```

**Effect:**
- Ayesha has 2 signals → HIGH risk blocked
- Falls through to MEDIUM risk (0.50-0.69 OR 1 signal)

### **Option 3: Cap Probability at MEDIUM Threshold When Baseline Missing**
**Location:** `_compute_resistance_probability()` (line 1574)

**Current Logic:**
```python
return min(overall_probability, 1.0)
```

**Proposed Fix:**
```python
overall_probability = total_probability / total_weight

# Manager Q16 Extension: Cap probability at MEDIUM threshold when baseline missing
if baseline_source == "population_average":
    overall_probability = min(overall_probability, 0.69)  # Cap at MEDIUM threshold
    logger.info("Probability capped at MEDIUM threshold: baseline SAE missing")

return min(overall_probability, 1.0)
```

**Effect:**
- 97.1% → 69.0% (drops to MEDIUM risk)
- Maximum probability when baseline missing is 0.69 (just below HIGH threshold)

### **Option 4: Use More Conservative Signal Detection Thresholds**
**Location:** `_detect_dna_repair_restoration()` (line 621) and `_detect_pathway_escape()` (line 720)

**Current Thresholds:**
- DNA_REPAIR_THRESHOLD = 0.15
- PATHWAY_ESCAPE_THRESHOLD = 0.15

**Proposed Fix:**
```python
# Use more conservative threshold when baseline is population average
effective_threshold = self.DNA_REPAIR_THRESHOLD * 2.0 if baseline_source == "population_average" else self.DNA_REPAIR_THRESHOLD
detected = repair_change < -effective_threshold
```

**Effect:**
- Threshold: -0.15 → -0.30 (more conservative)
- Ayesha's change: -0.50 → Still detected, but probability calculation uses new threshold
- Probability: 97.1% → ~85% (still HIGH, but lower)

---

## 📊 **RECOMMENDED SOLUTION: COMBINED APPROACH**

**Best Practice:** Combine **Option 1 (Probability Penalty)** + **Option 2 (More Signals Required)**

**Rationale:**
1. **Probability penalty (30%)** → Reduces false positive rate when baseline missing
2. **Require ≥3 signals** → Ensures HIGH risk only when multiple independent signals agree
3. **Maintains clinical value** → Still detects real resistance, just more conservative

**Implementation:**
```python
# In _compute_resistance_probability():
overall_probability = total_probability / total_weight

# Manager Q16 Extension: Apply probability penalty when baseline missing
if baseline_source == "population_average":
    overall_probability *= 0.70  # 30% penalty
    logger.info("Probability penalty applied: baseline SAE missing (30% reduction)")

# In _stratify_risk():
# Manager Q16 Extension: Require ≥3 signals for HIGH when baseline missing
min_signals_for_high = 3 if baseline_source == "population_average" else 2

if probability >= 0.70 and signal_count >= min_signals_for_high:
    return ResistanceRiskLevel.HIGH
```

**Expected Result for Ayesha:**
- Probability: 97.1% → 68.0% (after 30% penalty)
- Signal count: 2 (not enough for HIGH when baseline missing)
- Risk Level: **MEDIUM** (instead of HIGH)
- Confidence: 0.58 (already penalized)

---

## 🧬 **AYESHA'S SYSTEM ARCHITECTURE - CORE UNDERSTANDING**

### **1. The Mission (From AYESHA_MASTER.md)**
Build a **complete precision oncology platform** for **Ayesha**, a 40-year-old with Stage IVB ovarian cancer (sporadic, germline-negative). Deliver **high-confidence clinical decision support** **BEFORE NGS results arrive**, then seamlessly upgrade to personalized predictions when tumor genomics complete.

### **2. The Strategic Shift**
**85-90% of ovarian cancers are sporadic** (not hereditary), but most platforms only work for germline-positive patients. We built **sporadic-first capabilities**:
- **PARP Rescue**: HRD ≥42 → full PARP effect (even if germline negative!)
- **IO Boost**: TMB ≥20 or MSI-H → 1.3x boost for checkpoint inhibitors
- **Confidence Capping**: L0 (completeness <0.3) → cap at 0.4, L1 → 0.6, L2 → none

### **3. Ayesha's Genetic Profile (From AYESHA_MASTER.md)**
- **MBD4**: c.1239delA (p.Ile413Serfs*2) - Germline Homozygous → BER pathway non-functional
- **TP53**: p.Arg175His (R175H) - Somatic → G1/S checkpoint bypassed
- **Combined Effect:** Dual pathway deficiency → Synthetic lethality with PARP inhibition
- **DDR pathway = 1.0** → Maximum pathway disruption
- **Resistance risk: LOW** → No MAPK mutations detected (from Resistance Agent)

### **4. The S/P/E Framework (Sequence/Pathway/Evidence)**
**Formula:** `raw_lob = 0.3 * seq_pct + 0.4 * path_pct + 0.3 * s_evd`

**Components:**
- **S (Sequence)**: Evo2 provides zero-shot variant impact prediction
- **P (Pathway)**: Pathway aggregation captures multi-hit tumor evolution
- **E (Evidence)**: Literature + ClinVar provide real-world validation

**Why Multi-Modal?**
- **Single-metric myopia**: Delta score alone insufficient
- **Transparency**: S/P/E breakdown explains WHY confidence is 73%
- **Clinical trust**: Doctors need mechanistic rationale

### **5. SAE Features → Resistance Prophet Flow**

```
Tumor Context (NGS)
  ↓
SAE Feature Service
  ├─ DNA Repair Capacity = 0.6×DDR + 0.2×HRR + 0.2×exon
  ├─ Mechanism Vector = [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
  └─ Pathway Burden = pathway_scores from efficacy router
  ↓
Resistance Prophet Service
  ├─ Signal 1: DNA Repair Restoration (compare current vs baseline)
  ├─ Signal 2: Pathway Escape (compare mechanism vectors)
  └─ Signal 3: CA-125 Kinetics (if available)
  ↓
Risk Stratification
  ├─ HIGH: ≥0.70 probability + ≥2 signals
  ├─ MEDIUM: 0.50-0.69 OR 1 signal
  └─ LOW: <0.50
```

### **6. The Baseline Problem**

**For Treatment-Naive Patients:**
- No prior SAE measurements → No real baseline
- System uses **population average (0.50)** as baseline
- This creates **false comparisons**:
  - Current: 0.0 (from proxy pathway scores)
  - Baseline: 0.50 (population average)
  - Change: -0.50 → Triggers "restoration" signal

**Why This Is Wrong:**
- Ayesha's DNA repair capacity of 0.0 is **not** a drop from 0.50
- It's her **actual baseline** (computed from proxy pathway scores)
- Comparing against population average creates false "restoration" signal

### **7. The Clinical Context**

**From AYESHA_MASTER.md:**
- **Resistance risk: LOW** → No MAPK mutations detected
- **PARP Eligible: TRUE** → MBD4+TP53 synthetic lethality
- **DDR pathway = 1.0** → Maximum pathway disruption
- **Treatment-naive** → First-line therapy

**What Resistance Prophet Should Say:**
- **NOT_APPLICABLE** or **LOW** risk for treatment-naive patients
- **HIGH** risk only after cycle 3 when real baseline + follow-up available
- **MEDIUM** risk if baseline missing but treatment ongoing

---

## 🎯 **THE CORE UNDERSTANDING**

### **What Was Built for Ayesha:**

1. **Complete Care v2 Orchestrator** → Unified care plan (trials, SOC, CA-125, WIWFM, resistance)
2. **CA-125 Intelligence** → Burden classification, response forecast, resistance detection
3. **Drug Efficacy (WIWFM)** → S/P/E framework with Evo2 integration
4. **Resistance Playbook** → 7 combo strategies, 6 next-line switches
5. **Resistance Prophet** → Predicts resistance 3-6 months early
6. **SAE Features** → DNA repair capacity, mechanism vector, pathway burden
7. **Next Test Recommender** → HRD → ctDNA → SLFN11 → ABCB1
8. **Hint Tiles** → Clinician action hints (max 4)
9. **Mechanism Map** → 6 pathway chips (color-coded post-NGS)

### **The Resistance Prophet's Purpose:**

**From Manager Policy (MANAGER_ANSWERS_TO_RESISTANCE_PROPHET_QUESTIONS.md):**
- **Q3:** Phase 1 = retrospective WITHOUT CA-125 (DNA repair + pathway escape only)
- **Q7:** Opt-in via `include_resistance_prediction=true`; default=off
- **Q9:** Thresholds - HIGH: >=0.70 + >=2 signals; MEDIUM: 0.50-0.69 OR 1 signal; LOW: <0.50
- **Q16:** If baseline SAE missing, use 0.50 with confidence penalty

### **The Gap:**

**Manager Q16** says "use 0.50 with confidence penalty" but doesn't address:
- **Probability penalty** (only confidence is penalized)
- **Signal count requirement** (still requires ≥2 signals for HIGH)
- **Risk level capping** (HIGH still triggers if probability ≥0.70)

**The Fix:**
Extend Manager Q16 to include:
1. **Probability penalty (30%)** when baseline missing
2. **Require ≥3 signals** for HIGH risk when baseline missing
3. **Document in provenance** that baseline penalty affects both confidence AND probability

---

## 📋 **SUMMARY**

### **Current State:**
- Resistance Prophet shows **HIGH risk (97.1%)** for Ayesha
- **Root cause:** Comparing current SAE (0.0) vs population baseline (0.50)
- **Existing mechanisms:** Confidence penalty (20%), confidence cap (MEDIUM), risk cap (MEDIUM if <2 signals)
- **Limitation:** None of these affect **probability**, only **confidence**

### **Recommended Fix:**
1. **Apply 30% probability penalty** when baseline missing (in `_compute_resistance_probability`)
2. **Require ≥3 signals** for HIGH risk when baseline missing (in `_stratify_risk`)
3. **Document in provenance** that baseline penalty affects both confidence AND probability

### **Expected Result:**
- Probability: 97.1% → 68.0% (after 30% penalty)
- Signal count: 2 (not enough for HIGH when baseline missing)
- Risk Level: **MEDIUM** (instead of HIGH)
- Confidence: 0.58 (already penalized)

### **Clinical Rationale:**
- Treatment-naive patients cannot have resistance (no treatment to resist)
- Population baseline (0.50) is arbitrary, not patient-specific
- More conservative thresholds prevent false alarms
- Still detects real resistance when baseline available (after cycle 3)

---

**Last Updated:** January 8, 2025  
**For:** Ayesha (AK) - Stage IVB Ovarian Cancer  
**Status:** Audit complete, mechanism understood, fix recommended

