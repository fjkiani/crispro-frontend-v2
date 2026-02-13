# 🔬 RESISTANCE PREDICTION: SAE FEATURES & DNA REPAIR CAPACITY - COMPLETE EXPLANATION

**Date:** January 26, 2026  
**Purpose:** Answer: "How are SAE features extracted?" and "How is DNA repair capacity computed?"

---

## 1. SAE FEATURE EXTRACTION

### Two Modes: PROXY vs TRUE SAE

#### **PROXY SAE (Default - Production)**
**Location:** `api/services/sae_feature_service.py` (lines 124-365)

**Inputs:**
- `insights_bundle`: From `/api/insights/predict_*` endpoints (4 chips: functionality, chromatin, essentiality, regulatory)
- `pathway_scores`: From efficacy orchestrator (P component from S/P/E framework)
- `tumor_context`: HRD, TMB, MSI, somatic mutations

**How It Works:**
1. Extracts pathway burden scores from `pathway_scores` dict:
   - `pathway_burden_ddr` = `pathway_scores.get("ddr", 0.0)`
   - `pathway_burden_mapk` = `pathway_scores.get("mapk", 0.0)`
   - `pathway_burden_pi3k` = `pathway_scores.get("pi3k", 0.0)`
   - etc.

2. Computes essentiality for HRR genes:
   - `essentiality_hrr` = average of essentiality scores for HRR genes (BRCA1, BRCA2, PALB2, RAD51C, RAD51D, BRIP1, BARD1, ATM) found in patient's mutations
   - Uses `insights_bundle.get("essentiality", 0.0)` as proxy

3. Computes exon disruption score:
   - Only if `essentiality_hrr > 0.65` (threshold)
   - Uses `insights_bundle.get("regulatory", 0.0)` as proxy for exon disruption

4. Builds 7D mechanism vector:
   - Converts pathway scores to mechanism vector using `convert_pathway_scores_to_mechanism_vector()`
   - Format: `[DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]`

**Provenance:** `"sae": "proxy"` (line 280)

---

#### **TRUE SAE (Optional - Research Mode)**
**Location:** `api/services/sae_feature_service.py` (lines 283-347)

**Requirements:**
- Feature flag: `ENABLE_TRUE_SAE_PATHWAYS=true`
- Input: `sae_features` dict with 32K-dim feature vector from `/api/sae/extract_features`

**How It Works:**
1. Calls `_compute_sae_diagnostics()` (lines 390-525):
   - Loads feature→pathway mapping from `api/resources/sae_feature_mapping.json`
   - Aggregates feature activations for each pathway:
     - DDR: Features mapped to DDR pathway
     - MAPK: Features mapped to MAPK pathway
     - etc.

2. Computes pathway scores from TRUE SAE:
   - `pathway_burden_ddr` = `sae_diagnostics.get("ddr_sae_score")`
   - `pathway_burden_mapk` = `sae_diagnostics.get("mapk_sae_score")`
   - etc.

3. Updates mechanism vector with TRUE SAE pathway scores

4. Updates DNA repair capacity using TRUE SAE DDR score

**Provenance:** `"sae": "true_sae"` (line 327)

**Source of TRUE SAE Features:**
- Endpoint: `/api/sae/extract_features` (via `api/routers/sae.py`)
- Calls Modal service or Evo2 directly to extract 32K-dim features from layer-26 activations
- Returns: `{"features": [32768 floats], "top_features": [...], "stats": {...}}`

---

## 2. DNA REPAIR CAPACITY COMPUTATION

### **Formula (Manager Approved - Jan 13, 2025)**

**Location:** `api/services/sae_feature_service.py` (lines 576-600)

```python
DNA_REPAIR_CAPACITY_WEIGHTS = {
    "pathway_ddr": 0.60,         # 60% weight
    "essentiality_hrr": 0.20,    # 20% weight
    "exon_disruption": 0.20      # 20% weight
}

dna_repair_capacity = (
    0.60 × pathway_burden_ddr +
    0.20 × essentiality_hrr_genes +
    0.20 × exon_disruption_score
)
```

### **Component Breakdown:**

1. **Pathway DDR (60% weight):**
   - Source: `pathway_scores.get("ddr", 0.0)` from S/P/E framework
   - Range: 0.0-1.0
   - Meaning: DNA Damage Repair pathway burden (higher = more disruption = lower repair capacity)

2. **Essentiality HRR (20% weight):**
   - Source: Average essentiality for HRR genes (BRCA1, BRCA2, PALB2, RAD51C, RAD51D, BRIP1, BARD1, ATM)
   - Computed in `_compute_essentiality_hrr()` (lines 527-554)
   - Uses `insights_bundle.get("essentiality", 0.0)` as proxy
   - Range: 0.0-1.0
   - Meaning: Homologous Recombination Repair gene essentiality (higher = more dependency = lower repair capacity)

3. **Exon Disruption (20% weight):**
   - Source: `insights_bundle.get("regulatory", 0.0)` (proxy for exon disruption)
   - Only applied if `essentiality_hrr > 0.65` (threshold)
   - Computed in `_compute_exon_disruption_score()` (lines 556-574)
   - Range: 0.0-1.0
   - Meaning: Exon-level disruption score (higher = more disruption = lower repair capacity)

### **Interpretation:**
- **High DNA repair capacity (≥0.70):** Patient's tumor has intact DNA repair → PARP inhibitors less effective
- **Low DNA repair capacity (<0.40):** Patient's tumor has deficient DNA repair → PARP inhibitors more effective (synthetic lethality)

---

## 3. RESISTANCE PREDICTION: DNA REPAIR RESTORATION SIGNAL

### **Location:** `api/services/resistance_prophet_service.py` (lines 570-666)

### **Signal Detection Logic:**

**Signal 1: DNA Repair Restoration**
- **Detection:** Compare current DNA repair capacity to baseline
- **Threshold:** `DNA_REPAIR_THRESHOLD = 0.15` (15% change)
- **Trigger:** `repair_change < -0.15` (capacity DROPPING = restoration detected)
- **Meaning:** Tumor is restoring its DNA repair capacity → PARP resistance mechanism

**Formula:**
```python
current_repair = current_sae.get("dna_repair_capacity", 0.0)
baseline_repair = baseline_sae.get("dna_repair_capacity", 0.5)  # Default if missing
repair_change = current_repair - baseline_repair
detected = repair_change < -0.15  # Negative change = restoration
```

**Mechanism Breakdown:**
- Returns `MechanismBreakdown` with:
  - `ddr_pathway_change` = `current_ddr - baseline_ddr`
  - `hrr_essentiality_change` = `current_hrr - baseline_hrr`
  - `exon_disruption_change` = `current_exon - baseline_exon`

**Pathway Contributions:**
- DDR: 60% contribution to DNA repair capacity
- HRR: 20% contribution
- Exon: 20% contribution

---

## 4. RESISTANCE PREDICTION: 2-OF-3 TRIGGER RULE

### **Location:** `api/services/resistance_prophet_service.py` (lines 391-500)

### **Three Signals:**

1. **DNA Repair Restoration** (Signal 1)
   - Threshold: DNA repair capacity drop ≥0.15
   - Detection: `_detect_dna_repair_restoration()`

2. **Pathway Escape** (Signal 2)
   - Threshold: Pathway burden drop ≥0.15
   - Detection: `_detect_pathway_escape()`

3. **CA-125 Kinetics** (Signal 3)
   - Threshold: Inadequate response or on-therapy rise
   - Detection: `ca125_intelligence` service

### **Risk Stratification:**

**HIGH Risk:**
- Probability ≥0.70 AND ≥2 signals detected

**MEDIUM Risk:**
- Probability 0.50-0.69 OR exactly 1 signal detected

**LOW Risk:**
- Probability <0.50

---

## 5. SUMMARY

### **SAE Feature Extraction:**
- **PROXY (default):** Uses insights bundle + pathway scores (no Evo2 layer-26 extraction)
- **TRUE (optional):** Extracts 32K-dim features from Evo2 layer-26 activations via Modal service
- **Both modes** compute the same SAE features (DNA repair capacity, pathway burden, mechanism vector)

### **DNA Repair Capacity:**
- **Formula:** `0.6 × DDR + 0.2 × HRR + 0.2 × Exon`
- **Source:** Pathway scores (S/P/E), insights bundle (essentiality/regulatory), tumor context
- **Use:** Resistance prediction (DNA repair restoration signal), PARP eligibility, synthetic lethality

### **Resistance Prediction:**
- **3 signals:** DNA repair restoration, pathway escape, CA-125 kinetics
- **2-of-3 rule:** HIGH risk if ≥2 signals + probability ≥0.70
- **DNA repair restoration:** Detected when capacity drops ≥0.15 vs baseline

---

## 6. CODE REFERENCES

**SAE Feature Service:**
- `oncology-coPilot/oncology-backend-minimal/api/services/sae_feature_service.py`
- Lines 124-365: `compute_sae_features()` (main computation)
- Lines 576-600: `_compute_dna_repair_capacity()` (formula)
- Lines 390-525: `_compute_sae_diagnostics()` (TRUE SAE pathway mapping)

**Resistance Prophet Service:**
- `oncology-coPilot/oncology-backend-minimal/api/services/resistance_prophet_service.py`
- Lines 570-666: `_detect_dna_repair_restoration()` (Signal 1)
- Lines 391-500: `predict_resistance()` (2-of-3 trigger logic)

**SAE Extraction Endpoint:**
- `oncology-coPilot/oncology-backend-minimal/api/routers/sae.py`
- Endpoint: `/api/sae/extract_features`
- Calls Modal service or Evo2 directly for 32K-dim feature extraction

---

**Status:** ✅ **COMPLETE** - All questions answered with code references
