# 🔧 PLUMBER INVESTIGATION TASKS

**Date:** December 24, 2025  
**Purpose:** Concrete tasks Plumber can execute NOW to answer hard questions  
**Owner:** Plumber  
**Status:** ⏳ **READY FOR EXECUTION**

---

## 🎯 PRIORITY 1: VALIDATE DDR_bin BIOLOGY (Week 1)

### **Task P1.1: Test DDR_bin on Known Resistance Mutations**

**Goal:** Confirm DDR_bin actually drops when HR restoration happens.

**Method:**
1. Find patients in Tier-3 with known RAD51C or BRCA1 reversions (if any)
2. Check their DDR_bin scores
3. Compare to patients without reversions

**Files:**
- Input: `data/validation/sae_cohort/checkpoints/Tier3_validation_cohort.json`
- Output: `data/validation/resistance_mutation_ddr_bin_analysis.json`

**Expected Result:**
- Patients with HR restoration mutations should have LOWER DDR_bin
- If not → our DDR_bin interpretation is wrong

```python
# Pseudocode
for patient in tier3_cohort:
    if has_reversion(patient.mutations, ['RAD51C', 'RAD51D', 'BRCA1', 'BRCA2']):
        reversion_group.append(patient.ddr_bin)
    else:
        no_reversion_group.append(patient.ddr_bin)

# Statistical test
t_stat, p_value = ttest_ind(reversion_group, no_reversion_group)
```

---

### **Task P1.2: Test DDR_bin Stability with Low Variant Counts**

**Goal:** Confirm DDR_bin is stable when computed from ctDNA (fewer variants).

**Method:**
1. For each patient in Tier-3, randomly subsample to 5, 10, 20 variants
2. Recompute DDR_bin with subsampled variants
3. Measure variance (coefficient of variation)

**Acceptance Criteria:**
- CV < 0.15 (15% variation) at 10 variants → acceptable for ctDNA
- If CV > 0.30 → DDR_bin too unstable for clinical use

**Files:**
- Input: `data/validation/sae_cohort/checkpoints/Tier3_validation_cohort.json`
- Output: `data/validation/ddr_bin_stability_analysis.json`

```python
# Pseudocode
for patient in tier3_cohort:
    for n_variants in [5, 10, 20, 50]:
        subsampled_variants = random.sample(patient.variants, n_variants)
        ddr_bin = compute_ddr_bin(subsampled_variants)
        results[patient.id][n_variants].append(ddr_bin)

# Compute CV for each n_variants level
for n in [5, 10, 20, 50]:
    cv = std(results[:, n]) / mean(results[:, n])
```

---

### **Task P1.3: Build MAPK_bin and PI3K_bin**

**Goal:** Extend beyond DDR_bin to other pathway bins.

**Method:**
1. Define MAPK pathway genes: KRAS, NRAS, BRAF, NF1, MAP2K1, MAP2K2
2. Define PI3K pathway genes: PIK3CA, PTEN, AKT1, MTOR
3. For each patient, compute MAPK_bin (sum of SAE features for MAPK variants)
4. Compute PI3K_bin (sum of SAE features for PI3K variants)
5. Test if MAPK_bin/PI3K_bin correlate with resistance (like DDR_bin)

**Files:**
- Input: `data/validation/sae_cohort/checkpoints/Tier3_validation_cohort.json`
- Output: `api/resources/sae_feature_mapping.mapk_bin.v1.json`
- Output: `api/resources/sae_feature_mapping.pi3k_bin.v1.json`

**Acceptance Criteria:**
- At least one pathway bin (MAPK or PI3K) shows p < 0.05 for resistance
- If neither → we may not have enough signal in Tier-3 for these pathways

---

## 🎯 PRIORITY 2: LONGITUDINAL DATA HUNT (Week 1-2)

### **Task P2.1: Check TCGA-OV for Multiple Timepoints**

**Goal:** Find if TCGA-OV has serial samples (unlikely, but check).

**Method:**
1. Query TCGA-OV metadata for "sample_type" and "collection_date"
2. Identify patients with >1 sample
3. If found → compute DDR_bin trajectory

**Expected Result:**
- TCGA-OV is mostly single-timepoint (primary tumor only)
- May find a few recurrence samples

---

### **Task P2.2: Research Public Longitudinal Datasets**

**Goal:** Find datasets with serial ctDNA + outcomes.

**Datasets to Investigate:**

| Dataset | Cancer | Longitudinal? | ctDNA? | Outcomes? |
|---------|--------|---------------|--------|-----------|
| MMRF CoMMpass | Myeloma | ✅ Yes | ❓ Check | ✅ Yes |
| GRAIL CCGA | Multi-cancer | ✅ Yes | ✅ Yes | ✅ Yes |
| POG570 | Multi-cancer | ✅ Yes | ❓ Check | ✅ Yes |
| TRACERx | Lung | ✅ Yes | ✅ Yes | ✅ Yes |
| AACR GENIE | Multi-cancer | ❌ No | ❌ No | Partial |

**Method:**
1. For each dataset, check data access requirements
2. Check if raw sequence data available (needed for TRUE SAE)
3. Document findings

**Output:** `data/validation/longitudinal_data_sources.md`

---

### **Task P2.3: Compute DDR_bin Trajectories (If Data Found)**

**Goal:** If longitudinal data exists, compute DDR_bin over time.

**Method:**
1. For each patient with multiple timepoints
2. Compute DDR_bin at each timepoint
3. Correlate DDR_bin trajectory with time-to-progression
4. Identify threshold that predicts progression

**Output:** `data/validation/ddr_bin_trajectory_analysis.json`

---

## 🎯 PRIORITY 3: TECHNICAL BENCHMARKING (Week 2)

### **Task P3.1: Benchmark DDR_bin Computation Time**

**Goal:** Measure time and cost to compute DDR_bin for one patient.

**Method:**
1. Time the full pipeline: variant input → Evo2 → SAE → DDR_bin
2. Record Modal GPU cost
3. Test batch processing (10, 50, 100 patients at once)

**Metrics to Report:**
- Time per patient (single)
- Time per patient (batched)
- GPU cost per patient
- Total latency (variant input → DDR_bin output)

**Output:** `data/validation/ddr_bin_benchmark.json`

---

### **Task P3.2: Design Patient State Schema**

**Goal:** Define how to store longitudinal DDR_bin in patient state.

**Schema:**

```json
{
  "patient_id": "AK_001",
  "disease": "ovarian",
  "baseline_ddr_bin": 0.88,
  "ddr_bin_history": [
    {
      "date": "2025-01-15",
      "value": 0.88,
      "variant_count": 52,
      "model_version": "true_sae_diamonds.v1",
      "sample_type": "tumor",
      "alert_level": "NONE"
    },
    {
      "date": "2025-04-15",
      "value": 0.87,
      "variant_count": 48,
      "model_version": "true_sae_diamonds.v1",
      "sample_type": "ctdna",
      "alert_level": "NONE"
    },
    {
      "date": "2025-07-15",
      "value": 0.82,
      "variant_count": 55,
      "model_version": "true_sae_diamonds.v1",
      "sample_type": "ctdna",
      "alert_level": "WARNING",
      "delta_from_baseline": -0.06
    }
  ],
  "current_ddr_bin": 0.82,
  "trend": "DECLINING",
  "alert_triggered": true,
  "intervention_recommended": "COMBINATION_THERAPY"
}
```

**Files to Modify:**
- `api/models/patient_state.py` - Add DDR_bin fields
- `api/services/state_management.py` - Add DDR_bin tracking

---

### **Task P3.3: Implement Alert Thresholds**

**Goal:** Create alert logic based on DDR_bin.

**Logic:**

```python
def evaluate_ddr_bin_alert(current, baseline, history):
    delta = current - baseline
    
    if delta >= -0.03:
        return {"alert": "NONE", "message": "DDR_bin stable"}
    
    elif delta >= -0.08:
        return {
            "alert": "WARNING",
            "message": f"DDR_bin dropped {abs(delta)*100:.0f}% from baseline",
            "recommendation": "Consider combination therapy",
            "urgency": "HIGH"
        }
    
    elif delta >= -0.15:
        return {
            "alert": "ALERT",
            "message": f"DDR_bin dropped {abs(delta)*100:.0f}% - significant HR restoration",
            "recommendation": "Escalate therapy or enroll in trial",
            "urgency": "CRITICAL"
        }
    
    else:  # delta < -0.15
        return {
            "alert": "CRITICAL",
            "message": "DDR_bin severely depressed - resistant clone likely dominant",
            "recommendation": "Switch to platinum-based therapy",
            "urgency": "URGENT"
        }
```

**Files to Create:**
- `api/services/prevention_alert_service.py`

---

## 🎯 PRIORITY 4: SIMULATION & MODELING (Week 3)

### **Task P4.1: Model DDR_bin → Resistant Clone Relationship**

**Goal:** Create a mathematical model linking DDR_bin to resistant clone fraction.

**Assumptions:**
- Sensitive cells have DDR_bin ≈ 0.90 (high, HR-deficient)
- Resistant cells have DDR_bin ≈ 0.50 (low, HR-proficient)
- Observed DDR_bin = weighted average

**Model:**

```python
def estimate_resistant_clone(observed_ddr_bin, 
                              sensitive_ddr_bin=0.90, 
                              resistant_ddr_bin=0.50):
    """
    observed = (1 - f) * sensitive + f * resistant
    f = (sensitive - observed) / (sensitive - resistant)
    """
    f = (sensitive_ddr_bin - observed_ddr_bin) / (sensitive_ddr_bin - resistant_ddr_bin)
    return max(0, min(1, f))  # Clamp to [0, 1]

# Example:
# observed_ddr_bin = 0.82, sensitive = 0.90, resistant = 0.50
# f = (0.90 - 0.82) / (0.90 - 0.50) = 0.08 / 0.40 = 0.20 (20% resistant)
```

**Validation:**
- Apply to Tier-3 cohort
- Check if estimated resistant fraction correlates with resistance outcome

---

### **Task P4.2: Simulate Intervention Scenarios**

**Goal:** Model tumor dynamics under different intervention strategies.

**Model (Lotka-Volterra with treatment):**

```python
def simulate_tumor(initial_sensitive, initial_resistant, 
                   treatment_effect, fitness_cost, days):
    """
    dS/dt = r_s * S * (1 - (S+R)/K) - d_s * treatment * S
    dR/dt = r_r * R * (1 - (S+R)/K) - d_r * treatment * R
    
    Where:
    - r_s, r_r = growth rates (r_r < r_s due to fitness cost)
    - K = carrying capacity
    - d_s = treatment effect on sensitive (high for PARP)
    - d_r = treatment effect on resistant (low for PARP)
    """
    # Implement ODE solver
    pass
```

**Scenarios to Simulate:**
1. Continuous PARP (standard care)
2. PARP + carboplatin at DDR_bin drop
3. Adaptive therapy (dose reduction at DDR_bin drop)
4. Evolutionary steering (PARP → ATR → PARP cycling)

**Output:** `data/simulation/intervention_scenarios.json`

---

### **Task P4.3: Build Multi-Bin Decision Engine**

**Goal:** Create decision logic that considers DDR_bin + MAPK_bin + PI3K_bin.

```python
def recommend_intervention(ddr_bin, mapk_bin, pi3k_bin, efflux_bin):
    """
    Multi-pathway decision logic.
    """
    if ddr_bin < 0.80 and mapk_bin < 0.20 and pi3k_bin < 0.20:
        # Pure DDR escape
        return {
            "primary": "COMBINATION_THERAPY",
            "mechanism": "HR restoration",
            "alternative": "ATR inhibitor"
        }
    
    elif ddr_bin >= 0.80 and mapk_bin > 0.30:
        # MAPK bypass
        return {
            "primary": "ADD_MEK_INHIBITOR",
            "mechanism": "MAPK pathway activation",
            "alternative": "Trametinib + platinum"
        }
    
    elif ddr_bin >= 0.80 and pi3k_bin > 0.30:
        # PI3K bypass
        return {
            "primary": "ADD_PI3K_INHIBITOR",
            "mechanism": "PI3K pathway activation",
            "alternative": "Alpelisib + platinum"
        }
    
    elif efflux_bin > 0.30:
        # Drug efflux
        return {
            "primary": "SWITCH_NON_PARP",
            "mechanism": "Drug efflux upregulation",
            "alternative": "Immunotherapy or anti-VEGF"
        }
    
    else:
        return {"primary": "CONTINUE_CURRENT", "mechanism": "No escape detected"}
```

**Files to Create:**
- `api/services/prevention_recommendation_service.py`

---

## 📊 DELIVERABLES CHECKLIST

### **Week 1:**
- [ ] P1.1: Resistance mutation → DDR_bin analysis
- [ ] P1.2: DDR_bin stability with low variant counts
- [ ] P2.1: TCGA-OV longitudinal check
- [ ] P2.2: Longitudinal data sources research

### **Week 2:**
- [ ] P1.3: MAPK_bin and PI3K_bin extraction
- [ ] P3.1: DDR_bin computation benchmark
- [ ] P3.2: Patient state schema for DDR_bin
- [ ] P3.3: Alert threshold implementation

### **Week 3:**
- [ ] P4.1: DDR_bin → resistant clone model
- [ ] P4.2: Intervention scenario simulation
- [ ] P4.3: Multi-bin decision engine

---

## 🚨 BLOCKERS FOR MANAGER

Before Plumber can proceed with certain tasks, Manager needs to answer:

1. **Q1.1:** What DDR_bin thresholds should we use? (Plumber can propose, but need validation)
2. **Q2.1:** Is there ANY clinical data on DDR_bin response to combination therapy?
3. **Q5.1:** Where can we get longitudinal data? (Plumber researched options, need decision)

---

*Document Owner: Plumber*  
*Last Updated: December 24, 2025*  
*Status: ⏳ READY FOR EXECUTION*

