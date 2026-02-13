# RESISTANCE PREDICTION: SAE & DNA REPAIR LOGIC

**Topic**: SAE Feature Extraction & DNA Repair Capacity Formulas.

## 1. SAE FEATURE EXTRACTIONMODES

### PROXY SAE (Production Default)
**Location**: `api/services/sae_feature_service.py`
**Inputs**: Insights Bundle (Chips), Pathway Scores (Efficacy).
**Logic**:
1. **Pathway Burden**: Directly from `pathway_scores` (e.g., DDR, MAPK).
2. **Essentiality**: Average of HRR genes (BRCA1, etc.) from `insights_bundle`.
3. **Exon Disruption**: Proxy via `regulatory` score.
4. **Vector**: Maps scores to 7D mechanism vector.
**Provenance**: `"sae": "proxy"`

### TRUE SAE (Research/Future)
**Location**: `api/services/sae_feature_service.py`
**Inputs**: 32K-dim feature vector from Evo2 (Layer 26).
**Logic**:
1. **Extraction**: Call Modal service/Evo2.
2. **Mapping**: Map 32K features to pathways via `sae_feature_mapping.json`.
3. **Scoring**: Aggregate activations per pathway.
**Provenance**: `"sae": "true_sae"`

---

## 2. DNA REPAIR CAPACITY FORMULA

**Weights** (Manager Approved):
- **DDR Pathway Burden**: 60%
- **HRR Essentiality**: 20%
- **Exon Disruption**: 20%

**Formula**:
```python
dna_repair_capacity = (
    0.60 * pathway_burden_ddr +
    0.20 * essentiality_hrr_genes +
    0.20 * exon_disruption_score
)
```

**Interpretation**:
- **High (≥0.70)**: Intact Repair -> PARP Resistance.
- **Low (<0.40)**: Deficient Repair -> PARP Sensitivity (Synthetic Lethality).

---

## 3. RESISTANCE SIGNAL: 2-OF-3 RULE

**Trigger Condition**: High Risk if **≥2 signals** detected AND **Prob ≥ 0.70**.

### Signal 1: DNA Repair Restoration
- **Logic**: `current_repair - baseline_repair < -0.15` (Capacity restoring/dropping deficiency).
- **Meaning**: Tumor fixing its DNA repair to survive PARP.

### Signal 2: Pathway Escape
- **Logic**: Target pathway burden drops ≥ 0.15.
- **Meaning**: Tumor switching dependence away from drug target.

### Signal 3: CA-125 Kinetics
- **Logic**: Rustin criteria (Doubling, Rising, Monotonic).
- **Meaning**: Systemic tumor burden increase.
