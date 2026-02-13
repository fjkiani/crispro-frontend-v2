# SAE Proxy Architecture Deep Dive

**Status**: CURRENT PRODUCTION IMPLEMENTATION
**Reference**: Derived from `SAE_CONTEXT_MASTER_REFERENCE.md`

## 1. The Critical Distinction

While the "True SAE" extracts features from Evo2's neural activations, the **Production Proxy SAE** (`sae_feature_service.py`) operates by synthesizing known metadata into a compatible format.

**Gap Analysis**:
*   **True SAE**: `Evo2 Activations -> Sparse Features -> Pathway Score`
*   **Proxy SAE**: `Pathway Aggregation + Insights Bundle -> Proxy Features -> Pathway Score`

## 2. Component Implementation

### A. DNA Repair Capacity
computed via Manager's "C1 Formula":
```python
dna_repair_capacity = (
    0.6 * pathway_burden_ddr +           # From Pathway Aggregation
    0.2 * essentiality_hrr_genes +       # From Insights Bundle (Essentiality Chip)
    0.2 * exon_disruption_score          # From Insights Bundle (Regulatory Chip)
)
```
*   **Thresholds**:
    *   High (≥0.70): PARP Sensitive.
    *   Low (<0.40): PARP Resistant.

### B. Mechanism Vector (7D)
Constructed from the `Pathways` module outputs + Tumor Context:
```python
mechanism_vector = [
    pathway_burden_ddr,       # Index 0
    pathway_burden_mapk,      # Index 1
    pathway_burden_pi3k,      # Index 2
    pathway_burden_vegf,      # Index 3
    pathway_burden_her2,      # Index 4
    io_eligibility_score,     # Index 5 (Derived from TMB/MSI)
    cross_resistance_risk     # Index 6 (Derived from drug history)
]
```

### C. Resistance Triggers ("2-of-3" Rule)
The system flags "High Risk" if 2 of these 3 conditions are met:
1.  **HRD Drop**: ≥10 points vs baseline.
2.  **Capacity Drop**: DNA Repair Capacity decreases ≥0.15.
3.  **Kinetics**: CA-125 fails to drop 50% by Cycle 3.

## 3. Data Flow
1.  **Input**: Patient JSON (Variants, TMB, Treatment History).
2.  **Step 1**: `PathwaysModule` computes burdens (0.0 - 1.0).
3.  **Step 2**: `InsightsModule` computes essentiality/functionality scores.
4.  **Step 3**: `SaeFeatureService` ingests these outputs to assemble the `SaeBundle`.
5.  **Output**: UI receives `mechanism_chips`, `hint_tiles`, and `dna_repair_capacity`.
