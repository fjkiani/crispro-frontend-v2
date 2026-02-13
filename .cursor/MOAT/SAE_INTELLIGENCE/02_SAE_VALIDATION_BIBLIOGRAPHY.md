# SAE VALIDATION BIBLIOGRAPHY: VARIANT-LEVEL REPLACEMENT

**Mission**: Prove variant-level SAE features > gene-level markers (Oncotype, PD-L1, TMB).
**Status**: STRATEGIC ROADMAP

## 1. THE CORE DIFFERENTIATOR

**SAE = Variant-Level Biology**
| Traditional | SAE Approach | Improvement |
|---|---|---|
| BRCA mutated (Yes/No) | BRCA p.R1699Q SAE features | +15.5 pp AUROC (0.783 vs 0.628) |
| Gene-level markers | Evo2 Embeddings -> SAE (32K) | Platform Moat |

---

## 2. CURRENT STATUS

### ✅ Ovarian Cancer (TCGA-OV)
- **Status**: PRODUCTION READY (149 extracted).
- **Pathway**: DDR (AUROC 0.783).
- **Replaces**: HRD Test ($3,000, 60% acc).
- **Next**: Extract remaining 7 pathways (PI3K, MAPK, Immune, etc.).

### ❓ Breast Cancer (TCGA-BRCA)
- **Target**: **Oncotype DX** ($500M market).
- **Status**: Checkpoint exists (validation pending).
- **Priority**: HIGHEST.

---

## 3. VALIDATION ROADMAP

### PHASE 1: Breast (Oncotype DX Replacement)
- **Market**: $500M.
- **Baseline**: 21-gene recurrence score (AUROC ~0.65).
- **SAE Target**: AUROC 0.75 (+10pp).
- **Plan**:
    1. Validate BRCA checkpoint.
    2. Extract missing pathways (DDR, Proliferation, Immune).
    3. Train SAE Model vs Recurrence.

### PHASE 2: Lung (PD-L1/EGFR Replacement)
- **Market**: $350M.
- **Baseline**: PD-L1 IHC (55% acc).
- **SAE Target**: AUROC 0.72 (+14pp).
- **Plan**: Extract TCGA-LUAD (Immune, MAPK pathways).

### PHASE 3: Melanoma (TMB Replacement)
- **Market**: $300M.
- **Baseline**: TMB (Mutation Count).
- **SAE Target**: AUROC 0.74 (+11pp).
- **Plan**: Extract TCGA-SKCM (Immune Pathway).

---

## 4. AGENT TASK TEMPLATE

```markdown
MISSION: Prove SAE variant-level > Gene-level [TEST]
DATASET: [TCGA-ID]
1. Extract SAE [PATHWAYS].
2. Compute [BASELINE] score.
3. Train SAE -> [OUTCOME].
4. Compare AUROC.
```
