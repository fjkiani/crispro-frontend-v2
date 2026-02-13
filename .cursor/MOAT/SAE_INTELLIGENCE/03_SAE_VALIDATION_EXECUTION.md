# SAE VALIDATION EXECUTION PLAN

**Strategy**: Use `extract_true_sae_cohort_from_cbioportal.py` (Script 2).
**Protocol**: 72-hour sprints per target.

## PHASE 0: PREPARATION

### 0.1 Validate BRCA Checkpoint
- Locate `BRCA_TCGA_TRUE_SAE_cohort.json`.
- Count patients/variants.
- Verify recurrence labels.

### 0.2 Infrastructure Prep
- Verify Modal Service health (`/health`).
- Set Feature Flags (`ENABLE_EVO2_SAE=1`).
- Run Preflight Check.

---

## PHASE 1: BREAST CANCER (Oncotype DX)

### Day 1: Extraction
**Script**: `extract_true_sae_cohort_from_cbioportal.py`
**Config**:
- Study: `brca_tcga`
- Budget: 72,000s
- Cap: 100 patients, 50 variants/pat.
**Output**: `BRCA_TCGA_TRUE_SAE_cohort.json`

### Day 2: Training
1. **Baseline**: Compute Oncotype DX score (21-gene expression).
2. **SAE Model**: Logistic Regression (L1) on 32K SAE features.
3. **Compare**: AUROC SAE vs Oncotype.

### Day 3: Validation
- Generate Validation Receipt (JSON).
- Update Bibliography.

---

## PHASE 2: LUNG CANCER (PD-L1/EGFR)

### Day 1: Extraction
**Config**:
- Study: `luad_tcga`
- Outcomes: IO Response, EGFR Response.

### Day 2: Training
1. **Model 1**: SAE vs PD-L1 (IO Response).
2. **Model 2**: SAE vs EGFR Mutation (Inhibitor Response).

---

## INFRASTRUCTURE UPGRADES

**Standardized Script**: `scripts/sae/extract_validation_cohort.py`
- Preflight assembly detection.
- Circuit breaker (30% error).
- Atomic writes.
- JSONL Error logging.

**Monitoring**:
- Real-time extraction dashboard.
- Cost tracking.
