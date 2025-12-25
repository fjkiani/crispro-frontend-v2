# TRUE SAE "DIAMONDS" — EXTRACTION + MAPPING WORKPACK (PLUMBER)

**Date:** 2025-12-24  
**Owner:** Plumber  
**Status:** ✅ **COMPLETE** - Deliverables verified (mapping + baseline exist)  
**Verification:** See `AUDIT_REPORT.md` for systematic review

This is a **modular workpack** derived from `.cursor/MOAT/SAE_FAILURE_POSTMORTEM.mdc` and the Tier‑3 artifacts. The goal is to **mine** the existing Tier‑3 cohort for (a) defensible **feature→biology mappings** and (b) a **repeatable predictive baseline** that justifies Steerability V1 (hybrid bins) without overclaiming.

---

## Key artifacts (source of truth)

- **Tier‑3 cohort (patient + variant + top_features):** `data/validation/sae_cohort/checkpoints/Tier3_validation_cohort.json`
- **Tier‑3 stats + large/medium effect feature lists:** `data/validation/sae_cohort/checkpoints/sae_validation_results.json`
- **Tier‑3 report:** `data/validation/sae_cohort/checkpoints/sae_validation_report.txt`

---

## Critical label definition (do not get this wrong)

Tier‑3 cohort outcomes are:
- `sensitive`: 125
- `refractory`: 17
- `resistant`: 7

**For analysis:** treat **`refractory + resistant = resistant`** (pos class = 24), which matches `sae_validation_results.json`.

If you only use `resistant` (n=7), AUROC becomes unstable and you will falsely conclude "no signal".

---

## The "diamonds" (large-effect, higher-in-resistant)

From `sae_validation_results.json → analysis.large_effect_list` (direction=`higher_in_resistant`), prioritize:

- **27607** (d=0.635, p=0.0146)
- **16337** (d=0.634, p=0.0247)
- **26220** (d=0.609, p=0.0215)
- **12893** (d=0.597, p=0.0246)
- **6020** (d=0.574, p=0.0324)
- **22868** (d=0.544, p=0.0355)
- **1407** (d=0.537, p=0.0414)
- **9738** (d=0.530, p=0.0495)
- **31362** (d=0.517, p=0.0466)

---

## Quant check (already measured; reproduce it)

Using patient-level aggregation = **sum(feature_value across all variants where feature appears)**:

- **Single-feature AUROC** is typically ~0.62–0.65 for the "diamonds".
- **Multi-feature baseline** (logistic regression on 29 candidate features, 5-fold CV): **mean AUROC ≈ 0.78**

Your job is to **make this reproducible** (script + exact feature list + seed) so this becomes an internal "go/no-go" gate for Steerability V1.

---

## Deliverable A — feature→biology mapping (minimum viable)

### Goal
Produce a *defensible*, minimal mapping (even if coarse) that supports **bin-level steerability**, not feature-level claims.

### Output file (deliver separate)
Create:

`api/resources/sae_feature_mapping.true_sae_diamonds.v1.json`

Schema (start here; extend if needed):

```json
{
  "version": "true_sae_diamonds.v1",
  "source_artifacts": {
    "tier3_cohort": "data/validation/sae_cohort/checkpoints/Tier3_validation_cohort.json",
    "tier3_results": "data/validation/sae_cohort/checkpoints/sae_validation_results.json"
  },
  "labeling": { "positive_class": ["resistant", "refractory"] },
  "features": [
    {
      "feature_index": 27607,
      "direction": "higher_in_resistant",
      "effect_size_d": 0.635,
      "p_value": 0.01455,
      "mapping": {
        "hypothesis": "DDR_bin | MAPK_bin | PI3K_bin | OTHER",
        "confidence": "low|medium|high",
        "rationale": [
          "Top activating variants are enriched for gene X",
          "Correlation with PROXY pathway burden Y"
        ]
      },
      "evidence": {
        "top_patients": ["TCGA-..", "TCGA-.."],
        "top_variants": [
          { "patient_id": "TCGA-..", "gene": "NF1", "variant": "…", "value": 1.23 }
        ]
      }
    }
  ]
}
```

### Method (how to mine the mapping)

- **Step A1 — Identify top-activating variants**
  - For each diamond feature, scan Tier‑3 cohort and collect the **top K variant occurrences** by feature value.
  - Output a table: `feature → {gene counts, top variants}`.

- **Step A2 — Gene/pathway enrichment (coarse)**
  - Are the top-activating variants enriched for DDR genes (BRCA1/2, ATM/ATR, CHEK2), MAPK (NF1, KRAS, BRAF), PI3K (PIK3CA, PTEN)?
  - This can be **heuristic** at first; we want a mapping hypothesis, not a publication.

- **Step A3 — Proxy correlation sanity check**
  - For each patient, compute a simple proxy signal (e.g., DDR marker count, MAPK marker presence).
  - Correlate feature value (per patient) with the proxy score.
  - If correlation is directional and stable → upgrade confidence.

**Acceptance criteria for A:**
- ≥ **3** diamond features mapped into **DDR/MAPK/PI3K/OTHER** bins with explicit evidence (top variants list).

---

## Deliverable B — reproducible predictive baseline (go/no-go gate)

### Goal
Establish "TRUE SAE has signal" **without** relying on FDR significance.

### Output file
Create:

`data/validation/sae_cohort/checkpoints/true_sae_diamonds_baseline.v1.json`

Include:
- feature list used
- exact label definition
- aggregation choice
- CV splits seed
- AUROC per fold + mean

**Acceptance criteria for B:**
- Mean AUROC ≥ **0.70** with pos class = refractory+resistant  
- Report fold variance; if variance is huge, propose stabilization (bootstraps, repeated CV)

---

## "Do not do this" (waste)

- Do **not** re-run extraction.
- Do **not** chase FDR significance as a requirement to begin mapping.
- Do **not** claim monosemantic biology for a feature without evidence from the cohort.

---

## Open questions (ask Zo if needed)

1. Should "refractory" be treated as resistant **globally** across the pipeline, or only inside the Tier‑3 experiment?
2. Do we want the mapping bins to align to the 7D mechanism vector dims (DDR/MAPK/PI3K/VEGF/HER2/IO/Efflux), or start narrower (DDR/MAPK/PI3K/OTHER)?

