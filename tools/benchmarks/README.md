# Benchmarks

## Variant AUROC/AUPRC (ClinVar)

- Script: `tools/benchmarks/variant_auroc.py`
- Usage:
```bash
# Start backend first (and optional chromatin proxies)
make backend

# Run 200/200 balanced SNVs (includes likely_* and GRCh37 fallback)
make bench_variant

# Custom run
python3 tools/benchmarks/variant_auroc.py \
  --download --n_pos 1000 --n_neg 1000 --allow_grch37 \
  --api_base http://127.0.0.1:8000 --model_id evo2_7b \
  --out tools/benchmarks/variant_auroc_results.json
```
- Output: JSON with AUROC/AUPRC and run config

## Guidance tier evaluation (planned)
- Script to be added: `tools/benchmarks/guidance_tier_eval.py`
- Compares Tier I vs guideline/on‑label truth (FDA/DailyMed + curated rules)
- Metrics: Precision/Recall/F1; optional stratification by disease/therapy

## HRD/Platinum AUPRC (planned)
- Script to be added: `tools/benchmarks/hrd_platinum_auprc.py`
- Computes damage + dependency → therapy mapping vs platinum outcomes
- Metrics: AUPRC/AUROC; ablations ±Chromatin/±Essentiality lifts


