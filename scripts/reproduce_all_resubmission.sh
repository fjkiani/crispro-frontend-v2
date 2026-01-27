#!/bin/bash
# Resubmission-safe reproduction script
# - Stages canonical submission-bundle inputs into repo root `publication/`
# - Runs the same validation suite as `scripts/reproduce_all.sh`
# - Avoids any possibility of silently generating synthetic datasets
#
# Usage:
#   ./scripts/reproduce_all_resubmission.sh
#
set -e

echo "================================================================================"
echo "METASTASIS INTERCEPTION - RESUBMISSION-SAFE REPRODUCTION"
echo "================================================================================"
echo ""

START_TIME=$(date +%s)

echo "Step 1/6: Environment Setup..."
if [ ! -d "venv" ]; then
  echo "  Creating virtual environment..."
  python3 -m venv venv
fi
echo "  Installing dependencies..."
venv/bin/pip install -q --upgrade pip
venv/bin/pip install -q -r requirements.txt
echo ""

echo "Step 2/6: Configuration..."
export PYTHONPATH=$(pwd)
export SEED=42
export EVO_FORCE_MODEL=evo2_1b
export EVO_USE_DELTA_ONLY=1
echo ""

echo "Step 3/6: Stage canonical publication inputs..."
venv/bin/python scripts/metastasis/stage_publication_inputs.py
echo ""

echo "Step 3.5/6: Update chromatin with Enformer (Option A)..."
venv/bin/python scripts/metastasis/update_target_lock_chromatin_enformer.py
echo ""

echo "Step 4/6: Validation Metrics (Day 1)..."
venv/bin/python scripts/metastasis/compute_per_step_validation.py
venv/bin/python scripts/metastasis/compute_specificity_matrix.py
venv/bin/python scripts/metastasis/compute_precision_at_k.py
venv/bin/python scripts/metastasis/compute_ablation_study.py
venv/bin/python scripts/metastasis/compute_confounder_analysis.py
echo ""

echo "Step 5/6: Enhanced Validation (Day 2)..."
venv/bin/python scripts/metastasis/generate_calibration_curves.py
venv/bin/python scripts/metastasis/compute_effect_sizes.py
venv/bin/python scripts/metastasis/generate_table_s2.py
echo ""

echo "Step 6/6: Verification..."
EXPECTED_FILES=(
  "publication/figures/figure2a_per_step_roc.png"
  "publication/figures/figure2b_specificity_matrix.png"
  "publication/figures/figure2c_precision_at_k.png"
  "publication/figures/figure2d_ablation.png"
  "publication/figures/figure_s1_confounders.png"
  "publication/figures/figure_s2_calibration_curves.png"
  "publication/figures/figure_s3_effect_sizes.png"
  "publication/data/per_step_validation_metrics.csv"
  "publication/data/specificity_enrichment.csv"
  "publication/data/precision_at_k.csv"
  "publication/data/ablation_study.csv"
  "publication/data/chromatin_audit_enformer_38genes.csv"
  "publication/data/confounder_analysis.csv"
  "publication/data/effect_sizes.csv"
  "publication/tables/table_s2_validation_metrics.csv"
)

MISSING=0
for file in "${EXPECTED_FILES[@]}"; do
  if [ -f "$file" ]; then
    echo "  OK: $file"
  else
    echo "  MISSING: $file"
    MISSING=$((MISSING + 1))
  fi
done

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
MINUTES=$((DURATION / 60))
SECONDS=$((DURATION % 60))

echo ""
echo "================================================================================"
echo "DONE"
echo "================================================================================"
echo "Time elapsed: ${MINUTES}m ${SECONDS}s"
if [ $MISSING -ne 0 ]; then
  echo "WARNING: $MISSING expected output file(s) missing"
  exit 1
fi
echo "All expected files generated."

