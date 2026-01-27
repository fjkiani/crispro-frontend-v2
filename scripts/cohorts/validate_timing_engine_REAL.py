#!/usr/bin/env python3
"""
TIMING ENGINE PFI VALIDATION - TRULY NO DATA LEAKAGE

This version uses INDEPENDENT data sources:
- Last platinum date from clinical data
- Days to recurrence (NOT derived from PFI)
- Computes PFI using timing engine
- Compares to ground truth PFI

NO CIRCULAR DEPENDENCY.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
import sys

sys.path.insert(0, 'oncology-coPilot/oncology-backend-minimal')
from api.services.resistance.biomarkers.therapeutic.timing_chemo_features import build_timing_chemo_features

print("="*80)
print("TIMING ENGINE VALIDATION - TRULY INDEPENDENT DATA")
print("NO CIRCULAR DEPENDENCY")
print("="*80)
print()

# Load data
df = pd.read_excel('oncology-coPilot/oncology-backend-minimal/data/TCGA-OV/ds_cci.17.00096-1.xlsx', 
                   sheet_name='Master clinical dataset')

# Independent columns
platinum_col = 'Last day of platinum 1st line'  # Days from diagnosis to last platinum
progression_col = 'days_to_tumor_recurrence'     # Days from diagnosis to recurrence
pfi_gt_col = 'Days off platinum prior to recurrence 1st line'  # Ground truth

# Filter to patients with INDEPENDENT data
platinum_days = pd.to_numeric(df[platinum_col], errors='coerce')
recurrence_days = pd.to_numeric(df[progression_col], errors='coerce')
pfi_gt = pd.to_numeric(df[pfi_gt_col], errors='coerce')

valid_mask = platinum_days.notna() & recurrence_days.notna() & pfi_gt.notna()
valid_df = df[valid_mask].copy()

print(f"Patients with independent data: {len(valid_df)}")
print(f"  - Last platinum date: from {platinum_col}")
print(f"  - Recurrence date: from {progression_col}")
print(f"  - Ground truth PFI kept separate for comparison")
print()

# Create regimen table with INDEPENDENT dates
regimen_table = []
survival_table = []
clinical_table = []

base_date = datetime(2010, 1, 1)  # Arbitrary baseline

for idx, row in valid_df.iterrows():
    patient_id = row['bcr_patient_barcode']
    
    plat_day = int(pd.to_numeric(row[platinum_col], errors='coerce'))
    recur_day = int(pd.to_numeric(row[progression_col], errors='coerce'))
    
    # Dates derived from INDEPENDENT columns, NOT from PFI
    last_platinum_date = base_date + timedelta(days=plat_day)
    progression_date = base_date + timedelta(days=recur_day)
    
    # Estimate regimen start/end
    regimen_start = last_platinum_date - timedelta(days=168)  # ~6 cycles back
    regimen_end = last_platinum_date + timedelta(days=21)
    
    regimen_table.append({
        'patient_id': patient_id,
        'regimen_id': f'{patient_id}_R1',
        'regimen_start_date': regimen_start,
        'regimen_end_date': regimen_end,
        'regimen_type': 'platinum',
        'line_of_therapy': 1,
        'setting': 'frontline',
        'last_platinum_dose_date': last_platinum_date,
        'progression_date': progression_date,
    })
    
    survival_table.append({
        'patient_id': patient_id,
        'vital_status': 'Unknown',
        'last_followup_date': progression_date + timedelta(days=365),
    })
    
    clinical_table.append({
        'patient_id': patient_id,
        'disease_site': 'ovary',
    })

print(f"Created {len(regimen_table)} regimens using INDEPENDENT dates")
print()

# Run timing engine
print("="*80)
print("Running timing engine...")
print("="*80)
print()

results = build_timing_chemo_features(
    regimen_table=regimen_table,
    survival_table=survival_table,
    clinical_table=clinical_table,
)

print(f"Timing engine completed: {len(results)} regimens")
print()

# Extract computed PFI
computed = pd.DataFrame(results)
computed = computed[['patient_id', 'PFI_days', 'PFI_category']]

# Get ground truth (kept separate)
ground_truth = valid_df[['bcr_patient_barcode', pfi_gt_col]].copy()
ground_truth.columns = ['patient_id', 'ground_truth_pfi']

# Merge
comparison = computed.merge(ground_truth, on='patient_id')

# Calculate accuracy
comparison['diff'] = comparison['PFI_days'] - comparison['ground_truth_pfi']
comparison['abs_diff'] = comparison['diff'].abs()

# Category comparison
def categorize(days):
    if pd.isna(days) or days < 0:
        return None
    elif days < 180:
        return '<6m'
    elif days < 365:
        return '6-12m'
    else:
        return '>12m'

comparison['computed_cat'] = comparison['PFI_days'].apply(categorize)
comparison['gt_cat'] = comparison['ground_truth_pfi'].apply(categorize)
comparison['cat_match'] = comparison['computed_cat'] == comparison['gt_cat']

# Results
print("="*80)
print("VALIDATION RESULTS - INDEPENDENT DATA")
print("="*80)
print()

print(f"Total patients: {len(comparison)}")
print()

# PFI accuracy
exact_match = (comparison['abs_diff'] == 0).sum()
within_7 = (comparison['abs_diff'] <= 7).sum()
within_30 = (comparison['abs_diff'] <= 30).sum()

print("PFI DAYS ACCURACY:")
print(f"  Exact match (0 days):  {exact_match}/{len(comparison)} ({100*exact_match/len(comparison):.1f}%)")
print(f"  Within 7 days:         {within_7}/{len(comparison)} ({100*within_7/len(comparison):.1f}%)")
print(f"  Within 30 days:        {within_30}/{len(comparison)} ({100*within_30/len(comparison):.1f}%)")
print(f"  Mean absolute error:   {comparison['abs_diff'].mean():.1f} days")
print(f"  Median absolute error: {comparison['abs_diff'].median():.1f} days")
print()

# Category accuracy
cat_correct = comparison['cat_match'].sum()
print("CATEGORY ACCURACY:")
print(f"  Correct category: {cat_correct}/{len(comparison)} ({100*cat_correct/len(comparison):.1f}%)")
print()

# Confusion matrix
print("CONFUSION MATRIX:")
print(pd.crosstab(comparison['gt_cat'], comparison['computed_cat'], 
                  rownames=['Ground Truth'], colnames=['Computed'], margins=True))
print()

# Sample discrepancies
discrepancies = comparison[comparison['abs_diff'] > 0].head(10)
if len(discrepancies) > 0:
    print("SAMPLE DISCREPANCIES (first 10):")
    print()
    for idx, row in discrepancies.iterrows():
        print(f"  {row['patient_id']}")
        print(f"    Computed: {row['PFI_days']:.0f} days ({row['computed_cat']})")
        print(f"    Ground truth: {row['ground_truth_pfi']:.0f} days ({row['gt_cat']})")
        print(f"    Difference: {row['diff']:.0f} days")
        print()

# Save
output_dir = Path('oncology-coPilot/oncology-backend-minimal/api/services/resistance')
comparison.to_csv(output_dir / 'timing_engine_real_validation.csv', index=False)

print("="*80)
print("CONCLUSION")
print("="*80)
print()

if cat_correct / len(comparison) >= 0.95:
    print("✅ VALIDATION PASSED - Category accuracy ≥95%")
elif cat_correct / len(comparison) >= 0.90:
    print("⚠️ VALIDATION MARGINAL - Category accuracy 90-95%")
else:
    print("❌ NEEDS INVESTIGATION - Category accuracy <90%")
print()
print(f"Results saved to: {output_dir / 'timing_engine_real_validation.csv'}")
