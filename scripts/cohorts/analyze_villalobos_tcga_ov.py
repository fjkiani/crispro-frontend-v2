#!/usr/bin/env python3
"""
VILLALOBOS 2018 TCGA-OV REANNOTATED DATASET - TIMING ENGINE VALIDATION

This script extracts and validates the perfect dataset we already have:
- 603 patients with detailed treatment timelines
- 504 patients with platinum-free interval (PFI) data
- Up to 6 lines of therapy per patient
- Time-to-treatment failure (TTF) for each regimen
"""

import pandas as pd
from pathlib import Path

print("="*80)
print("🏆 VILLALOBOS 2018 - TCGA-OV REANNOTATED DATASET")
print("PERFECT DATA FOR TIMING ENGINE VALIDATION!")
print("="*80)
print()

# Load data
data_file = Path("oncology-coPilot/oncology-backend-minimal/data/TCGA-OV/ds_cci.17.00096-1.xlsx")
df = pd.read_excel(data_file, sheet_name='Master clinical dataset')

print(f"📊 Total Patients: {len(df)}")
print()

# PFI Data Analysis
print("="*80)
print("1. PLATINUM-FREE INTERVAL (PFI) DATA")
print("="*80)
print()

pfi_col = 'Days off platinum prior to recurrence 1st line'
pfi_data = pd.to_numeric(df[pfi_col], errors='coerce').dropna()

print(f"Patients with PFI data: {len(pfi_data)}/{len(df)} ({100*len(pfi_data)/len(df):.1f}%)")
print(f"PFI Range: {pfi_data.min():.0f} - {pfi_data.max():.0f} days")
print(f"PFI Median: {pfi_data.median():.0f} days")
print()

# PFI Categories
resistant = (pfi_data < 180).sum()
partial = ((pfi_data >= 180) & (pfi_data < 365)).sum()
sensitive = (pfi_data >= 365).sum()

print("PFI Categories:")
print(f"  Resistant (<6 months):   {resistant:3d} patients ({100*resistant/len(pfi_data):5.1f}%)")
print(f"  Partial (6-12 months):   {partial:3d} patients ({100*partial/len(pfi_data):5.1f}%)")
print(f"  Sensitive (>12 months):  {sensitive:3d} patients ({100*sensitive/len(pfi_data):5.1f}%)")
print()

# Regimen Coverage
print("="*80)
print("2. TREATMENT REGIMEN DATA (LINE-BY-LINE)")
print("="*80)
print()

regimen_coverage = []
for i in range(1, 7):
    ordinal = {1: 'st', 2: 'nd', 3: 'rd'}.get(i, 'th')
    regimen_col = f'{i}{ordinal}_regimen'
    ttf_col = f'{i}{ordinal}_chemo_regimen_TTF_days'
    
    if regimen_col in df.columns:
        regimen_count = df[regimen_col].notna().sum()
        ttf_count = pd.to_numeric(df[ttf_col], errors='coerce').notna().sum() if ttf_col in df.columns else 0
        
        regimen_coverage.append({
            'Line': i,
            'Patients with Regimen': regimen_count,
            'Patients with TTF': ttf_count,
            'Regimen %': 100 * regimen_count / len(df),
            'TTF %': 100 * ttf_count / len(df) if ttf_count > 0 else 0
        })
        
        print(f"Line {i}: {regimen_count:3d} regimens ({100*regimen_count/len(df):5.1f}%), "
              f"{ttf_count:3d} TTF values ({100*ttf_count/len(df):5.1f}%)")

print()

# Survival Data
print("="*80)
print("3. SURVIVAL OUTCOMES")
print("="*80)
print()

survival_cols = {
    'total_days_overall_survival': 'Overall Survival',
    'days_to_tumor_progression': 'Time to Progression',
    'days_to_death': 'Time to Death',
    'days_to_last_followup': 'Last Follow-up'
}

for col, label in survival_cols.items():
    if col in df.columns:
        count = pd.to_numeric(df[col], errors='coerce').notna().sum()
        print(f"{label:25s}: {count:3d}/{len(df)} ({100*count/len(df):5.1f}%)")

print()

# Sample Data
print("="*80)
print("4. SAMPLE PATIENT DATA (First 5 Patients)")
print("="*80)
print()

sample_cols = [
    'bcr_patient_barcode',
    '1st_regimen',
    '1st_chemo_regimen_TTF_days',
    '2nd_regimen',
    '2nd_chemo_regimen_TTF_days',
    'Days off platinum prior to recurrence 1st line',
    'total_days_overall_survival'
]

sample_df = df[sample_cols].head(5)

for idx, row in sample_df.iterrows():
    print(f"Patient: {row['bcr_patient_barcode']}")
    print(f"  1st Regimen: {row['1st_regimen']}")
    print(f"  1st TTF: {row['1st_chemo_regimen_TTF_days']} days")
    print(f"  2nd Regimen: {row['2nd_regimen']}")
    print(f"  2nd TTF: {row['2nd_chemo_regimen_TTF_days']} days")
    print(f"  PFI: {row['Days off platinum prior to recurrence 1st line']} days")
    print(f"  OS: {row['total_days_overall_survival']} days")
    print()

# Summary
print("="*80)
print("✅ VALIDATION CAPABILITY SUMMARY")
print("="*80)
print()
print("WHAT WE CAN VALIDATE WITH THIS DATA:")
print()
print("✅ PFI (Platinum-Free Interval):")
print("   - 504 patients with ground-truth PFI values")
print("   - Can validate PFI computation accuracy")
print("   - Can validate PFI categorization (<6m, 6-12m, >12m)")
print()
print("✅ TTF (Time-to-Treatment Failure):")
print("   - Per-regimen TTF for lines 1-6")
print("   - Can validate treatment sequence timing")
print()
print("✅ Regimen Sequencing:")
print("   - Up to 6 lines of therapy per patient")
print("   - Can validate treatment-free intervals (TFI)")
print("   - Can validate line-of-therapy logic")
print()
print("✅ Survival Outcomes:")
print("   - OS, PFS, time-to-progression")
print("   - Can validate PFS/OS computation from regimen start")
print()
print("⚠️ LIMITATION:")
print("   - Lacks EXACT regimen start/end dates (has TTF instead)")
print("   - Cannot validate PTPI (platinum-to-PARPi) - TCGA pre-dates PARPi era")
print()
print("="*80)
print("RECOMMENDATION: USE THIS FOR PFI VALIDATION NOW!")
print("="*80)
print()
print("Next Steps:")
print("1. ✅ Extract PFI ground truth (504 patients)")
print("2. ✅ Run timing engine with inferred regimen dates")
print("3. ✅ Compare computed PFI vs ground truth PFI")
print("4. ✅ Generate validation report with accuracy metrics")
print()

# Save summary
output_file = Path("data/validation/timing_engine/villalobos_2018_summary.txt")
output_file.parent.mkdir(parents=True, exist_ok=True)

with open(output_file, 'w') as f:
    f.write(f"Villalobos 2018 TCGA-OV Reannotated Dataset Summary\n")
    f.write(f"="*80 + "\n\n")
    f.write(f"Total Patients: {len(df)}\n")
    f.write(f"Patients with PFI: {len(pfi_data)} ({100*len(pfi_data)/len(df):.1f}%)\n")
    f.write(f"\nPFI Categories:\n")
    f.write(f"  Resistant (<6m):  {resistant} ({100*resistant/len(pfi_data):.1f}%)\n")
    f.write(f"  Partial (6-12m):  {partial} ({100*partial/len(pfi_data):.1f}%)\n")
    f.write(f"  Sensitive (>12m): {sensitive} ({100*sensitive/len(pfi_data):.1f}%)\n")
    
print(f"✅ Summary saved to: {output_file}")
