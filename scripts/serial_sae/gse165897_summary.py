#!/usr/bin/env python3
"""Generate summary report for GSE165897 analysis."""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime

# Load pathway kinetics
kinetics_path = Path("data/serial_sae/gse165897/results/pathway_kinetics.csv")
results_dir = Path("data/serial_sae/gse165897/results")

if not kinetics_path.exists():
    print("❌ Pathway kinetics file not found. Run pathway_kinetics_gse165897.py first.")
    exit(1)

df = pd.read_csv(kinetics_path)

print("="*80)
print("GSE165897 Pathway Kinetics Summary")
print("="*80)
print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

print(f"Patients analyzed: {len(df)}")
print(f"Pathways: DDR, MAPK, PI3K, VEGF\n")

# Summary statistics
delta_cols = [c for c in df.columns if '_delta' in c]
print("Pathway Kinetics (Δ = post - pre):")
print(df[delta_cols].describe().round(4))

print("\n" + "="*80)
print("Key Findings:")
print("="*80)

# Count significant changes
for col in delta_cols:
    pathway = col.replace('_delta', '').upper()
    increases = (df[col] > 0.05).sum()
    decreases = (df[col] < -0.05).sum()
    print(f"\n{pathway}:")
    print(f"  Increases (Δ > 0.05): {increases}/{len(df)} patients")
    print(f"  Decreases (Δ < -0.05): {decreases}/{len(df)} patients")
    print(f"  Mean Δ: {df[col].mean():.4f}")
    print(f"  Median Δ: {df[col].median():.4f}")

print("\n" + "="*80)
print("Status:")
print("="*80)
print("✅ Data downloaded and processed")
print("✅ Pathway kinetics computed (pre/post-NACT)")
print("⚠️  Treatment response labels not available in metadata")
print("\nNext Steps:")
print("1. Extract treatment response labels from publication (RECIST criteria)")
print("2. Correlate pathway kinetics with treatment response")
print("3. Identify predictive patterns in pre-treatment pathway scores")
print("4. Compare with TCGA-OV findings")

# Save summary
summary_path = results_dir / "gse165897_analysis_summary.txt"
with open(summary_path, "w") as f:
    f.write("="*80 + "\n")
    f.write("GSE165897 Analysis Summary\n")
    f.write("="*80 + "\n")
    f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    f.write(f"Patients: {len(df)}\n")
    f.write(f"Pathways: DDR, MAPK, PI3K, VEGF\n\n")
    f.write("Pathway Kinetics Summary:\n")
    f.write(df[delta_cols].describe().to_string())
    f.write("\n\n")
    f.write("Status: Data processed, pathway kinetics computed.\n")
    f.write("Missing: Treatment response labels (need to extract from publication)\n")

print(f"\n💾 Summary saved to: {summary_path}")
