#!/usr/bin/env python3
"""Inspect supplementary data files for resistance labels."""

import pandas as pd
from pathlib import Path

data_dir = Path("scripts/data_acquisition/sae/sciadv.abm1831_data_s1_and_s2 (1)")

print("="*80)
print("INSPECTING SUPPLEMENTARY DATA FILES")
print("="*80)

# Check Data S1
f1 = data_dir / "sciadv.abm1831_data_s1.xlsx"
if f1.exists():
    print(f"\n📄 {f1.name}")
    xl1 = pd.ExcelFile(f1)
    print(f"   Sheets: {xl1.sheet_names}")
    for sheet in xl1.sheet_names:
        try:
            df = pd.read_excel(f1, sheet_name=sheet, nrows=20)
            print(f"\n   Sheet '{sheet}': {df.shape}")
            print(f"   Columns: {df.columns.tolist()[:10]}")
            print(f"   First 3 rows:")
            print(df.head(3).to_string(max_cols=5))
        except Exception as e:
            print(f"   Error reading sheet '{sheet}': {str(e)[:100]}")

# Check Data S2
f2 = data_dir / "sciadv.abm1831_data_s2.xlsx"
if f2.exists():
    print(f"\n📄 {f2.name}")
    try:
        xl2 = pd.ExcelFile(f2)
        print(f"   Sheets: {xl2.sheet_names}")
        for sheet in xl2.sheet_names[:5]:  # Check first 5 sheets
            try:
                df = pd.read_excel(f2, sheet_name=sheet, header=None, nrows=20)
                print(f"\n   Sheet '{sheet}': {df.shape}")
                print(f"   First 10 rows (raw):")
                print(df.head(10).to_string(max_cols=5))
            except Exception as e:
                print(f"   Error reading sheet '{sheet}': {str(e)[:100]}")
    except Exception as e:
        print(f"   Error opening file: {str(e)[:200]}")

print("\n" + "="*80)
