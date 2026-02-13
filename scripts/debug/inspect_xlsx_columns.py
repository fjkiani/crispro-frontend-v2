import pandas as pd
import os

XLSX_PATH = "oncology-coPilot/oncology-backend-minimal/data/MSK-Spectrum/Tracking-Clonal/SupplementaryTables.xlsx"

try:
    xl = pd.ExcelFile(XLSX_PATH)
    print(f"Sheets: {xl.sheet_names}")
    
    # S13
    s13 = next(s for s in xl.sheet_names if "S13" in s and "Raw" not in s)
    df13 = pd.read_excel(XLSX_PATH, sheet_name=s13, nrows=5)
    print(f"\n--- {s13} Columns ---")
    print(df13.columns.tolist())
    print(df13.head(2).to_string())

    # S13 Raw
    s13r = next((s for s in xl.sheet_names if "S13 Raw" in s), None)
    if s13r:
        df13r = pd.read_excel(XLSX_PATH, sheet_name=s13r, nrows=5)
        print(f"\n--- {s13r} Columns ---")
        print(df13r.columns.tolist())

    # S3
    s3 = next(s for s in xl.sheet_names if "S3" in s)
    df3 = pd.read_excel(XLSX_PATH, sheet_name=s3, nrows=5)
    print(f"\n--- {s3} Columns ---")
    print(df3.columns.tolist())
    print(df3.head(2).to_string())

except Exception as e:
    print(e)
