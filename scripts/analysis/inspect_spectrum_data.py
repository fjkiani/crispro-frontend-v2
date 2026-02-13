import pandas as pd
import os

files = [
    "oncology-coPilot/oncology-backend-minimal/data/external/spectrum/patient_metadata.xlsx",
    "oncology-coPilot/oncology-backend-minimal/data/external/spectrum/sample_metadata.xlsx"
]

for f in files:
    print(f"\n--- Inspecting {os.path.basename(f)} ---")
    try:
        df = pd.read_excel(f)
        print("Columns:", list(df.columns))
        print("Shape:", df.shape)
        print("Head (2 rows):")
        print(df.head(2).to_string())
        
        # Check for CA125 keywords
        ca125_cols = [c for c in df.columns if "125" in c.lower() or "marker" in c.lower()]
        if ca125_cols:
            print(">>> CA-125 Columns Found:", ca125_cols)
        else:
            print(">>> No explicit CA-125 columns in headers.")
            
    except Exception as e:
        print(f"Error reading {f}: {e}")
