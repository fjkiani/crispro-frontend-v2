#!/usr/bin/env python3
"""Check Excel files for PFI data"""

import pandas as pd
from pathlib import Path

data_dir = Path('scripts/data_acquisition/sae/sciadv.abm1831_data_s1_and_s2 (1)')
print('Checking Excel files for PFI data...\n')

GSE165897_PATIENTS = [
    'EOC1005', 'EOC136', 'EOC153', 'EOC227', 'EOC349', 
    'EOC372', 'EOC3', 'EOC443', 'EOC540', 'EOC733', 'EOC87'
]

for excel_file in sorted(data_dir.glob('*.xlsx')):
    print(f'📄 {excel_file.name}')
    try:
        xls = pd.ExcelFile(excel_file)
        print(f'   Sheets: {xls.sheet_names}\n')
        
        for sheet in xls.sheet_names:
            df = pd.read_excel(excel_file, sheet_name=sheet, nrows=50)
            print(f'   Sheet "{sheet}": {df.shape[0]} rows, {df.shape[1]} cols')
            print(f'      Columns: {list(df.columns[:8])}')
            
            # Check for patient IDs
            cols_lower = [str(c).lower() for c in df.columns]
            has_patient = any('patient' in c or 'eoc' in c or 'id' in c for c in cols_lower)
            has_pfi = any('pfi' in c or 'platinum' in c or 'free' in c for c in cols_lower)
            
            if has_patient or has_pfi:
                print(f'      ✅ Contains patient/PFI data!')
                
                # Show first few rows if relevant
                if has_patient:
                    print(f'      First rows:')
                    print(df.head(3).to_string(max_cols=5))
                    print()
                    
                    # Check if any GSE165897 patient IDs are present
                    for col in df.columns:
                        if any(pid.lower() in str(col).lower() for pid in GSE165897_PATIENTS):
                            print(f'      ✅ Column "{col}" contains patient IDs!')
                        elif df[col].dtype == 'object':
                            sample_vals = df[col].dropna().astype(str).head(10).tolist()
                            matching = [v for v in sample_vals if any(pid in v for pid in GSE165897_PATIENTS)]
                            if matching:
                                print(f'      ✅ Found patient IDs in column "{col}": {matching[:3]}')
        print()
    except Exception as e:
        print(f'   ⚠️  Error: {e}\n')
