
import pandas as pd
import numpy as np
from pathlib import Path

# Paths to Raw Files in Cloned Repos
BASE_PATH = Path("oncology-coPilot/oncology-backend-minimal/data/britoc/source_repos/britroc-1-cn-analysis")
PATIENT_DATA = BASE_PATH / "britroc_cohort_patient_data.tsv"
TREATMENT_DATA = BASE_PATH / "britroc_patient_treatment_data.tsv"
RESISTANCE_LIST = BASE_PATH / "britroc_primary_platinum_resistant_patient_list.tsv"
OUTPUT_FILE = "BRITROC_MASTER.csv"

def load_and_merge():
    print("Loading BriTROC-1 Raw Data...")
    
    # 1. Patient Cohort Data (Demographics + pre-calc PFS)
    df_patient = pd.read_csv(PATIENT_DATA, sep='\t')
    print(f"Loaded Patient Data: {len(df_patient)} rows")
    
    # 2. Treatment Data (Detailed lines)
    df_treatment = pd.read_csv(TREATMENT_DATA, sep='\t')
    print(f"Loaded Treatment Data: {len(df_treatment)} rows")
    
    # 3. Resistance Labels (The 11 primary resistant cases)
    df_res_list = pd.read_csv(RESISTANCE_LIST, sep='\t')
    resistant_ids = df_res_list['PATIENT_ID'].unique() # e.g. BRITROC-1
    # Note: britroc_cohort_patient_data uses integer IDs in 'britroc_number' (e.g. 1)
    
    # standardize IDs
    # df_patient['britroc_number'] is int (1, 2...)
    # df_res_list['PATIENT_ID'] is string ("BRITROC-1")
    # Need to extract int from "BRITROC-X"
    resistant_nums = [int(p.split('-')[1]) for p in resistant_ids]
    
    # 4. Merge Logic
    # Start with patient level
    master_df = df_patient.copy()
    
    # Add Binary Resistance Label
    master_df['is_primary_resistant'] = master_df['britroc_number'].isin(resistant_nums).astype(int)
    
    # Add Treatment Summary (Optional: Count lines)
    # Group treatment by patient (fk_britroc_number)
    tx_summary = df_treatment.groupby('fk_britroc_number').agg({
        'chemotherapy_line': 'max',
        'line_start': 'min', # First treatment date
        'line_end': 'max'    # Last treatment date
    }).rename(columns={'chemotherapy_line': 'total_tx_lines'})
    
    master_df = master_df.merge(tx_summary, left_on='britroc_number', right_index=True, how='left')
    
    # 5. Check for Signatures (Placeholder)
    # Signatures are in .rds. For now, we note they are missing from this CSV build.
    master_df['sig7_exposure'] = np.nan # To be filled via R script or audit data
    
    print("-" * 40)
    print("Master Dataset Summary")
    print("-" * 40)
    print(f"Total Patients: {len(master_df)}")
    print(f"Primary Resistant: {master_df['is_primary_resistant'].sum()}")
    print(f"PFS Available: {master_df['pfs'].notnull().sum()}")
    
    # Save
    master_df.to_csv(OUTPUT_FILE, index=False)
    print(f"\n✅ Master CSV built: {OUTPUT_FILE}")
    
    # Comparison check with resistance list
    # Verify the 11 patients are marked
    print("\nVerifying Resistance Labels:")
    res_in_master = master_df[master_df['is_primary_resistant'] == 1]['britroc_number'].tolist()
    print(f"Found {len(res_in_master)} resistant patients in master: {sorted(res_in_master)}")
    print(f"Expected: {sorted(resistant_nums)}")
    
    if set(res_in_master) == set(resistant_nums):
        print("✅ Labels match perfectly.")
    else:
        print("❌ Label mismatch!")

if __name__ == "__main__":
    load_and_merge()
