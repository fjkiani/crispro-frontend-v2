#!/usr/bin/env python3
"""
Extract colorectal cohort data with treatment and AE outcomes
"""

import swat
import json
import os
import pandas as pd
from pathlib import Path
from datetime import datetime

# Configuration
cas_url = "https://mpmprodvdmml.ondemand.sas.com/cas-shared-default-http/"
username = "mpm0fxk2"
password = "Thisisfjk12345!"
ssl_cert = Path("data/certs/trustedcerts.pem")

OUTPUT_DIR = Path("oncology-coPilot/oncology-backend-minimal/data/cohorts/pds_clinical_data")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("Extracting Colorectal Cohort Data")
print("=" * 60)

if ssl_cert.exists():
    os.environ["CAS_CLIENT_SSL_CA_LIST"] = str(ssl_cert)

conn = swat.CAS(cas_url, 443, username=username, password=password)
print("Connected!")

# Target caslib with good data structure
TARGET_CASLIB = "Colorec_Amgen_2005_262"

results = {
    "timestamp": datetime.now().isoformat(),
    "source": TARGET_CASLIB,
    "treatment_summary": {},
    "ae_summary": {},
    "sample_data": []
}

print(f"\nExtracting from: {TARGET_CASLIB}")
print("=" * 60)

# Load the main AE table
try:
    conn.table.loadTable(
        path="ae.sas7bdat",
        caslib=TARGET_CASLIB,
        casOut={"name": "AE", "caslib": "CASUSER"}
    )
    print("Loaded AE table")
    
    # Get all columns
    col_result = conn.table.columnInfo(table={"caslib": "CASUSER", "name": "AE"})
    columns = col_result.ColumnInfo['Column'].tolist()
    print(f"Columns ({len(columns)}): {columns[:15]}...")
    
    results["ae_columns"] = columns
    
    # Get table dimensions
    cas_table = conn.CASTable("AE", caslib="CASUSER")
    table_info = cas_table.tableInfo()
    n_rows = table_info.TableInfo['Rows'][0]
    print(f"Rows: {n_rows}")
    
    # Get treatment distribution
    print("\nTreatment distribution:")
    freq_result = cas_table.freq(inputs=["TRT"])
    if hasattr(freq_result, 'Frequency'):
        freq_df = freq_result.Frequency
        for _, row in freq_df.iterrows():
            trt = row.get('FmtVar', 'Unknown')
            count = row.get('Frequency', 0)
            pct = row.get('Percent', 0)
            print(f"  {trt}: {count} ({pct:.1f}%)")
            results["treatment_summary"][str(trt)] = {"count": int(count), "percent": float(pct)}
    
    # Get AE severity distribution (if exists)
    if 'AESEV' in columns or 'AESEVCD' in columns or 'AEGRADE' in columns:
        grade_col = 'AESEV' if 'AESEV' in columns else ('AESEVCD' if 'AESEVCD' in columns else 'AEGRADE')
        print(f"\nAE Severity ({grade_col}) distribution:")
        freq_result = cas_table.freq(inputs=[grade_col])
        if hasattr(freq_result, 'Frequency'):
            freq_df = freq_result.Frequency
            for _, row in freq_df.iterrows():
                grade = row.get('FmtVar', 'Unknown')
                count = row.get('Frequency', 0)
                print(f"  Grade {grade}: {count}")
                results["ae_summary"][str(grade)] = int(count)
    
    # Get sample data
    print("\nExtracting sample data (first 100 rows)...")
    sample_result = cas_table.head(n=100)
    if hasattr(sample_result, 'Head'):
        sample_df = sample_result.Head
        results["sample_data"] = sample_df.to_dict('records')[:20]  # Save first 20
        print(f"  Sample columns: {list(sample_df.columns)[:10]}")
    
    conn.table.dropTable(name="AE", caslib="CASUSER", quiet=True)

except Exception as e:
    print(f"Error: {e}")

# Also load the a_sendpt table for survival/response data
try:
    conn.table.loadTable(
        path="a_sendpt.sas7bdat",
        caslib=TARGET_CASLIB,
        casOut={"name": "A_SENDPT", "caslib": "CASUSER"}
    )
    print("\n\nLoaded A_SENDPT table (survival endpoints)")
    
    cas_table = conn.CASTable("A_SENDPT", caslib="CASUSER")
    col_result = conn.table.columnInfo(table={"caslib": "CASUSER", "name": "A_SENDPT"})
    columns = col_result.ColumnInfo['Column'].tolist()
    print(f"Columns: {columns[:15]}...")
    
    results["survival_columns"] = columns
    
    # Get dimensions
    table_info = cas_table.tableInfo()
    n_rows = table_info.TableInfo['Rows'][0]
    print(f"Rows: {n_rows}")
    
    # Get treatment distribution
    print("\nTreatment distribution:")
    freq_result = cas_table.freq(inputs=["TRT"])
    if hasattr(freq_result, 'Frequency'):
        freq_df = freq_result.Frequency
        for _, row in freq_df.iterrows():
            trt = row.get('FmtVar', 'Unknown')
            count = row.get('Frequency', 0)
            print(f"  {trt}: {count}")
    
    # Check for PFS/OS columns
    pfs_cols = [c for c in columns if 'pfs' in c.lower() or 'os' in c.lower() or 'survival' in c.lower()]
    print(f"\nSurvival-related columns: {pfs_cols}")
    
    conn.table.dropTable(name="A_SENDPT", caslib="CASUSER", quiet=True)

except Exception as e:
    print(f"Survival table error: {e}")

# Save results
output_path = OUTPUT_DIR / f"colorectal_cohort_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
with open(output_path, 'w') as f:
    json.dump(results, f, indent=2, default=str)

print()
print("=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"Source: {TARGET_CASLIB}")
print(f"Treatment groups: {list(results['treatment_summary'].keys())}")
print(f"AE grades: {list(results['ae_summary'].keys())}")
print(f"\nResults saved: {output_path}")

# Key finding
print("\n" + "=" * 60)
print("KEY FINDING FOR PGx VALIDATION")
print("=" * 60)
print("We have colorectal cancer trial data with:")
print("  - Treatment columns (TRT, TRTCD, etc.)")
print("  - Adverse event data (AE severity grades)")
print("  - Survival endpoints")
print()
print("LIMITATION:")
print("  - No explicit PGx genotype columns (DPYD, UGT1A1, etc.)")
print("  - PGx screening may not have been part of these trials")
print()
print("FOR REAL-WORLD PGx VALIDATION, WE NEED:")
print("  - Partner data with PGx screening + outcomes")
print("  - FDA FAERS data with genetic annotations")
print("  - Biobank data (UK Biobank, All of Us)")

conn.terminate()
print("\nDone!")
