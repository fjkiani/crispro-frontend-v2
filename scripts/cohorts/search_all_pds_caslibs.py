#!/usr/bin/env python3
"""
Search ALL Project Data Sphere caslibs for PGx-relevant data
"""

import swat
import json
import os
from pathlib import Path
from datetime import datetime

# Configuration
cas_url = "https://mpmprodvdmml.ondemand.sas.com/cas-shared-default-http/"
username = "mpm0fxk2"
password = "Thisisfjk12345!"
ssl_cert = Path("data/certs/trustedcerts.pem")

# Keywords
PGX_KEYWORDS = ["dpyd", "tpmt", "ugt1a1", "pgx", "pharmacogen", "genotype", "metabolizer", "cyp2d6"]
TREATMENT_KEYWORDS = ["trt", "drug", "dose", "treatment", "5fu", "fluorouracil", "capecitabine", "irinotecan", "mercaptopurine", "thiopurine"]

# Output directory
OUTPUT_DIR = Path("oncology-coPilot/oncology-backend-minimal/data/cohorts/pds_clinical_data")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("Searching ALL PDS Caslibs for PGx Data")
print("=" * 60)
print()

# Set SSL cert
if ssl_cert.exists():
    os.environ["CAS_CLIENT_SSL_CA_LIST"] = str(ssl_cert)

# Connect
print("Connecting...")
try:
    conn = swat.CAS(cas_url, 443, username=username, password=password)
    print("Connected!")
except Exception as e:
    print(f"Connection failed: {e}")
    exit(1)

# Get all caslibs
print("\nListing all caslibs...")
caslib_result = conn.table.caslibInfo()
caslibs = []
if hasattr(caslib_result, 'CASLibInfo'):
    caslib_df = caslib_result.CASLibInfo
    caslibs = caslib_df['Name'].tolist()
    print(f"Found {len(caslibs)} caslibs")

# Filter for potentially relevant caslibs (colorectal, gastric, breast for fluoropyrimidines)
relevant_keywords = ["colorectal", "colon", "gastric", "breast", "pancrea", "head", "neck"]
relevant_caslibs = []
for c in caslibs:
    c_lower = c.lower()
    for kw in relevant_keywords:
        if kw in c_lower:
            relevant_caslibs.append(c)
            break

print(f"\nFiltered to {len(relevant_caslibs)} potentially relevant caslibs:")
for c in relevant_caslibs[:20]:
    print(f"  - {c}")
if len(relevant_caslibs) > 20:
    print(f"  ... and {len(relevant_caslibs) - 20} more")

results = {
    "timestamp": datetime.now().isoformat(),
    "caslibs_searched": [],
    "sas_files_found": [],
    "tables_with_treatment": [],
    "tables_with_pgx": []
}

# Search each relevant caslib for SAS files
print("\n" + "=" * 60)
print("Searching caslibs for clinical data files...")
print("=" * 60)

for caslib in relevant_caslibs[:30]:  # Limit to first 30
    print(f"\n{caslib}:")
    
    try:
        file_result = conn.table.fileInfo(caslib=caslib, allFiles=True)
        if hasattr(file_result, 'FileInfo'):
            files = file_result.FileInfo.to_dict('records')
            sas_files = [f for f in files if str(f.get('Name', '')).endswith('.sas7bdat')]
            
            if sas_files:
                print(f"  Found {len(sas_files)} SAS files")
                for sf in sas_files[:5]:
                    print(f"    - {sf.get('Name')}")
                    results["sas_files_found"].append({
                        "caslib": caslib,
                        "file": sf.get('Name'),
                        "size": sf.get('Size')
                    })
                
                # Load and check first SAS file for columns
                first_file = sas_files[0].get('Name', '')
                table_name = first_file.replace('.sas7bdat', '')
                
                try:
                    conn.table.loadTable(
                        path=first_file,
                        caslib=caslib,
                        casOut={"name": table_name, "caslib": "CASUSER"}
                    )
                    
                    col_result = conn.table.columnInfo(table={"caslib": "CASUSER", "name": table_name})
                    if hasattr(col_result, 'ColumnInfo'):
                        columns = col_result.ColumnInfo['Column'].tolist()
                        
                        # Check for PGx columns
                        for col in columns:
                            col_lower = col.lower()
                            for kw in PGX_KEYWORDS:
                                if kw in col_lower:
                                    print(f"    *** PGx COLUMN: {col} ***")
                                    results["tables_with_pgx"].append({
                                        "caslib": caslib,
                                        "table": table_name,
                                        "column": col
                                    })
                        
                        # Check for treatment columns
                        for col in columns:
                            col_lower = col.lower()
                            for kw in TREATMENT_KEYWORDS:
                                if kw in col_lower:
                                    print(f"    Treatment: {col}")
                                    results["tables_with_treatment"].append({
                                        "caslib": caslib,
                                        "table": table_name,
                                        "column": col
                                    })
                    
                    conn.table.dropTable(name=table_name, caslib="CASUSER", quiet=True)
                except Exception as e:
                    print(f"    Load error: {str(e)[:40]}")
            else:
                print(f"  No SAS files")
        
        results["caslibs_searched"].append(caslib)
    
    except Exception as e:
        print(f"  Error: {str(e)[:40]}")

# Save results
output_path = OUTPUT_DIR / f"pds_full_search_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
with open(output_path, 'w') as f:
    json.dump(results, f, indent=2)

print()
print("=" * 60)
print("FINAL SUMMARY")
print("=" * 60)
print(f"Caslibs searched: {len(results['caslibs_searched'])}")
print(f"SAS files found: {len(results['sas_files_found'])}")
print(f"Tables with treatment columns: {len(results['tables_with_treatment'])}")
print(f"Tables with PGx columns: {len(results['tables_with_pgx'])}")

if results["tables_with_pgx"]:
    print("\n*** PGx DATA FOUND! ***")
    for item in results["tables_with_pgx"]:
        print(f"  {item['caslib']}/{item['table']}: {item['column']}")
else:
    print("\nNo explicit PGx columns found in searched caslibs.")
    print("Note: PGx data may be in:")
    print("  - ZIP files (need extraction)")
    print("  - Separate biomarker caslibs")
    print("  - Not available in this dataset")

print(f"\nResults saved: {output_path}")

conn.terminate()
print("\nDone!")
