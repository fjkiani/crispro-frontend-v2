#!/usr/bin/env python3
"""
Deep search for colorectal caslibs and biomarker data
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

OUTPUT_DIR = Path("oncology-coPilot/oncology-backend-minimal/data/cohorts/pds_clinical_data")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("Deep Search: Colorectal & Biomarker Caslibs")
print("=" * 60)

if ssl_cert.exists():
    os.environ["CAS_CLIENT_SSL_CA_LIST"] = str(ssl_cert)

conn = swat.CAS(cas_url, 443, username=username, password=password)
print("Connected!")

# Get all caslibs
caslib_result = conn.table.caslibInfo()
caslibs = caslib_result.CASLibInfo['Name'].tolist()

# Find colorectal and biomarker caslibs
colorectal = [c for c in caslibs if 'color' in c.lower() or 'colon' in c.lower() or 'crc' in c.lower()]
biomarker = [c for c in caslibs if 'biomarker' in c.lower() or 'genetic' in c.lower() or 'genom' in c.lower()]

print(f"\nColorectal caslibs: {len(colorectal)}")
for c in colorectal:
    print(f"  - {c}")

print(f"\nBiomarker/genetic caslibs: {len(biomarker)}")
for c in biomarker:
    print(f"  - {c}")

# Search colorectal caslibs for clinical data with treatment
results = {
    "timestamp": datetime.now().isoformat(),
    "colorectal_caslibs": colorectal,
    "tables_with_treatment": [],
    "tables_with_ae": [],
    "all_columns_found": {}
}

print("\n" + "=" * 60)
print("Examining colorectal caslibs in detail...")
print("=" * 60)

for caslib in colorectal[:15]:
    print(f"\n{caslib}:")
    
    try:
        file_result = conn.table.fileInfo(caslib=caslib, allFiles=True)
        if hasattr(file_result, 'FileInfo'):
            files = file_result.FileInfo.to_dict('records')
            sas_files = [f for f in files if str(f.get('Name', '')).endswith('.sas7bdat')]
            
            if sas_files:
                print(f"  SAS files: {len(sas_files)}")
                
                for sf in sas_files[:3]:  # Check first 3 files
                    file_name = sf.get('Name', '')
                    table_name = file_name.replace('.sas7bdat', '')
                    
                    try:
                        conn.table.loadTable(
                            path=file_name,
                            caslib=caslib,
                            casOut={"name": table_name, "caslib": "CASUSER"}
                        )
                        
                        col_result = conn.table.columnInfo(table={"caslib": "CASUSER", "name": table_name})
                        if hasattr(col_result, 'ColumnInfo'):
                            columns = col_result.ColumnInfo['Column'].tolist()
                            results["all_columns_found"][f"{caslib}/{table_name}"] = columns
                            
                            # Look for key columns
                            cols_lower = [c.lower() for c in columns]
                            has_trt = any('trt' in c or 'drug' in c or 'dose' in c for c in cols_lower)
                            has_ae = any('ae' in c or 'toxicity' in c or 'grade' in c for c in cols_lower)
                            has_dpyd = any('dpyd' in c or 'dpy' in c for c in cols_lower)
                            
                            if has_trt:
                                trt_cols = [c for c in columns if any(x in c.lower() for x in ['trt', 'drug', 'dose'])]
                                print(f"    {table_name}: TRT columns = {trt_cols[:5]}")
                                results["tables_with_treatment"].append({
                                    "caslib": caslib,
                                    "table": table_name,
                                    "columns": trt_cols
                                })
                            
                            if has_ae:
                                ae_cols = [c for c in columns if any(x in c.lower() for x in ['ae', 'toxicity', 'grade'])]
                                print(f"    {table_name}: AE columns = {ae_cols[:5]}")
                                results["tables_with_ae"].append({
                                    "caslib": caslib,
                                    "table": table_name,
                                    "columns": ae_cols
                                })
                            
                            if has_dpyd:
                                print(f"    *** DPYD FOUND! ***")
                        
                        conn.table.dropTable(name=table_name, caslib="CASUSER", quiet=True)
                    except Exception as e:
                        pass
            else:
                # Check for ZIP files
                zip_files = [f for f in files if str(f.get('Name', '')).endswith('.zip')]
                if zip_files:
                    print(f"  ZIP files: {len(zip_files)}")
                    for zf in zip_files:
                        print(f"    - {zf.get('Name')}")
    except Exception as e:
        print(f"  Error: {str(e)[:40]}")

# Save results
output_path = OUTPUT_DIR / f"colorectal_deep_search_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
with open(output_path, 'w') as f:
    json.dump(results, f, indent=2)

print()
print("=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"Tables with treatment: {len(results['tables_with_treatment'])}")
print(f"Tables with AE data: {len(results['tables_with_ae'])}")
print(f"\nResults saved: {output_path}")

# Check for FOLFOX/FOLFIRI treatments specifically
print("\nSearching for fluoropyrimidine treatments (FOLFOX, FOLFIRI, 5-FU, capecitabine)...")

for item in results["tables_with_treatment"]:
    caslib = item["caslib"]
    table = item["table"]
    
    try:
        # Load table
        file_name = table + ".sas7bdat"
        conn.table.loadTable(
            path=file_name,
            caslib=caslib,
            casOut={"name": table, "caslib": "CASUSER"}
        )
        
        # Get unique treatment values
        cas_table = conn.CASTable(table, caslib="CASUSER")
        
        for col in item["columns"][:2]:
            try:
                freq_result = cas_table.freq(inputs=[col])
                if hasattr(freq_result, 'Frequency'):
                    vals = freq_result.Frequency['FmtVar'].tolist()[:20]
                    vals_str = [str(v).lower() for v in vals]
                    
                    if any('folfox' in v or 'folfiri' in v or '5-fu' in v or 'fluorouracil' in v or 'capecitabine' in v for v in vals_str):
                        print(f"\n*** FLUOROPYRIMIDINE TREATMENT FOUND! ***")
                        print(f"  Caslib: {caslib}")
                        print(f"  Table: {table}")
                        print(f"  Column: {col}")
                        print(f"  Values: {vals[:10]}")
            except:
                pass
        
        conn.table.dropTable(name=table, caslib="CASUSER", quiet=True)
    except:
        pass

conn.terminate()
print("\nDone!")
