#!/usr/bin/env python3
"""
Download and extract clinical data files from Project Data Sphere gastric caslibs
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

# Output directory
OUTPUT_DIR = Path("oncology-coPilot/oncology-backend-minimal/data/cohorts/pds_clinical_data")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Target caslibs with clinical data
GASTRIC_CASLIBS = [
    "Gastric_SanofiU_1999_143",
    "Gastric_MerckKG_2008_130",
    "Gastric_Multipl_1999_416",
    "Gastric_Multipl_2008_415"
]

# Keywords to search for
PGX_KEYWORDS = ["dpyd", "tpmt", "ugt1a1", "pgx", "pharmacogen", "genotype", "variant", "metabolizer"]
TREATMENT_KEYWORDS = ["trt", "drug", "dose", "treatment", "therapy", "regimen", "arm", "5fu", "fluorouracil", "capecitabine"]
OUTCOME_KEYWORDS = ["ae", "adverse", "toxicity", "event", "grade", "response", "death", "survival", "os", "pfs"]

print("=" * 60)
print("Searching Clinical Data in Project Data Sphere")
print("=" * 60)
print()

# Set SSL cert
if ssl_cert.exists():
    os.environ["CAS_CLIENT_SSL_CA_LIST"] = str(ssl_cert)
    print(f"SSL cert set")

# Connect
print(f"Connecting to CAS server...")
try:
    conn = swat.CAS(cas_url, 443, username=username, password=password)
    print("Connected!")
except Exception as e:
    print(f"Connection failed: {e}")
    exit(1)

results = {
    "timestamp": datetime.now().isoformat(),
    "tables_analyzed": [],
    "treatment_columns_found": [],
    "outcome_columns_found": [],
    "pgx_columns_found": []
}

def analyze_columns(columns, table_name, caslib):
    """Analyze columns for PGx, treatment, and outcome relevance"""
    findings = {"treatment": [], "outcome": [], "pgx": []}
    
    for col in columns:
        col_lower = col.lower()
        for kw in PGX_KEYWORDS:
            if kw in col_lower:
                findings["pgx"].append(col)
                break
        for kw in TREATMENT_KEYWORDS:
            if kw in col_lower:
                findings["treatment"].append(col)
                break
        for kw in OUTCOME_KEYWORDS:
            if kw in col_lower:
                findings["outcome"].append(col)
                break
    
    return findings

# Explore each caslib
for caslib in GASTRIC_CASLIBS:
    print(f"\n{'='*60}")
    print(f"CASLIB: {caslib}")
    print(f"{'='*60}")
    
    try:
        # Get file listing
        file_result = conn.table.fileInfo(caslib=caslib, allFiles=True)
        if hasattr(file_result, 'FileInfo'):
            files = file_result.FileInfo.to_dict('records')
            sas_files = [f for f in files if str(f.get('Name', '')).endswith('.sas7bdat')]
            
            print(f"Found {len(sas_files)} SAS data files")
            
            for sas_file in sas_files:
                file_name = sas_file.get('Name', '')
                table_name = file_name.replace('.sas7bdat', '')
                
                print(f"\n  Loading: {table_name}...")
                
                try:
                    # Load the table
                    conn.table.loadTable(
                        path=file_name,
                        caslib=caslib,
                        casOut={"name": table_name, "caslib": "CASUSER"}
                    )
                    
                    # Get column info
                    col_result = conn.table.columnInfo(table={"caslib": "CASUSER", "name": table_name})
                    
                    if hasattr(col_result, 'ColumnInfo'):
                        col_df = col_result.ColumnInfo
                        columns = col_df['Column'].tolist() if 'Column' in col_df.columns else []
                        
                        findings = analyze_columns(columns, table_name, caslib)
                        
                        results["tables_analyzed"].append({
                            "caslib": caslib,
                            "table": table_name,
                            "total_columns": len(columns),
                            "columns": columns[:30]  # First 30 for reference
                        })
                        
                        if findings["treatment"]:
                            print(f"    TREATMENT: {findings['treatment']}")
                            results["treatment_columns_found"].extend([
                                {"caslib": caslib, "table": table_name, "column": c} 
                                for c in findings["treatment"]
                            ])
                        
                        if findings["outcome"]:
                            print(f"    OUTCOME: {findings['outcome']}")
                            results["outcome_columns_found"].extend([
                                {"caslib": caslib, "table": table_name, "column": c}
                                for c in findings["outcome"]
                            ])
                        
                        if findings["pgx"]:
                            print(f"    PGx: {findings['pgx']}")
                            results["pgx_columns_found"].extend([
                                {"caslib": caslib, "table": table_name, "column": c}
                                for c in findings["pgx"]
                            ])
                        
                        if not any([findings["treatment"], findings["outcome"], findings["pgx"]]):
                            print(f"    Columns: {columns[:10]}...")
                    
                    # Drop the table to free memory
                    conn.table.dropTable(name=table_name, caslib="CASUSER", quiet=True)
                
                except Exception as e:
                    print(f"    Error: {str(e)[:50]}")
    
    except Exception as e:
        print(f"  Error: {e}")

# Save results
output_path = OUTPUT_DIR / f"clinical_data_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
with open(output_path, 'w') as f:
    json.dump(results, f, indent=2)

print()
print("=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"Tables analyzed: {len(results['tables_analyzed'])}")
print(f"Treatment columns found: {len(results['treatment_columns_found'])}")
print(f"Outcome columns found: {len(results['outcome_columns_found'])}")
print(f"PGx columns found: {len(results['pgx_columns_found'])}")
print(f"\nResults saved: {output_path}")

if results["pgx_columns_found"]:
    print("\n*** FOUND PGx DATA! ***")
    for item in results["pgx_columns_found"]:
        print(f"   {item['caslib']}/{item['table']}: {item['column']}")

# Disconnect
try:
    conn.terminate()
except:
    pass

print("\nDone!")
