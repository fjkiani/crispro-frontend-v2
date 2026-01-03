#!/usr/bin/env python3
"""
Examine gastric caslib SAS data files for PGx-relevant columns
"""

import swat
import json
from pathlib import Path
from datetime import datetime

# Configuration
cas_url = "https://mpmprodvdmml.ondemand.sas.com/cas-shared-default-http/"
username = "mpm0fxk2"
password = "Thisisfjk12345!"
ssl_cert = Path("data/certs/trustedcerts.pem")

# Data files to examine
DATA_FILES = [
    {"caslib": "Gastric_Multipl_1999_416", "file": "linked_gastric143.sas7bdat"},
    {"caslib": "Gastric_Multipl_2008_415", "file": "linked_gastric130.sas7bdat"}
]

# Output directory
OUTPUT_DIR = Path("oncology-coPilot/oncology-backend-minimal/data/cohorts")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("Examining Gastric Data Files for PGx-Relevant Columns")
print("=" * 60)
print()

# Set SSL cert
if ssl_cert.exists():
    import os
    os.environ["CAS_CLIENT_SSL_CA_LIST"] = str(ssl_cert)
    print(f"✅ SSL cert set: {ssl_cert}")

# Connect
print(f"🔌 Connecting to CAS server...")
try:
    conn = swat.CAS(cas_url, 443, username=username, password=password)
    print("✅ Connected!")
except Exception as e:
    print(f"❌ Connection failed: {e}")
    exit(1)

results = {
    "timestamp": datetime.now().isoformat(),
    "files_examined": [],
    "summary": {
        "total_files": 0,
        "files_with_treatment": 0,
        "files_with_outcome": 0,
        "files_with_pgx": 0
    }
}

# Keywords to search for in column names
TREATMENT_KEYWORDS = [
    "drug", "treatment", "therapy", "medication", "dose", "dosing",
    "fluorouracil", "5-fu", "capecitabine", "5fu", "folfox", "folfiri"
]

OUTCOME_KEYWORDS = [
    "ae", "adverse", "toxicity", "event", "response", "outcome",
    "grade", "severity", "death", "survival", "os", "pfs", "dfs"
]

PGX_KEYWORDS = [
    "dpyd", "tpmt", "ugt1a1", "pgx", "pharmacogen", "genotype",
    "variant", "mutation", "metabolizer", "phenotype"
]

def examine_file(caslib_name, file_name):
    """Examine a single SAS data file"""
    print(f"\n📂 Examining: {caslib_name}/{file_name}")
    print("-" * 60)
    
    file_result = {
        "caslib": caslib_name,
        "file": file_name,
        "columns": [],
        "n_rows": 0,
        "n_columns": 0,
        "treatment_columns": [],
        "outcome_columns": [],
        "pgx_columns": [],
        "error": None
    }
    
    try:
        # Load table - CAS creates table automatically
        table_name = file_name.replace(".sas7bdat", "").upper().replace(".", "_")
        
        # Load the file
        cas_table = conn.CASTable(table_name, caslib='CASUSER')
        conn.table.loadTable(
            sourceCaslib=caslib_name,
            casOut=cas_table,
            path=file_name
        )
        
        print(f"  ✅ Table loaded: {table_name}")
        
        # Get column info - reference the loaded table
        cas_table_ref = conn.CASTable(table_name, caslib='CASUSER')
        result = conn.table.columnInfo(table=cas_table_ref)
        if hasattr(result, 'ColumnInfo'):
            col_df = result.ColumnInfo
            columns = col_df.to_dict('records')
            
            file_result["n_columns"] = len(columns)
            file_result["columns"] = [col.get('Column', '') for col in columns]
            
            print(f"  Found {len(columns)} columns")
            
            # Search for relevant columns
            for col_info in columns:
                col_name = str(col_info.get('Column', '')).lower()
                
                # Check treatment
                if any(kw in col_name for kw in TREATMENT_KEYWORDS):
                    file_result["treatment_columns"].append(col_info.get('Column', ''))
                
                # Check outcome
                if any(kw in col_name for kw in OUTCOME_KEYWORDS):
                    file_result["outcome_columns"].append(col_info.get('Column', ''))
                
                # Check PGx
                if any(kw in col_name for kw in PGX_KEYWORDS):
                    file_result["pgx_columns"].append(col_info.get('Column', ''))
            
            # Get row count
            try:
                result = conn.table.tableInfo(table=cas_table)
                if hasattr(result, 'TableInfo'):
                    info_df = result.TableInfo
                    if len(info_df) > 0:
                        file_result["n_rows"] = int(info_df.iloc[0].get('Rows', 0))
            except:
                pass
            
            print(f"  Rows: {file_result['n_rows']}")
            print(f"  Treatment columns: {len(file_result['treatment_columns'])}")
            if file_result['treatment_columns']:
                print(f"    - {', '.join(file_result['treatment_columns'][:5])}")
            print(f"  Outcome columns: {len(file_result['outcome_columns'])}")
            if file_result['outcome_columns']:
                print(f"    - {', '.join(file_result['outcome_columns'][:5])}")
            print(f"  PGx columns: {len(file_result['pgx_columns'])}")
            if file_result['pgx_columns']:
                print(f"    - {', '.join(file_result['pgx_columns'])}")
            
        else:
            file_result["error"] = "No column info returned"
            print(f"  ⚠️  Could not get column info")
            
    except Exception as e:
        file_result["error"] = str(e)
        print(f"  ⚠️  Error: {e}")
    
    return file_result

# Examine each file
for data_file in DATA_FILES:
    result = examine_file(data_file["caslib"], data_file["file"])
    results["files_examined"].append(result)
    
    results["summary"]["total_files"] += 1
    if result["treatment_columns"]:
        results["summary"]["files_with_treatment"] += 1
    if result["outcome_columns"]:
        results["summary"]["files_with_outcome"] += 1
    if result["pgx_columns"]:
        results["summary"]["files_with_pgx"] += 1

# Save results
output_path = OUTPUT_DIR / f"pds_gastric_column_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
with open(output_path, 'w') as f:
    json.dump(results, f, indent=2)

print()
print("=" * 60)
print("Summary")
print("=" * 60)
print(f"  Files examined: {results['summary']['total_files']}")
print(f"  Files with treatment columns: {results['summary']['files_with_treatment']}")
print(f"  Files with outcome columns: {results['summary']['files_with_outcome']}")
print(f"  Files with PGx columns: {results['summary']['files_with_pgx']}")
print(f"  Results saved: {output_path}")

# Disconnect
try:
    conn.terminate()
    print("\n✅ Disconnected")
except:
    pass

print("\n✅ Analysis complete!")

