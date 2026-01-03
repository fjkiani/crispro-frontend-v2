#!/usr/bin/env python3
"""
Direct exploration of gastric caslibs using existing connection pattern
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

# Gastric caslibs
GASTRIC_CASLIBS = [
    "Gastric_MerckKG_2008_130",
    "Gastric_Multipl_1999_416",
    "Gastric_Multipl_2008_415",
    "Gastric_SanofiU_1999_143"
]

# Output directory
OUTPUT_DIR = Path("oncology-coPilot/oncology-backend-minimal/data/cohorts")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
RECEIPTS_DIR = OUTPUT_DIR / "receipts"
RECEIPTS_DIR.mkdir(exist_ok=True)

print("=" * 60)
print("Project Data Sphere: Gastric Caslib Exploration")
print("=" * 60)
print()

# Set SSL cert
if ssl_cert.exists():
    import os
    os.environ["CAS_CLIENT_SSL_CA_LIST"] = str(ssl_cert)
    print(f"✅ SSL cert set: {ssl_cert}")
else:
    print(f"⚠️  SSL cert not found: {ssl_cert}")

# Connect
print(f"🔌 Connecting to CAS server...")
try:
    conn = swat.CAS(cas_url, 443, username=username, password=password)
    print("✅ Connected!")
except Exception as e:
    print(f"❌ Connection failed: {e}")
    exit(1)

# Explore each caslib
results = {
    "timestamp": datetime.now().isoformat(),
    "caslibs_explored": [],
    "summary": {
        "total_files": 0,
        "total_clinical_files": 0
    }
}

print()
print("=" * 60)
print("Exploring Gastric Caslibs")
print("=" * 60)

for caslib_name in GASTRIC_CASLIBS:
    print(f"\n📂 {caslib_name}")
    print("-" * 60)
    
    caslib_result = {
        "caslib_name": caslib_name,
        "files": [],
        "file_count": 0,
        "clinical_data_files": [],
        "error": None
    }
    
    try:
        # List files
        params = {'allFiles': True, 'caslib': caslib_name}
        result = conn.table.fileInfo(**params)
        
        files = []
        if hasattr(result, 'FileInfo'):
            file_df = result.FileInfo
            files = file_df.to_dict('records')
        
        caslib_result["file_count"] = len(files)
        caslib_result["files"] = files
        
        print(f"  Found {len(files)} files")
        
        # Look for clinical data files
        clinical_keywords = [
            "clinical", "patient", "treatment", "outcome", 
            "ae", "adverse", "toxicity", "response",
            "csv", "sas7bdat", "xlsx", "data"
        ]
        
        for file_info in files:
            file_name = str(file_info.get("Name", "")).lower()
            file_path = str(file_info.get("Path", "")).lower()
            full_text = f"{file_name} {file_path}"
            
            if any(keyword in full_text for keyword in clinical_keywords):
                caslib_result["clinical_data_files"].append(file_info)
                print(f"    📋 {file_info.get('Name', 'Unknown')}")
        
        print(f"  ✅ Clinical data files: {len(caslib_result['clinical_data_files'])}")
        
    except Exception as e:
        caslib_result["error"] = str(e)
        print(f"  ⚠️  Error: {e}")
    
    results["caslibs_explored"].append(caslib_result)
    results["summary"]["total_files"] += caslib_result["file_count"]
    results["summary"]["total_clinical_files"] += len(caslib_result["clinical_data_files"])

# Save results
output_path = OUTPUT_DIR / f"pds_gastric_exploration_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
with open(output_path, 'w') as f:
    json.dump(results, f, indent=2)

print()
print("=" * 60)
print("Summary")
print("=" * 60)
print(f"  Total files found: {results['summary']['total_files']}")
print(f"  Clinical data files: {results['summary']['total_clinical_files']}")
print(f"  Results saved: {output_path}")

# Save receipt
receipt = {
    "timestamp": datetime.now().isoformat(),
    "output_file": str(output_path),
    "caslibs_explored": GASTRIC_CASLIBS,
    "method": "Direct CAS API connection"
}

receipt_path = RECEIPTS_DIR / f"pds_gastric_exploration_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
with open(receipt_path, 'w') as f:
    json.dump(receipt, f, indent=2)

print(f"  Receipt saved: {receipt_path}")

# Disconnect
try:
    conn.terminate()
    print("\n✅ Disconnected")
except:
    pass

print("\n✅ Exploration complete!")

