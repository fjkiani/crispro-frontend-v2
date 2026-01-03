#!/usr/bin/env python3
"""
Connect to Project Data Sphere and explore gastric caslibs for DPYD/fluoropyrimidine data
"""

import sys
import os
from pathlib import Path
import json
from datetime import datetime

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "oncology-coPilot" / "oncology-backend-minimal"))

try:
    from scripts.data_acquisition.pds.project_data_sphere_client import ProjectDataSphereClient
except ImportError as e:
    print(f"❌ Could not import ProjectDataSphereClient: {e}")
    print("   Make sure you're in the correct directory")
    sys.exit(1)

# Gastric caslibs to explore
GASTRIC_CASLIBS = [
    "Gastric_MerckKG_2008_130",
    "Gastric_Multipl_1999_416",
    "Gastric_Multipl_2008_415",
    "Gastric_SanofiU_1999_143"
]

# Configuration
CAS_URL = "https://mpmprodvdmml.ondemand.sas.com/cas-shared-default-http/"
SSL_CERT = Path("data/certs/trustedcerts.pem")
USERNAME = "mpm0fxk2"

# Output directory
OUTPUT_DIR = Path("oncology-coPilot/oncology-backend-minimal/data/cohorts")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
RECEIPTS_DIR = OUTPUT_DIR / "receipts"
RECEIPTS_DIR.mkdir(exist_ok=True)


def explore_caslib(client: ProjectDataSphereClient, caslib_name: str) -> dict:
    """Explore a single caslib"""
    print(f"\n📂 Exploring caslib: {caslib_name}")
    print("-" * 60)
    
    result = {
        "caslib_name": caslib_name,
        "files": [],
        "file_count": 0,
        "clinical_data_files": [],
        "error": None
    }
    
    try:
        # List files
        files = client.list_files_in_caslib(caslib_name)
        result["file_count"] = len(files)
        result["files"] = files
        
        print(f"  Found {len(files)} files")
        
        # Look for clinical data files
        clinical_keywords = [
            "clinical", "patient", "treatment", "outcome", 
            "ae", "adverse", "toxicity", "response",
            "csv", "sas7bdat", "xlsx"
        ]
        
        for file_info in files:
            file_name = file_info.get("name", "").lower()
            file_path = file_info.get("path", "").lower()
            full_text = f"{file_name} {file_path}"
            
            if any(keyword in full_text for keyword in clinical_keywords):
                result["clinical_data_files"].append(file_info)
                print(f"    📋 {file_info.get('name', 'Unknown')}")
        
        print(f"  Clinical data files: {len(result['clinical_data_files'])}")
        
    except Exception as e:
        result["error"] = str(e)
        print(f"  ⚠️  Error: {e}")
    
    return result


def main():
    print("=" * 60)
    print("Project Data Sphere: Gastric Caslib Exploration")
    print("=" * 60)
    print()
    
    # Check SSL cert
    ssl_cert_path = None
    if SSL_CERT.exists():
        ssl_cert_path = str(SSL_CERT)
        print(f"✅ SSL cert found: {SSL_CERT}")
    else:
        print(f"⚠️  SSL cert not found: {SSL_CERT}")
        print("   Will try without explicit cert")
    
    # Initialize client
    try:
        client = ProjectDataSphereClient(
            cas_url=CAS_URL,
            ssl_cert_path=ssl_cert_path
        )
        print(f"✅ Client initialized")
    except Exception as e:
        print(f"❌ Error initializing client: {e}")
        return
    
    # Get password
    print()
    print("🔐 Connection requires password")
    password = os.getenv("PDS_PASSWORD")
    if not password:
        print("   Set PDS_PASSWORD environment variable or enter password:")
        import getpass
        password = getpass.getpass("   Password: ")
    
    # Connect
    print()
    print("🔌 Connecting to Project Data Sphere...")
    if not client.connect(username=USERNAME, password=password):
        print("❌ Connection failed")
        return
    
    # Explore each caslib
    print()
    print("=" * 60)
    print("Exploring Gastric Caslibs")
    print("=" * 60)
    
    all_results = {
        "timestamp": datetime.now().isoformat(),
        "caslibs_explored": [],
        "summary": {
            "total_files": 0,
            "total_clinical_files": 0
        }
    }
    
    for caslib in GASTRIC_CASLIBS:
        result = explore_caslib(client, caslib)
        all_results["caslibs_explored"].append(result)
        all_results["summary"]["total_files"] += result["file_count"]
        all_results["summary"]["total_clinical_files"] += len(result["clinical_data_files"])
    
    # Save results
    output_path = OUTPUT_DIR / f"pds_gastric_exploration_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(output_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    
    print()
    print("=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"  Total files found: {all_results['summary']['total_files']}")
    print(f"  Clinical data files: {all_results['summary']['total_clinical_files']}")
    print(f"  Results saved: {output_path}")
    
    # Save receipt
    receipt = {
        "timestamp": datetime.now().isoformat(),
        "output_file": str(output_path),
        "caslibs_explored": GASTRIC_CASLIBS,
        "method": "Project Data Sphere CAS API"
    }
    
    receipt_path = RECEIPTS_DIR / f"pds_gastric_exploration_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(receipt_path, 'w') as f:
        json.dump(receipt, f, indent=2)
    
    print(f"  Receipt saved: {receipt_path}")
    
    # Disconnect
    try:
        client.disconnect()
        print("\n✅ Disconnected from Project Data Sphere")
    except:
        pass
    
    return all_results


if __name__ == "__main__":
    try:
        results = main()
        sys.exit(0)
    except KeyboardInterrupt:
        print("\n\n⚠️  Interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

