#!/usr/bin/env python3
"""
Search Project Data Sphere for PARP inhibitor + IO combination trials
Focus: Breast, Ovarian, or "Multiple" caslibs
Goal: Find trials for Holistic Score validation
"""

import sys
import os
from pathlib import Path
import json
from datetime import datetime

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from data_acquisition.utils.project_data_sphere_client import ProjectDataSphereClient

# Configuration
CAS_URL = "https://mpmprodvdmml.ondemand.sas.com/cas-shared-default-http/"
SSL_CERT = Path("data/certs/trustedcerts.pem")
USERNAME = "mpm0fxk2"

# Target caslibs - focus on breast, multiple, and others
TARGET_CASLIB_PREFIXES = ["Breast", "Multiple", "LungNo", "HeadNe"]

# Keywords for PARP/IO trials
TRIAL_KEYWORDS = {
    "parp": ["niraparib", "olaparib", "rucaparib", "talazoparib", "parp", "zejula", "lynparza", "rubraca"],
    "io": ["pembrolizumab", "nivolumab", "atezolizumab", "durvalumab", "keytruda", "opdivo", "tecentriq", "imfinzi", "pd-1", "pd-l1", "checkpoint"],
    "ddr": ["brca", "hrd", "homologous", "recombination", "dna repair", "platinum"],
    "trial_names": ["topacio", "keynote", "paola", "prima", "nova", "ariel"]
}

# Output directory
OUTPUT_DIR = Path("data/cohorts/pds_parp_io_search")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def search_caslib_for_keywords(client: ProjectDataSphereClient, caslib_name: str) -> dict:
    """Search a caslib for PARP/IO trial keywords"""
    print(f"\n📂 Searching {caslib_name}...")
    
    result = {
        "caslib_name": caslib_name,
        "files_found": [],
        "matched_keywords": set(),
        "clinical_files": []
    }
    
    try:
        files = client.list_files_in_caslib(caslib_name)
        
        for file_info in files:
            file_name_lower = file_info.get("name", "").lower()
            file_path_lower = file_info.get("path", "").lower()
            full_text = f"{file_name_lower} {file_path_lower}"
            
            # Check all keyword categories
            matched = False
            for category, keywords in TRIAL_KEYWORDS.items():
                for keyword in keywords:
                    if keyword in full_text:
                        result["files_found"].append({
                            "file": file_info,
                            "keyword": keyword,
                            "category": category
                        })
                        result["matched_keywords"].add(keyword)
                        matched = True
                        print(f"  ✅ {file_info.get('name', 'Unknown')} - matched '{keyword}' ({category})")
                        break
                if matched:
                    break
            
            # Also track general clinical files
            if any(kw in full_text for kw in ["clinical", "patient", "treatment", "outcome", "response"]):
                result["clinical_files"].append(file_info)
        
        print(f"  Found {len(result['files_found'])} PARP/IO files, {len(result['clinical_files'])} clinical files")
        
    except Exception as e:
        print(f"  ⚠️ Error: {e}")
        result["error"] = str(e)
    
    result["matched_keywords"] = list(result["matched_keywords"])
    return result


def main():
    print("=" * 80)
    print("Project Data Sphere: PARP/IO Trial Search")
    print("=" * 80)
    print("\nGoal: Find trials with PARP inhibitors + IO for Holistic Score validation")
    print()
    
    # Check SSL cert
    ssl_cert_path = str(SSL_CERT) if SSL_CERT.exists() else None
    if ssl_cert_path:
        print(f"✅ SSL cert: {SSL_CERT}")
    else:
        print(f"⚠️ SSL cert not found: {SSL_CERT}")
    
    # Initialize client
    try:
        client = ProjectDataSphereClient(cas_url=CAS_URL, ssl_cert_path=ssl_cert_path)
        print("✅ Client initialized")
    except Exception as e:
        print(f"❌ Error: {e}")
        return
    
    # Get password
    password = os.getenv("PDS_PASSWORD")
    if not password:
        print("\n🔐 Enter PDS password (same as SAS password):")
        import getpass
        password = getpass.getpass("Password: ")
    else:
        print("\n✅ Using password from PDS_PASSWORD environment variable")

    
    # Connect
    print("\n🔌 Connecting...")
    if not client.connect(username=USERNAME, password=password):
        print("❌ Connection failed")
        return
    
    # Get all caslibs
    print("\n📋 Fetching caslib list...")
    all_caslibs = client.list_caslibs()
    print(f"Total caslibs available: {len(all_caslibs)}")
    
    # Filter to target caslibs
    target_caslibs = [
        c for c in all_caslibs 
        if any(c.get("name", "").startswith(prefix) for prefix in TARGET_CASLIB_PREFIXES)
    ]
    
    print(f"\n🎯 Targeting {len(target_caslibs)} caslibs:")
    for c in target_caslibs:
        print(f"  - {c.get('name', 'Unknown')}")
    
    # Search each caslib
    print("\n" + "=" * 80)
    print("Searching Caslibs")
    print("=" * 80)
    
    all_results = {
        "timestamp": datetime.now().isoformat(),
        "caslibs_searched": len(target_caslibs),
        "results": [],
        "summary": {
            "total_parp_io_files": 0,
            "total_clinical_files": 0,
            "matched_keywords": set()
        }
    }
    
    for caslib in target_caslibs:
        caslib_name = caslib.get("name", "")
        result = search_caslib_for_keywords(client, caslib_name)
        all_results["results"].append(result)
        all_results["summary"]["total_parp_io_files"] += len(result["files_found"])
        all_results["summary"]["total_clinical_files"] += len(result["clinical_files"])
        all_results["summary"]["matched_keywords"].update(result["matched_keywords"])
    
    # Convert set to list for JSON serialization
    all_results["summary"]["matched_keywords"] = list(all_results["summary"]["matched_keywords"])
    
    # Save results
    output_file = OUTPUT_DIR / f"pds_parp_io_search_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    
    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"  Caslibs searched: {len(target_caslibs)}")
    print(f"  PARP/IO files found: {all_results['summary']['total_parp_io_files']}")
    print(f"  Clinical files found: {all_results['summary']['total_clinical_files']}")
    print(f"  Matched keywords: {', '.join(all_results['summary']['matched_keywords']) or 'None'}")
    print(f"\n  Results saved: {output_file}")
    
    # Show top candidates
    print("\n📊 Top Candidate Caslibs:")
    sorted_results = sorted(
        all_results["results"], 
        key=lambda x: len(x["files_found"]), 
        reverse=True
    )
    for i, result in enumerate(sorted_results[:5], 1):
        if result["files_found"]:
            print(f"  {i}. {result['caslib_name']}: {len(result['files_found'])} PARP/IO files")
            print(f"     Keywords: {', '.join(result['matched_keywords'][:5])}")
    
    # Disconnect
    client.disconnect()
    print("\n✅ Disconnected")
    
    return all_results


if __name__ == "__main__":
    try:
        results = main()
        sys.exit(0)
    except KeyboardInterrupt:
        print("\n⚠️ Interrupted")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
