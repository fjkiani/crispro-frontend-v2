#!/usr/bin/env python3
"""
Explore Project Data Sphere for PGx-relevant clinical trial data

Purpose: Find caslibs with data for:
- Fluoropyrimidines (colorectal, gastric cancer) - DPYD
- Thiopurines (ALL, leukemia) - TPMT  
- Irinotecan (colorectal, pancreatic) - UGT1A1
"""

import json
import sys
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any

# Load existing caslibs data
CASLIBS_FILE = Path("oncology-coPilot/oncology-backend-minimal/data/kelim_resurrection/s4_deliverables/project_data_sphere_caslibs.json")

# PGx-relevant cancer types
PGX_CANCER_TYPES = {
    "DPYD": ["colorectal", "colon", "rectal", "gastric", "stomach", "esoph"],
    "TPMT": ["leukemia", "all", "lymphoma", "childhood"],
    "UGT1A1": ["colorectal", "colon", "rectal", "pancreatic", "pancreas"],
}

# Output directory
OUTPUT_DIR = Path("oncology-coPilot/oncology-backend-minimal/data/cohorts")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
RECEIPTS_DIR = OUTPUT_DIR / "receipts"
RECEIPTS_DIR.mkdir(exist_ok=True)


def load_caslibs() -> List[Dict[str, Any]]:
    """Load caslibs from existing file"""
    if not CASLIBS_FILE.exists():
        print(f"❌ Caslibs file not found: {CASLIBS_FILE}")
        return []
    
    with open(CASLIBS_FILE, 'r') as f:
        caslibs = json.load(f)
    
    return caslibs


def search_pgx_relevant_caslibs(caslibs: List[Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
    """
    Search for caslibs relevant to PGx drugs
    """
    print("🔍 Searching for PGx-relevant caslibs...\n")
    
    results = {
        "DPYD": [],  # Fluoropyrimidines (colorectal, gastric)
        "TPMT": [],  # Thiopurines (ALL, leukemia)
        "UGT1A1": [],  # Irinotecan (colorectal, pancreatic)
        "multiple": []  # Caslibs matching multiple categories
    }
    
    for caslib in caslibs:
        name = caslib.get("Name", "").lower()
        description = caslib.get("Description", "").lower() if caslib.get("Description") else ""
        full_text = f"{name} {description}"
        
        matches = []
        
        # Check DPYD (colorectal, gastric)
        if any(term in full_text for term in PGX_CANCER_TYPES["DPYD"]):
            matches.append("DPYD")
            results["DPYD"].append(caslib)
        
        # Check TPMT (ALL, leukemia)
        if any(term in full_text for term in PGX_CANCER_TYPES["TPMT"]):
            matches.append("TPMT")
            results["TPMT"].append(caslib)
        
        # Check UGT1A1 (colorectal, pancreatic)
        if any(term in full_text for term in PGX_CANCER_TYPES["UGT1A1"]):
            matches.append("UGT1A1")
            results["UGT1A1"].append(caslib)
        
        # Multiple matches
        if len(matches) > 1:
            results["multiple"].append({
                "caslib": caslib,
                "matches": matches
            })
    
    return results


def print_results(results: Dict[str, List[Dict[str, Any]]]):
    """Print search results"""
    print("=" * 60)
    print("PGx-Relevant Caslibs in Project Data Sphere")
    print("=" * 60)
    print()
    
    for gene, caslibs in results.items():
        if gene == "multiple":
            continue
        
        print(f"📊 {gene} (Relevant Cancer Types):")
        print(f"   Found {len(caslibs)} caslibs")
        
        if caslibs:
            print("   Sample caslibs:")
            for caslib in caslibs[:5]:
                name = caslib.get("Name", "Unknown")
                path = caslib.get("Path", "Unknown")
                print(f"     - {name}")
                print(f"       Path: {path}")
            if len(caslibs) > 5:
                print(f"     ... and {len(caslibs) - 5} more")
        print()
    
    if results["multiple"]:
        print(f"🔄 Caslibs matching multiple categories: {len(results['multiple'])}")
        for item in results["multiple"][:3]:
            name = item["caslib"].get("Name", "Unknown")
            matches = ", ".join(item["matches"])
            print(f"   - {name} ({matches})")
        print()


def save_results(results: Dict[str, List[Dict[str, Any]]]):
    """Save results to JSON"""
    output = {
        "timestamp": datetime.now().isoformat(),
        "search_criteria": {
            "DPYD": "colorectal, gastric (fluoropyrimidines)",
            "TPMT": "ALL, leukemia (thiopurines)",
            "UGT1A1": "colorectal, pancreatic (irinotecan)"
        },
        "results": {
            "DPYD": results["DPYD"],
            "TPMT": results["TPMT"],
            "UGT1A1": results["UGT1A1"],
            "multiple": results["multiple"]
        },
        "summary": {
            "DPYD_count": len(results["DPYD"]),
            "TPMT_count": len(results["TPMT"]),
            "UGT1A1_count": len(results["UGT1A1"]),
            "multiple_count": len(results["multiple"])
        }
    }
    
    output_path = OUTPUT_DIR / f"pds_pgx_caslibs_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"💾 Results saved: {output_path}")
    
    # Save receipt
    receipt = {
        "timestamp": datetime.now().isoformat(),
        "output_file": str(output_path),
        "method": "Caslib name/description search for PGx-relevant cancer types"
    }
    
    receipt_path = RECEIPTS_DIR / f"pds_pgx_search_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(receipt_path, 'w') as f:
        json.dump(receipt, f, indent=2)
    
    return output_path


def main():
    """Main execution"""
    print("=" * 60)
    print("Project Data Sphere: PGx-Relevant Caslib Search")
    print("=" * 60)
    print()
    
    # Load caslibs
    caslibs = load_caslibs()
    if not caslibs:
        print("❌ No caslibs found. Exiting.")
        sys.exit(1)
    
    print(f"📂 Loaded {len(caslibs)} caslibs from Project Data Sphere\n")
    
    # Search for PGx-relevant caslibs
    results = search_pgx_relevant_caslibs(caslibs)
    
    # Print results
    print_results(results)
    
    # Save results
    output_path = save_results(results)
    
    # Summary
    total_relevant = len(results["DPYD"]) + len(results["TPMT"]) + len(results["UGT1A1"])
    print(f"\n✅ Summary:")
    print(f"   Total PGx-relevant caslibs: {total_relevant}")
    print(f"   DPYD (colorectal/gastric): {len(results['DPYD'])}")
    print(f"   TPMT (ALL/leukemia): {len(results['TPMT'])}")
    print(f"   UGT1A1 (colorectal/pancreatic): {len(results['UGT1A1'])}")
    
    if total_relevant > 0:
        print(f"\n🎯 Next Steps:")
        print(f"   1. Review caslibs in: {output_path}")
        print(f"   2. Connect to Project Data Sphere to explore files in these caslibs")
        print(f"   3. Look for clinical data tables with PGx, treatment, and outcome columns")
    else:
        print(f"\n⚠️  No PGx-relevant caslibs found by name search.")
        print(f"   May need to:")
        print(f"   1. Search 'Multiple' caslibs")
        print(f"   2. Explore clinical data tables directly")
        print(f"   3. Search for specific drug names in caslib contents")
    
    return results


if __name__ == "__main__":
    try:
        results = main()
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

