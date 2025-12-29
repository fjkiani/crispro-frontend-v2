#!/usr/bin/env python3
"""
Investigate External Ovarian Cancer Data Sources

Searches for publicly available non-TCGA ovarian cancer cohorts with:
- Platinum response labels (sensitive/resistant/refractory)
- Pre-treatment somatic variants with coordinates
- N ≥ 50 patients

Sources to investigate:
1. cBioPortal public studies
2. GEO/SRA datasets
3. ProteomeXchange
4. TCIA
"""

import json
import requests
import sys
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime

# Add project root to path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root))

OUTPUT_DIR = project_root / "data" / "external" / "ov_platinum_non_tcga" / "raw"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def search_cbioportal_ovarian_studies() -> List[Dict]:
    """
    Search cBioPortal for ovarian cancer studies.
    
    Note: cBioPortal API requires authentication for some endpoints.
    This is a placeholder for manual investigation.
    """
    print("🔍 Searching cBioPortal for ovarian cancer studies...")
    
    # cBioPortal API endpoint
    # Note: May need API key for full access
    api_base = "https://www.cbioportal.org/api"
    
    results = []
    
    try:
        # Search for ovarian cancer studies
        response = requests.get(
            f"{api_base}/studies",
            params={"keyword": "ovarian"},
            timeout=30
        )
        
        if response.status_code == 200:
            studies = response.json()
            print(f"  Found {len(studies)} studies matching 'ovarian'")
            
            for study in studies[:20]:  # Limit to first 20
                study_id = study.get("studyId", "")
                name = study.get("name", "")
                
                # Filter out TCGA studies
                if "tcga" in study_id.lower() or "tcga" in name.lower():
                    continue
                
                results.append({
                    "study_id": study_id,
                    "name": name,
                    "source": "cbioportal",
                    "url": f"https://www.cbioportal.org/study/summary?id={study_id}"
                })
        else:
            print(f"  ⚠️  cBioPortal API returned {response.status_code}")
            print("  Note: May require manual investigation via web portal")
            
    except Exception as e:
        print(f"  ❌ Error accessing cBioPortal API: {e}")
        print("  Note: Manual investigation recommended")
    
    return results


def check_proteomexchange_pxd028225() -> Dict:
    """
    Check ProteomeXchange dataset PXD028225.
    
    Returns metadata about the dataset.
    """
    print("🔍 Checking ProteomeXchange PXD028225...")
    
    # ProteomeXchange API
    px_id = "PXD028225"
    
    result = {
        "px_id": px_id,
        "source": "proteomexchange",
        "url": f"https://www.ebi.ac.uk/pride/archive/projects/{px_id}",
        "status": "needs_verification",
        "notes": []
    }
    
    try:
        # Try to access PRIDE API
        response = requests.get(
            f"https://www.ebi.ac.uk/pride/archive/projects/{px_id}",
            timeout=30
        )
        
        if response.status_code == 200:
            result["status"] = "accessible"
            result["notes"].append("Dataset exists, need to verify contents")
        else:
            result["status"] = "access_error"
            result["notes"].append(f"HTTP {response.status_code}")
            
    except Exception as e:
        result["status"] = "error"
        result["notes"].append(f"Error: {e}")
    
    result["notes"].append("Critical: Verify if WES/WGS variants included or only proteomics")
    result["notes"].append("Critical: Verify if platinum response labels present")
    
    return result


def search_geo_ovarian_platinum() -> List[Dict]:
    """
    Search GEO for ovarian cancer + platinum + sequencing datasets.
    
    Note: This is a placeholder - actual GEO search requires NCBI E-utilities API.
    """
    print("🔍 Searching GEO for ovarian cancer + platinum datasets...")
    
    results = []
    
    # Example search terms
    search_terms = [
        "ovarian cancer platinum chemotherapy whole exome sequencing",
        "ovarian cancer platinum response WES",
        "HGSOC platinum sensitive resistant sequencing"
    ]
    
    print("  Note: GEO search requires NCBI E-utilities API")
    print("  Recommended: Manual search at https://www.ncbi.nlm.nih.gov/geo/")
    print(f"  Search terms: {', '.join(search_terms)}")
    
    # Placeholder for actual API call
    results.append({
        "source": "geo",
        "search_terms": search_terms,
        "url": "https://www.ncbi.nlm.nih.gov/geo/",
        "status": "manual_search_required"
    })
    
    return results


def document_icgc_argo_access() -> Dict:
    """
    Document ICGC-ARGO access requirements.
    """
    print("📋 Documenting ICGC-ARGO access requirements...")
    
    return {
        "source": "icgc_argo",
        "portal": "https://platform.icgc-argo.org",
        "daco": "https://daco.icgc-argo.org",
        "access_type": "controlled",
        "requirements": [
            "DACO (Data Access Compliance Office) application",
            "Institutional approval",
            "Data use agreement"
        ],
        "timeline": "2-4 weeks for approval",
        "status": "requires_application",
        "notes": [
            "Best shot for real clinical annotations",
            "May have platinum response metadata",
            "Start application if no public alternatives found"
        ]
    }


def main():
    """Run investigation of all sources."""
    print("="*80)
    print("EXTERNAL OVARIAN CANCER DATA SOURCE INVESTIGATION")
    print("="*80)
    print()
    
    all_results = {
        "investigation_date": datetime.now().isoformat(),
        "sources": {}
    }
    
    # 1. cBioPortal
    print("\n1. cBioPortal Investigation")
    print("-" * 80)
    cbioportal_results = search_cbioportal_ovarian_studies()
    all_results["sources"]["cbioportal"] = {
        "results": cbioportal_results,
        "count": len(cbioportal_results)
    }
    
    # 2. ProteomeXchange
    print("\n2. ProteomeXchange Investigation")
    print("-" * 80)
    px_result = check_proteomexchange_pxd028225()
    all_results["sources"]["proteomexchange"] = px_result
    
    # 3. GEO
    print("\n3. GEO Investigation")
    print("-" * 80)
    geo_results = search_geo_ovarian_platinum()
    all_results["sources"]["geo"] = geo_results
    
    # 4. ICGC-ARGO
    print("\n4. ICGC-ARGO Documentation")
    print("-" * 80)
    icgc_info = document_icgc_argo_access()
    all_results["sources"]["icgc_argo"] = icgc_info
    
    # Save results
    output_file = OUTPUT_DIR / "source_investigation_results.json"
    with open(output_file, "w") as f:
        json.dump(all_results, f, indent=2)
    
    print(f"\n✅ Results saved to: {output_file}")
    
    # Summary
    print("\n" + "="*80)
    print("INVESTIGATION SUMMARY")
    print("="*80)
    print(f"cBioPortal: {len(cbioportal_results)} non-TCGA studies found")
    print(f"ProteomeXchange: {px_result['status']}")
    print(f"GEO: Manual search required")
    print(f"ICGC-ARGO: Controlled access (requires application)")
    print()
    print("Next Steps:")
    print("1. Manually verify cBioPortal studies for platinum response labels")
    print("2. Check ProteomeXchange PXD028225 for variant data")
    print("3. Search GEO manually for ovarian + platinum datasets")
    print("4. If no public sources found, start ICGC-ARGO DACO application")
    print()


if __name__ == "__main__":
    main()

