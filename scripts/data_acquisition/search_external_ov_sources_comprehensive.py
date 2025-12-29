#!/usr/bin/env python3
"""
Comprehensive Search for External Ovarian Cancer Datasets

Searches multiple sources for non-TCGA ovarian cancer cohorts with:
- Platinum response labels (sensitive/resistant/refractory)
- N ≥ 50 patients
- Pre-treatment somatic variants with coordinates
"""

import json
import requests
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime

# Add project root to path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root))

OUTPUT_DIR = project_root / "data" / "external" / "ov_platinum_non_tcga" / "raw"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SEARCH_RESULTS_FILE = OUTPUT_DIR / "comprehensive_search_results.json"


def search_cbioportal_all_ovarian_studies() -> List[Dict]:
    """Search cBioPortal for ALL ovarian cancer studies (not just MSK)."""
    print("\n" + "="*80)
    print("1. cBioPortal - Comprehensive Ovarian Study Search")
    print("="*80)
    
    api_base = "https://www.cbioportal.org/api"
    results = []
    
    try:
        # Get all studies
        print("📥 Fetching all studies from cBioPortal...")
        response = requests.get(f"{api_base}/studies", timeout=60)
        
        if response.status_code == 200:
            all_studies = response.json()
            print(f"   ✅ Found {len(all_studies)} total studies")
            
            # Filter for ovarian cancer studies
            ovarian_keywords = ["ovarian", "ov", "hgsoc", "serous"]
            
            for study in all_studies:
                study_id = study.get("studyId", "").lower()
                name = study.get("name", "").lower()
                
                # Check if ovarian-related
                is_ovarian = any(keyword in study_id or keyword in name for keyword in ovarian_keywords)
                
                # Exclude TCGA
                is_tcga = "tcga" in study_id or "tcga" in name
                
                if is_ovarian and not is_tcga:
                    results.append({
                        "study_id": study.get("studyId"),
                        "name": study.get("name"),
                        "description": study.get("description", ""),
                        "source": "cbioportal",
                        "url": f"https://www.cbioportal.org/study/summary?id={study.get('studyId')}"
                    })
            
            print(f"   ✅ Found {len(results)} non-TCGA ovarian studies")
            
            # Try to get patient counts for each
            for result in results[:10]:  # Limit to first 10 to avoid rate limiting
                try:
                    study_id = result["study_id"]
                    time.sleep(1)  # Rate limiting
                    
                    patients_response = requests.get(
                        f"{api_base}/studies/{study_id}/patients",
                        params={"pageSize": 10000},
                        timeout=30
                    )
                    
                    if patients_response.status_code == 200:
                        patients = patients_response.json()
                        result["patient_count"] = len(patients)
                        result["meets_minimum"] = len(patients) >= 50
                    else:
                        result["patient_count"] = "unknown"
                        result["meets_minimum"] = "unknown"
                except Exception as e:
                    result["patient_count"] = f"error: {e}"
                    result["meets_minimum"] = "unknown"
        
        else:
            print(f"   ⚠️  API returned {response.status_code}")
    
    except Exception as e:
        print(f"   ❌ Error: {e}")
    
    return results


def search_geo_ovarian_platinum() -> List[Dict]:
    """Search GEO for ovarian cancer + platinum datasets."""
    print("\n" + "="*80)
    print("2. GEO (Gene Expression Omnibus) - Ovarian + Platinum Search")
    print("="*80)
    
    print("   Note: GEO search requires NCBI E-utilities API or manual search")
    print("   Recommended manual search URLs:")
    
    search_terms = [
        "ovarian cancer platinum chemotherapy",
        "ovarian cancer platinum response",
        "HGSOC platinum sensitive resistant",
        "ovarian cancer whole exome sequencing platinum"
    ]
    
    results = []
    
    for term in search_terms:
        encoded_term = term.replace(" ", "+")
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE&term={encoded_term}"
        
        results.append({
            "source": "geo",
            "search_term": term,
            "url": url,
            "status": "manual_search_required",
            "notes": "Check for datasets with: WES/WGS variants + platinum response labels"
        })
    
    print(f"   ✅ Prepared {len(results)} search queries")
    print("   Action: Manual search recommended at https://www.ncbi.nlm.nih.gov/geo/")
    
    return results


def search_sra_ovarian_platinum() -> List[Dict]:
    """Search SRA for ovarian cancer + platinum sequencing datasets."""
    print("\n" + "="*80)
    print("3. SRA (Sequence Read Archive) - Ovarian + Platinum Search")
    print("="*80)
    
    print("   Note: SRA search requires NCBI E-utilities API or manual search")
    
    search_terms = [
        "ovarian cancer platinum whole exome",
        "ovarian cancer platinum whole genome",
        "HGSOC platinum chemotherapy sequencing"
    ]
    
    results = []
    
    for term in search_terms:
        encoded_term = term.replace(" ", "+")
        url = f"https://www.ncbi.nlm.nih.gov/sra/?term={encoded_term}"
        
        results.append({
            "source": "sra",
            "search_term": term,
            "url": url,
            "status": "manual_search_required",
            "notes": "Check for: raw sequencing data + associated clinical metadata with platinum response"
        })
    
    print(f"   ✅ Prepared {len(results)} search queries")
    print("   Action: Manual search recommended at https://www.ncbi.nlm.nih.gov/sra/")
    
    return results


def investigate_proteomexchange_pxd028225() -> Dict:
    """Investigate ProteomeXchange dataset PXD028225."""
    print("\n" + "="*80)
    print("4. ProteomeXchange PXD028225 Investigation")
    print("="*80)
    
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
        print(f"📥 Checking ProteomeXchange dataset {px_id}...")
        
        # PRIDE API endpoint
        response = requests.get(
            f"https://www.ebi.ac.uk/pride/archive/projects/{px_id}",
            timeout=30
        )
        
        if response.status_code == 200:
            result["status"] = "accessible"
            result["notes"].append("Dataset exists and is accessible")
            
            # Try to get project details via PRIDE API
            try:
                api_response = requests.get(
                    f"https://www.ebi.ac.uk/pride/archive/api/projects/{px_id}",
                    timeout=30
                )
                if api_response.status_code == 200:
                    project_data = api_response.json()
                    result["project_title"] = project_data.get("title", "")
                    result["project_description"] = project_data.get("projectDescription", "")
                    result["notes"].append(f"Title: {result.get('project_title', 'N/A')}")
            except:
                pass
        else:
            result["status"] = f"access_error_{response.status_code}"
            result["notes"].append(f"HTTP {response.status_code}")
    
    except Exception as e:
        result["status"] = "error"
        result["notes"].append(f"Error: {e}")
    
    result["notes"].append("CRITICAL: Verify if WES/WGS variants included or only proteomics")
    result["notes"].append("CRITICAL: Verify if platinum response labels present")
    result["notes"].append("CRITICAL: Check patient count (need ≥50)")
    
    print(f"   Status: {result['status']}")
    if "project_title" in result:
        print(f"   Title: {result['project_title']}")
    
    return result


def search_pubmed_ovarian_platinum_studies() -> List[Dict]:
    """Search PubMed for studies that might have published datasets."""
    print("\n" + "="*80)
    print("5. PubMed - Ovarian + Platinum + Sequencing Studies")
    print("="*80)
    
    print("   Searching for published studies that might have deposited data...")
    
    # PubMed E-utilities search
    search_queries = [
        "ovarian cancer[Title] AND platinum[Title] AND (whole exome[Title/Abstract] OR WES[Title/Abstract] OR sequencing[Title/Abstract])",
        "HGSOC[Title] AND platinum response[Title/Abstract] AND sequencing[Title/Abstract]",
        "ovarian cancer[Title] AND platinum sensitive[Title] AND genomic[Title/Abstract]"
    ]
    
    results = []
    
    for query in search_queries:
        encoded_query = query.replace(" ", "+").replace("[", "%5B").replace("]", "%5D")
        url = f"https://pubmed.ncbi.nlm.nih.gov/?term={encoded_query}"
        
        results.append({
            "source": "pubmed",
            "search_query": query,
            "url": url,
            "status": "manual_search_required",
            "notes": "Check publications for: data availability statements, GEO/SRA accession numbers, supplementary data"
        })
    
    print(f"   ✅ Prepared {len(results)} PubMed search queries")
    print("   Action: Manual search recommended - look for data availability statements")
    
    return results


def document_icgc_argo_application() -> Dict:
    """Document ICGC-ARGO application process."""
    print("\n" + "="*80)
    print("6. ICGC-ARGO - Controlled Access Documentation")
    print("="*80)
    
    return {
        "source": "icgc_argo",
        "portal": "https://platform.icgc-argo.org",
        "daco": "https://daco.icgc-argo.org",
        "access_type": "controlled",
        "status": "requires_application",
        "requirements": [
            "DACO (Data Access Compliance Office) application",
            "Institutional approval",
            "Data use agreement",
            "Research proposal"
        ],
        "timeline": "2-4 weeks for approval (typical)",
        "best_for": "High-quality clinical annotations + platinum response metadata",
        "notes": [
            "Best shot for real clinical annotations",
            "May have platinum response metadata",
            "Start application if no public alternatives found",
            "Check available ovarian programs before applying"
        ],
        "action_items": [
            "Identify candidate ovarian projects/programs",
            "Confirm available fields include platinum response proxy",
            "Start DACO application",
            "Document requirements + expected timeline"
        ]
    }


def main():
    """Run comprehensive search across all sources."""
    print("="*80)
    print("COMPREHENSIVE EXTERNAL OVARIAN CANCER DATA SEARCH")
    print("="*80)
    print(f"Date: {datetime.now().isoformat()}")
    print()
    
    all_results = {
        "search_date": datetime.now().isoformat(),
        "sources": {}
    }
    
    # 1. cBioPortal comprehensive search
    cbioportal_results = search_cbioportal_all_ovarian_studies()
    all_results["sources"]["cbioportal"] = {
        "results": cbioportal_results,
        "count": len(cbioportal_results),
        "meets_minimum": [r for r in cbioportal_results if r.get("meets_minimum") == True]
    }
    
    # 2. GEO search
    geo_results = search_geo_ovarian_platinum()
    all_results["sources"]["geo"] = geo_results
    
    # 3. SRA search
    sra_results = search_sra_ovarian_platinum()
    all_results["sources"]["sra"] = sra_results
    
    # 4. ProteomeXchange
    px_result = investigate_proteomexchange_pxd028225()
    all_results["sources"]["proteomexchange"] = px_result
    
    # 5. PubMed
    pubmed_results = search_pubmed_ovarian_platinum_studies()
    all_results["sources"]["pubmed"] = pubmed_results
    
    # 6. ICGC-ARGO
    icgc_info = document_icgc_argo_application()
    all_results["sources"]["icgc_argo"] = icgc_info
    
    # Save results
    with open(SEARCH_RESULTS_FILE, "w") as f:
        json.dump(all_results, f, indent=2)
    
    print("\n" + "="*80)
    print("SEARCH SUMMARY")
    print("="*80)
    print(f"✅ cBioPortal: {len(cbioportal_results)} non-TCGA ovarian studies found")
    print(f"   Studies meeting minimum (≥50 patients): {len([r for r in cbioportal_results if r.get('meets_minimum') == True])}")
    print(f"✅ GEO: {len(geo_results)} search queries prepared")
    print(f"✅ SRA: {len(sra_results)} search queries prepared")
    print(f"✅ ProteomeXchange: {px_result['status']}")
    print(f"✅ PubMed: {len(pubmed_results)} search queries prepared")
    print(f"✅ ICGC-ARGO: Documentation complete")
    print()
    print(f"📁 Results saved to: {SEARCH_RESULTS_FILE}")
    print()
    print("Next Steps:")
    print("1. Review cBioPortal studies for platinum response labels")
    print("2. Manually search GEO/SRA using provided queries")
    print("3. Verify ProteomeXchange PXD028225 contents")
    print("4. Check PubMed publications for data availability")
    print("5. If no public sources found, start ICGC-ARGO DACO application")
    print()


if __name__ == "__main__":
    main()

