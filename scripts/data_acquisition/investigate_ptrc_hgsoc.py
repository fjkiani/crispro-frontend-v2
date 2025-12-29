#!/usr/bin/env python3
"""
Investigate PTRC-HGSOC Collection

PTRC-HGSOC (Proteogenomic Tumor Research Consortium - High-Grade Serous Ovarian Cancer)
Collection on The Cancer Imaging Archive (TCIA)

Reported to have:
- Proteogenomic data from HGSOC tumor biopsies
- Platinum sensitivity annotations (refractory or sensitive)
- Need to verify: WES/WGS variant data, patient count
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

TCIA_BASE = "https://www.cancerimagingarchive.net"
PTRC_HGSOC_URL = "https://www.cancerimagingarchive.net/collection/ptrc-hgsoc/"


def investigate_ptrc_hgsoc() -> Dict:
    """Investigate PTRC-HGSOC collection."""
    print("="*80)
    print("PTRC-HGSOC COLLECTION INVESTIGATION")
    print("="*80)
    print()
    
    result = {
        "collection_name": "PTRC-HGSOC",
        "source": "TCIA (The Cancer Imaging Archive)",
        "url": PTRC_HGSOC_URL,
        "investigation_date": datetime.now().isoformat(),
        "status": "needs_manual_investigation",
        "findings": {}
    }
    
    print("📋 Collection Information:")
    print(f"   Name: PTRC-HGSOC (Proteogenomic Tumor Research Consortium - HGSOC)")
    print(f"   Source: The Cancer Imaging Archive (TCIA)")
    print(f"   URL: {PTRC_HGSOC_URL}")
    print()
    
    print("🔍 Reported Features:")
    print("   ✅ Proteogenomic data (proteomics + genomics)")
    print("   ✅ Platinum sensitivity annotations (refractory/sensitive)")
    print("   ✅ High-grade serous ovarian cancer (HGSOC)")
    print()
    
    print("❓ Critical Verification Needed:")
    print("   1. Does it include WES/WGS variant data (VCF/MAF)?")
    print("   2. Patient count (need ≥50)")
    print("   3. Variant coordinates (chrom, pos, ref, alt)")
    print("   4. Platinum response label format")
    print()
    
    print("📥 Access Information:")
    print("   TCIA requires registration for data download")
    print("   Registration: https://www.cancerimagingarchive.net/registration/")
    print("   Data Use Agreement required")
    print()
    
    # Try to access collection page
    try:
        print("📥 Attempting to access collection page...")
        response = requests.get(PTRC_HGSOC_URL, timeout=30)
        
        if response.status_code == 200:
            result["findings"]["page_accessible"] = True
            result["findings"]["page_content_length"] = len(response.text)
            
            # Check for keywords in page content
            content_lower = response.text.lower()
            keywords = {
                "whole_exome": "whole exome" in content_lower or "wes" in content_lower,
                "whole_genome": "whole genome" in content_lower or "wgs" in content_lower,
                "variant": "variant" in content_lower or "mutation" in content_lower,
                "vcf": "vcf" in content_lower,
                "maf": "maf" in content_lower,
                "platinum": "platinum" in content_lower,
                "sensitive": "sensitive" in content_lower,
                "refractory": "refractory" in content_lower
            }
            
            result["findings"]["keywords_found"] = keywords
            result["findings"]["has_variant_keywords"] = keywords.get("whole_exome") or keywords.get("whole_genome") or keywords.get("variant") or keywords.get("vcf") or keywords.get("maf")
            
            print(f"   ✅ Page accessible (HTTP {response.status_code})")
            print(f"   Keywords found:")
            for key, found in keywords.items():
                status = "✅" if found else "❌"
                print(f"      {status} {key}: {found}")
        else:
            result["findings"]["page_accessible"] = False
            result["findings"]["http_status"] = response.status_code
            print(f"   ⚠️  Page returned HTTP {response.status_code}")
    
    except Exception as e:
        result["findings"]["error"] = str(e)
        print(f"   ❌ Error accessing page: {e}")
    
    result["findings"]["recommendation"] = "Manual investigation required - check TCIA portal for dataset contents"
    result["findings"]["action_items"] = [
        "Register for TCIA access (if not already registered)",
        "Access PTRC-HGSOC collection page",
        "Check for WES/WGS variant files (VCF/MAF)",
        "Verify patient count",
        "Check platinum response label format",
        "Download if suitable"
    ]
    
    return result


def investigate_proteomexchange_pxd028225() -> Dict:
    """Investigate ProteomeXchange PXD028225 in detail."""
    print("\n" + "="*80)
    print("PROTEOMEXCHANGE PXD028225 DETAILED INVESTIGATION")
    print("="*80)
    print()
    
    px_id = "PXD028225"
    
    result = {
        "px_id": px_id,
        "source": "ProteomeXchange / PRIDE",
        "url": f"https://www.ebi.ac.uk/pride/archive/projects/{px_id}",
        "investigation_date": datetime.now().isoformat(),
        "status": "needs_verification"
    }
    
    print("📋 Dataset Information:")
    print(f"   ID: {px_id}")
    print(f"   Source: ProteomeXchange / PRIDE Archive")
    print(f"   URL: {result['url']}")
    print()
    
    print("🔍 Reported Features:")
    print("   ✅ Proteomic profiles")
    print("   ✅ Platinum-sensitive vs platinum-resistant HGSOC patients")
    print("   ✅ Has platinum response labels")
    print()
    
    print("❓ Critical Verification Needed:")
    print("   ⚠️  CRITICAL: Does it include WES/WGS variant data or only proteomics?")
    print("   ❓ Patient count (need ≥50)")
    print("   ❓ Variant coordinates if variants included")
    print()
    
    # Try to access PRIDE API
    try:
        print("📥 Attempting to access PRIDE API...")
        api_url = f"https://www.ebi.ac.uk/pride/archive/api/projects/{px_id}"
        response = requests.get(api_url, timeout=30)
        
        if response.status_code == 200:
            project_data = response.json()
            result["project_title"] = project_data.get("title", "")
            result["project_description"] = project_data.get("projectDescription", "")
            result["publication_date"] = project_data.get("publicationDate", "")
            result["submission_date"] = project_data.get("submissionDate", "")
            
            print(f"   ✅ API accessible")
            print(f"   Title: {result.get('project_title', 'N/A')}")
            print(f"   Description: {result.get('project_description', 'N/A')[:200]}...")
            
            # Check for variant-related keywords
            desc_lower = result.get("project_description", "").lower()
            result["has_variant_keywords"] = any(kw in desc_lower for kw in ["variant", "mutation", "wes", "wgs", "exome", "genome", "sequencing", "vcf", "maf"])
            
            if result["has_variant_keywords"]:
                print(f"   ✅ Description mentions variant/genomic keywords")
            else:
                print(f"   ⚠️  Description does not mention variant/genomic keywords")
                print(f"      May be proteomics-only dataset")
        else:
            result["api_status"] = f"error_{response.status_code}"
            print(f"   ⚠️  API returned {response.status_code}")
    
    except Exception as e:
        result["api_error"] = str(e)
        print(f"   ❌ Error accessing API: {e}")
    
    result["recommendation"] = "Manual verification required - check if WES/WGS variants included"
    result["action_items"] = [
        "Access PRIDE/ProteomeXchange portal",
        "Check dataset file list for VCF/MAF files",
        "Verify patient count",
        "Check if variants are in supplementary files",
        "Download if suitable"
    ]
    
    return result


def main():
    """Investigate promising leads."""
    print("="*80)
    print("PROMISING LEADS INVESTIGATION")
    print("="*80)
    print()
    
    all_results = {}
    
    # 1. PTRC-HGSOC
    ptrc_result = investigate_ptrc_hgsoc()
    all_results["ptrc_hgsoc"] = ptrc_result
    
    # 2. ProteomeXchange
    px_result = investigate_proteomexchange_pxd028225()
    all_results["proteomexchange_pxd028225"] = px_result
    
    # Save results
    output_file = OUTPUT_DIR / "promising_leads_investigation.json"
    with open(output_file, "w") as f:
        json.dump(all_results, f, indent=2)
    
    print("\n" + "="*80)
    print("INVESTIGATION SUMMARY")
    print("="*80)
    
    print(f"\n✅ PTRC-HGSOC:")
    print(f"   Status: {ptrc_result['status']}")
    if ptrc_result.get("findings", {}).get("has_variant_keywords"):
        print(f"   ✅ Page mentions variant/genomic keywords")
    else:
        print(f"   ⚠️  Need to verify variant data availability")
    
    print(f"\n✅ ProteomeXchange PXD028225:")
    print(f"   Status: {px_result['status']}")
    if px_result.get("has_variant_keywords"):
        print(f"   ✅ Description mentions variant/genomic keywords")
    else:
        print(f"   ⚠️  May be proteomics-only (needs verification)")
    
    print(f"\n📁 Results saved to: {output_file}")
    print()
    print("Next Steps:")
    print("1. Manually access PTRC-HGSOC collection on TCIA")
    print("2. Check for WES/WGS variant files (VCF/MAF)")
    print("3. Verify ProteomeXchange PXD028225 file contents")
    print("4. Extract suitable datasets")
    print()


if __name__ == "__main__":
    main()

