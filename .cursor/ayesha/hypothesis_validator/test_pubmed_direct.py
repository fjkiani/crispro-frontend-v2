#!/usr/bin/env python3
"""
Direct PubMed Test - Diagnose the issue
"""
import asyncio
import sys
from pathlib import Path

# Add backend to path
backend_path = Path(__file__).parent.parent.parent.parent / "oncology-coPilot" / "oncology-backend-minimal"
sys.path.insert(0, str(backend_path))

from api.services.enhanced_evidence_service import EnhancedEvidenceService

async def test_pubmed():
    print("🔍 Testing PubMed API Integration\n")
    
    service = EnhancedEvidenceService()
    
    # Test 1: Direct search_pubmed
    print("TEST 1: search_pubmed()")
    query = "curcumin ovarian cancer"
    print(f"Query: {query}")
    
    papers = await service.search_pubmed(query, max_results=5)
    print(f"Papers returned: {len(papers)}")
    
    if papers:
        print(f"✅ SUCCESS - Sample paper:")
        print(f"  PMID: {papers[0].get('pmid', 'N/A')}")
        print(f"  Title: {papers[0].get('title', 'N/A')[:100]}")
    else:
        print(f"❌ FAILED - No papers returned")
    
    print("\n" + "="*80 + "\n")
    
    # Test 2: get_complete_evidence
    print("TEST 2: get_complete_evidence()")
    print(f"Compound: Curcumin")
    print(f"Disease: ovarian_cancer_hgs")
    
    evidence = await service.get_complete_evidence(
        compound="Curcumin",
        disease="ovarian_cancer_hgs"
    )
    
    print(f"Total papers: {evidence.get('total_papers', 0)}")
    print(f"Evidence grade: {evidence.get('evidence_grade', 'N/A')}")
    print(f"Query used: {evidence.get('query_used', 'N/A')}")
    print(f"Mechanisms: {evidence.get('mechanisms', [])}")
    
    if evidence.get('total_papers', 0) > 0:
        print("✅ SUCCESS - Evidence mining working")
    else:
        print("❌ FAILED - No papers in evidence")

if __name__ == "__main__":
    asyncio.run(test_pubmed())


