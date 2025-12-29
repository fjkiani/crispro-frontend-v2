# 🔬 FRAMEWORKS INTEGRATION - IMPLEMENTATION GUIDE

**Date**: January 15, 2025  
**Status**: ✅ **AUDITED, TESTED, READY TO BUILD**  
**Frameworks**: pubmearch ✅, pubmed_parser ✅

---

## 📦 FRAMEWORKS SUMMARY

### **Framework 1: pubmearch** ✅

**Location**: `.github/frameworks/pubmearch-main/`

**What It Does**:
- Advanced PubMed search with Bio.Entrez
- Keyword hotspot analysis
- Trend tracking over time
- Publication count analysis
- Batch retrieval (1000+ results)

**Key Classes**:
- `PubMedSearcher` - Search PubMed, export results
- `PubMedAnalyzer` - Analyze keywords, trends, counts

**Status**: ✅ **TESTED & WORKING**

---

### **Framework 2: pubmed_parser** ✅

**Location**: `.github/frameworks/pubmearch-main/pubmed_parser-master/`

**What It Does**:
- Parse PubMed Open Access XML (full-text articles)
- Parse MEDLINE XML (abstracts + metadata)
- Extract citations, references
- Extract tables, figures, captions
- Parse paragraphs with citation context

**Key Functions**:
- `parse_pubmed_xml()` - Parse PubMed OA XML
- `parse_medline_xml()` - Parse MEDLINE XML
- `parse_pubmed_references()` - Extract citations
- `parse_pubmed_paragraph()` - Extract paragraphs with citations
- `parse_pubmed_table()` - Extract tables

**Status**: ✅ **CORE FUNCTIONS WORK** (needs install for version check)

---

## 🎯 MOAT EXPLOITATION STRATEGY

### **1. Enhanced PubMed Search** (pubmearch)

**Current Problem**: Limited to 100 results, no keyword analysis  
**Solution**: Use pubmearch for advanced search + analysis

**MOAT Benefits**:
- **Keyword hotspot analysis** → Identify emerging mechanisms
- **Trend tracking** → Track mechanism evolution over time
- **Batch retrieval** → Get 1000+ papers for comprehensive analysis

**Integration**:
```python
# Replace current PubMed search
from pubmearch.pubmed_searcher import PubMedSearcher
from pubmearch.analyzer import PubMedAnalyzer

# Search with advanced syntax
searcher = PubMedSearcher(email=..., api_key=...)
articles = searcher.search(
    advanced_search="(purple potatoes OR anthocyanins) AND ovarian cancer",
    date_range=("2020/01/01", "2024/12/31"),
    max_results=1000
)

# Analyze keywords (identify mechanisms)
analyzer = PubMedAnalyzer()
keyword_analysis = analyzer.analyze_research_keywords(articles, top_n=20)

# Top keywords = potential mechanisms!
mechanisms = [kw['keyword'] for kw in keyword_analysis['top_keywords']]
```

---

### **2. Deep Full-Text Parsing** (pubmed_parser)

**Current Problem**: Only abstracts, missing full-text mechanisms  
**Solution**: Parse full-text PubMed OA articles

**MOAT Benefits**:
- **Full-text parsing** → Extract mechanisms from Methods/Results (not just abstracts)
- **Table extraction** → Get quantitative data (IC50, EC50, dosages)
- **Citation context** → Assess evidence strength from citation networks
- **Paragraph parsing** → Extract mechanism mentions with context

**Integration**:
```python
from pubmed_parser.pubmed_oa_parser import parse_pubmed_xml, parse_pubmed_paragraph

# Parse full-text article
full_text = parse_pubmed_xml(pmc_xml_path)

# Extract paragraphs with citations
paragraphs = parse_pubmed_paragraph(pmc_xml_path, all_paragraph=True)

# Filter for mechanism-relevant sections
methods_results = [
    p for p in paragraphs 
    if p['section'] in ['Methods', 'Results', 'Discussion']
]

# LLM extracts mechanisms from full context (not just abstract!)
mechanisms = await llm.extract_mechanisms(methods_results)
```

---

### **3. Research Intelligence Orchestrator**

**Combines**: pubmearch + pubmed_parser + LLM + MOAT

**Flow**:
```
Question: "How do purple potatoes help with ovarian cancer?"
    ↓
[1] LLM Question Formulator
    → Sub-questions, entities, portal queries
    ↓
[2] Multi-Portal Query (Parallel)
    → pubmearch: Search + keyword analysis
    → cBioPortal: Genomics data
    → ClinVar: Variant data
    ↓
[3] Deep Parsing (pubmed_parser)
    → Parse full-text for top papers
    → Extract mechanisms from Methods/Results
    → Extract tables for quantitative data
    ↓
[4] LLM Synthesis
    → Read all parsed content
    → Extract mechanisms, evidence, confidence
    ↓
[5] MOAT Analysis
    → Pathway alignment
    → Treatment line appropriateness
    → Biomarker matching
```

---

## 🔧 MODULAR INTEGRATION ARCHITECTURE

### **Directory Structure**

```
api/services/research_intelligence/
├── __init__.py
├── orchestrator.py                    # Main orchestrator
├── question_formulator.py              # LLM question decomposition
├── synthesis_engine.py                 # LLM synthesis
├── moat_integrator.py                 # MOAT analysis
│
├── portals/                            # Portal clients (modular)
│   ├── __init__.py
│   ├── pubmed_enhanced.py            # pubmearch wrapper ✅
│   ├── cbioportal_client.py          # pyBioPortal wrapper
│   ├── clinvar_client.py             # Existing
│   └── chembl_client.py              # Existing
│
└── parsers/                            # Deep parsers (modular)
    ├── __init__.py
    ├── pubmed_deep_parser.py          # pubmed_parser wrapper ✅
    └── medline_batch_parser.py        # Batch MEDLINE parsing
```

---

## 📝 IMPLEMENTATION FILES

### **File 1: `api/services/research_intelligence/portals/pubmed_enhanced.py`**

**Purpose**: Wrapper around pubmearch for enhanced PubMed search

```python
"""
Enhanced PubMed Portal using pubmearch framework.

Provides:
- Advanced PubMed search with batch retrieval
- Keyword hotspot analysis
- Trend tracking
- Publication count analysis
"""

import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import asyncio

# Add pubmearch to path
pubmearch_path = Path(__file__).resolve().parent.parent.parent.parent.parent / ".github" / "frameworks" / "pubmearch-main"
sys.path.insert(0, str(pubmearch_path))

from pubmearch.pubmed_searcher import PubMedSearcher
from pubmearch.analyzer import PubMedAnalyzer


class EnhancedPubMedPortal:
    """
    Enhanced PubMed portal using pubmearch framework.
    
    Advantages over current implementation:
    - Advanced search syntax support
    - Batch retrieval (handles 1000+ results)
    - Keyword hotspot analysis
    - Trend tracking
    - Publication count analysis
    """
    
    def __init__(self, email: Optional[str] = None, api_key: Optional[str] = None):
        """
        Initialize enhanced PubMed portal.
        
        Args:
            email: NCBI email (defaults to NCBI_USER_EMAIL env var)
            api_key: NCBI API key (defaults to NCBI_USER_API_KEY env var)
        """
        self.email = email or os.getenv('NCBI_USER_EMAIL')
        self.api_key = api_key or os.getenv('NCBI_USER_API_KEY')
        
        if not self.email:
            raise ValueError("NCBI_USER_EMAIL must be set")
        
        self.searcher = PubMedSearcher(email=self.email, api_key=self.api_key)
        self.analyzer = PubMedAnalyzer()
    
    async def search_with_analysis(
        self,
        query: str,
        date_range: Optional[Tuple[str, str]] = None,
        max_results: int = 1000,
        analyze_keywords: bool = True,
        include_trends: bool = True
    ) -> Dict[str, Any]:
        """
        Search PubMed and return results + keyword analysis.
        
        Args:
            query: PubMed search query (advanced syntax supported)
            date_range: Optional tuple (start_date, end_date) in YYYY/MM/DD format
            max_results: Maximum results to retrieve (default: 1000)
            analyze_keywords: Whether to analyze keywords (default: True)
            include_trends: Whether to include trend analysis (default: True)
        
        Returns:
        {
            "articles": List[Dict],  # Article data
            "keyword_analysis": {
                "top_keywords": [...],  # Research hotspots
                "trends": {...}  # Trend data over time
            },
            "publication_counts": {...},  # Publication count analysis
            "query_used": str,
            "article_count": int
        }
        """
        # Run search in executor (Bio.Entrez is blocking)
        loop = asyncio.get_event_loop()
        articles = await loop.run_in_executor(
            None,
            lambda: self.searcher.search(
                advanced_search=query,
                date_range=date_range,
                max_results=max_results
            )
        )
        
        result = {
            "articles": articles,
            "query_used": query,
            "article_count": len(articles)
        }
        
        # Analyze keywords if requested
        if analyze_keywords and articles:
            # Export to temporary JSON for analyzer
            import tempfile
            import json
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
                json.dump({
                    "metadata": {"export_time": "temp"},
                    "articles": articles
                }, f)
                temp_file = f.name
            
            try:
                keyword_analysis = await loop.run_in_executor(
                    None,
                    lambda: self.analyzer.analyze_research_keywords(
                        articles,
                        top_n=20,
                        include_trends=include_trends
                    )
                )
                
                pub_counts = await loop.run_in_executor(
                    None,
                    lambda: self.analyzer.analyze_publication_count(articles, months_per_period=3)
                )
                
                result["keyword_analysis"] = keyword_analysis
                result["publication_counts"] = pub_counts
            finally:
                # Clean up temp file
                try:
                    os.unlink(temp_file)
                except:
                    pass
        
        return result
    
    def get_top_keywords(self, analysis_result: Dict[str, Any], top_n: int = 10) -> List[str]:
        """
        Extract top keywords from analysis result.
        
        Useful for mechanism discovery.
        """
        if "keyword_analysis" not in analysis_result:
            return []
        
        top_keywords = analysis_result["keyword_analysis"].get("top_keywords", [])
        return [kw["keyword"] for kw in top_keywords[:top_n]]
```

---

### **File 2: `api/services/research_intelligence/parsers/pubmed_deep_parser.py`**

**Purpose**: Wrapper around pubmed_parser for deep full-text parsing

```python
"""
Deep PubMed Parser using pubmed_parser framework.

Provides:
- Full-text article parsing (not just abstracts)
- Citation/reference extraction
- Table and figure extraction
- Paragraph parsing with citation context
"""

import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Any
import httpx
import asyncio

# Add pubmed_parser to path
pubmed_parser_path = Path(__file__).resolve().parent.parent.parent.parent.parent / ".github" / "frameworks" / "pubmearch-main" / "pubmed_parser-master"
sys.path.insert(0, str(pubmed_parser_path))

try:
    from pubmed_parser.pubmed_oa_parser import (
        parse_pubmed_xml,
        parse_pubmed_references,
        parse_pubmed_paragraph,
        parse_pubmed_table,
        parse_pubmed_caption
    )
    from pubmed_parser.medline_parser import parse_medline_xml
    from pubmed_parser.pubmed_web_parser import parse_xml_web
    PARSER_AVAILABLE = True
except ImportError as e:
    print(f"⚠️ pubmed_parser not available: {e}")
    PARSER_AVAILABLE = False


class DeepPubMedParser:
    """
    Deep parsing of PubMed papers using pubmed_parser.
    
    Advantages:
    - Parse full-text PubMed OA XML (not just abstracts!)
    - Extract citations/references
    - Extract tables, figures, captions
    - Parse paragraphs with citation context
    """
    
    def __init__(self):
        if not PARSER_AVAILABLE:
            raise RuntimeError("pubmed_parser not available")
    
    async def parse_full_text_from_pmc(
        self,
        pmc_id: str,
        download_if_needed: bool = True
    ) -> Dict[str, Any]:
        """
        Parse full-text article from PubMed Central.
        
        Args:
            pmc_id: PMC ID (e.g., "PMC1234567" or "1234567")
            download_if_needed: Download XML if not cached
        
        Returns:
        {
            "pmid": str,
            "pmc": str,
            "title": str,
            "abstract": str,
            "full_text": str,
            "paragraphs": [
                {
                    "text": "...",
                    "section": "Methods",
                    "reference_ids": ["ref1", "ref2"]
                }
            ],
            "tables": [...],
            "figures": [...],
            "citations": [...]
        }
        """
        # Clean PMC ID
        pmc_id = pmc_id.replace("PMC", "").strip()
        
        # Download XML if needed
        xml_path = await self._download_pmc_xml(pmc_id) if download_if_needed else None
        
        if not xml_path:
            # Fallback to web parsing
            return await self._parse_from_web(pmc_id)
        
        # Parse XML
        loop = asyncio.get_event_loop()
        parsed = await loop.run_in_executor(
            None,
            lambda: parse_pubmed_xml(str(xml_path))
        )
        
        # Extract additional data
        paragraphs = await loop.run_in_executor(
            None,
            lambda: parse_pubmed_paragraph(str(xml_path), all_paragraph=True)
        )
        
        tables = await loop.run_in_executor(
            None,
            lambda: parse_pubmed_table(str(xml_path), return_xml=False)
        )
        
        citations = await loop.run_in_executor(
            None,
            lambda: parse_pubmed_references(str(xml_path))
        )
        
        return {
            **parsed,
            "paragraphs": paragraphs,
            "tables": tables,
            "citations": citations,
            "full_text": self._combine_paragraphs(paragraphs)
        }
    
    async def parse_medline_batch(
        self,
        pmids: List[str]
    ) -> List[Dict[str, Any]]:
        """
        Parse multiple MEDLINE records efficiently.
        
        Uses E-utils web API for batch parsing.
        """
        results = []
        
        # Parse in batches of 100 (E-utils limit)
        batch_size = 100
        for i in range(0, len(pmids), batch_size):
            batch = pmids[i:i+batch_size]
            
            # Parse each in executor
            loop = asyncio.get_event_loop()
            batch_results = await asyncio.gather(*[
                loop.run_in_executor(None, lambda pmid=pmid: parse_xml_web(pmid, save_xml=False))
                for pmid in batch
            ])
            
            results.extend([r for r in batch_results if r])
        
        return results
    
    def _combine_paragraphs(self, paragraphs: List[Dict]) -> str:
        """Combine paragraphs into full text."""
        return "\n\n".join([
            f"{p.get('section', 'Unknown')}: {p.get('text', '')}"
            for p in paragraphs
        ])
    
    async def _download_pmc_xml(self, pmc_id: str) -> Optional[Path]:
        """Download PMC XML if available."""
        # Implementation: Download from PMC FTP or API
        # For now, return None (use web parsing fallback)
        return None
    
    async def _parse_from_web(self, pmc_id: str) -> Dict[str, Any]:
        """Fallback to web parsing."""
        loop = asyncio.get_event_loop()
        parsed = await loop.run_in_executor(
            None,
            lambda: parse_xml_web(pmc_id, save_xml=False)
        )
        return parsed or {}
```

---

### **File 3: `api/services/research_intelligence/orchestrator.py`**

**Purpose**: Orchestrates all portals + parsers + LLM + MOAT

```python
"""
Research Intelligence Orchestrator

Combines:
- pubmearch (PubMed search + analysis)
- pubmed_parser (deep parsing)
- LLM (synthesis)
- MOAT (mechanism analysis)
"""

from typing import Dict, List, Optional, Any
from .portals.pubmed_enhanced import EnhancedPubMedPortal
from .parsers.pubmed_deep_parser import DeepPubMedParser
from api.services.gpt_service import get_gpt_service


class ResearchIntelligenceOrchestrator:
    """
    Orchestrates multi-portal research with deep parsing and LLM synthesis.
    """
    
    def __init__(self):
        # Portals
        self.pubmed = EnhancedPubMedPortal()
        self.pubmed_parser = DeepPubMedParser()
        
        # LLM
        self.llm = get_gpt_service()
    
    async def research_question(
        self,
        question: str,
        context: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Full research pipeline for a question.
        
        Args:
            question: Natural language question (e.g., "How do purple potatoes help with ovarian cancer?")
            context: Patient context (disease, treatment_line, biomarkers)
        
        Returns:
        {
            "research_plan": {...},
            "portal_results": {
                "pubmed": {...},
                "keyword_analysis": {...}
            },
            "parsed_content": {
                "full_text_articles": [...],
                "mechanisms_extracted": [...]
            },
            "synthesized_findings": {...},
            "moat_analysis": {...}
        }
        """
        # [1] Formulate research plan (LLM)
        research_plan = await self._formulate_research_plan(question, context)
        
        # [2] Query portals (parallel)
        portal_results = await self._query_portals(research_plan)
        
        # [3] Deep parse top papers
        parsed_content = await self._deep_parse_top_papers(portal_results)
        
        # [4] LLM synthesis
        synthesized_findings = await self._synthesize_findings(
            portal_results,
            parsed_content,
            research_plan
        )
        
        # [5] MOAT analysis
        moat_analysis = await self._run_moat_analysis(
            synthesized_findings,
            context
        )
        
        return {
            "research_plan": research_plan,
            "portal_results": portal_results,
            "parsed_content": parsed_content,
            "synthesized_findings": synthesized_findings,
            "moat_analysis": moat_analysis
        }
    
    async def _formulate_research_plan(
        self,
        question: str,
        context: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Use LLM to decompose question into research plan."""
        prompt = f"""Decompose this research question into a structured research plan:

Question: {question}
Context: {context}

Return JSON:
{{
    "primary_question": "...",
    "entities": {{
        "compound": "...",
        "disease": "...",
        "mechanisms_of_interest": [...]
    }},
    "sub_questions": [...],
    "portal_queries": {{
        "pubmed": [...]
    }}
}}"""
        
        response = await self.llm.chat(prompt, temperature=0.3)
        # Parse JSON response
        import json
        return json.loads(response)
    
    async def _query_portals(
        self,
        research_plan: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Query all portals in parallel."""
        pubmed_queries = research_plan.get("portal_queries", {}).get("pubmed", [])
        
        # Query PubMed with analysis
        pubmed_results = await self.pubmed.search_with_analysis(
            query=pubmed_queries[0] if pubmed_queries else "",
            max_results=1000,
            analyze_keywords=True,
            include_trends=True
        )
        
        return {
            "pubmed": pubmed_results,
            "keyword_analysis": pubmed_results.get("keyword_analysis", {}),
            "top_keywords": self.pubmed.get_top_keywords(pubmed_results, top_n=20)
        }
    
    async def _deep_parse_top_papers(
        self,
        portal_results: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Parse full-text for top papers."""
        articles = portal_results.get("pubmed", {}).get("articles", [])
        
        # Get top 10 papers (by keyword relevance or date)
        top_papers = articles[:10]
        
        # Extract PMC IDs
        pmc_ids = [
            article.get("pmc", "").replace("PMC", "")
            for article in top_papers
            if article.get("pmc")
        ]
        
        # Parse full-text
        parsed_articles = []
        for pmc_id in pmc_ids[:5]:  # Limit to 5 for now (expensive)
            try:
                full_text = await self.pubmed_parser.parse_full_text_from_pmc(pmc_id)
                parsed_articles.append(full_text)
            except Exception as e:
                print(f"⚠️ Failed to parse PMC{pmc_id}: {e}")
                continue
        
        return {
            "full_text_articles": parsed_articles,
            "parsed_count": len(parsed_articles)
        }
    
    async def _synthesize_findings(
        self,
        portal_results: Dict[str, Any],
        parsed_content: Dict[str, Any],
        research_plan: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Use LLM to synthesize findings from all sources."""
        # Combine all content
        abstracts = "\n\n".join([
            f"Title: {a.get('title', '')}\nAbstract: {a.get('abstract', '')}"
            for a in portal_results.get("pubmed", {}).get("articles", [])[:20]
        ])
        
        full_texts = "\n\n".join([
            f"Title: {ft.get('full_title', '')}\nFull Text: {ft.get('full_text', '')[:5000]}"
            for ft in parsed_content.get("full_text_articles", [])
        ])
        
        keywords = portal_results.get("top_keywords", [])
        
        prompt = f"""Synthesize research findings from these sources:

Research Question: {research_plan.get('primary_question', '')}

Abstracts ({len(portal_results.get('pubmed', {}).get('articles', []))} papers):
{abstracts[:10000]}

Full-Text Articles ({len(parsed_content.get('full_text_articles', []))} papers):
{full_texts[:10000]}

Top Keywords (Research Hotspots):
{', '.join(keywords[:20])}

Extract and synthesize:
1. Mechanisms of action (how it works)
2. Evidence strength (RCTs, in vitro, etc.)
3. Confidence scores
4. Knowledge gaps

Return JSON:
{{
    "mechanisms": [
        {{
            "mechanism": "...",
            "target": "...",
            "evidence": "...",
            "confidence": 0.0-1.0,
            "sources": ["pmid1", "pmid2"]
        }}
    ],
    "evidence_summary": "...",
    "knowledge_gaps": [...],
    "overall_confidence": 0.0-1.0
}}"""
        
        response = await self.llm.chat(prompt, temperature=0.3, max_tokens=2000)
        import json
        return json.loads(response)
    
    async def _run_moat_analysis(
        self,
        synthesized_findings: Dict[str, Any],
        context: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Run MOAT mechanism analysis on synthesized findings."""
        # Import MOAT services
        from api.services.food_spe_integration import FoodSPEIntegrationService
        from api.services.food_treatment_line_service import compute_food_treatment_line_features
        
        mechanisms = synthesized_findings.get("mechanisms", [])
        
        # Map mechanisms to pathways
        pathways = []
        for mech in mechanisms:
            # Map mechanism to cancer pathway
            pathway = self._map_mechanism_to_pathway(mech.get("mechanism", ""))
            if pathway:
                pathways.append(pathway)
        
        # Run MOAT analysis
        # (Use existing MOAT framework)
        
        return {
            "pathways": pathways,
            "mechanisms": mechanisms,
            "treatment_line_analysis": {},  # From MOAT
            "biomarker_analysis": {}  # From MOAT
        }
    
    def _map_mechanism_to_pathway(self, mechanism: str) -> Optional[str]:
        """Map mechanism name to cancer pathway."""
        mechanism_lower = mechanism.lower()
        
        pathway_map = {
            "angiogenesis": "angiogenesis",
            "vegf": "angiogenesis",
            "inflammation": "inflammation",
            "nf-κb": "inflammation",
            "dna repair": "dna_repair",
            "brca": "dna_repair",
            "apoptosis": "apoptosis",
            "cell cycle": "cell_cycle"
        }
        
        for key, pathway in pathway_map.items():
            if key in mechanism_lower:
                return pathway
        
        return None
```

---

### **File 4: `api/routers/research_intelligence.py`**

**Purpose**: API endpoint for research intelligence

```python
"""
Research Intelligence API Router

Endpoint: POST /api/research/intelligence
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import Dict, List, Optional, Any
from api.services.research_intelligence.orchestrator import ResearchIntelligenceOrchestrator

router = APIRouter(prefix="/api/research", tags=["Research Intelligence"])

class ResearchIntelligenceRequest(BaseModel):
    """Request for research intelligence."""
    question: str
    context: Dict[str, Any]  # disease, treatment_line, biomarkers
    portals: Optional[List[str]] = ["pubmed"]  # Which portals to query
    synthesize: bool = True
    run_moat_analysis: bool = True

@router.post("/intelligence")
async def research_intelligence(request: ResearchIntelligenceRequest):
    """
    Research intelligence endpoint.
    
    Uses pubmearch + pubmed_parser + LLM + MOAT to answer research questions.
    
    Example:
    {
        "question": "How do purple potatoes help with ovarian cancer?",
        "context": {
            "disease": "ovarian_cancer_hgs",
            "treatment_line": "L2",
            "biomarkers": {"HRD": "POSITIVE"}
        }
    }
    """
    try:
        orchestrator = ResearchIntelligenceOrchestrator()
        
        result = await orchestrator.research_question(
            question=request.question,
            context=request.context
        )
        
        return result
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
```

---

## 🧪 TESTING PLAN

### **Test 1: pubmearch Integration**

```python
from api.services.research_intelligence.portals.pubmed_enhanced import EnhancedPubMedPortal

portal = EnhancedPubMedPortal()
results = await portal.search_with_analysis(
    query="purple potatoes AND ovarian cancer",
    max_results=100,
    analyze_keywords=True
)

assert len(results['articles']) > 0
assert 'keyword_analysis' in results
assert len(results['keyword_analysis']['top_keywords']) > 0
```

---

### **Test 2: pubmed_parser Integration**

```python
from api.services.research_intelligence.parsers.pubmed_deep_parser import DeepPubMedParser

parser = DeepPubMedParser()
full_text = await parser.parse_full_text_from_pmc("PMC1234567")

assert 'full_text' in full_text or 'abstract' in full_text
assert 'paragraphs' in full_text or 'title' in full_text
```

---

### **Test 3: End-to-End Research**

```python
from api.services.research_intelligence.orchestrator import ResearchIntelligenceOrchestrator

orchestrator = ResearchIntelligenceOrchestrator()
result = await orchestrator.research_question(
    question="How do purple potatoes help with ovarian cancer?",
    context={"disease": "ovarian_cancer_hgs", "treatment_line": "L2"}
)

assert 'research_plan' in result
assert 'portal_results' in result
assert 'synthesized_findings' in result
assert 'moat_analysis' in result
```

---

## ✅ DEPENDENCIES

### **Install pubmearch**:

```bash
cd .github/frameworks/pubmearch-main
pip install -e .
```

**Dependencies**: `biopython`, `mcp` (optional)

---

### **Install pubmed_parser**:

```bash
cd .github/frameworks/pubmearch-main/pubmed_parser-master
pip install -e .
```

**Dependencies**: `lxml`, `unidecode`, `requests`

---

## 🚀 READY TO BUILD

**Frameworks**: ✅ **AUDITED & TESTED**  
**Integration Plan**: ✅ **DESIGNED**  
**Implementation Files**: ✅ **SPECIFIED**

**Next**: Build the modular research intelligence layer!




