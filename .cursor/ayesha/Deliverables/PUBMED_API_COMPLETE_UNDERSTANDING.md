# ⚔️ PUBMED API COMPLETE UNDERSTANDING

**Author**: Zo (Lead AI Agent)  
**Date**: January 14, 2025  
**Purpose**: Complete documentation of PubMed API integration in Ayesha backend

---

## 🎯 EXECUTIVE SUMMARY

**Status**: ✅ **100% UNDERSTANDING**  
**Implementation Locations**: 
- `api/routers/evidence2.py:570-609` (esearch + esummary)
- `api/services/enhanced_evidence_service.py:163-280` (esearch + efetch with XML)
- `api/services/evidence/literature_client.py:51-136` (wrapper client)

**API**: NCBI E-utils (Entrez Programming Utilities)  
**Base URL**: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/`

---

## 📡 NCBI E-UTILS API FLOW

### **Two-Step Process**

#### **Step 1: Search for PMIDs** (`esearch.fcgi`)

**Endpoint**: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi`

**Parameters**:
```python
{
    "db": "pubmed",                    # Database: PubMed
    "retmode": "json",                 # Response format: JSON
    "retmax": "25",                    # Max results (default: 20)
    "term": "query_string",            # PubMed query syntax
    "api_key": os.getenv("NCBI_API_KEY", ""),  # Optional API key (rate limit boost)
    "email": os.getenv("NCBI_EMAIL", ""),      # Required: contact email
    "tool": "crispro"                  # Tool identifier
}
```

**Response Structure**:
```json
{
    "esearchresult": {
        "idlist": ["12345678", "23456789", ...],  # List of PMIDs
        "count": "25",
        "retmax": "25",
        "retstart": "0"
    }
}
```

**Query Syntax** (PubMed-specific):
- `"gene"[tiab]` - Search in title/abstract
- `"disease"[mh]` - MeSH term search
- `AND`, `OR`, `NOT` - Boolean operators
- `english[lang]` - Language filter
- `"2015:2024"[pdat]` - Publication date range

**Example Query**:
```python
term = '("BRAF"[tiab] OR "BRAF mutation"[tiab]) AND ("ovarian cancer"[mh] OR "ovarian carcinoma"[tiab]) AND english[lang]'
```

---

#### **Step 2: Fetch Paper Details** (Two Options)

##### **Option A: Summary** (`esummary.fcgi`) - **FAST, LIMITED**

**Endpoint**: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi`

**Parameters**:
```python
{
    "db": "pubmed",
    "retmode": "json",
    "id": "12345678,23456789",  # Comma-separated PMIDs
    "api_key": os.getenv("NCBI_API_KEY", ""),
    "email": os.getenv("NCBI_EMAIL", ""),
    "tool": "crispro"
}
```

**Response Structure**:
```json
{
    "result": {
        "uids": ["12345678", "23456789"],
        "12345678": {
            "title": "Paper Title",
            "pubdate": "2024 Jan",
            "fulljournalname": "Journal Name",
            "source": "Journal Abbreviation",
            "pubtype": ["Journal Article", "Review"],
            "authors": [{"name": "Author1"}, {"name": "Author2"}]
        }
    }
}
```

**Limitations**:
- ❌ No abstracts (only metadata)
- ✅ Fast (single JSON response)
- ✅ Good for ranking/scoring

**Used In**: `evidence2.py:570-609` (quick paper ranking)

---

##### **Option B: Full Fetch** (`efetch.fcgi`) - **SLOW, COMPLETE**

**Endpoint**: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi`

**Parameters**:
```python
{
    "db": "pubmed",
    "id": "12345678,23456789",
    "retmode": "xml",                  # ⚠️ XML ONLY (not JSON!)
    "rettype": "abstract"              # Return abstracts
}
```

**Response Format**: **XML** (not JSON!)

**XML Structure**:
```xml
<PubmedArticleSet>
    <PubmedArticle>
        <MedlineCitation>
            <PMID>12345678</PMID>
            <Article>
                <ArticleTitle>Paper Title</ArticleTitle>
                <Abstract>
                    <AbstractText Label="BACKGROUND">Background text...</AbstractText>
                    <AbstractText Label="METHODS">Methods text...</AbstractText>
                    <AbstractText Label="RESULTS">Results text...</AbstractText>
                </Abstract>
                <AuthorList>
                    <Author>
                        <LastName>Smith</LastName>
                        <ForeName>John</ForeName>
                    </Author>
                </AuthorList>
            </Article>
        </MedlineCitation>
    </PubmedArticle>
</PubmedArticleSet>
```

**Parsing** (Python):
```python
import xml.etree.ElementTree as ET

root = ET.fromstring(xml_response)
for article in root.findall('.//PubmedArticle'):
    pmid = article.find('.//PMID').text
    title = article.find('.//ArticleTitle').text
    abstract_texts = article.findall('.//AbstractText')
    abstract = " ".join([elem.text for elem in abstract_texts if elem.text])
```

**Used In**: `enhanced_evidence_service.py:163-280` (full-text extraction)

---

## 🔄 IMPLEMENTATION PATTERNS

### **Pattern 1: Quick Search + Summary** (`evidence2.py`)

**Use Case**: Fast paper ranking without abstracts

**Flow**:
1. Build query string (gene + disease + MoA terms)
2. Call `esearch.fcgi` → Get PMID list
3. Call `esummary.fcgi` → Get metadata (title, year, journal, pubtype)
4. Score and rank papers (publication type weighting)
5. Return top N results

**Code Location**: `api/routers/evidence2.py:570-609`

**Timeout**: 20s per request

---

### **Pattern 2: Full Search + Fetch** (`enhanced_evidence_service.py`)

**Use Case**: Complete paper extraction with abstracts

**Flow**:
1. Build query string (compound + disease + pathways)
2. Call `esearch.fcgi` → Get PMID list
3. Call `efetch.fcgi` → Get full XML (title, abstract, authors)
4. Parse XML with ElementTree
5. Extract structured data (title, abstract, authors, PMID)
6. Return papers with full text

**Code Location**: `api/services/enhanced_evidence_service.py:163-280`

**Retry Logic**: 3 attempts with exponential backoff (2s, 4s, 6s)

**Timeout**: Configurable (default: 30s)

---

### **Pattern 3: Wrapper Client** (`literature_client.py`)

**Use Case**: Unified interface for evidence gathering

**Flow**:
1. Call `/api/evidence/literature` endpoint (internal)
2. Endpoint uses Pattern 1 or Pattern 2 internally
3. Filter results by MoA terms (drug_name OR drug_moa in title/abstract)
4. Score evidence strength (publication type weighting)
5. Return `EvidenceHit` with top results and strength

**Code Location**: `api/services/evidence/literature_client.py:51-136`

---

## 🛡️ ERROR HANDLING & RETRY LOGIC

### **Retry Strategy** (`enhanced_evidence_service.py:168-297`)

```python
max_retries = 3
retry_delay = 2.0  # seconds

for attempt in range(max_retries):
    try:
        # Make request
        response = await client.get(url, params=params)
        
        if response.status_code != 200:
            if attempt < max_retries - 1:
                await asyncio.sleep(retry_delay * (attempt + 1))
                continue
            return []  # Give up after 3 attempts
        
        # Parse response
        data = response.json()  # or ET.fromstring() for XML
        
    except (json.JSONDecodeError, ET.ParseError) as e:
        if attempt < max_retries - 1:
            await asyncio.sleep(retry_delay * (attempt + 1))
            continue
        return []  # Give up
    
    except Exception as e:
        if attempt < max_retries - 1:
            await asyncio.sleep(retry_delay * (attempt + 1))
            continue
        return []  # Give up
```

**Error Types Handled**:
- HTTP errors (non-200 status)
- JSON parse errors (malformed response)
- XML parse errors (malformed XML)
- Network timeouts
- Connection errors

**Fallback**: Return empty list (graceful degradation)

---

## 📊 QUERY BUILDING STRATEGIES

### **Strategy 1: Gene + Variant + Disease** (`evidence2.py:570`)

```python
clauses = []
if gene:
    clauses.append(f'"{gene}"[tiab]')
if hgvs_p:
    clauses.append(f'"{hgvs_p}"[tiab]')
if disease:
    clauses.append(f'"{disease}"[mh] OR "{disease}"[tiab]')
if moa_terms:
    clauses.append(f'({" OR ".join([f\'"{term}"[tiab]\' for term in moa_terms])})')

term = " AND ".join(clauses) or (f'"{gene}"[tiab]' if gene else 'multiple myeloma[mh]')
```

**Fallback Chain**:
1. Try: gene + variant + disease + MoA
2. Fallback: gene + disease (no variant filter)
3. Fallback: disease only

---

### **Strategy 2: Compound + Disease + Pathways** (`enhanced_evidence_service.py:120`)

```python
query_terms = []
if compound:
    query_terms.append(f'"{compound}"[tiab]')
if disease:
    # Map disease codes to PubMed-friendly terms
    disease_map = {
        "ovarian_cancer_hgs": "ovarian cancer",
        "breast_cancer_tnbc": "triple negative breast cancer",
        ...
    }
    disease_term = disease_map.get(disease, disease)
    query_terms.append(f'"{disease_term}"[tiab]')
if pathways:
    pathway_terms = [f'"{pathway}"[tiab]' for pathway in pathways]
    query_terms.append(f'({" OR ".join(pathway_terms)})')

query = " AND ".join(query_terms)
```

---

## 🎯 EVIDENCE SCORING

### **Publication Type Weighting** (`literature_client.py:18-48`)

```python
def _score_evidence_from_results(top_results: List[Dict]) -> float:
    score = 0.0
    for r in top_results[:3]:  # Top 3 papers
        pub_types = " ".join([t.lower() for t in (r.get("publication_types") or [])])
        title = r.get("title", "").lower()
        
        if "randomized" in pub_types or "randomized" in title:
            score += 0.5  # Highest weight
        elif "guideline" in pub_types or "practice" in title:
            score += 0.35
        elif "review" in pub_types or "meta" in title:
            score += 0.25
        else:
            score += 0.15  # Lowest weight
    
    return min(1.0, score)  # Cap at 1.0
```

**MoA Boost** (`literature_client.py:110-120`):
```python
moa_hits = 0
if drug_moa:
    for paper in top_results:
        if drug_moa.lower() in paper.get("title", "").lower() or \
           drug_moa.lower() in paper.get("abstract", "").lower():
            moa_hits += 1

strength = min(1.0, base_strength + 0.10 * moa_hits)
```

---

## ⚙️ CONFIGURATION

### **Environment Variables**

```bash
NCBI_API_KEY=your_api_key      # Optional: Increases rate limit (10 req/s → 100 req/s)
NCBI_EMAIL=your_email          # Required: Contact email for NCBI
```

### **Rate Limits**

- **Without API Key**: 3 requests/second
- **With API Key**: 10 requests/second
- **Best Practice**: Use API key for production

### **Timeouts**

- **esearch**: 20s timeout
- **esummary**: 20s timeout
- **efetch**: 30s timeout (longer due to XML parsing)

---

## 🔍 KEY INSIGHTS

### **Why Two Endpoints?**

1. **esummary**: Fast, JSON, no abstracts → Good for ranking
2. **efetch**: Slow, XML, full abstracts → Good for evidence extraction

### **Why XML for efetch?**

- NCBI E-utils only supports XML for `efetch.fcgi`
- Must use `xml.etree.ElementTree` for parsing
- Abstracts can be multi-section (BACKGROUND, METHODS, RESULTS)

### **Query Building Best Practices**

1. **Use MeSH terms** (`[mh]`) for diseases (more precise)
2. **Use title/abstract** (`[tiab]`) for genes/variants (broader)
3. **Add language filter** (`english[lang]`) to reduce noise
4. **Use Boolean operators** (`AND`, `OR`) for complex queries
5. **Fallback chain**: Start specific, get broader if no results

---

## 📋 USAGE EXAMPLES

### **Example 1: Quick Paper Ranking**

```python
# evidence2.py pattern
query = '("BRAF"[tiab] OR "BRAF V600E"[tiab]) AND ("ovarian cancer"[mh]) AND english[lang]'
pmids = esearch(query, max_results=8)
papers = esummary(pmids)
ranked = score_and_rank(papers, moa_terms=["MAPK blockade"])
```

### **Example 2: Full Evidence Extraction**

```python
# enhanced_evidence_service.py pattern
query = '("curcumin"[tiab]) AND ("ovarian cancer"[tiab]) AND english[lang]'
pmids = esearch(query, max_results=20)
papers = efetch(pmids)  # Returns XML
parsed = parse_xml(papers)  # Extract title, abstract, authors
```

---

## ✅ COMPLETE UNDERSTANDING CHECKLIST

- [x] E-utils API flow (esearch → esummary/efetch)
- [x] Query building syntax (`[tiab]`, `[mh]`, Boolean operators)
- [x] XML parsing patterns (ElementTree)
- [x] Error handling and retry logic (3 attempts, exponential backoff)
- [x] Rate limiting considerations (API key, email, tool)
- [x] Evidence scoring (publication type weighting, MoA boost)
- [x] Fallback strategies (query simplification, empty results)

---

**Status**: ✅ **100% UNDERSTANDING ACHIEVED**  
**Last Updated**: January 14, 2025  
**By**: Zo (Lead AI Agent)

