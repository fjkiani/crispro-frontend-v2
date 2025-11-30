# ⚔️ JR2 IMPLEMENTATION GUIDE - CODE & SCHEMAS ⚔️

**Purpose**: Complete code examples, schemas, and implementation details for JR2

---

## 📋 **ASTRA DB SCHEMA**

```python
# AstraDB document structure (from clinical_trials_eligibility2)
trial_doc = {
    "_id": "https://clinicaltrials.gov/study/NCT06819007",
    "nct_id": "NCT06819007",
    "title": "Study to Evaluate INCB123667...",
    "status": "RECRUITING",  # RECRUITING|NOT_YET_RECRUITING|COMPLETED|TERMINATED
    "phase": "PHASE3",  # PHASE1|PHASE2|PHASE3|PHASE4|N/A
    "disease_category": "gynecologic_oncology",
    "disease_subcategory": "ovarian_cancer",
    "biomarker_requirements": None,  # Free text (parse with NLP)
    "locations_data": [  # Array of location objects
        {
            "facility": "Memorial Sloan Kettering Cancer Center",
            "city": "New York",
            "state": "NY",
            "country": "United States",
            "status": "RECRUITING",
            "contact": {...}
        }
    ],
    "eligibility_text": "Inclusion Criteria: Stage III/IV...",  # Truncated to 7500 bytes
    "description_text": "This study will compare...",  # Truncated to 7500 bytes
    "source_url": "https://clinicaltrials.gov/study/NCT06819007",
    "sponsor_name": "Incyte Corporation",
    "principal_investigator_name": None,
    "pi_contact_email": None,
    "study_coordinator_email": None,
    "primary_endpoint": None,  # Not in AstraDB (scrape from full page)
    "site_count": 50,
    "estimated_enrollment": 300,
    "mechanism_tags": None,
    "biomarker_requirements_gtm": None,
    "$vector": [0.123, 0.456, ...]  # 768-dim embedding
}
```

---

## 🕷️ **TRIAL SCRAPING (Diffbot - Already Integrated)**

**✅ Diffbot is already set up!** Use the existing endpoint.

```python
import httpx
from api.config import DIFFBOT_TOKEN

async def scrape_trial_page(nct_id: str) -> dict:
    """
    Scrape full ClinicalTrials.gov page using Diffbot.
    
    Uses existing Diffbot integration at api/routers/evidence/extraction.py
    
    Returns:
        {
            'inclusion_criteria_full': str,
            'exclusion_criteria_full': str,
            'interventions': List[str],
            'primary_endpoint': str,
            'study_start_date': str,
            'primary_completion_date': str,
            'locations_full': List[dict],
            'full_html': str,
            'full_text': str
        }
    """
    if not DIFFBOT_TOKEN:
        raise ValueError("DIFFBOT_TOKEN not configured")
    
    url = f"https://clinicaltrials.gov/study/{nct_id}"
    
    # Use Diffbot API directly (same as existing extraction endpoint)
    api_url = "https://api.diffbot.com/v3/article"
    params = {
        "token": DIFFBOT_TOKEN,
        "url": url,
        "fields": "title,author,date,siteName,tags,images,html,text",
    }
    
    async with httpx.AsyncClient(timeout=30) as client:
        response = await client.get(api_url, params=params)
        response.raise_for_status()
        js = response.json()
    
    obj = (js.get("objects") or [None])[0]
    if not obj:
        raise ValueError(f"Diffbot failed to extract data from {url}")
    
    # Parse HTML to extract structured data
    from bs4 import BeautifulSoup
    soup = BeautifulSoup(obj.get("html", ""), 'html.parser')
    
    # Extract study start date
    start_date_elem = soup.find('td', string='Study Start')
    study_start_date = start_date_elem.find_next_sibling('td').text.strip() if start_date_elem else None
    
    # Extract primary completion date
    completion_elem = soup.find('td', string='Primary Completion')
    primary_completion_date = completion_elem.find_next_sibling('td').text.strip() if completion_elem else None
    
    return {
        'inclusion_criteria_full': extract_inclusion(soup),
        'exclusion_criteria_full': extract_exclusion(soup),
        'interventions': extract_interventions(soup),
        'primary_endpoint': extract_primary_endpoint(soup),
        'study_start_date': study_start_date,
        'primary_completion_date': primary_completion_date,
        'locations_full': extract_locations(soup),
        'full_html': obj.get("html"),
        'full_text': obj.get("text")
    }

# Alternative: Use existing endpoint via HTTP
async def scrape_trial_page_via_endpoint(nct_id: str) -> dict:
    """Use existing /api/evidence/extract endpoint."""
    url = f"https://clinicaltrials.gov/study/{nct_id}"
    
    async with httpx.AsyncClient() as client:
        response = await client.post(
            "http://localhost:8000/api/evidence/extract",  # Adjust base URL
            json={"url": url}
        )
        response.raise_for_status()
        data = response.json()
    
    # Parse HTML from Diffbot response
    from bs4 import BeautifulSoup
    soup = BeautifulSoup(data.get("html", ""), 'html.parser')
    
    # Extract structured data (same as above)
    return {
        'inclusion_criteria_full': extract_inclusion(soup),
        'exclusion_criteria_full': extract_exclusion(soup),
        'interventions': extract_interventions(soup),
        'primary_endpoint': extract_primary_endpoint(soup),
        'study_start_date': extract_start_date(soup),
        'primary_completion_date': extract_completion_date(soup),
        'locations_full': extract_locations(soup),
        'full_html': data.get("html"),
        'full_text': data.get("text")
    }
```

**Existing Integration**:
- ✅ Endpoint: `POST /api/evidence/extract`
- ✅ File: `api/routers/evidence/extraction.py`
- ✅ Config: `DIFFBOT_TOKEN` from `api/config.py`
- ✅ No rate limits (handled by Diffbot)

---

## 💾 **CACHING STRATEGY**

```python
import json
from pathlib import Path
from datetime import datetime, timedelta

CACHE_DIR = Path(".cursor/ayesha/cache/")
CACHE_TTL_HOURS = 24  # 24 hours for scraped trial data

def cache_trial_data(nct_id: str, data: dict):
    """Cache scraped trial data to file."""
    cache_file = CACHE_DIR / f"trial_{nct_id}.json"
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(cache_file, 'w') as f:
        json.dump({
            'cached_at': datetime.now().isoformat(),
            'nct_id': nct_id,
            'data': data
        }, f, indent=2)

def get_cached_trial_data(nct_id: str) -> dict | None:
    """Get cached trial data if not expired."""
    cache_file = CACHE_DIR / f"trial_{nct_id}.json"
    
    if not cache_file.exists():
        return None
    
    with open(cache_file, 'r') as f:
        cached = json.load(f)
    
    # Check if expired
    cached_at = datetime.fromisoformat(cached['cached_at'])
    if datetime.now() - cached_at > timedelta(hours=CACHE_TTL_HOURS):
        return None  # Expired
    
    return cached['data']
```

---

## 💊 **DRUG MECHANISM DATABASE**

See complete 20-drug database in separate file: `DRUG_MECHANISM_DB.json`

**Structure**:
```python
DRUG_MECHANISM_DB = {
    'Trastuzumab Deruxtecan': {
        'class': 'Antibody-Drug Conjugate (ADC)',
        'target': 'HER2',
        'moa': 'HER2-targeted antibody delivers topoisomerase I inhibitor payload',
        'layman': 'Trojan horse targeting HER2, works even with low expression (IHC 1+)',
        'breakthrough': 'First ADC for low HER2 expression, bystander effect kills nearby cells'
    },
    # ... 19 more drugs
}
```

**Fallback Strategy**:
```python
def get_drug_mechanism(drug_name: str) -> dict:
    """Get drug mechanism with 3-tier fallback."""
    
    # Tier 1: Check DRUG_MECHANISM_DB
    if drug_name in DRUG_MECHANISM_DB:
        return DRUG_MECHANISM_DB[drug_name]
    
    # Tier 2: Check EnhancedEvidenceService
    evidence = evidence_service.query_drug_info(drug_name)
    if evidence:
        return {
            'class': evidence.get('drug_class', '[UNKNOWN CLASS]'),
            'moa': evidence.get('mechanism', '[UNKNOWN MECHANISM]'),
            'layman': f"[INFERRED FROM EVIDENCE] {evidence.get('summary', '')}",
            'breakthrough': '[NEEDS MANUAL RESEARCH]'
        }
    
    # Tier 3: Flag for Zo review
    return {
        'class': '[UNKNOWN - NEEDS ZO REVIEW]',
        'moa': f'[UNKNOWN DRUG: {drug_name}] - Add to DRUG_MECHANISM_DB',
        'layman': '[CANNOT GENERATE - DRUG NOT IN DATABASE]',
        'breakthrough': '[MANUAL RESEARCH REQUIRED]',
        'confidence': 'LOW - FLAG FOR ZO REVIEW'
    }
```

---

## 🔗 **EVIDENCE SYNTHESIS INTEGRATION**

```python
from api.services.enhanced_evidence_service import EnhancedEvidenceService

evidence_service = EnhancedEvidenceService()

def generate_clinical_evidence(drugs: List[str], cancer_type: str) -> str:
    """Generate evidence section using existing service."""
    markdown = ""
    
    for drug in drugs:
        # Query existing service
        evidence = evidence_service.query_rct_data(drug, cancer_type)
        
        if evidence:
            markdown += f"### **{drug} Track Record**:\n"
            for trial in evidence.get('rct_citations', []):
                markdown += f"- **{trial['cancer_type']} ({trial['trial_name']})**: {trial['result']} ([{trial['journal']} {trial['year']}](https://pubmed.ncbi.nlm.nih.gov/{trial['pubmed_id']}))\n"
            
            # Add key insights for ovarian
            if 'key_insights' in evidence and cancer_type in evidence['key_insights']:
                markdown += f"\n### **Key Insight for {cancer_type.replace('_', ' ').title()}**:\n"
                for insight in evidence['key_insights'][cancer_type]:
                    markdown += f"- {insight}\n"
            markdown += "\n"
        else:
            # Fallback: Use DRUG_MECHANISM_DB
            if drug in DRUG_MECHANISM_DB:
                mech = DRUG_MECHANISM_DB[drug]
                markdown += f"### **{drug}**:\n"
                markdown += f"- **Mechanism**: {mech['moa']}\n"
                markdown += f"- **Breakthrough**: {mech['breakthrough']}\n\n"
            else:
                markdown += f"### **{drug}**:\n"
                markdown += f"- **[NEEDS EVIDENCE RESEARCH]** - Unknown drug, no evidence in database\n\n"
    
    return markdown
```

---

## 🏗️ **DEPENDENCIES & UTILITIES**

```python
# Logging (ALREADY EXISTS)
from api.utils.logger import setup_logger
logger = setup_logger(__name__)

# Database connections (ALREADY EXISTS)
from api.services.database_connections import get_db_connections
db = get_db_connections()

# Error handling (CREATE NEW)
# File: api/utils/dossier_errors.py
class DossierGenerationError(Exception):
    """Base exception for dossier generation errors."""
    pass

class TrialScrapingError(DossierGenerationError):
    """Trial scraping failed."""
    pass

class EligibilityMatchingError(DossierGenerationError):
    """Eligibility matching failed."""
    pass
```

**Dependencies**:
- ✅ Diffbot: Already integrated (uses `DIFFBOT_TOKEN` from env)
- ✅ httpx: Already in requirements.txt
- ✅ BeautifulSoup: Needed for parsing Diffbot HTML (add if not present):
  ```
  beautifulsoup4==4.12.2
  lxml==4.9.3
  ```

---

## 📝 **CONFIDENCE FLAGS**

```python
# HIGH CONFIDENCE (backed by data)
"Patient is BRCA wildtype (germline test confirmed)"  # ✅ No flag needed

# MEDIUM CONFIDENCE (inferred from data)
"Trial likely requires HER2 IHC 1+ based on eligibility text [INFERRED]"

# LOW CONFIDENCE (needs verification)
"Drug mechanism appears to be PARP inhibition [NEEDS VERIFICATION - Zo review required]"

# MISSING DATA (critical gap)
"Primary endpoint not available [MISSING - scrape full page or contact trial site]"

# CONFLICTING DATA (needs resolution)
"AstraDB shows RECRUITING but ClinicalTrials.gov shows COMPLETED [CONFLICT - verify status]"
```

---

## 📁 **FILE NAMING & STORAGE**

```python
from datetime import datetime

def generate_dossier_filename(nct_id: str, patient_id: str) -> str:
    """Generate dossier filename with timestamp."""
    timestamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    return f"dossier_{nct_id}_{patient_id}_{timestamp}.md"

# Storage locations:
# - Generated: .cursor/ayesha/dossiers/{nct_id}/
# - Approved: .cursor/ayesha/dossiers/approved/{nct_id}/
# - Rejected: .cursor/ayesha/dossiers/rejected/{nct_id}/
```

---

**Next**: See [09_API_SPECIFICATIONS.md](./09_API_SPECIFICATIONS.md) for API endpoint code

