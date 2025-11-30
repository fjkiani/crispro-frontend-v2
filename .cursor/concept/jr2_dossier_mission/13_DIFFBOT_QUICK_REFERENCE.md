# ⚔️ DIFFBOT QUICK REFERENCE - TRIAL SCRAPING ⚔️

**Status**: ✅ **ALREADY INTEGRATED** - No setup needed!  
**Location**: `api/routers/evidence/extraction.py`

---

## 🎯 **USING DIFFBOT FOR TRIAL SCRAPING**

### **Option 1: Use Existing Endpoint (Recommended)**

```python
import httpx

async def scrape_trial_via_endpoint(nct_id: str) -> dict:
    """Use existing /api/evidence/extract endpoint."""
    url = f"https://clinicaltrials.gov/study/{nct_id}"
    
    async with httpx.AsyncClient() as client:
        response = await client.post(
            "http://localhost:8000/api/evidence/extract",  # Your backend URL
            json={"url": url}
        )
        response.raise_for_status()
        data = response.json()
    
    # Diffbot returns: {title, author, date, site_name, text, html, tags}
    # Parse HTML to extract structured trial data
    from bs4 import BeautifulSoup
    soup = BeautifulSoup(data.get("html", ""), 'html.parser')
    
    return {
        'full_html': data.get("html"),
        'full_text': data.get("text"),
        'inclusion_criteria_full': extract_inclusion(soup),
        'exclusion_criteria_full': extract_exclusion(soup),
        # ... other fields
    }
```

---

### **Option 2: Call Diffbot API Directly**

```python
import httpx
from api.config import DIFFBOT_TOKEN

async def scrape_trial_direct(nct_id: str) -> dict:
    """Call Diffbot API directly."""
    if not DIFFBOT_TOKEN:
        raise ValueError("DIFFBOT_TOKEN not configured")
    
    url = f"https://clinicaltrials.gov/study/{nct_id}"
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
    
    # Parse HTML
    from bs4 import BeautifulSoup
    soup = BeautifulSoup(obj.get("html", ""), 'html.parser')
    
    return {
        'full_html': obj.get("html"),
        'full_text': obj.get("text"),
        'inclusion_criteria_full': extract_inclusion(soup),
        # ... other fields
    }
```

---

## 📋 **DIFFBOT RESPONSE STRUCTURE**

```python
{
    "title": "Study to Evaluate INCB123667...",
    "author": None,  # Usually None for ClinicalTrials.gov
    "date": None,  # Usually None
    "site_name": "ClinicalTrials.gov",
    "text": "Full text content...",
    "html": "<html>...</html>",  # Full HTML page
    "tags": []  # Usually empty for ClinicalTrials.gov
}
```

**Key Fields**:
- `html`: Full HTML page (parse with BeautifulSoup)
- `text`: Clean text extraction (useful for NLP)
- `title`: Page title

---

## 🔧 **EXTRACTING STRUCTURED DATA FROM HTML**

```python
from bs4 import BeautifulSoup

def extract_inclusion(soup: BeautifulSoup) -> str:
    """Extract inclusion criteria from ClinicalTrials.gov HTML."""
    # Look for "Eligibility Criteria" section
    criteria_section = soup.find('div', {'id': 'eligibility-criteria'})
    if criteria_section:
        inclusion = criteria_section.find('div', string='Inclusion Criteria')
        if inclusion:
            return inclusion.find_next_sibling('div').get_text(strip=True)
    return ""

def extract_exclusion(soup: BeautifulSoup) -> str:
    """Extract exclusion criteria."""
    criteria_section = soup.find('div', {'id': 'eligibility-criteria'})
    if criteria_section:
        exclusion = criteria_section.find('div', string='Exclusion Criteria')
        if exclusion:
            return exclusion.find_next_sibling('div').get_text(strip=True)
    return ""

def extract_interventions(soup: BeautifulSoup) -> list:
    """Extract intervention names."""
    interventions = []
    interv_section = soup.find('div', {'id': 'interventions'})
    if interv_section:
        for li in interv_section.find_all('li'):
            interventions.append(li.get_text(strip=True))
    return interventions

def extract_start_date(soup: BeautifulSoup) -> str:
    """Extract study start date."""
    start_elem = soup.find('td', string='Study Start')
    if start_elem:
        return start_elem.find_next_sibling('td').text.strip()
    return None
```

---

## ⚙️ **CONFIGURATION**

**Environment Variable**:
```bash
DIFFBOT_TOKEN=your_token_here
```

**Check if Configured**:
```python
from api.config import DIFFBOT_TOKEN

if not DIFFBOT_TOKEN:
    raise ValueError("DIFFBOT_TOKEN not set in environment")
```

---

## 🚨 **ERROR HANDLING**

```python
try:
    data = await scrape_trial_page(nct_id)
except ValueError as e:
    # Diffbot failed or token missing
    logger.error(f"Diffbot error: {e}")
    # Fallback: Use AstraDB record only
    data = astradb_record
except httpx.HTTPError as e:
    # Network error
    logger.error(f"HTTP error: {e}")
    # Retry or fallback
```

---

## 💡 **BEST PRACTICES**

1. **Cache Results**: Diffbot responses are expensive - cache for 24 hours
2. **Parse HTML**: Diffbot gives you HTML - use BeautifulSoup to extract structured data
3. **Fallback Strategy**: If Diffbot fails, use `astradb_record` (truncated but available)
4. **Rate Limits**: Diffbot handles rate limiting - no manual delays needed

---

**Reference**: See [05_IMPLEMENTATION_GUIDE.md](./05_IMPLEMENTATION_GUIDE.md) for complete implementation

