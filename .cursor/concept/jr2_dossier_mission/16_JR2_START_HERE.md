# ⚔️ JR2 - START HERE - IMMEDIATE ACTION PLAN ⚔️

**Date**: January 14, 2025  
**Status**: 🎯 **READY TO BUILD** - All data located, all tools ready

---

## ✅ **WHAT WE KNOW NOW**

1. ✅ **1,000 trials** in SQLite: `oncology-coPilot/oncology-backend-minimal/data/clinical_trials.db`
2. ✅ **Table name**: `clinical_trials` (not `trials`)
3. ✅ **Diffbot** already integrated: `api/routers/evidence/extraction.py`
4. ✅ **50 existing candidates**: `.cursor/ayesha/50_vector_candidates_for_jr2.json` (ready to use)

---

## 🎯 **IMMEDIATE NEXT STEPS (DO NOW)**

### **STEP 1: Verify You Can Access the Data** (5 min)

```python
# Test script: test_data_access.py
import sqlite3
import json

db_path = "oncology-coPilot/oncology-backend-minimal/data/clinical_trials.db"
conn = sqlite3.connect(db_path)
conn.row_factory = sqlite3.Row

cursor = conn.cursor()
cursor.execute("SELECT COUNT(*) as total FROM clinical_trials")
total = cursor.fetchone()['total']
print(f"✅ Found {total} trials in SQLite")

# Get sample trial
cursor.execute("SELECT nct_id, title, status FROM clinical_trials LIMIT 1")
sample = dict(cursor.fetchone())
print(f"Sample trial: {sample}")

conn.close()
```

**Run this**: `cd oncology-coPilot/oncology-backend-minimal && python test_data_access.py`

---

### **STEP 2: Build Trial Query Function** (15 min)

Create: `api/services/client_dossier/trial_querier.py`

```python
import sqlite3
from typing import List, Dict, Any
from pathlib import Path

def get_trials_from_sqlite(limit: int = 100) -> List[Dict[str, Any]]:
    """
    Get trials from SQLite for JR2's filtering pipeline.
    
    Args:
        limit: Max trials to return (0 = all)
    
    Returns:
        List of trial dictionaries with all fields
    """
    # Get absolute path to database
    project_root = Path(__file__).parent.parent.parent.parent.parent
    db_path = project_root / "oncology-coPilot" / "oncology-backend-minimal" / "data" / "clinical_trials.db"
    
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row
    
    cursor = conn.cursor()
    
    if limit > 0:
        cursor.execute("SELECT * FROM clinical_trials LIMIT ?", (limit,))
    else:
        cursor.execute("SELECT * FROM clinical_trials")
    
    trials = [dict(row) for row in cursor.fetchall()]
    
    # Parse JSON fields
    for trial in trials:
        for json_field in ['pis_json', 'orgs_json', 'sites_json']:
            if trial.get(json_field):
                try:
                    trial[json_field] = json.loads(trial[json_field])
                except:
                    trial[json_field] = []
    
    conn.close()
    return trials
```

---

### **STEP 3: Build Diffbot Scraper** (20 min)

Create: `api/services/client_dossier/trial_scraper.py`

```python
import httpx
import os
from typing import Dict, Any, Optional
from bs4 import BeautifulSoup
import json
from pathlib import Path
from datetime import datetime, timedelta

# Cache directory
CACHE_DIR = Path(".cursor/ayesha/cache/")
CACHE_TTL_HOURS = 24

def scrape_trial_via_diffbot(nct_id: str, use_cache: bool = True) -> Dict[str, Any]:
    """
    Scrape ClinicalTrials.gov page using Diffbot.
    
    Uses existing /api/evidence/extract endpoint or Diffbot API directly.
    """
    # Check cache first
    if use_cache:
        cached = get_cached_trial_data(nct_id)
        if cached:
            return cached
    
    # Use Diffbot via existing endpoint
    url = f"https://clinicaltrials.gov/study/{nct_id}"
    
    # Option 1: Use existing backend endpoint
    try:
        async with httpx.AsyncClient() as client:
            response = await client.post(
                "http://localhost:8000/api/evidence/extract",
                json={"url": url},
                timeout=30.0
            )
            response.raise_for_status()
            diffbot_data = response.json()
    except Exception as e:
        logger.warning(f"Diffbot endpoint failed: {e}, trying direct API...")
        # Option 2: Call Diffbot API directly
        diffbot_token = os.getenv("DIFFBOT_TOKEN")
        if not diffbot_token:
            raise ValueError("DIFFBOT_TOKEN not set")
        
        diffbot_url = f"https://api.diffbot.com/v3/article?token={diffbot_token}&url={url}"
        response = httpx.get(diffbot_url, timeout=30.0)
        response.raise_for_status()
        diffbot_data = response.json()
    
    # Parse Diffbot HTML to extract structured data
    html = diffbot_data.get('html', '')
    soup = BeautifulSoup(html, 'html.parser')
    
    # Extract structured fields
    trial_data = {
        'nct_id': nct_id,
        'inclusion_criteria_full': extract_inclusion(soup),
        'exclusion_criteria_full': extract_exclusion(soup),
        'interventions': extract_interventions(soup),
        'primary_endpoint': extract_primary_endpoint(soup),
        'study_start_date': extract_start_date(soup),
        'primary_completion_date': extract_completion_date(soup),
        'locations_full': extract_locations(soup),
        'scraped_at': datetime.now().isoformat()
    }
    
    # Cache result
    cache_trial_data(nct_id, trial_data)
    
    return trial_data

def extract_inclusion(soup) -> str:
    """Extract inclusion criteria from ClinicalTrials.gov HTML."""
    # Implementation: Find inclusion section in HTML
    # ... (use BeautifulSoup to parse)
    pass

def extract_exclusion(soup) -> str:
    """Extract exclusion criteria."""
    pass

def extract_interventions(soup) -> List[str]:
    """Extract intervention list."""
    pass

def extract_primary_endpoint(soup) -> str:
    """Extract primary endpoint."""
    pass

def extract_start_date(soup) -> Optional[str]:
    """Extract study start date."""
    pass

def extract_completion_date(soup) -> Optional[str]:
    """Extract primary completion date."""
    pass

def extract_locations(soup) -> List[Dict]:
    """Extract location data."""
    pass

# Caching functions
def cache_trial_data(nct_id: str, data: dict):
    """Cache scraped trial data."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file = CACHE_DIR / f"trial_{nct_id}.json"
    
    with open(cache_file, 'w') as f:
        json.dump({
            'cached_at': datetime.now().isoformat(),
            'nct_id': nct_id,
            'data': data
        }, f, indent=2)

def get_cached_trial_data(nct_id: str) -> Optional[dict]:
    """Get cached trial data if not expired."""
    cache_file = CACHE_DIR / f"trial_{nct_id}.json"
    
    if not cache_file.exists():
        return None
    
    with open(cache_file, 'r') as f:
        cached = json.load(f)
    
    # Check expiration
    cached_at = datetime.fromisoformat(cached['cached_at'])
    if datetime.now() - cached_at > timedelta(hours=CACHE_TTL_HOURS):
        return None
    
    return cached['data']
```

---

### **STEP 4: Start Filtering Pipeline** (30 min)

Create: `api/services/client_dossier/trial_filter.py`

```python
from typing import List, Dict, Any
from .trial_querier import get_trials_from_sqlite

def filter_50_candidates(patient_profile: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Replicate Zo's "1 in 700" filtering strategy.
    
    Hard filters (must pass):
    - Stage IV ✅
    - First-line ✅
    - Recruiting ✅
    - NYC metro (NY/NJ/CT) ✅
    
    Returns top 50 candidates.
    """
    # Get all trials from SQLite
    all_trials = get_trials_from_sqlite(limit=0)  # Get all 1000
    
    filtered = []
    
    for trial in all_trials:
        # Hard filter 1: Stage IV
        if not is_stage_iv(trial, patient_profile):
            continue
        
        # Hard filter 2: First-line
        if not is_first_line(trial):
            continue
        
        # Hard filter 3: Recruiting
        if trial.get('status') != 'RECRUITING':
            continue
        
        # Hard filter 4: NYC metro
        if not is_nyc_metro(trial):
            continue
        
        filtered.append(trial)
    
    # Sort by relevance (add scoring logic)
    # Return top 50
    return filtered[:50]
```

---

### **STEP 5: Test End-to-End** (15 min)

Create: `test_jr2_pipeline.py`

```python
from api.services.client_dossier.trial_querier import get_trials_from_sqlite
from api.services.client_dossier.trial_filter import filter_50_candidates
from api.services.client_dossier.trial_scraper import scrape_trial_via_diffbot

# Test 1: Query SQLite
print("Test 1: Querying SQLite...")
trials = get_trials_from_sqlite(limit=10)
print(f"✅ Got {len(trials)} trials")

# Test 2: Filter candidates
print("\nTest 2: Filtering candidates...")
patient = {
    "stage": "IVB",
    "treatment_line": "1",
    "location_state": "NY"
}
candidates = filter_50_candidates(patient)
print(f"✅ Filtered to {len(candidates)} candidates")

# Test 3: Scrape one trial
print("\nTest 3: Scraping trial with Diffbot...")
if candidates:
    nct_id = candidates[0]['nct_id']
    scraped = scrape_trial_via_diffbot(nct_id)
    print(f"✅ Scraped {nct_id}: {len(scraped)} fields")
```

---

## 📋 **BUILD ORDER (PRIORITY)**

1. ✅ **Verify data access** (5 min) - DO THIS FIRST
2. ✅ **Build trial querier** (15 min) - Read from SQLite
3. ✅ **Build Diffbot scraper** (20 min) - Scrape full eligibility
4. ✅ **Build filter logic** (30 min) - Replicate Zo's filters
5. ✅ **Test end-to-end** (15 min) - Verify pipeline works

**Total**: ~1.5 hours to working pipeline

---

## 🎯 **AFTER PIPELINE WORKS**

1. Filter 50 candidates → Get top 5-10
2. Scrape top 10 with Diffbot
3. Generate eligibility assessments
4. Generate 5-10 dossiers
5. Submit to Zo for review

---

## 📁 **FILES TO CREATE**

```
oncology-coPilot/oncology-backend-minimal/
└── api/
    └── services/
        └── client_dossier/
            ├── __init__.py
            ├── trial_querier.py      # Read from SQLite
            ├── trial_scraper.py      # Diffbot scraping
            ├── trial_filter.py       # Filtering logic
            └── dossier_generator.py # Generate dossiers (later)
```

---

## ⚔️ **START BUILDING NOW!**

**First Command**:
```bash
cd oncology-coPilot/oncology-backend-minimal
python -c "import sqlite3; conn = sqlite3.connect('data/clinical_trials.db'); print(f'Trials: {conn.execute(\"SELECT COUNT(*) FROM clinical_trials\").fetchone()[0]}')"
```

**Expected Output**: `Trials: 1000` ✅

**Then**: Start building the pipeline files above!

---

**Last Updated**: January 14, 2025  
**Status**: 🎯 **READY TO BUILD** - All data located, all tools ready

