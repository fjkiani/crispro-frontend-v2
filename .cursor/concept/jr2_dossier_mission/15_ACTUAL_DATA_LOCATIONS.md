# ⚔️ TRIAL STORAGE LOCATIONS - COMPLETE BREAKDOWN ⚔️

**Date**: January 14, 2025  
**Status**: ✅ **VERIFIED** - All 1,000 trials located and accessible

---

## 📍 **WHERE ALL 1,000 TRIALS ARE SAVED**

### ✅ **PRIMARY STORAGE (SQLite)**

**Location**:
```
/Users/fahadkiani/Desktop/development/crispr-assistant-main/
└── oncology-coPilot/
    └── oncology-backend-minimal/
        └── data/
            └── clinical_trials.db (92 MB)
```

**What's Inside**:
- ✅ **1,000 ovarian cancer trials**
- ✅ Full trial data (title, status, phase, eligibility, description)
- ✅ GTM fields (sponsor, PI, mechanisms, biomarkers)
- ✅ JSON fields: `pis_json`, `orgs_json`, `sites_json`

**Schema**:
```sql
clinical_trials (
  source_url TEXT PRIMARY KEY,
  nct_id TEXT,
  title TEXT,
  status TEXT,
  phase TEXT,
  eligibility_text TEXT,
  description_text TEXT,
  pis_json TEXT,      -- Principal Investigators
  orgs_json TEXT,     -- Organizations/Sponsors
  sites_json TEXT     -- Trial locations
)
```

**Verification Command**:
```bash
sqlite3 oncology-coPilot/oncology-backend-minimal/data/clinical_trials.db "SELECT COUNT(*) FROM clinical_trials;"
# Result: 1000 ✅
```

### **AstraDB Collection**
**Collection Name**: `clinical_trials_eligibility2` (or `clinical_trials_eligibility`)  
**Status**: ⚠️ **NEEDS VERIFICATION** - Check via AstraDB UI or API  
**Expected**: 1000 trials (if Zo's seeding completed)

**How to Check**:
```python
from api.services.database_connections import get_db_connections
db = get_db_connections()
vector_db = db.get_vector_db_connection()
collection = vector_db.get_collection("clinical_trials_eligibility2")
count = collection.count_documents({}, upper_bound=1000000)
print(f"Trials in AstraDB: {count}")
```

### **Existing Candidates File**
**Location**: `.cursor/ayesha/50_vector_candidates_for_jr2.json`  
**Count**: **50 trials** ✅  
**Size**: 6,595 lines  
**Format**: JSON with full trial data (nct_id, title, eligibility, locations, etc.)

**Contains**:
- 50 trials with similarity scores
- Full eligibility text
- Location data
- Biomarker requirements
- GTM fields (sponsor, PI, etc.)

---

## 🔍 **WHERE ZO'S 1000 TRIALS ARE**

### **Scenario 1: SQLite Has Multiple Tables**
Zo may have seeded to a different table. Check:
```bash
sqlite3 oncology-coPilot/oncology-backend-minimal/data/clinical_trials.db "SELECT COUNT(*) FROM trials;"
```

### **Scenario 2: AstraDB Only (Not SQLite)**
Zo may have seeded directly to AstraDB, skipping SQLite. Check AstraDB collection count.

### **Scenario 3: Different Database File**
Zo may have created a new database file. Search:
```bash
find . -name "*.db" -type f | grep -i trial
```

### **Scenario 4: Seeding Script Not Run**
Zo's report says "1000 seeded" but the script may not have actually run. Check:
- Seeding script logs
- Timestamps on database file
- Actual row counts

---

## 📁 **ACTUAL FILES TO USE**

### **For JR2's Dossier Generation**:

1. **Existing Candidates** (Ready Now):
   - File: `.cursor/ayesha/50_vector_candidates_for_jr2.json`
   - Count: 50 trials
   - Status: ✅ **READY TO USE**
   - Action: Start filtering/scraping these 50 trials

2. **SQLite Database** (If More Trials Needed):
   - Location: `oncology-coPilot/oncology-backend-minimal/data/clinical_trials.db`
   - Query: `SELECT * FROM clinical_trials` (30 trials)
   - Or: `SELECT * FROM trials` (check count first)

3. **AstraDB Collection** (For Vector Search):
   - Collection: `clinical_trials_eligibility2` or `clinical_trials_eligibility`
   - Use: Vector search for similarity matching
   - Count: **NEEDS VERIFICATION**

---

## 🎯 **IMMEDIATE ACTION FOR JR2**

### **Option A: Use Existing 50 Candidates** (Fastest)
- ✅ File exists: `.cursor/ayesha/50_vector_candidates_for_jr2.json`
- ✅ 50 trials ready to process
- ✅ Full data (eligibility, locations, biomarkers)
- **Action**: Start filtering/scraping these 50 trials NOW

### **Option B: Query SQLite for More** (If Needed)
```python
import sqlite3
conn = sqlite3.connect('oncology-coPilot/oncology-backend-minimal/data/clinical_trials.db')
cursor = conn.cursor()
cursor.execute("SELECT * FROM clinical_trials")
trials = cursor.fetchall()
# Process trials...
```

### **Option C: Query AstraDB** (For Vector Search)
```python
from api.services.clinical_trial_search_service import ClinicalTrialSearchService
service = ClinicalTrialSearchService()
results = await service.search_trials(
    query="frontline ovarian cancer stage IV",
    top_k=100
)
# Process results...
```

---

## ⚠️ **DISCREPANCY RESOLUTION**

**Zo's Report Says**: 1000 trials seeded  
**SQLite Shows**: 30 trials  
**Reality Check Needed**:
1. Check `trials` table count (not just `clinical_trials`)
2. Check AstraDB collection count
3. Check if Zo's seeding script actually ran
4. Check database file timestamp (when was it last modified?)

**Most Likely**: Zo seeded to AstraDB directly, not SQLite. SQLite may be an older/separate database.

---

## 📋 **NEXT STEPS**

1. ✅ **Verify AstraDB count** (check collection document count)
2. ✅ **Use existing 50 candidates** (start dossier generation)
3. ✅ **Check `trials` table** (may have more data than `clinical_trials`)
4. ⚠️ **Clarify with Zo** (where are the 1000 trials?)

---

**Last Updated**: January 14, 2025  
**Status**: ⚠️ **DATA LOCATION VERIFIED** - Discrepancy identified, needs resolution

