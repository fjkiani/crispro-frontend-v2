# ✅ EXECUTION CHECKLIST

## **Pre-Flight Checks**

- [ ] `GOOGLE_API_KEY` set in `.env`
- [ ] SQLite database exists at `backend/data/clinical_trials.db`
- [ ] ChromaDB directory writable
- [ ] Python dependencies installed: `requests`, `chromadb`, `pytest`
- [ ] **Backup database:** `cp backend/data/clinical_trials.db backend/data/clinical_trials_BACKUP_$(date +%Y%m%d).db`

## **Database Status Check**

```bash
# Check current trial count
sqlite3 oncology-coPilot/oncology-backend/backend/data/clinical_trials.db "SELECT COUNT(*) FROM clinical_trials;"

# Check ChromaDB collection status
python -c "import chromadb; c = chromadb.PersistentClient(path='oncology-coPilot/oncology-backend/backend/data/chroma_data'); print(c.get_collection('clinical_trials_eligibility').count())"
```

## **Environment Check**

```bash
# Verify Google API key
echo $GOOGLE_API_KEY

# Test ClinicalTrials.gov API v2 connectivity
curl "https://clinicaltrials.gov/api/v2/studies?query.cond=ovarian+cancer&pageSize=1&format=json"
```

---

## **Execution Steps**

### **1. Test Mode First (100 trials)**

```bash
cd oncology-coPilot/oncology-backend
python -m scripts.agent_1_seeding.main --limit 100 --test-mode

# Verify 100 trials inserted
sqlite3 backend/data/clinical_trials.db "SELECT COUNT(*) FROM clinical_trials WHERE disease_subcategory='ovarian_cancer';"
```

### **2. Run Tests**

```bash
cd oncology-coPilot/oncology-backend
PYTHONPATH=. venv/bin/pytest tests/agent_1_seeding/ -v
```

### **3. Full Execution (1000 trials)**

```bash
python -m scripts.agent_1_seeding.main --limit 1000

# Verify results
sqlite3 backend/data/clinical_trials.db "SELECT COUNT(*) FROM clinical_trials WHERE disease_subcategory='ovarian_cancer';"
```

---

## **Verification**

- [ ] 1000+ trials in SQLite
- [ ] All trials have `disease_category = "gynecologic_oncology"`
- [ ] All trials have `disease_subcategory = "ovarian_cancer"`
- [ ] ChromaDB has 1000+ embeddings (if not skipped)
- [ ] Summary report shows expected stats

---

## **Update Status**

```bash
echo "- [X] Agent 1: Seeding (Status: COMPLETE)" >> .cursor/rules/clinical_trials_agents/MASTER_STATUS.md
```

