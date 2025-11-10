# 🖥️ EXECUTION COMMANDS

## **Setup Commands**

### **Create Folder Structure**
```bash
cd oncology-coPilot/oncology-backend
mkdir -p scripts/agent_1_seeding/{api,parsers,database,utils}
mkdir -p tests/agent_1_seeding
touch scripts/agent_1_seeding/{__init__,main,config}.py
touch scripts/agent_1_seeding/api/{__init__,ctgov_client}.py
touch scripts/agent_1_seeding/parsers/{__init__,study_parser,biomarker_extractor,locations_parser}.py
touch scripts/agent_1_seeding/database/{__init__,migration,sqlite_client,chromadb_client}.py
touch scripts/agent_1_seeding/utils/{__init__,error_handler,logger}.py
touch scripts/migrate_schema_v2.sql
touch tests/agent_1_seeding/{__init__,test_api_client,test_parsers,test_database,test_integration,conftest}.py
```

---

## **Pre-Flight Commands**

### **Backup Database**
```bash
cp oncology-coPilot/oncology-backend/backend/data/clinical_trials.db \
   oncology-coPilot/oncology-backend/backend/data/clinical_trials_BACKUP_$(date +%Y%m%d).db
```

### **Check Database Status**
```bash
sqlite3 oncology-coPilot/oncology-backend/backend/data/clinical_trials.db \
  "SELECT COUNT(*) FROM clinical_trials;"
```

### **Check ChromaDB Status**
```bash
python -c "
import chromadb
c = chromadb.PersistentClient(path='oncology-coPilot/oncology-backend/backend/data/chroma_data')
print(c.get_collection('clinical_trials_eligibility').count())
"
```

### **Test API Connectivity**
```bash
curl "https://clinicaltrials.gov/api/v2/studies?query.cond=ovarian+cancer&pageSize=1&format=json"
```

---

## **Execution Commands**

### **Test Mode (100 trials, no embeddings)**
```bash
cd oncology-coPilot/oncology-backend
python -m scripts.agent_1_seeding.main --limit 100 --test-mode
```

### **Full Execution (1000 trials)**
```bash
cd oncology-coPilot/oncology-backend
python -m scripts.agent_1_seeding.main --limit 1000
```

### **Skip Migration (if already run)**
```bash
python -m scripts.agent_1_seeding.main --limit 1000 --skip-migration
```

### **Skip Embeddings (faster testing)**
```bash
python -m scripts.agent_1_seeding.main --limit 1000 --skip-embeddings
```

---

## **Verification Commands**

### **Count Ovarian Trials**
```bash
sqlite3 oncology-coPilot/oncology-backend/backend/data/clinical_trials.db \
  "SELECT COUNT(*) FROM clinical_trials WHERE disease_subcategory='ovarian_cancer';"
```

### **Verify Disease Tags**
```bash
sqlite3 oncology-coPilot/oncology-backend/backend/data/clinical_trials.db \
  "SELECT COUNT(*) FROM clinical_trials WHERE disease_category='gynecologic_oncology' AND disease_subcategory='ovarian_cancer';"
```

### **Check Biomarker Coverage**
```bash
sqlite3 oncology-coPilot/oncology-backend/backend/data/clinical_trials.db \
  "SELECT COUNT(*) FROM clinical_trials WHERE biomarker_requirements IS NOT NULL AND disease_subcategory='ovarian_cancer';"
```

### **Check Locations Coverage**
```bash
sqlite3 oncology-coPilot/oncology-backend/backend/data/clinical_trials.db \
  "SELECT COUNT(*) FROM clinical_trials WHERE locations_data IS NOT NULL AND disease_subcategory='ovarian_cancer';"
```

### **Check NY Locations**
```bash
sqlite3 oncology-coPilot/oncology-backend/backend/data/clinical_trials.db \
  "SELECT COUNT(*) FROM clinical_trials WHERE locations_data LIKE '%\"state\": \"NY\"%' AND disease_subcategory='ovarian_cancer';"
```

---

## **Testing Commands**

### **Run All Tests**
```bash
cd oncology-coPilot/oncology-backend
PYTHONPATH=. venv/bin/pytest tests/agent_1_seeding/ -v
```

### **Run Specific Test Module**
```bash
PYTHONPATH=. venv/bin/pytest tests/agent_1_seeding/test_api_client.py -v
```

### **Run With Test Mode Flag**
```bash
PYTHONPATH=. venv/bin/pytest tests/agent_1_seeding/ -v --test-mode
```

---

## **Status Update Command**

```bash
echo "- [X] Agent 1: Seeding (Status: COMPLETE)" >> \
  .cursor/rules/clinical_trials_agents/MASTER_STATUS.md
```

