import sqlite3
import pandas as pd
import os

db_path = 'data/clinical_trials.db'
if not os.path.exists(db_path):
    # Try looking in oncology-backend-minimal context
    db_path = 'oncology-coPilot/oncology-backend-minimal/data/clinical_trials.db'

print(f"Connecting to {db_path}...")
conn = sqlite3.connect(db_path)

nct_ids = ['NCT04284969', 'NCT02655016']
query = f"SELECT id, status, conditions, title, phases FROM trials WHERE id IN {tuple(nct_ids)}"

try:
    df = pd.read_sql_query(query, conn)
    print(df.to_string())
except Exception as e:
    print(f"Error: {e}")

conn.close()
