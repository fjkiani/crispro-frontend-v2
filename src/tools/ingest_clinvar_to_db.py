import pandas as pd
import sqlite3
import gzip
import os

def ingest_clinvar_to_db(file_path, db_path):
    """
    Parses the ClinVar variant_summary.txt.gz file and ingests it into a SQLite database.

    Args:
        file_path (str): The path to the variant_summary.txt.gz file.
        db_path (str): The path to the SQLite database file.
    """
    print(f"Starting ingestion of {file_path} into {db_path}...")

    # Define chunk size for memory-efficient processing
    chunk_size = 100000

    # Connect to SQLite database
    conn = sqlite3.connect(db_path)
    
    # Read the gzipped file in chunks and write to the database
    with gzip.open(file_path, 'rt') as f:
        # Read the header to get column names
        header = f.readline().strip().split('\t')
        
        # Create a pandas DataFrame iterator
        chunks = pd.read_csv(f, sep='\t', chunksize=chunk_size, names=header, low_memory=False)

        for i, chunk in enumerate(chunks):
            print(f"Processing chunk {i+1}...")
            # Clean column names to be valid SQL identifiers
            chunk.columns = chunk.columns.str.replace('[^A-Za-z0-9_]+', '', regex=True)
            chunk.to_sql('variants', conn, if_exists='append', index=False)

    print("Ingestion complete.")
    conn.close()

if __name__ == "__main__":
    file_path = "data/variant_summary.txt.gz"
    db_path = "data/databases/clinvar.db"

    # Ensure the database directory exists
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    
    ingest_clinvar_to_db(file_path, db_path) 