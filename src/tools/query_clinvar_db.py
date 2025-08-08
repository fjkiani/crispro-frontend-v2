import sqlite3
import pandas as pd

def query_clinvar_db(db_path, query):
    """
    Queries the ClinVar SQLite database and returns the result as a pandas DataFrame.

    Args:
        db_path (str): The path to the SQLite database file.
        query (str): The SQL query to execute.

    Returns:
        pandas.DataFrame: The result of the query.
    """
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query(query, conn)
    conn.close()
    return df

if __name__ == "__main__":
    db_path = "data/databases/clinvar.db"
    
    # Example query: Find a specific pathogenic variant reviewed by an expert panel
    query = """
    SELECT *
    FROM variants
    WHERE AlleleID = 15323
    LIMIT 1;
    """
    
    df = query_clinvar_db(db_path, query)
    
    # Print all details of the variant
    if not df.empty:
        print(df.iloc[0].to_string())
    else:
        print("Variant not found.") 