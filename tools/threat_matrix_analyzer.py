import sqlite3
import pandas as pd
from pathlib import Path
import argparse

DB_PATH = Path("data/databases/threat_matrix.db")

class ThreatMatrixAnalyzer:
    """
    Analyzes the local Threat Matrix database to extract strategic insights.
    """
    def __init__(self):
        if not DB_PATH.exists():
            raise FileNotFoundError(f"Database not found at {DB_PATH}. Please run cosmic_importer.py first.")
        self.conn = sqlite3.connect(DB_PATH)
        print("✅ Connected to Threat Matrix database.")

    def get_tissue_distribution_for_gene(self, gene_symbol: str):
        """
        Calculates and displays the tissue distribution of mutations for a given gene
        by querying the local database.
        """
        print(f"\n🔬 Analyzing tissue distribution for gene: {gene_symbol}")
        try:
            query = "SELECT primary_site FROM variants WHERE gene_symbol = ? AND primary_site IS NOT NULL"
            df = pd.read_sql_query(query, self.conn, params=(gene_symbol,))

            if df.empty:
                print(f"   - ⚠️ No variants with primary site information found for {gene_symbol} in the database.")
                return

            print(f"--- Tissue Distribution for {gene_symbol} (from local DB) ---")
            tissue_counts = df['primary_site'].value_counts()
            
            # Filter for top 15 to keep it clean
            print(tissue_counts.head(15).to_string())
            print("----------------------------------------------------------\n")

        except Exception as e:
            print(f"   - ❌ An error occurred during database analysis: {e}")

    def get_histology_distribution(self, gene_symbol: str, tissue_site: str):
        """
        Calculates and displays the histology distribution for a given gene within a specific primary tissue site.
        """
        print(f"\n🔬 Analyzing histology distribution for {gene_symbol} in tissue: {tissue_site}")
        try:
            query = """
                SELECT primary_histology 
                FROM variants 
                WHERE gene_symbol = ? AND primary_site = ? AND primary_histology IS NOT NULL
            """
            df = pd.read_sql_query(query, self.conn, params=(gene_symbol, tissue_site))

            if df.empty:
                print(f"   - ⚠️ No variants with histology information found for {gene_symbol} in {tissue_site}.")
                return

            print(f"--- Histology Distribution for {gene_symbol} in {tissue_site} ---")
            histology_counts = df['primary_histology'].value_counts()
            
            print(histology_counts.head(15).to_string())
            print("----------------------------------------------------------------" + "-" * len(gene_symbol) + "-" * len(tissue_site) + "\n")

        except Exception as e:
            print(f"   - ❌ An error occurred during database analysis: {e}")

    def close_connection(self):
        self.conn.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze the Threat Matrix database.")
    parser.add_argument(
        "--tissue-distribution",
        metavar="GENE",
        help="Get the tissue distribution of mutations for a specific gene (e.g., ASXL1)."
    )
    parser.add_argument(
        "--histology-distribution",
        nargs=2,
        metavar=("GENE", "TISSUE"),
        help="Get histology distribution for a gene within a primary tissue site (e.g., ASXL1 haematopoietic_and_lymphoid_tissue)."
    )
    args = parser.parse_args()

    analyzer = ThreatMatrixAnalyzer()
    if args.tissue_distribution:
        analyzer.get_tissue_distribution_for_gene(args.tissue_distribution)
    elif args.histology_distribution:
        gene, tissue = args.histology_distribution
        analyzer.get_histology_distribution(gene, tissue)
    else:
        parser.print_help()

    analyzer.close_connection() 