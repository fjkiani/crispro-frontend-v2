import os
import requests
import sqlite3
from pathlib import Path
from dotenv import load_dotenv
from typing import Optional
import pandas as pd
import argparse
import json

# --- Configuration ---
load_dotenv()
DB_PATH = Path("data/databases/threat_matrix.db")

# List of key driver genes for our use case
TARGET_GENES = ["ASXL1", "RUNX1", "ABL1", "BRAF"]

class CosmicApiImporter:
    def __init__(self):
        DB_PATH.parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(DB_PATH)
        self.cursor = self.conn.cursor()
        self.base_url = "https://clinicaltables.nlm.nih.gov/api/cosmic/v3/search"

    def run_import(self, genes: list):
        print("🚀 Starting COSMIC API import process...")
        self._setup_database()

        for gene in genes:
            print(f"🧬 Processing gene: {gene}")
            api_data = self._fetch_data_for_gene(gene)
            if not api_data:
                print(f"   - ⚠️ No data returned from API for {gene}. Skipping.")
                continue

            self._parse_and_insert_gene_data(gene, api_data)
            self._parse_and_insert_variant_data(gene, api_data)

        self.conn.close()
        print("✅ Import process complete.")

    def _fetch_data_for_gene(self, gene_symbol: str) -> Optional[list]:
        params = {
            "terms": gene_symbol,
            "maxList": 500,
            "ef": "MutationDescription,MutationZygosity,LOH,GRCh,MutationCDS,MutationAA,MutationGenomePosition,MutationStrand,SNP,ResistanceMutation,FATHMMprediction,FATHMMscore,CGCgene,PrimarySite,SiteSubtype1,SiteSubtype2,SiteSubtype3,PrimaryHistology,HistologySubtype1,HistologySubtype2,HistologySubtype3,PubmedPMID,ID_STUDY,SampleType"
        }
        try:
            print(f"   - 📞 Calling COSMIC API for {gene_symbol}...")
            response = requests.get(self.base_url, params=params)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            print(f"   - ❌ API request failed for {gene_symbol}: {e}")
            return None
        except requests.exceptions.JSONDecodeError:
            print(f"   - ❌ Failed to decode JSON response for {gene_symbol}.")
            return None

    def _parse_and_insert_gene_data(self, gene_symbol: str, api_data: list):
        total_mutations = api_data[0]
        self.cursor.execute("INSERT OR REPLACE INTO genes (gene_symbol, total_mutations) VALUES (?, ?)", (gene_symbol, total_mutations))
        self.conn.commit()
        print(f"     - 📝 Upserted gene record for {gene_symbol} with {total_mutations} total mutations.")

    def _parse_and_insert_variant_data(self, gene_symbol: str, api_data: list):
        if len(api_data) < 4:
            print("     - ⚠️ Variant data not found in the expected format. Skipping variant insertion.")
            return

        mutation_ids = api_data[1]
        extra_fields = api_data[2]
        display_fields = api_data[3]

        variants_to_insert = []
        for i, mut_id in enumerate(mutation_ids):
            display_row = next((row for row in display_fields if row[0] == mut_id), None)
            if not display_row:
                continue

            pmid_raw = extra_fields.get("PubmedPMID", [])[i]
            pmid = int(pmid_raw) if str(pmid_raw).isdigit() else None

            variant = {
                "mutation_id": mut_id,
                "gene_symbol": gene_symbol,
                "mutation_cds": display_row[2],
                "mutation_aa": display_row[3],
                "mutation_description": extra_fields.get("MutationDescription", [])[i],
                "genome_position": extra_fields.get("MutationGenomePosition", [])[i],
                "primary_site": extra_fields.get("PrimarySite", [])[i],
                "primary_histology": extra_fields.get("PrimaryHistology", [])[i],
                "pubmed_pmid": pmid
            }
            variants_to_insert.append(tuple(variant.values()))

        if variants_to_insert:
            self.cursor.executemany("""
                INSERT OR REPLACE INTO variants 
                (mutation_id, gene_symbol, mutation_cds, mutation_aa, mutation_description, genome_position, primary_site, primary_histology, pubmed_pmid)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, variants_to_insert)
            self.conn.commit()
            print(f"     - ➕ Inserted {len(variants_to_insert)} variant records.")

    def _setup_database(self):
        print("  - Setting up database schema...")
        self.cursor.execute("DROP TABLE IF EXISTS genes")
        self.cursor.execute("DROP TABLE IF EXISTS variants")

        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS genes (
                gene_symbol TEXT PRIMARY KEY,
                total_mutations INTEGER
            )
        """)

        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS variants (
                mutation_id TEXT PRIMARY KEY,
                gene_symbol TEXT,
                mutation_cds TEXT,
                mutation_aa TEXT,
                mutation_description TEXT,
                genome_position TEXT,
                primary_site TEXT,
                primary_histology TEXT,
                pubmed_pmid INTEGER,
                FOREIGN KEY (gene_symbol) REFERENCES genes (gene_symbol)
            )
        """)
        self.conn.commit()
        print("✅ Database setup complete.")

    def get_tissue_distribution(self, gene_symbol: str):
        """
        Queries the COSMIC API to get the tissue distribution of mutations for a specific gene.
        """
        print(f"🔬 Querying COSMIC for tissue distribution of '{gene_symbol}' mutations...")
        
        params = {
            "terms": gene_symbol,
            "fields": "sample.tumor_site"
        }

        try:
            response = requests.get(self.base_url, params=params)
            response.raise_for_status()
            
            # The API returns a list where the first element is a list of headers,
            # and subsequent elements are the data rows.
            data = response.json()
            print("--- RAW API RESPONSE ---")
            print(json.dumps(data, indent=4))
            print("------------------------")
            
            if len(data) < 3: # Expecting [count, [headers], [rows...]]
                print(f"   - ⚠️ No tissue data found for {gene_symbol}. API response format may have changed.")
                return

            total_count = data[0]
            headers = data[1]
            rows = data[2:]
            
            print(f"   - ℹ️  Found {total_count} total mutation records.")

            # Use pandas for easy analysis and aggregation
            if len(headers) == 1:
                # If only one column is requested, the data is a flat list.
                # We need to wrap it in a list of lists for pandas.
                df = pd.DataFrame([[row] for row in rows], columns=headers)
            else:
                df = pd.DataFrame(rows, columns=headers)
            
            if 'sample.tumor_site' not in df.columns:
                print(f"   - ❌ 'sample.tumor_site' field not found in API response for {gene_symbol}.")
                return

            print(f"\n--- Tissue Distribution for {gene_symbol} ---")
            tissue_counts = df['sample.tumor_site'].value_counts()
            print(tissue_counts.to_string())
            print("------------------------------------------\n")

        except requests.RequestException as e:
            print(f"   - ❌ API request failed for {gene_symbol}: {e}")
        except (ValueError, KeyError) as e:
            print(f"   - ❌ Failed to parse API response for {gene_symbol}: {e}")

    def close_connection(self):
        self.conn.close()

# --- Main Execution ---
if __name__ == "__main__":
    importer = CosmicApiImporter()
    importer.run_import(['ASXL1', 'RUNX1'])
    importer.close_connection() 