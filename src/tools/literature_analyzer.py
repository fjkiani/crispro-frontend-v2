# \"\"\"
# Literature Intelligence Analyzer (Diffbot Edition)

# This tool enriches the "Threat Matrix" database by analyzing scientific
# literature associated with variants. It uses the Diffbot Article API to ensure
# clean, reliable text extraction for high-quality LLM summarization.

# Workflow:
# 1.  Reads all unique PubMed IDs (PMIDs) from the `variants` table.
# 2.  Converts each PMID to a PubMed Central (PMC) URL.
# 3.  Calls the Diffbot Article API to extract the clean text of the article.
# 4.  Sends the text to an LLM for a deep, strategic summary.
# 5.  Stores the summary and source URL in the `literature_intelligence` table.
# \"\"\"
import os
import sqlite3
import requests
from pathlib import Path
from dotenv import load_dotenv
from llm_api import query_llm
import json
from sentence_transformers import SentenceTransformer
import textwrap

# --- Configuration ---
load_dotenv()
DB_PATH = Path("data/databases/threat_matrix.db")
DIFFBOT_TOKEN = os.getenv("DIFFBOT_TOKEN")
DIFFBOT_API_URL = "https://api.diffbot.com/v3/article"
ID_CONVERTER_API_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"

class LiteratureAnalyzer:
    def __init__(self):
        if not DB_PATH.exists():
            raise FileNotFoundError(f"Database not found at {DB_PATH}. Please run cosmic_importer.py first.")
        if not DIFFBOT_TOKEN:
            raise ValueError("DIFFBOT_TOKEN not found in environment variables. Please set it in your .env file.")
        self.conn = sqlite3.connect(DB_PATH)
        self.conn.enable_load_extension(True)
        self.cursor = self.conn.cursor()
        self._setup_database()
        print("🧠 Loading sentence-transformer model...")
        self.embedding_model = SentenceTransformer('all-MiniLM-L6-v2')
        print("✅ Model loaded.")

    def run_ingestion(self):
        print("🚀 Starting Literature Intelligence Ingestion (Phase 1: Metadata & Embeddings)...")
        pmids = self._get_unique_pmids()
        print(f"  - Found {len(pmids)} unique PubMed IDs to process.")

        for pmid in pmids[:5]: # Let's process a few for demonstration
            if self._is_already_processed(pmid):
                print(f"  - ⏭️ PMID: {pmid} already processed. Skipping.")
                continue
            
            print(f"📄 Processing PMID: {pmid}")
            article_url = self._convert_pmid_to_pmc_url(pmid)
            if not article_url:
                print(f"   - ❌ Could not find a full-text PMC URL for PMID {pmid}. Skipping.")
                self._mark_as_processed(pmid, success=False)
                continue
            
            print(f"   - ✅ Found PMC URL: {article_url}")
            article_data = self._get_article_data_from_diffbot(article_url)
            if not article_data:
                print(f"   - ❌ Diffbot failed to extract data from {article_url}. Skipping.")
                self._mark_as_processed(pmid, success=False)
                continue

            self._store_structured_intelligence(pmid, article_url, article_data)
            self._generate_and_store_embeddings(pmid, article_data.get("text", ""))
            print(f"   - ✅ Successfully processed and stored structured data for PMID {pmid}.")

        # self.conn.close()
        print("✅ Literature ingestion complete.")

    def run_strategic_summary(self, pmid: int, prompt_template: str, provider: str = "gemini"):
        """
        Runs a targeted, strategic analysis on a single pre-ingested article.
        """
        print(f"🕵️  Running strategic summary for PMID: {pmid}")
        self.cursor.execute("SELECT full_text, title, author FROM literature_intelligence WHERE pubmed_id = ? AND processed_successfully = TRUE", (pmid,))
        result = self.cursor.fetchone()

        if not result:
            print(f"  - ❌ Could not find pre-processed data for PMID {pmid}. Run ingestion first.")
            return

        full_text, title, author = result
        if not full_text:
            print(f"  - ❌ No full text available for PMID {pmid}.")
            return

        print(f"  - Title: {title}")

        prompt = prompt_template.format(
            pmid=pmid,
            title=title,
            full_text=full_text[:25000] # Stay within reasonable token limits
        )

        print("  - 🧠 Sending request to LLM for targeted analysis...")
        summary = query_llm(prompt, provider=provider)

        self._store_strategic_summary(pmid, prompt_template, summary, provider)
        print(f"  - ✅ Stored strategic summary for PMID {pmid}.")

    def force_reingest(self, pmid: int):
        """Forces the re-ingestion of a specific PMID by deleting its existing records."""
        print(f"🔥 Forcing re-ingestion for PMID: {pmid}")
        self.cursor.execute("DELETE FROM strategic_summaries WHERE pubmed_id = ?", (pmid,))
        self.cursor.execute("DELETE FROM article_tags WHERE pubmed_id = ?", (pmid,))
        self.cursor.execute("DELETE FROM literature_intelligence WHERE pubmed_id = ?", (pmid,))
        self.conn.commit()
        print(f"  - ✅ Cleared all previous records for PMID: {pmid}")

    def close_connection(self):
        """Closes the database connection."""
        if self.conn:
            self.conn.close()
            print("🔒 Database connection closed.")

    def _convert_pmid_to_pmc_url(self, pmid: int) -> str | None:
        params = {"ids": pmid, "format": "json"}
        try:
            response = requests.get(ID_CONVERTER_API_URL, params=params)
            response.raise_for_status()
            data = response.json()
            if data.get("status") == "ok" and data.get("records"):
                record = data["records"][0]
                if "pmcid" in record:
                    return f"https://www.ncbi.nlm.nih.gov/pmc/articles/{record['pmcid']}/"
        except requests.RequestException as e:
            print(f"   - ⚠️ ID Converter API request failed for PMID {pmid}: {e}")
        return None

    def _get_article_data_from_diffbot(self, article_url: str) -> dict | None:
        params = {"token": DIFFBOT_TOKEN, "url": article_url, "fields": "title,author,date,siteName,tags,images,html,text"}
        try:
            response = requests.get(DIFFBOT_API_URL, params=params, timeout=60)
            response.raise_for_status()
            data = response.json()
            if "objects" in data and data["objects"]:
                return data["objects"][0]
        except requests.RequestException as e:
            print(f"   - ⚠️ Diffbot API request failed for URL {article_url}: {e}")
        return None

    def _store_structured_intelligence(self, pmid: int, url: str, data: dict):
        self.cursor.execute(
            """
            INSERT OR REPLACE INTO literature_intelligence (
                pubmed_id, full_text_url, title, author, publication_date, 
                site_name, full_text, full_html, diffbot_response_json, processed_successfully
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                pmid, url, data.get("title"), data.get("author"), data.get("date"),
                data.get("siteName"), data.get("text"), data.get("html"),
                json.dumps(data), True
            )
        )

        # Clear existing tags for this pmid before inserting new ones
        # self.cursor.execute("DELETE FROM article_tags WHERE pubmed_id = ?", (pmid,))
        
        tags = data.get("tags", [])
        for tag in tags:
            self.cursor.execute(
                """
                INSERT INTO article_tags (pubmed_id, label, uri, score, sentiment) 
                VALUES (?, ?, ?, ?, ?)
                """,
                (pmid, tag.get("label"), tag.get("uri"), tag.get("score"), tag.get("sentiment"))
            )
        
        self.conn.commit()

    def _mark_as_processed(self, pmid: int, success: bool):
        self.cursor.execute(
            "INSERT OR REPLACE INTO literature_intelligence (pubmed_id, processed_successfully) VALUES (?, ?)",
            (pmid, success)
        )
        self.conn.commit()

    def _is_already_processed(self, pmid: int) -> bool:
        self.cursor.execute("SELECT processed_successfully FROM literature_intelligence WHERE pubmed_id = ?", (pmid,))
        result = self.cursor.fetchone()
        return result is not None

    def _store_strategic_summary(self, pmid: int, prompt: str, summary: str, provider: str):
        self.cursor.execute(
            """
            INSERT INTO strategic_summaries (pubmed_id, prompt, summary, model_used)
            VALUES (?, ?, ?, ?)
            """,
            (pmid, prompt, summary, provider)
        )
        self.conn.commit()

    def _generate_and_store_embeddings(self, pmid: int, text: str):
        if not text:
            return
            
        print(f"   - 🔮 Generating embeddings for PMID: {pmid}")
        
        # More robust chunking strategy
        chunks = textwrap.wrap(text, width=1500, break_long_words=False, replace_whitespace=False)
        
        if not chunks:
            print("   - ⚠️ No suitable text chunks found for embedding.")
            return

        embeddings = self.embedding_model.encode(chunks, convert_to_tensor=False)
        
        self.cursor.execute("DELETE FROM vector_store WHERE pubmed_id = ?", (pmid,))
        
        for i, chunk in enumerate(chunks):
            embedding_blob = json.dumps(embeddings[i].tolist())
            self.cursor.execute(
                "INSERT INTO vector_store (pubmed_id, chunk_id, chunk_text, embedding) VALUES (?, ?, ?, ?)",
                (pmid, i, chunk, embedding_blob)
            )
        self.conn.commit()
        print(f"   - ✅ Stored {len(chunks)} text chunks and embeddings.")

    def _setup_database(self):
        # self.cursor.execute("DROP TABLE IF EXISTS literature_intelligence")
        # self.cursor.execute("DROP TABLE IF EXISTS article_tags")
        # self.cursor.execute("DROP TABLE IF EXISTS strategic_summaries")
        # self.cursor.execute("DROP TABLE IF EXISTS vector_store")

        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS literature_intelligence (
                pubmed_id INTEGER PRIMARY KEY,
                full_text_url TEXT,
                title TEXT,
                author TEXT,
                publication_date TEXT,
                site_name TEXT,
                summary TEXT, -- For future use in Phase 2
                full_text TEXT,
                full_html TEXT,
                diffbot_response_json TEXT,
                processed_successfully BOOLEAN DEFAULT FALSE,
                ingested_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        
        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS article_tags (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                pubmed_id INTEGER,
                label TEXT,
                uri TEXT,
                score REAL,
                sentiment REAL,
                FOREIGN KEY (pubmed_id) REFERENCES literature_intelligence (pubmed_id)
            )
        """)
        
        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS strategic_summaries (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                pubmed_id INTEGER,
                prompt TEXT,
                summary TEXT,
                model_used TEXT,
                analyzed_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (pubmed_id) REFERENCES literature_intelligence (pubmed_id)
            )
        """)

        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS vector_store (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                pubmed_id INTEGER,
                chunk_id INTEGER,
                chunk_text TEXT,
                embedding TEXT, -- Storing as JSON string
                FOREIGN KEY (pubmed_id) REFERENCES literature_intelligence (pubmed_id)
            )
        """)
        self.conn.commit()
    
    def _get_unique_pmids(self) -> list[int]:
        self.cursor.execute("SELECT DISTINCT pubmed_pmid FROM variants WHERE pubmed_pmid IS NOT NULL")
        return [row[0] for row in self.cursor.fetchall()]

if __name__ == "__main__":
    analyzer = LiteratureAnalyzer()
    try:
        target_pmid = 21828135

        # Force a clean re-ingestion for our target to ensure we have the full text
        analyzer.force_reingest(target_pmid)

        # == PHASE 1: Run once to ingest metadata from all new papers ==
        print("--- Running Phase 1: Metadata Ingestion ---")
        analyzer.run_ingestion() 
        print("--- Phase 1 Complete ---\n")

        # == PHASE 2: Run targeted analysis on specific papers of interest ==
        print("--- Running Phase 2: Strategic Analysis ---")
        # Using one of the PMIDs we know was successfully ingested in the last run.
        
        # This is a much more focused prompt than a generic summary.
        strategic_prompt = (
            "This article is being analyzed for its relevance to leukemogenesis, specifically concerning the genes ASXL1 and RUNX1, and its potential metastatic sites (or 'soils'). "
            "Based on the article text provided, answer the following questions with direct quotes and citations where possible:\n\n"
            "1. What is the role of epigenetic dysregulation in Chronic Myelomonocytic Leukemia (CMML)?\n"
            "2. Does this paper mention ASXL1 mutations? If so, what is their frequency and what is their presumed functional consequence (e.g., loss-of-function, truncating)?\n"
            "3. Does the paper discuss the co-occurrence of ASXL1 mutations with other known driver mutations like TET2 or CBL?\n"
            "4. **Crucially, what tissues or organs (e.g., spleen, liver, CNS, skin, lymph nodes) are mentioned as sites of extramedullary disease or leukemic infiltration? List them directly.**\n"
            "5. Is there any mention of RUNX1 in the context of CMML in this paper?\n\n"
            "Article Title: {title}\n"
            "Article Text:\n---\n{full_text}"
        )

        analyzer.run_strategic_summary(target_pmid, strategic_prompt, provider="gemini")
        
        print("\n✅ Strategic analysis demonstration complete.")
    finally:
        analyzer.close_connection() 