import asyncio
import json
import logging
from typing import Dict, Any

# --- DOCTRINE: DYNAMIC PATH CORRECTION ---
import os
import sys
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from src.tools.search_engine import search
from src.tools.web_scraper import run_scrape
from src.tools.llm_api import query_llm
from src.services.command_center.schemas import IntelligenceDossier

logger = logging.getLogger(__name__)

class IntelligenceAnalyst:
    """A dedicated tool for generating high-fidelity intelligence dossiers."""

    async def _conduct_reconnaissance(self, query: str) -> str:
        """Performs web search and scraping to gather raw intelligence."""
        raw_intelligence = ""
        try:
            logger.info(f" analyst conducting web search for '{query}'...")
            search_results = await asyncio.to_thread(search, query)
            
            urls_to_scrape = [result['URL'] for result in search_results[:2]]
            logger.info(f" analyst scrapping URLs: {urls_to_scrape}")
            
            scraped_content_list = await asyncio.wait_for(run_scrape(urls_to_scrape), timeout=45.0)

            for item in scraped_content_list:
                raw_intelligence += f"\\n\\n--- START OF SOURCE: {item['url']} ---\\n{item['content']}\\n--- END OF SOURCE ---\\n\\n"
            
            if not raw_intelligence:
                raw_intelligence = "No content gathered from web reconnaissance."
                logger.warning("Web reconnaissance yielded no content.")
            return raw_intelligence

        except asyncio.TimeoutError:
            logger.error("Web scraping timed out.")
            return "Web scraping timed out. The primary analyst will have to proceed with its internal knowledge base."
        except Exception as e:
            logger.error(f"An error occurred during web reconnaissance: {e}")
            return f"An error occurred during reconnaissance: {e}. The primary analyst will proceed with its internal knowledge base."

    async def _run_synthesis(self, raw_intelligence: str, gene_symbol: str, mutation_details: str) -> Dict[str, Any]:
        """Runs the LLM synthesis process and parses the response."""
        synthesis_prompt = f\"\"\"
        **Role:** You are a senior biomedical researcher and intelligence analyst.

        **Source Material:**
        {raw_intelligence}

        **Task:**
        Based *only* on the source material provided above, analyze the user's query and generate a comprehensive intelligence dossier.
        User Query: "Provide a detailed analysis of the biological function and role in cancer of {gene_symbol} with mutation {mutation_details}."

        **Output Format (JSON ONLY):**
        You MUST respond with ONLY a valid JSON object matching this structure. The 'relevance' field for each entity is mandatory and must be a float between 0.0 and 1.0.
        {{
          "summary": "<Your 2-3 sentence summary of the most critical findings based on the source material.>",
          "entities": [
            {{"name": "<Entity Name>", "type": "<Gene/Protein/Compound>", "relevance": 1.0}}
          ],
          "mechanisms": [
            "<Detailed mechanism 1 based on the source material.>",
            "<Detailed mechanism 2 based on the source material.>"
          ],
          "conclusions": [
            "<Actionable conclusion 1 based on the source material.>",
            "<Actionable conclusion 2 based on the source material.>"
          ]
        }}
        \"\"\"
        try:
            logger.info("Engaging primary AI analyst for synthesis...")
            synthesis_result_str = await asyncio.to_thread(
                query_llm,
                prompt=synthesis_prompt,
                provider="gemini"
            )
            
            clean_response = synthesis_result_str.strip().replace("```json", "").replace("```", "").strip()
            logger.info(f"Primary analyst response: \\n{clean_response}")
            return json.loads(clean_response)

        except Exception as e:
            logger.error(f"Failed to generate high-fidelity dossier: {e}. Falling back to mock data.")
            return {
                "summary": "Primary analyst failed to generate dossier. This is mock data.",
                "entities": [],
                "mechanisms": [],
                "conclusions": []
            }

    async def generate_dossier(self, gene_symbol: str, mutation_details: str) -> IntelligenceDossier:
        """The main entry point for the tool."""
        recon_query = f"{gene_symbol} {mutation_details} role in cancer metastasis and therapeutic targeting"
        
        # Step 1: Reconnaissance
        raw_intelligence = await self._conduct_reconnaissance(recon_query)
        
        # Step 2: Synthesis
        synthesis_data = await self._run_synthesis(raw_intelligence, gene_symbol, mutation_details)
        
        # Step 3: Pydantic Validation
        return IntelligenceDossier(**synthesis_data) 