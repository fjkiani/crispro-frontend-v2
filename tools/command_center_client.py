# THIS IS THE V9 CLIENT. IF YOU SEE THE OLD LOGS, THE CACHE IS FUCKED.
import os
import httpx
import json
import asyncio
import logging
import modal
from xml.etree import ElementTree as ET

# --- App Setup ---
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# --- Service URLs ---
# These should match your live Modal deployments.
COMMAND_CENTER_URL = os.environ.get("COMMAND_CENTER_URL", "https://crispro--crispr-assistant-command-center-v3-commandcenter-api.modal.run")
ZETA_FORGE_URL = "https://crispro--evo-service-evoservice-api.modal.run" # Corrected from deployment logs
ZETA_ORACLE_URL = "https://crispro--zetascorer-api.modal.run" # Placeholder - Needs to be deployed
BLAST_SERVICE_URL = "https://crispro--blast-service-human-blastservice-web-app.modal.run" # Corrected to use the live web endpoint
IMMUNOGENICITY_SERVICE_URL = "mock" # Placeholder

class CommandCenterClient:
    """
    An asynchronous client to orchestrate guide design campaigns 
    through the CRISPR Assistant's microservices.
    """

    async def _make_request(self, client, method, url, **kwargs):
        """A helper to make async HTTP requests with error handling."""
        try:
            response = await client.request(method, url, **kwargs)
            response.raise_for_status()
            return response.json()
        except httpx.HTTPStatusError as e:
            logger.error(f"HTTP error for {e.request.url}: {e.response.status_code} - {e.response.text}")
            return {"error": f"HTTP {e.response.status_code}", "details": e.response.text}
        except httpx.RequestError as e:
            logger.error(f"Request error for {e.request.url}: {e}")
            return {"error": "Request failed", "details": str(e)}

    async def forge_guides(self, target_sequence: str, num_guides: int = 10):
        """
        Calls ZetaForge (evo-service) to generate guide RNA candidates.
        """
        async with httpx.AsyncClient(timeout=600.0, follow_redirects=True) as client:
            logger.info(f"Calling ZetaForge to generate {num_guides} guides for sequence.")
            
            # Step 1: Start the generation job
            start_payload = {"target_sequence": target_sequence, "num_guides": num_guides}
            start_url = f"{ZETA_FORGE_URL}/generate"
            start_response = await self._make_request(client, "POST", start_url, json=start_payload)
            
            if "error" in start_response:
                return start_response
            
            job_id = start_response.get("job_id")
            if not job_id:
                return {"error": "Failed to get job_id from ZetaForge"}

            logger.info(f"ZetaForge job initiated with ID: {job_id}. Polling for results...")

            # Step 2: Poll for the results
            status_url = f"{ZETA_FORGE_URL}/status/{job_id}"
            for _ in range(20): # Poll for up to 10 minutes
                await asyncio.sleep(30)
                status_response = await self._make_request(client, "GET", status_url)
                if "error" in status_response:
                    # Log the polling error but continue, as it might be transient
                    logger.warning(f"Polling ZetaForge status failed: {status_response}. Retrying...")
                    continue
                
                status = status_response.get("status")
                logger.info(f"Polling ZetaForge... current status: '{status}'")
                
                if status == "completed":
                    logger.info("ZetaForge job complete.")
                    return status_response.get("guides", [])
                elif status == "failed":
                    logger.error(f"ZetaForge job failed: {status_response.get('error')}")
                    return {"error": "ZetaForge generation failed", "details": status_response.get('error')}
            
            return {"error": "ZetaForge polling timed out."}

    async def score_guides_with_oracle(self, guides: list[str], context_sequence: str):
        """
        Calls ZetaOracle to score guides for on-target efficacy.
        (Placeholder - returns mock scores)
        """
        logger.warning("ZetaOracle client is a placeholder. Returning mock scores.")
        await asyncio.sleep(1) # Simulate network latency
        mock_scores = []
        for i, guide in enumerate(guides):
            mock_scores.append({
                "guide_sequence": guide,
                "zeta_score": 0.85 + (i * 0.01), # Mock score
                "confidence": 0.95
            })
        return mock_scores

    def _parse_blast_xml(self, xml_string: str, query_len: int, mismatch_threshold: int = 2) -> int:
        """
        Parses BLAST XML output to count significant off-target hits.
        A 'significant' hit is one that is not the query itself and has few mismatches.
        """
        if not xml_string:
            return 999 # Return a high number if BLAST fails
        
        root = ET.fromstring(xml_string)
        off_target_count = 0
        
        # We look for 'Hit' elements in the XML
        for hit in root.findall('.//Hit'):
            # Each hit has one or more 'Hsp' (High-scoring Pair)
            for hsp in hit.findall('.//Hsp'):
                # We are interested in the number of mismatches
                mismatch_node = hsp.find('Hsp_mismatch')
                mismatches = int(mismatch_node.text) if mismatch_node is not None else 0
                
                identity_node = hsp.find('Hsp_identity')
                identity = int(identity_node.text) if identity_node is not None else 0
                
                # The on-target hit will have perfect identity and 0 mismatches
                if identity == query_len and mismatches == 0:
                    continue # This is the on-target sequence, ignore it
                
                if mismatches <= mismatch_threshold:
                    off_target_count += 1
                    
        return off_target_count

    async def check_off_targets_with_blast(self, guides: list[str]):
        """
        Calls the live BLAST service via its web endpoint to check for potential off-targets.
        """
        logger.info(f"Connecting via HTTP to BLAST service at: {BLAST_SERVICE_URL}")

        async def run_search(client, guide):
            logger.info(f"  BLASTing guide (V9 HTTP METHOD): {guide}")
            payload = {"query_sequence": guide}
            blast_url = f"{BLAST_SERVICE_URL}/search"
            logger.info(f"    -> Calling BLAST URL: {blast_url}")
            
            # Use the existing helper for the HTTP request
            result = await self._make_request(client, "POST", blast_url, json=payload, timeout=300.0)
            
            if "error" in result:
                logger.error(f"BLAST search for {guide} failed: {result.get('details', result['error'])}")
                off_targets = 999
            else:
                off_targets = self._parse_blast_xml(result.get("raw_blast_xml", ""), len(guide))
            
            return {
                "guide_sequence": guide,
                "off_target_count": off_targets,
                "confidence": 0.99 # High confidence as it's from a direct BLAST query
            }

        # Create a single client for all requests
        async with httpx.AsyncClient(follow_redirects=True) as client:
            tasks = [run_search(client, guide) for guide in guides]
            safety_data = await asyncio.gather(*tasks)
        
        logger.info("BLAST off-target analysis complete.")
        return safety_data

    async def check_immunogenicity(self, guides: list[str]):
        """
        Calls the Immunogenicity service to predict immune response.
        (Placeholder - returns mock data)
        """
        logger.warning("Immunogenicity client is a placeholder. Returning mock immunogenicity data.")
        await asyncio.sleep(1) # Simulate network latency
        mock_immuno_data = []
        for guide in guides:
            mock_immuno_data.append({
                "guide_sequence": guide,
                "immunogenicity_score": 0.05, # Low is good
                "confidence": 0.90
            })
        return mock_immuno_data 