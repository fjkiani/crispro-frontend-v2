import httpx
from loguru import logger
from typing import List, Dict, Any

ZETA_ORACLE_URL = "https://crispro--fusion-engine-v1-fusionengine-api.modal.run"

class ZetaOracleClient:
    """
    A client for interacting with the Zeta Oracle (fusion-engine-v1) service.
    """

    async def get_zeta_score(self, hgvs: str, alphamissense_variant: str, protein_sequence: str) -> float:
        """
        Gets the Zeta Score for a single variant.
        """
        logger.info(f"Querying Zeta Oracle for variant: {hgvs}")
        
        variant_payload = {
            "variant_id": hgvs,
            "hgvs": hgvs,
            "alphamissense_variant_str": alphamissense_variant
        }
        
        request_payload = {
            "protein_sequence": protein_sequence,
            "variants": [variant_payload]
        }
        
        async with httpx.AsyncClient(timeout=300.0) as client:
            try:
                response = await client.post(f"{ZETA_ORACLE_URL}/score_variants", json=request_payload)
                response.raise_for_status()
                
                api_response = response.json()
                scored_variants = api_response.get("scored_variants", [])
                
                if not scored_variants:
                    logger.error("Zeta Oracle returned no scored variants.")
                    return -999.0
                
                # The Zeta Score is the primary output of the oracle
                zeta_score = scored_variants[0].get("zeta_score", -999.0)
                logger.success(f"Zeta Oracle returned score: {zeta_score} for {hgvs}")
                return zeta_score
                
            except httpx.HTTPStatusError as e:
                logger.error(f"Zeta Oracle HTTP Error: {e.response.status_code} - {e.response.text}")
                return -999.0
            except Exception as e:
                logger.error(f"An unexpected error occurred while calling Zeta Oracle: {e}")
                return -999.0 
from loguru import logger
from typing import List, Dict, Any

ZETA_ORACLE_URL = "https://crispro--fusion-engine-v1-fusionengine-api.modal.run"

class ZetaOracleClient:
    """
    A client for interacting with the Zeta Oracle (fusion-engine-v1) service.
    """

    async def get_zeta_score(self, hgvs: str, alphamissense_variant: str, protein_sequence: str) -> float:
        """
        Gets the Zeta Score for a single variant.
        """
        logger.info(f"Querying Zeta Oracle for variant: {hgvs}")
        
        variant_payload = {
            "variant_id": hgvs,
            "hgvs": hgvs,
            "alphamissense_variant_str": alphamissense_variant
        }
        
        request_payload = {
            "protein_sequence": protein_sequence,
            "variants": [variant_payload]
        }
        
        async with httpx.AsyncClient(timeout=300.0) as client:
            try:
                response = await client.post(f"{ZETA_ORACLE_URL}/score_variants", json=request_payload)
                response.raise_for_status()
                
                api_response = response.json()
                scored_variants = api_response.get("scored_variants", [])
                
                if not scored_variants:
                    logger.error("Zeta Oracle returned no scored variants.")
                    return -999.0
                
                # The Zeta Score is the primary output of the oracle
                zeta_score = scored_variants[0].get("zeta_score", -999.0)
                logger.success(f"Zeta Oracle returned score: {zeta_score} for {hgvs}")
                return zeta_score
                
            except httpx.HTTPStatusError as e:
                logger.error(f"Zeta Oracle HTTP Error: {e.response.status_code} - {e.response.text}")
                return -999.0
            except Exception as e:
                logger.error(f"An unexpected error occurred while calling Zeta Oracle: {e}")
                return -999.0 
from loguru import logger
from typing import List, Dict, Any

ZETA_ORACLE_URL = "https://crispro--fusion-engine-v1-fusionengine-api.modal.run"

class ZetaOracleClient:
    """
    A client for interacting with the Zeta Oracle (fusion-engine-v1) service.
    """

    async def get_zeta_score(self, hgvs: str, alphamissense_variant: str, protein_sequence: str) -> float:
        """
        Gets the Zeta Score for a single variant.
        """
        logger.info(f"Querying Zeta Oracle for variant: {hgvs}")
        
        variant_payload = {
            "variant_id": hgvs,
            "hgvs": hgvs,
            "alphamissense_variant_str": alphamissense_variant
        }
        
        request_payload = {
            "protein_sequence": protein_sequence,
            "variants": [variant_payload]
        }
        
        async with httpx.AsyncClient(timeout=300.0) as client:
            try:
                response = await client.post(f"{ZETA_ORACLE_URL}/score_variants", json=request_payload)
                response.raise_for_status()
                
                api_response = response.json()
                scored_variants = api_response.get("scored_variants", [])
                
                if not scored_variants:
                    logger.error("Zeta Oracle returned no scored variants.")
                    return -999.0
                
                # The Zeta Score is the primary output of the oracle
                zeta_score = scored_variants[0].get("zeta_score", -999.0)
                logger.success(f"Zeta Oracle returned score: {zeta_score} for {hgvs}")
                return zeta_score
                
            except httpx.HTTPStatusError as e:
                logger.error(f"Zeta Oracle HTTP Error: {e.response.status_code} - {e.response.text}")
                return -999.0
            except Exception as e:
                logger.error(f"An unexpected error occurred while calling Zeta Oracle: {e}")
                return -999.0 
from loguru import logger
from typing import List, Dict, Any

ZETA_ORACLE_URL = "https://crispro--fusion-engine-v1-fusionengine-api.modal.run"

class ZetaOracleClient:
    """
    A client for interacting with the Zeta Oracle (fusion-engine-v1) service.
    """

    async def get_zeta_score(self, hgvs: str, alphamissense_variant: str, protein_sequence: str) -> float:
        """
        Gets the Zeta Score for a single variant.
        """
        logger.info(f"Querying Zeta Oracle for variant: {hgvs}")
        
        variant_payload = {
            "variant_id": hgvs,
            "hgvs": hgvs,
            "alphamissense_variant_str": alphamissense_variant
        }
        
        request_payload = {
            "protein_sequence": protein_sequence,
            "variants": [variant_payload]
        }
        
        async with httpx.AsyncClient(timeout=300.0) as client:
            try:
                response = await client.post(f"{ZETA_ORACLE_URL}/score_variants", json=request_payload)
                response.raise_for_status()
                
                api_response = response.json()
                scored_variants = api_response.get("scored_variants", [])
                
                if not scored_variants:
                    logger.error("Zeta Oracle returned no scored variants.")
                    return -999.0
                
                # The Zeta Score is the primary output of the oracle
                zeta_score = scored_variants[0].get("zeta_score", -999.0)
                logger.success(f"Zeta Oracle returned score: {zeta_score} for {hgvs}")
                return zeta_score
                
            except httpx.HTTPStatusError as e:
                logger.error(f"Zeta Oracle HTTP Error: {e.response.status_code} - {e.response.text}")
                return -999.0
            except Exception as e:
                logger.error(f"An unexpected error occurred while calling Zeta Oracle: {e}")
                return -999.0 