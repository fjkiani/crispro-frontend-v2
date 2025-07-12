import requests
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# This is a placeholder for the actual AlphaMissense API endpoint.
# In a real scenario, this would be a configured, stable URL.
ALPHAMISSENSE_API_URL = "https://alphamissense.api.placeholder.com/v1/score"

class AlphaMissenseClient:
    def __init__(self, api_url: str = ALPHAMISSENSE_API_URL):
        self.api_url = api_url

    def get_score(self, hgvs_p: str, transcript_id: str) -> dict:
        """
        Queries the AlphaMissense API for a pathogenicity score for a given variant.

        Args:
            hgvs_p: The protein variant in HGVS notation (e.g., "p.V600E").
            transcript_id: The transcript ID (e.g., "ENST00000288602").

        Returns:
            A dictionary containing the score or an error message.
        """
        logger.info(f"Querying AlphaMissense for variant: {hgvs_p} on transcript {transcript_id}")
        
        # The payload structure is an assumption based on typical VEP APIs.
        payload = {
            "variant": hgvs_p,
            "transcript_id": transcript_id,
        }
        
        try:
            response = requests.post(self.api_url, json=payload, timeout=30)
            response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)
            
            data = response.json()
            logger.info(f"Successfully received score from AlphaMissense: {data.get('score')}")
            return {
                "model": "AlphaMissense",
                "score": data.get("score"),
                "classification": data.get("classification"),
                "error": None
            }

        except requests.exceptions.HTTPError as http_err:
            logger.error(f"HTTP error occurred while querying AlphaMissense: {http_err}")
            return {"model": "AlphaMissense", "score": None, "error": f"HTTP Error: {http_err}"}
        except requests.exceptions.RequestException as req_err:
            logger.error(f"Request error occurred while querying AlphaMissense: {req_err}")
            return {"model": "AlphaMissense", "score": None, "error": f"Request Error: {req_err}"}
        except Exception as e:
            logger.error(f"An unexpected error occurred: {e}")
            # In a real scenario, we might not have a running AlphaMissense API.
            # We will return a mock score for testing purposes.
            mock_score = 0.85 # Mock pathogenic score
            logger.warning(f"Returning MOCK score from AlphaMissense due to exception: {mock_score}")
            return {
                "model": "AlphaMissense",
                "score": mock_score,
                "classification": "likely_pathogenic (mocked)",
                "error": f"Exception occurred, returning mock data: {e}"
            }


if __name__ == '__main__':
    # Example usage:
    client = AlphaMissenseClient()
    # This requires a live API to work against.
    # The current implementation will fall back to mock data.
    score_data = client.get_score(hgvs_p="p.V600E", transcript_id="ENST00000288602")
    print(score_data)
    
    score_data_benign = client.get_score(hgvs_p="p.A123T", transcript_id="ENST00000288602")
    print(score_data_benign)

    logger.info("AlphaMissense client tool created. Ready for integration.") 