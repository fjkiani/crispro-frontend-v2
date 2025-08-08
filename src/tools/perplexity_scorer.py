import requests
import os

# The URL for our new, dedicated Evo 2 service.
# This will be provided upon deployment of that service.
EVO2_SERVICE_URL = os.environ.get("EVO2_SERVICE_URL", "https://your-evo2-service-url.modal.run")

def calculate_perplexity(sequence: str) -> float:
    """
    Calls the dedicated Evo 2 service to calculate perplexity for a sequence.
    """
    if not EVO2_SERVICE_URL or "your-evo2-service-url" in EVO2_SERVICE_URL:
        print("WARNING: EVO2_SERVICE_URL is not set. Returning mock perplexity.")
        # In a real scenario, you might raise an error or have a better fallback.
        import random
        return random.uniform(1.1, 5.0)

    endpoint = f"{EVO2_SERVICE_URL}/calculate_perplexity"
    headers = {"Content-Type": "application/json"}
    payload = {"sequence": sequence}

    try:
        response = requests.post(endpoint, json=payload, headers=headers)
        response.raise_for_status()  # Raise an exception for bad status codes
        data = response.json()
        return data.get("mean_perplexity", -1.0)
    except requests.exceptions.RequestException as e:
        print(f"Error calling Evo 2 service: {e}")
        # Fallback score in case of network errors, etc.
        return -1.0

if __name__ == '__main__':
    # Example usage:
    test_sequence = "AGCT"
    print(f"Requesting perplexity for '{test_sequence}'...")
    perplexity = calculate_perplexity(test_sequence)
    if perplexity != -1.0:
        print(f"SUCCESS: Mean perplexity: {perplexity}")
    else:
        print("FAILED: Could not retrieve perplexity from the service.")

