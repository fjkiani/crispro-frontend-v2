import requests
import os
from dotenv import load_dotenv
import urllib3

# --- Configuration ---
load_dotenv()

# Suppress the InsecureRequestWarning from urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

COMMAND_CENTER_URL = os.getenv("COMMAND_CENTER_URL")
if not COMMAND_CENTER_URL:
    print("⚠️ WARNING: COMMAND_CENTER_URL not found in .env file.")
    print("Using default URL: https://crispro--command-center-commandcenter-api.modal.run")
    COMMAND_CENTER_URL = "https://crispro--command-center-commandcenter-api.modal.run"

def main():
    """Tests only the health check endpoint of the Command Center."""
    print("🚀 Running Minimal Health Check Test 🚀")
    
    health_endpoint = "/health"
    url = f"{COMMAND_CENTER_URL}{health_endpoint}"
    
    print(f"  - Testing endpoint: {url}")
    
    try:
        response = requests.post(url, timeout=60, verify=False)
        response.raise_for_status()
        
        print("\n✅ SUCCESS: Health check successful.")
        print("--- RESPONSE ---")
        print(response.json())
        print("----------------")
        
    except requests.exceptions.RequestException as e:
        print(f"\n❌ FAILED: Health check failed.")
        print(f"  - Error: {e}")
        if e.response is not None:
            print(f"  - Status Code: {e.response.status_code}")
            try:
                print(f"  - Response Body: {e.response.json()}")
            except requests.exceptions.JSONDecodeError:
                print(f"  - Response Body: {e.response.text}")

if __name__ == "__main__":
    main() 