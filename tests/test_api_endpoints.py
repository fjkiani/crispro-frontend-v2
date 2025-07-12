import requests
import unittest
import json
import os

# Base URL for the deployed Modal Command Center service
BASE_URL = os.getenv("COMMAND_CENTER_API_URL", "https://crispro--command-center-commandcenter-api.modal.run")

class TestCommandCenterAPI(unittest.TestCase):
    # BASE_URL = "http://127.0.0.1:8000" # For local testing
    # BASE_URL = os.getenv("COMMAND_CENTER_API_URL", "YOUR_DEPLOYED_URL_HERE") # For deployed service

    # API Endpoints
    HEALTH_CHECK_ENDPOINT = f"{BASE_URL}/health"
    ASSESS_THREAT_ENDPOINT = f"{BASE_URL}/workflow/assess_threat"

    def test_01_health_check(self):
        """Tests if the CommandCenter service is alive."""
        print(f"\n🔬 Testing: API Health Check on {self.HEALTH_CHECK_ENDPOINT}")
        try:
            response = requests.post(self.HEALTH_CHECK_ENDPOINT, timeout=10)
            self.assertEqual(response.status_code, 200, f"Expected 200 OK, but got {response.status_code}. Response: {response.text}")
            self.assertEqual(response.json(), {"status": "healthy", "service": "command-center"})
            print("  - Passed: Service is healthy.")
        except requests.exceptions.RequestException as e:
            self.fail(f"API request failed for health check. Error: {e}")

    def test_02_assess_threat_with_specific_variant(self):
        """Tests the /workflow/assess_threat endpoint with a specific, known variant."""
        print("\n🔬 Testing: /workflow/assess_threat (Specific Variant)")
        payload = {"gene_symbol": "ASXL1", "protein_change": "p.Gly646fs"}
        print(f"  - Payload: {payload}")

        response = requests.post(f"{BASE_URL}/workflow/assess_threat", json=payload)
        self.assertEqual(response.status_code, 200)
        dossier = response.json()

        # Assert key values
        self.assertEqual(dossier["gene_symbol"], "ASXL1")
        self.assertIn("fs", dossier["protein_change"])
        self.assertEqual(dossier["variant_profile"]["mutation_id"], "COSM7413279") # Updated to match actual DB content
        self.assertIn("Frameshift", dossier["variant_profile"]["mutation_description"]) # Check for 'Frameshift' within the description

        print("  - Passed: Received a valid, structured Threat Dossier.")

    def test_03_assess_threat_gene_only(self):
        """Tests the endpoint with only a gene symbol, letting the backend pick a variant."""
        print("\n🔬 Testing: /workflow/assess_threat (Gene Only)")
        payload = {"gene_symbol": "BRAF"}
        print(f"  - Payload: {payload}")
        response = requests.post(f"{BASE_URL}/workflow/assess_threat", json=payload)
        self.assertEqual(response.status_code, 200)
        dossier = response.json()
        self.assertEqual(dossier["gene_symbol"], "BRAF")
        self.assertIn("p.V600E", dossier["protein_change"]) # Should find the famous V600E mutation if prioritized

    def test_04_assess_threat_invalid_gene(self):
        """Tests the endpoint's error handling for a non-existent gene."""
        print("\n🔬 Testing: /workflow/assess_threat (Invalid Gene)")
        payload = {"gene_symbol": "NONEXISTENTGENE123"}
        print(f"  - Payload: {payload}")
        response = requests.post(f"{BASE_URL}/workflow/assess_threat", json=payload)
        self.assertEqual(response.status_code, 200) # Still 200 OK even for errors, as per current design
        error_response = response.json()
        self.assertIn("error", error_response)
        self.assertIn("No variant profile found", error_response["error"])
        print("  - Passed: Correctly handled non-existent gene with an error message.")

if __name__ == '__main__':
    # Discover and run tests
    unittest.main() 