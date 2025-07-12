#!/usr/bin/env python3
"""
Test runner for the AI Research Assistant streamlit application.
Run this script to verify that the application is working correctly.
"""

import os
import sys
import unittest
import requests
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# --- Test Configuration ---
# DEPLOYMENT OVERRIDE: Forcing the test to use the new v2 endpoint.
COMMAND_CENTER_URL = "https://crispro--command-center-v2-commandcenter-api.modal.run"
# if not COMMAND_CENTER_URL:
#     print("❌ CRITICAL: COMMAND_CENTER_URL environment variable is not set and no fallback is available. Cannot run tests.")
#     sys.exit(1)

BASE_URL = COMMAND_CENTER_URL.rstrip('/')

class TestCommandCenterAPI(unittest.TestCase):

    def test_health_check(self):
        """Tests if the CommandCenter health check endpoint is responsive."""
        print(f"\n--- Running Test: Health Check Endpoint ---")
        response = requests.post(f"{BASE_URL}/health", timeout=30)
        self.assertEqual(response.status_code, 200, f"Expected status 200, got {response.status_code}")
        self.assertEqual(response.json(), {"status": "healthy", "service": "command-center"}, "Health check response body mismatch")
        print(f"✅ PASS: Health Check Endpoint")

    def test_assess_threat_endpoint_frameshift(self):
        """Tests the /workflow/assess_threat endpoint with a known frameshift mutation."""
        print(f"\n--- Running Test: Threat Assessor (Frameshift) ---")
        endpoint = "/workflow/assess_threat"
        payload = {"gene_symbol": "ASXL1", "protein_change": "p.Gly646fs"}
        response = requests.post(f"{BASE_URL}{endpoint}", json=payload, timeout=120)
        
        self.assertEqual(response.status_code, 200, f"Expected status 200, got {response.status_code}")
        data = response.json()
        print(f"DEBUG: Full API response for frameshift test: {data}")
        
        self.assertEqual(data['gene_symbol'], "ASXL1", "Returned gene symbol mismatch")
        self.assertEqual(data['protein_change'], "p.Gly646fs", "Returned protein change mismatch")
        self.assertEqual(data['assessment']['assessment_source'], "Truncation Sieve", "Assessment source should be Truncation Sieve")
        self.assertTrue(data['assessment']['is_pathogenic'], "Frameshift mutation should be marked as pathogenic")
        
        print(f"✅ PASS: Threat Assessor (Frameshift)")

    def test_assess_threat_endpoint_missense(self):
        """
        Tests the /workflow/assess_threat endpoint with a known missense mutation.
        This should trigger the 'Zeta Oracle' path. Since the oracle may not be configured,
        we'll check for either a proper oracle response or an error condition.
        """
        print(f"\n--- Running Test: Threat Assessor (Missense) ---")
        endpoint = "/workflow/assess_threat"
        payload = {"gene_symbol": "BRAF", "protein_change": "p.Val600Glu"} # V600E
        response = requests.post(f"{BASE_URL}{endpoint}", json=payload, timeout=120)
        
        self.assertEqual(response.status_code, 200, f"Expected status 200, got {response.status_code}")
        data = response.json()
        print(f"DEBUG: Full API response for missense test: {data}")
        
        self.assertEqual(data['gene_symbol'], "BRAF", "Returned gene symbol mismatch")
        self.assertEqual(data['protein_change'], "p.Val600Glu", "Returned protein change mismatch")
        
        # Assertions for a successful missense assessment
        self.assertIsNone(data.get('error'), "Error field should be null in a successful response")
        self.assertIsNotNone(data['assessment'], "Assessment dictionary should be present")
        
        # DEFINITIVE FIX: The Triumvirate Protocol correctly uses the Adjudicator for missense mutations.
        # The test must be updated to reflect this new, correct logic.
        self.assertEqual(data['assessment']['assessment_source'], "Adjudicator", "Assessment source should be Adjudicator")
        self.assertTrue(data['assessment']['is_pathogenic'], "BRAF V600E should be classified as pathogenic")
        self.assertIn(data['assessment']['prediction'], ["Pathogenic", 1], "Prediction should be 'Pathogenic' or 1")
        print("✅ PASS: Threat Assessor (Missense)")

if __name__ == "__main__":
    print("🚀 Starting CommandCenter API Validation Suite 🚀\n")
    unittest.main() 