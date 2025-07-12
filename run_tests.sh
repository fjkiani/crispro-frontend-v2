#!/bin/bash
# This script runs the entire test suite for the CasPro platform.

echo "--- 🚀 Starting CasPro Test Suite 🚀 ---"

# Stop script if any command fails
set -e

echo -e "\n--- 🧪 Running Database Client Unit Tests ---"
python3 tests/test_threat_matrix_client.py

echo -e "\n--- 🛰️  Running CommandCenter API Integration Tests ---"
python3 tests/test_api_endpoints.py

echo -e "\n--- ✅ All tests passed successfully! ---" 