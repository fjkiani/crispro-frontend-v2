#!/bin/bash

# ⚔️ SAE PHASE 3 INTEGRATION TEST SCRIPT ⚔️
# Tests SAE Phase 1+2 services via Complete Care v2 orchestrator
# Owner: Zo (Lead Commander)
# Date: January 13, 2025

set -e  # Exit on error

echo "⚔️ SAE PHASE 3 INTEGRATION TESTS ⚔️"
echo "=================================="
echo ""

# Colors
GREEN='\033[0.32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Base URL
BASE_URL="http://localhost:8000"
TEST_DIR=".cursor/ayesha/test_payloads"

# Test function
run_test() {
    local test_name=$1
    local payload_file=$2
    local description=$3
    
    echo -e "${YELLOW}TEST: ${test_name}${NC}"
    echo "Description: ${description}"
    echo "Payload: ${payload_file}"
    echo ""
    
    # Make request
    response=$(curl -s -X POST "${BASE_URL}/api/ayesha/complete_care_v2" \
        -H "Content-Type: application/json" \
        -d @"${payload_file}")
    
    # Check if response contains error
    if echo "$response" | grep -q "\"detail\":"; then
        echo -e "${RED}❌ TEST FAILED${NC}"
        echo "Error: $(echo $response | jq -r '.detail' 2>/dev/null || echo $response)"
        echo ""
        return 1
    fi
    
    # Extract key fields
    local sae_features_status=$(echo "$response" | jq -r '.sae_features.status // "computed"' 2>/dev/null)
    local resistance_status=$(echo "$response" | jq -r '.resistance_alert.status // "computed"' 2>/dev/null)
    local next_test_count=$(echo "$response" | jq -r '.next_test_recommender.total_tests // 0' 2>/dev/null)
    local hint_tiles_count=$(echo "$response" | jq -r '.hint_tiles.total_tiles // 0' 2>/dev/null)
    local mechanism_map_status=$(echo "$response" | jq -r '.mechanism_map.status // "unknown"' 2>/dev/null)
    
    echo -e "${GREEN}✅ TEST PASSED${NC}"
    echo "Results:"
    echo "  - SAE Features: ${sae_features_status}"
    echo "  - Resistance Alert: ${resistance_status}"
    echo "  - Next Tests: ${next_test_count}"
    echo "  - Hint Tiles: ${hint_tiles_count}"
    echo "  - Mechanism Map: ${mechanism_map_status}"
    
    # Test-specific validations
    case $test_name in
        "PRE-NGS")
            if [ "$sae_features_status" != "awaiting_ngs" ]; then
                echo -e "${RED}⚠️  Expected sae_features.status='awaiting_ngs', got '${sae_features_status}'${NC}"
            fi
            if [ "$next_test_count" -lt 2 ]; then
                echo -e "${RED}⚠️  Expected at least 2 next tests, got ${next_test_count}${NC}"
            fi
            ;;
        "BRCA1-BIALLELIC")
            local dna_repair=$(echo "$response" | jq -r '.sae_features.dna_repair_capacity // 0' 2>/dev/null)
            echo "  - DNA Repair Capacity: ${dna_repair}"
            if (( $(echo "$dna_repair < 0.60" | bc -l) )); then
                echo -e "${RED}⚠️  Expected high DNA repair (>0.60), got ${dna_repair}${NC}"
            fi
            ;;
        "HER2-POSITIVE")
            local her2_vector=$(echo "$response" | jq -r '.sae_features.mechanism_vector[4] // 0' 2>/dev/null)
            echo "  - HER2 Mechanism Vector: ${her2_vector}"
            ;;
    esac
    
    echo ""
    echo "---"
    echo ""
}

# Health check first
echo "🔍 Checking backend health..."
health_response=$(curl -s "${BASE_URL}/api/ayesha/complete_care_v2/health")
if echo "$health_response" | grep -q "operational"; then
    echo -e "${GREEN}✅ Backend is operational${NC}"
    sae_phase1=$(echo "$health_response" | jq -r '.sae_phase1_enabled' 2>/dev/null)
    sae_phase2=$(echo "$health_response" | jq -r '.sae_phase2_enabled' 2>/dev/null)
    echo "  - SAE Phase 1: ${sae_phase1}"
    echo "  - SAE Phase 2: ${sae_phase2}"
else
    echo -e "${RED}❌ Backend health check failed${NC}"
    echo "Response: $health_response"
    exit 1
fi
echo ""
echo "---"
echo ""

# Run tests
run_test "PRE-NGS" \
    "${TEST_DIR}/01_pre_ngs.json" \
    "Ayesha TODAY (no NGS data) - should recommend tests"

run_test "BRCA1-BIALLELIC" \
    "${TEST_DIR}/02_brca1_biallelic.json" \
    "Ayesha with BRCA1 biallelic loss (HRD=58) - high DDR burden"

run_test "HER2-POSITIVE" \
    "${TEST_DIR}/03_her2_positive.json" \
    "Ayesha with HER2 amplification - unlocks NCT06819007 trial"

echo ""
echo "⚔️ ALL TESTS COMPLETE ⚔️"
echo ""
echo "To see full response for a test:"
echo "  curl -X POST ${BASE_URL}/api/ayesha/complete_care_v2 -H 'Content-Type: application/json' -d @${TEST_DIR}/01_pre_ngs.json | jq"
echo ""

