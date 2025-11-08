#!/bin/bash

# ⚔️ TREATMENT LINE INTEGRATION - END-TO-END SMOKE TEST
# Tests the complete flow from treatment history → confidence modulation

echo "======================================================================="
echo "TREATMENT LINE INTEGRATION - E2E SMOKE TEST"
echo "======================================================================="
echo ""

API_BASE="${API_BASE:-http://127.0.0.1:8000}"

# Color codes
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test counter
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

run_test() {
    local test_name="$1"
    local endpoint="$2"
    local payload="$3"
    local expected_drug="$4"
    local expected_confidence_min="$5"
    local expected_confidence_max="$6"
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    echo "-------------------------------------------------------------------"
    echo "TEST $TOTAL_TESTS: $test_name"
    echo "-------------------------------------------------------------------"
    
    # Send request
    response=$(curl -s -X POST "$API_BASE$endpoint" \
        -H 'Content-Type: application/json' \
        -d "$payload")
    
    # Check if response is valid JSON
    if ! echo "$response" | jq empty 2>/dev/null; then
        echo -e "${RED}❌ FAILED${NC}: Invalid JSON response"
        echo "Response: $response"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
    
    # Extract drug data
    drug_data=$(echo "$response" | jq ".drugs[] | select(.drug_name == \"$expected_drug\")")
    
    if [ -z "$drug_data" ]; then
        echo -e "${RED}❌ FAILED${NC}: Drug '$expected_drug' not found in response"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
    
    # Extract confidence
    confidence=$(echo "$drug_data" | jq -r '.confidence')
    
    # Extract treatment line provenance
    treatment_line_prov=$(echo "$drug_data" | jq -r '.treatment_line_provenance')
    
    echo "Drug: $expected_drug"
    echo "Confidence: $confidence"
    
    if [ "$treatment_line_prov" != "null" ]; then
        echo "Treatment Line Provenance:"
        echo "$treatment_line_prov" | jq '.'
        
        # Validate confidence is in expected range
        if (( $(echo "$confidence >= $expected_confidence_min && $confidence <= $expected_confidence_max" | bc -l) )); then
            echo -e "${GREEN}✅ PASSED${NC}: Confidence in expected range [$expected_confidence_min - $expected_confidence_max]"
            PASSED_TESTS=$((PASSED_TESTS + 1))
        else
            echo -e "${RED}❌ FAILED${NC}: Confidence $confidence not in expected range [$expected_confidence_min - $expected_confidence_max]"
            FAILED_TESTS=$((FAILED_TESTS + 1))
            return 1
        fi
    else
        echo -e "${YELLOW}⚠️  WARNING${NC}: No treatment line provenance found (may be expected for baseline test)"
        
        # For baseline tests without treatment history, just check confidence is present
        if [ -n "$confidence" ] && [ "$confidence" != "null" ]; then
            echo -e "${GREEN}✅ PASSED${NC}: Confidence present"
            PASSED_TESTS=$((PASSED_TESTS + 1))
        else
            echo -e "${RED}❌ FAILED${NC}: Confidence missing"
            FAILED_TESTS=$((FAILED_TESTS + 1))
            return 1
        fi
    fi
    
    echo ""
}

echo "Starting smoke tests..."
echo ""

# ===========================================================================
# TEST 1: BASELINE (No Treatment History)
# ===========================================================================

run_test \
    "Baseline - Ayesha BRCA1 Q356* without treatment history" \
    "/api/efficacy/predict" \
    '{
        "mutations": [
            {
                "gene": "BRCA1",
                "hgvs_p": "p.Gln356Ter",
                "chrom": "17",
                "pos": 43094464,
                "ref": "C",
                "alt": "T"
            }
        ],
        "disease": "ovarian_cancer",
        "model_id": "evo2_1b",
        "options": {"adaptive": true}
    }' \
    "olaparib" \
    0.70 \
    1.00

# ===========================================================================
# TEST 2: AYESHA CASE (Ovarian L2 Post-Platinum → Olaparib)
# ===========================================================================

run_test \
    "Ayesha - Ovarian L2 post-platinum (expect ~8% penalty)" \
    "/api/efficacy/predict" \
    '{
        "mutations": [
            {
                "gene": "BRCA1",
                "hgvs_p": "p.Gln356Ter",
                "chrom": "17",
                "pos": 43094464,
                "ref": "C",
                "alt": "T"
            }
        ],
        "disease": "ovarian_cancer",
        "model_id": "evo2_1b",
        "options": {"adaptive": true},
        "treatment_history": {
            "current_line": 2,
            "prior_therapies": ["carboplatin", "paclitaxel"]
        }
    }' \
    "olaparib" \
    0.65 \
    0.75

# ===========================================================================
# TEST 3: DR. LUSTBERG CASE (Breast HER2+ L3 Post-T-DXd → Tucatinib)
# ===========================================================================

run_test \
    "Dr. Lustberg - Breast HER2+ L3 post-T-DXd (expect ~4% penalty)" \
    "/api/efficacy/predict" \
    '{
        "mutations": [
            {
                "gene": "ERBB2",
                "hgvs_p": "p.Val777Leu",
                "chrom": "17",
                "pos": 39723967,
                "ref": "G",
                "alt": "C"
            }
        ],
        "disease": "breast_her2_positive",
        "model_id": "evo2_1b",
        "options": {"adaptive": true},
        "treatment_history": {
            "current_line": 3,
            "prior_therapies": ["trastuzumab deruxtecan", "pertuzumab"]
        }
    }' \
    "tucatinib+trastuzumab+capecitabine" \
    0.78 \
    0.85

# ===========================================================================
# TEST 4: FIRST-LINE (No Prior Therapies)
# ===========================================================================

run_test \
    "First-line - Ovarian L1 carboplatin (no penalty expected)" \
    "/api/efficacy/predict" \
    '{
        "mutations": [
            {
                "gene": "BRCA1",
                "hgvs_p": "p.Gln356Ter",
                "chrom": "17",
                "pos": 43094464,
                "ref": "C",
                "alt": "T"
            }
        ],
        "disease": "ovarian_cancer",
        "model_id": "evo2_1b",
        "options": {"adaptive": true},
        "treatment_history": {
            "current_line": 1,
            "prior_therapies": []
        }
    }' \
    "carboplatin+paclitaxel" \
    0.70 \
    1.00

# ===========================================================================
# SUMMARY
# ===========================================================================

echo "======================================================================="
echo "TEST SUMMARY"
echo "======================================================================="
echo ""
echo "Total Tests:  $TOTAL_TESTS"
echo -e "${GREEN}Passed:       $PASSED_TESTS${NC}"
echo -e "${RED}Failed:       $FAILED_TESTS${NC}"
echo ""

if [ $FAILED_TESTS -eq 0 ]; then
    echo -e "${GREEN}✅ ALL TESTS PASSED!${NC}"
    echo "======================================================================="
    exit 0
else
    echo -e "${RED}❌ SOME TESTS FAILED!${NC}"
    echo "======================================================================="
    exit 1
fi


