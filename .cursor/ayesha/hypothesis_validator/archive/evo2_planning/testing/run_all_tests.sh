#!/bin/bash
# Run all proof-of-concept tests

echo "⚔️ EVO2 FOOD VALIDATOR - PROOF OF CONCEPT TESTS"
echo "================================================"
echo ""

# Set API base if needed
export API_BASE=${API_BASE:-"http://127.0.0.1:8000"}

# Phase 1: Quick Wins
echo "📋 PHASE 1: Quick Wins (30 min)"
echo "--------------------------------"
echo ""

echo "🧪 Test 1: Evo2 API Format..."
python3 test_01_evo2_api_format.py
echo ""

echo "🧪 Test 2: Ensembl Gene Fetching..."
python3 test_02_ensembl_gene_fetch.py
echo ""

echo "🧪 Test 4: Knowledge Base Extraction..."
python3 test_04_knowledge_base.py
echo ""

# Phase 2: Core Hypothesis
echo ""
echo "📋 PHASE 2: Core Hypothesis (1 hour)"
echo "------------------------------------"
echo ""

echo "🧪 Test 3: Baseline vs Intervention (CRITICAL)..."
python3 test_03_baseline_vs_intervention.py
echo ""

# If Test 3 fails, run alternatives
if [ $? -ne 0 ]; then
    echo "⚠️  Test 3 failed - trying alternatives..."
    echo ""
    echo "🧪 Test 6: Variant Scoring Alternative..."
    python3 test_06_variant_scoring_alternative.py
    echo ""
    echo "🧪 Test 7: Generation Alternative..."
    python3 test_07_generate_alternative.py
    echo ""
fi

# Phase 3: Integration
echo ""
echo "📋 PHASE 3: Integration (30 min)"
echo "--------------------------------"
echo ""

echo "🧪 Test 8: End-to-End Minimal..."
python3 test_08_end_to_end_minimal.py
echo ""

echo "🧪 Test 9: Multi-Compound Validation..."
python3 test_09_multi_compound.py
echo ""

# Phase 4: Validation
echo ""
echo "📋 PHASE 4: Validation (30 min)"
echo "-------------------------------"
echo ""

echo "🧪 Test 10: Biological Sanity Check..."
python3 test_10_sanity_check.py
echo ""

echo ""
echo "✅ ALL TESTS COMPLETE"
echo "===================="
echo ""
echo "📊 Review results above to make GO/NO-GO decision"

