#!/bin/bash

###############################################################################
# ⚔️ REAL CO-PILOT BACKEND ENDPOINT TESTS ⚔️
# Tests all 10 capabilities with REALISTIC payloads
# Shows ACTUAL input/output, not just pass/fail
###############################################################################

set -e

API_BASE="http://127.0.0.1:8000"
TEST_OUTPUT_DIR=".cursor/tests/copilot/output"
mkdir -p "$TEST_OUTPUT_DIR"

echo "⚔️ CO-PILOT BACKEND ENDPOINT TESTS ⚔️"
echo "======================================"
echo ""
echo "Testing ${API_BASE}"
echo "Output directory: ${TEST_OUTPUT_DIR}"
echo ""

###############################################################################
# TEST 1: DRUG EFFICACY (WIWFM)
###############################################################################
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "TEST 1: DRUG EFFICACY (WIWFM)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Query: 'Will Olaparib work for me?'"
echo "Context: BRCA1 p.Cys61Gly, Ovarian Cancer, HRD+, L3 post-platinum"
echo ""

curl -sS -X POST "${API_BASE}/api/efficacy/predict" \
  -H 'Content-Type: application/json' \
  -d '{
    "mutations": [{"gene":"BRCA1","hgvs_p":"p.Cys61Gly","chrom":"17","pos":43115779,"ref":"T","alt":"G","build":"GRCh38"}],
    "disease": "ovarian_cancer",
    "treatment_history": [
      {"line": 1, "drugs": ["carboplatin", "paclitaxel"], "outcome": "progression"},
      {"line": 2, "drugs": ["olaparib"], "outcome": "progression"}
    ],
    "biomarkers": {"HRD": "POSITIVE", "TMB": 8.2},
    "include_sae_features": true,
    "options": {"profile": "baseline"}
  }' | python3 -m json.tool > "${TEST_OUTPUT_DIR}/1_drug_efficacy.json" 2>&1

if [ $? -eq 0 ]; then
  echo "✅ SUCCESS"
  echo ""
  echo "📊 KEY RESULTS:"
  python3 -c "
import json, sys
try:
    with open('${TEST_OUTPUT_DIR}/1_drug_efficacy.json', 'r') as f:
        data = json.load(f)
    
    if 'drugs' in data and len(data['drugs']) > 0:
        top_drug = data['drugs'][0]
        print(f\"  Top Drug: {top_drug.get('name', 'N/A')}\")
        print(f\"  Efficacy Score: {top_drug.get('efficacy_score', 'N/A')}\")
        print(f\"  Confidence: {top_drug.get('confidence', 'N/A')}\")
        print(f\"  Evidence Tier: {top_drug.get('evidence_tier', 'N/A')}\")
        
        if 'sae_features' in top_drug:
            sae = top_drug['sae_features']
            print(f\"  SAE Line Appropriateness: {sae.get('line_appropriateness', 'N/A')}\")
            print(f\"  SAE Cross Resistance: {sae.get('cross_resistance', 'N/A')}\")
    else:
        print('  No drugs returned')
except Exception as e:
    print(f'  Error parsing response: {e}')
"
else
  echo "❌ FAILED"
fi
echo ""

###############################################################################
# TEST 2: FOOD/SUPPLEMENT VALIDATOR
###############################################################################
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "TEST 2: FOOD/SUPPLEMENT VALIDATOR"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Query: 'Should I take Vitamin D?'"
echo "Context: Ovarian Cancer, HRD+, L3 post-platinum, taking warfarin"
echo ""

curl -sS -X POST "${API_BASE}/api/hypothesis/validate_food_dynamic" \
  -H 'Content-Type: application/json' \
  -d '{
    "compound": "Vitamin D",
    "disease_context": {
      "disease": "ovarian_cancer_hgs",
      "biomarkers": {"HRD": "POSITIVE", "TMB": 8.2},
      "pathways_disrupted": ["DNA repair", "Cell cycle"]
    },
    "treatment_history": {
      "current_line": "L3",
      "prior_therapies": ["carboplatin", "paclitaxel", "olaparib"]
    },
    "patient_medications": ["warfarin"],
    "use_llm": true
  }' | python3 -m json.tool > "${TEST_OUTPUT_DIR}/2_food_validator.json" 2>&1

if [ $? -eq 0 ]; then
  echo "✅ SUCCESS"
  echo ""
  echo "📊 KEY RESULTS:"
  python3 -c "
import json, sys
try:
    with open('${TEST_OUTPUT_DIR}/2_food_validator.json', 'r') as f:
        data = json.load(f)
    
    print(f\"  Verdict: {data.get('verdict', 'N/A')}\")
    print(f\"  Overall Score: {data.get('overall_score', 'N/A')}\")
    
    if 'evidence' in data:
        print(f\"  Evidence Grade: {data['evidence'].get('grade', 'N/A')}\")
    
    if 'recommendations' in data:
        rec = data['recommendations']
        print(f\"  Dosage: {rec.get('dosage', 'N/A')}\")
        print(f\"  Timing: {rec.get('timing', 'N/A')}\")
        
        if rec.get('interactions'):
            print(f\"  Interactions: {rec['interactions']}\")
except Exception as e:
    print(f'  Error parsing response: {e}')
"
else
  echo "❌ FAILED"
fi
echo ""

###############################################################################
# TEST 3: COMPLETE CARE PLAN
###############################################################################
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "TEST 3: COMPLETE CARE PLAN (UNIFIED)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Query: 'Give me a complete care plan'"
echo "Context: BRCA1, Ovarian Cancer, HRD+, L3"
echo ""

curl -sS -X POST "${API_BASE}/api/ayesha/complete_care_plan" \
  -H 'Content-Type: application/json' \
  -d '{
    "patient_context": {
      "disease": "ovarian_cancer",
      "mutations": [{"gene":"BRCA1","hgvs_p":"p.Cys61Gly"}],
      "biomarkers": {"HRD": "POSITIVE", "TMB": 8.2},
      "treatment_history": {"line": 3, "prior_therapies": ["carboplatin", "olaparib"]}
    }
  }' | python3 -m json.tool > "${TEST_OUTPUT_DIR}/3_complete_care.json" 2>&1

if [ $? -eq 0 ]; then
  echo "✅ SUCCESS"
  echo ""
  echo "📊 KEY RESULTS:"
  python3 -c "
import json, sys
try:
    with open('${TEST_OUTPUT_DIR}/3_complete_care.json', 'r') as f:
        data = json.load(f)
    
    if 'drug_recommendations' in data:
        print(f\"  Drug Recommendations: {len(data['drug_recommendations'])} drugs\")
        if data['drug_recommendations']:
            top = data['drug_recommendations'][0]
            print(f\"    Top: {top.get('name', 'N/A')} (score: {top.get('efficacy_score', 'N/A')})\")
    
    if 'food_recommendations' in data:
        print(f\"  Food Recommendations: {len(data['food_recommendations'])} compounds\")
        if data['food_recommendations']:
            top = data['food_recommendations'][0]
            print(f\"    Top: {top.get('compound', 'N/A')} (score: {top.get('overall_score', 'N/A')})\")
    
    print(f\"  Integrated Confidence: {data.get('integrated_confidence', 'N/A')}\")
except Exception as e:
    print(f'  Error parsing response: {e}')
"
else
  echo "❌ FAILED"
fi
echo ""

###############################################################################
# TEST 4: CLINICAL TRIALS
###############################################################################
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "TEST 4: CLINICAL TRIALS"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Query: 'Find clinical trials for me'"
echo "Context: 55yo female, ovarian cancer, HRD+, post-platinum"
echo ""

curl -sS -X POST "${API_BASE}/api/trials/agent/search" \
  -H 'Content-Type: application/json' \
  -d '{
    "patient_summary": "55yo female, ovarian cancer, BRCA1 p.Cys61Gly, HRD+, TMB 8.2, Line 3, post-carboplatin, paclitaxel"
  }' | python3 -m json.tool > "${TEST_OUTPUT_DIR}/4_trials.json" 2>&1

if [ $? -eq 0 ]; then
  echo "✅ SUCCESS"
  echo ""
  echo "📊 KEY RESULTS:"
  python3 -c "
import json, sys
try:
    with open('${TEST_OUTPUT_DIR}/4_trials.json', 'r') as f:
        data = json.load(f)
    
    if 'trials' in data:
        print(f\"  Trials Found: {len(data['trials'])}\")
        if data['trials']:
            top = data['trials'][0]
            print(f\"    NCT: {top.get('nct_id', 'N/A')}\")
            print(f\"    Title: {top.get('title', 'N/A')[:60]}...\")
            print(f\"    Relevance: {top.get('relevance', 'N/A')}\")
            print(f\"    Status: {top.get('status', 'N/A')}\")
    else:
        print('  No trials field in response')
except Exception as e:
    print(f'  Error parsing response: {e}')
"
else
  echo "❌ FAILED"
fi
echo ""

###############################################################################
# TEST 5: TOXICITY RISK (PGx)
###############################################################################
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "TEST 5: TOXICITY RISK (PGx)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Query: 'Will platinum be toxic for me?'"
echo "Context: DPYD variant (germline)"
echo ""

curl -sS -X POST "${API_BASE}/api/safety/toxicity_risk" \
  -H 'Content-Type: application/json' \
  -d '{
    "patient": {
      "germlineVariants": [
        {"gene":"DPYD","chrom":"1","pos":97450058,"ref":"C","alt":"T"}
      ]
    },
    "candidate": {
      "type": "drug",
      "name": "carboplatin",
      "moa": "platinum_agent"
    }
  }' | python3 -m json.tool > "${TEST_OUTPUT_DIR}/5_toxicity.json" 2>&1

if [ $? -eq 0 ]; then
  echo "✅ SUCCESS"
  echo ""
  echo "📊 KEY RESULTS:"
  python3 -c "
import json, sys
try:
    with open('${TEST_OUTPUT_DIR}/5_toxicity.json', 'r') as f:
        data = json.load(f)
    
    print(f\"  Risk Score: {data.get('risk_score', 'N/A')}\")
    print(f\"  Risk Level: {data.get('risk_level', 'N/A')}\")
    
    if 'factors' in data:
        print(f\"  Factors: {data['factors']}\")
    
    if 'warnings' in data:
        print(f\"  Warnings: {data['warnings']}\")
except Exception as e:
    print(f'  Error parsing response: {e}')
"
else
  echo "❌ FAILED"
fi
echo ""

###############################################################################
# TEST 6: SYNTHETIC LETHALITY
###############################################################################
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "TEST 6: SYNTHETIC LETHALITY"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Query: 'What's my cancer's weakness?'"
echo "Context: BRCA1 mutation"
echo ""

curl -sS -X POST "${API_BASE}/api/guidance/synthetic_lethality" \
  -H 'Content-Type: application/json' \
  -d '{
    "disease": "ovarian_cancer",
    "mutations": [{"gene":"BRCA1","hgvs_p":"p.Cys61Gly","chrom":"17","pos":43115779,"ref":"T","alt":"G"}],
    "api_base": "http://127.0.0.1:8000"
  }' | python3 -m json.tool > "${TEST_OUTPUT_DIR}/6_synthetic_lethality.json" 2>&1

if [ $? -eq 0 ]; then
  echo "✅ SUCCESS"
  echo ""
  echo "📊 KEY RESULTS:"
  python3 -c "
import json, sys
try:
    with open('${TEST_OUTPUT_DIR}/6_synthetic_lethality.json', 'r') as f:
        data = json.load(f)
    
    print(f\"  Suggested Therapy: {data.get('suggested_therapy', 'N/A')}\")
    
    if 'damage_report' in data:
        print(f\"  Damage Report: {len(data['damage_report'])} genes\")
    
    if 'essentiality_report' in data:
        print(f\"  Essentiality Report: {len(data['essentiality_report'])} genes\")
    
    if 'guidance' in data:
        g = data['guidance']
        print(f\"  Guidance Tier: {g.get('tier', 'N/A')}\")
        print(f\"  Guidance Confidence: {g.get('confidence', 'N/A')}\")
except Exception as e:
    print(f'  Error parsing response: {e}')
"
else
  echo "❌ FAILED"
fi
echo ""

###############################################################################
# TEST 7: VARIANT IMPACT (DEEP ANALYSIS)
###############################################################################
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "TEST 7: VARIANT IMPACT (DEEP ANALYSIS)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Query: 'What is the functional impact of BRAF V600E?'"
echo ""

curl -sS -X POST "${API_BASE}/api/evidence/deep_analysis" \
  -H 'Content-Type: application/json' \
  -d '{
    "gene": "BRAF",
    "hgvs_p": "V600E",
    "disease": "melanoma"
  }' | python3 -m json.tool > "${TEST_OUTPUT_DIR}/7_variant_impact.json" 2>&1

if [ $? -eq 0 ]; then
  echo "✅ SUCCESS"
  echo ""
  echo "📊 KEY RESULTS:"
  python3 -c "
import json, sys
try:
    with open('${TEST_OUTPUT_DIR}/7_variant_impact.json', 'r') as f:
        data = json.load(f)
    
    if 'clinvar' in data:
        print(f\"  ClinVar Classification: {data['clinvar'].get('classification', 'N/A')}\")
        print(f\"  ClinVar Review Status: {data['clinvar'].get('review_status', 'N/A')}\")
    
    if 'functional_impact' in data:
        print(f\"  Functional Impact: {data['functional_impact']}\")
    
    if 'evidence_level' in data:
        print(f\"  Evidence Level: {data['evidence_level']}\")
except Exception as e:
    print(f'  Error parsing response: {e}')
"
else
  echo "❌ FAILED"
fi
echo ""

###############################################################################
# TEST 8: RADIATION GUIDANCE
###############################################################################
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "TEST 8: RADIATION GUIDANCE (PrecisionRad)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Query: 'Should I get radiation therapy?'"
echo "Context: TP53 R175H"
echo ""

curl -sS -X POST "${API_BASE}/api/guidance/radonc" \
  -H 'Content-Type: application/json' \
  -d '{
    "disease": "breast_cancer",
    "mutations": [{"gene":"TP53","hgvs_p":"R175H","chrom":"17","pos":7673803,"ref":"C","alt":"T","build":"GRCh38"}],
    "options": {"adaptive": true, "ensemble": true}
  }' | python3 -m json.tool > "${TEST_OUTPUT_DIR}/8_radiation_guidance.json" 2>&1

if [ $? -eq 0 ]; then
  echo "✅ SUCCESS"
  echo ""
  echo "📊 KEY RESULTS:"
  python3 -c "
import json, sys
try:
    with open('${TEST_OUTPUT_DIR}/8_radiation_guidance.json', 'r') as f:
        data = json.load(f)
    
    print(f\"  Modality: {data.get('modality', 'N/A')}\")
    print(f\"  Tier: {data.get('tier', 'N/A')}\")
    print(f\"  Radiosensitivity Score: {data.get('radiosensitivity_score', 'N/A')}\")
    print(f\"  Confidence: {data.get('confidence', 'N/A')}\")
    print(f\"  On-Label: {data.get('on_label', 'N/A')}\")
except Exception as e:
    print(f'  Error parsing response: {e}')
"
else
  echo "❌ FAILED"
fi
echo ""

###############################################################################
# TEST 9: CHEMO GUIDANCE
###############################################################################
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "TEST 9: CHEMO GUIDANCE"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Query: 'Should I use platinum chemotherapy?'"
echo "Context: BRCA1, HRD+"
echo ""

curl -sS -X POST "${API_BASE}/api/guidance/chemo" \
  -H 'Content-Type: application/json' \
  -d '{
    "disease": "ovarian_cancer",
    "drug_or_class": "platinum",
    "mutations": [{"gene":"BRCA1","hgvs_p":"p.Cys61Gly"}],
    "options": {"adaptive": true, "ensemble": true}
  }' | python3 -m json.tool > "${TEST_OUTPUT_DIR}/9_chemo_guidance.json" 2>&1

if [ $? -eq 0 ]; then
  echo "✅ SUCCESS"
  echo ""
  echo "📊 KEY RESULTS:"
  python3 -c "
import json, sys
try:
    with open('${TEST_OUTPUT_DIR}/9_chemo_guidance.json', 'r') as f:
        data = json.load(f)
    
    print(f\"  Tier: {data.get('tier', 'N/A')}\")
    print(f\"  Strength: {data.get('strength', 'N/A')}\")
    print(f\"  On-Label: {data.get('on_label', 'N/A')}\")
    print(f\"  Efficacy Score: {data.get('efficacy_score', 'N/A')}\")
    print(f\"  Confidence: {data.get('confidence', 'N/A')}\")
    print(f\"  Evidence Tier: {data.get('evidence_tier', 'N/A')}\")
except Exception as e:
    print(f'  Error parsing response: {e}')
"
else
  echo "❌ FAILED"
fi
echo ""

###############################################################################
# TEST 10: LITERATURE RETRIEVAL (RAG FALLBACK)
###############################################################################
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "TEST 10: LITERATURE RETRIEVAL (RAG)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Query: 'Find papers on BRCA1 mutations'"
echo ""

curl -sS -X POST "${API_BASE}/api/evidence/rag-query" \
  -H 'Content-Type: application/json' \
  -d '{
    "query": "Find papers on BRCA1 mutations in ovarian cancer",
    "gene": "BRCA1",
    "disease": "ovarian_cancer"
  }' | python3 -m json.tool > "${TEST_OUTPUT_DIR}/10_rag.json" 2>&1

if [ $? -eq 0 ]; then
  echo "✅ SUCCESS"
  echo ""
  echo "📊 KEY RESULTS:"
  python3 -c "
import json, sys
try:
    with open('${TEST_OUTPUT_DIR}/10_rag.json', 'r') as f:
        data = json.load(f)
    
    print(f\"  Answer Length: {len(data.get('answer', ''))} chars\")
    print(f\"  Evidence Level: {data.get('evidence_level', 'N/A')}\")
    print(f\"  Confidence Score: {data.get('confidence_score', 'N/A')}\")
    
    if 'supporting_papers' in data:
        print(f\"  Supporting Papers: {len(data['supporting_papers'])}\")
except Exception as e:
    print(f'  Error parsing response: {e}')
"
else
  echo "❌ FAILED"
fi
echo ""

###############################################################################
# SUMMARY
###############################################################################
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "SUMMARY"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "All output files saved to: ${TEST_OUTPUT_DIR}/"
echo ""
echo "Files:"
ls -lh "${TEST_OUTPUT_DIR}/"
echo ""
echo "⚔️ REAL BACKEND TESTS COMPLETE ⚔️"






