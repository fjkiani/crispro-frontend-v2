/**
 * ⚔️ COMPREHENSIVE Q2C ROUTER TEST ⚔️
 * Tests all 13 intents for the Clinical Co-Pilot
 * 
 * Run: node .cursor/tests/copilot/test_q2c_router_comprehensive.js
 */

// Import Q2C Router (adjust path if needed)
import { Q2C_ROUTER } from '../../../oncology-coPilot/oncology-frontend/src/components/CoPilot/Q2CRouter/intents.js';

const { classifyIntent } = Q2C_ROUTER;

// Test cases for all 13 intents
const TEST_CASES = [
  // 1. variant_impact
  {
    query: "What is the functional impact of BRAF V600E?",
    expected_intent: "variant_impact",
    expected_endpoint: "/api/evidence/deep_analysis"
  },
  
  // 2. drug_efficacy
  {
    query: "Will Olaparib work for me?",
    expected_intent: "drug_efficacy",
    expected_endpoint: "/api/efficacy/predict"
  },
  
  // 3. radonc_guidance
  {
    query: "Should I get radiation therapy?",
    expected_intent: "radonc_guidance",
    expected_endpoint: "/api/guidance/radonc"
  },
  
  // 4. chemo_guidance
  {
    query: "Should I use platinum chemotherapy?",
    expected_intent: "chemo_guidance",
    expected_endpoint: "/api/guidance/chemo"
  },
  
  // 5. literature_retrieval
  {
    query: "Find papers on BRCA1 mutations",
    expected_intent: "literature_retrieval",
    expected_endpoint: "/api/evidence/literature"
  },
  
  // 6. clinvar_context
  {
    query: "What does ClinVar say about this variant?",
    expected_intent: "clinvar_context",
    expected_endpoint: "/api/evidence/deep_analysis"
  },
  
  // 7. design_request
  {
    query: "Design a CRISPR guide RNA for this mutation",
    expected_intent: "design_request",
    expected_endpoint: "/api/design/guide_rna"
  },
  
  // 8. explain_result
  {
    query: "Why did I get this result?",
    expected_intent: "explain_result",
    expected_endpoint: "/api/evidence/explain"
  },
  
  // 9. food_validator ⚔️ NEW
  {
    query: "Should I take Vitamin D?",
    expected_intent: "food_validator",
    expected_endpoint: "/api/hypothesis/validate_food_dynamic"
  },
  
  // 10. trials ⚔️ NEW
  {
    query: "Find clinical trials for me",
    expected_intent: "trials",
    expected_endpoint: "/api/trials/agent/search"
  },
  
  // 11. complete_care ⚔️ NEW
  {
    query: "Give me a complete care plan",
    expected_intent: "complete_care",
    expected_endpoint: "/api/ayesha/complete_care_plan"
  },
  
  // 12. synthetic_lethality ⚔️ NEW
  {
    query: "What's my cancer's weakness?",
    expected_intent: "synthetic_lethality",
    expected_endpoint: "/api/guidance/synthetic_lethality"
  },
  
  // 13. toxicity_risk ⚔️ NEW
  {
    query: "Will this drug be toxic for me?",
    expected_intent: "toxicity_risk",
    expected_endpoint: "/api/safety/toxicity_risk"
  },
  
  // RAG Fallback (no intent match)
  {
    query: "Explain how cancer metastasis works",
    expected_intent: null,
    expected_endpoint: null,
    note: "Should fallback to RAG"
  }
];

// Run tests
console.log("⚔️ Q2C ROUTER COMPREHENSIVE TEST ⚔️\n");
console.log(`Testing ${TEST_CASES.length} test cases...\n`);

let passed = 0;
let failed = 0;

TEST_CASES.forEach((test, index) => {
  const result = classifyIntent(test.query);
  
  const intentMatch = result?.intent === test.expected_intent;
  const endpointMatch = result?.endpoint === test.expected_endpoint;
  
  if (intentMatch && endpointMatch) {
    console.log(`✅ Test ${index + 1}: PASS`);
    console.log(`   Query: "${test.query}"`);
    console.log(`   Intent: ${result?.intent || 'null'}`);
    console.log(`   Endpoint: ${result?.endpoint || 'null'}\n`);
    passed++;
  } else {
    console.log(`❌ Test ${index + 1}: FAIL`);
    console.log(`   Query: "${test.query}"`);
    console.log(`   Expected Intent: ${test.expected_intent}`);
    console.log(`   Got Intent: ${result?.intent || 'null'}`);
    console.log(`   Expected Endpoint: ${test.expected_endpoint}`);
    console.log(`   Got Endpoint: ${result?.endpoint || 'null'}`);
    if (test.note) console.log(`   Note: ${test.note}`);
    console.log();
    failed++;
  }
});

// Summary
console.log("═══════════════════════════════════════");
console.log(`SUMMARY: ${passed}/${TEST_CASES.length} tests passed`);
console.log(`Passed: ${passed}`);
console.log(`Failed: ${failed}`);
console.log("═══════════════════════════════════════");

if (failed === 0) {
  console.log("\n⚔️ ALL TESTS PASSED! Q2C ROUTER IS READY! ⚔️\n");
  process.exit(0);
} else {
  console.log("\n🔴 SOME TESTS FAILED. FIX INTENTS BEFORE PROCEEDING! 🔴\n");
  process.exit(1);
}






