# ⚔️ CO-PILOT REAL TEST PLAN ⚔️

**Date**: November 4, 2025  
**Status**: ⚠️ **PREVIOUS TESTS WERE INADEQUATE**  
**Commander's Feedback**: "none of these tests look like actual copilot tests"

---

## 💥 **WHAT WAS WRONG WITH PREVIOUS TESTS**

### **Previous "Tests" Were Just API Smoke Tests**
❌ Called `/api/evidence/rag-query` with curl  
❌ Got back "I couldn't find literature" (FAILURE, not success)  
❌ Knowledge base empty (0 papers)  
❌ Never tested frontend UI  
❌ Never tested Q2C Router  
❌ Never tested treatment history context  
❌ Never tested multi-turn conversations  

**Verdict**: I tested that the endpoint **exists**, not that the Co-Pilot **works**.

---

## 🎯 **REAL CO-PILOT TEST PLAN**

### **PREREQUISITE: Fix Blocking Issues**

#### **P0 BLOCKER: Populate Knowledge Base**
**Issue**: Empty KB means LLM can't answer questions  
**Fix**:
```bash
# Manual population (bypass async bug)
cd oncology-coPilot/oncology-backend-minimal/Pubmed-LLM-Agent-main
python3 rag_agent.py --add-variant BRAF V600E melanoma
python3 rag_agent.py --add-variant KRAS G12C colorectal_cancer
python3 rag_agent.py --add-variant TP53 R175H breast_cancer
```

**Expected**: 100+ papers added to KB

#### **P1 BLOCKER: Fix Async Event Loop**
**Issue**: `/rag-add-variant` endpoint fails  
**Impact**: Can't populate KB via API  
**Fix**: Refactor `PubMedClientEnhanced` to handle nested async contexts

---

## 🧪 **TEST SUITE (REAL CO-PILOT TESTING)**

### **PHASE 1: Knowledge Base Validation**

#### **Test 1.1: Populate KB Manually**
```bash
cd oncology-coPilot/oncology-backend-minimal/Pubmed-LLM-Agent-main
python3 rag_agent.py --add-variant BRAF V600E melanoma

# Expected output:
# ✅ Found 50+ papers
# ✅ Added to knowledge base
# ✅ Embeddings generated
```

#### **Test 1.2: Verify KB Stats**
```bash
curl -sS http://127.0.0.1:8000/api/evidence/rag-stats

# Expected:
{
  "total_papers": 50,
  "genes": ["BRAF"],
  "diseases": ["melanoma"]
}
```

#### **Test 1.3: RAG Query with Populated KB**
```bash
curl -X POST http://127.0.0.1:8000/api/evidence/rag-query \
  -d '{"query": "What is BRAF V600E?", "gene": "BRAF", "hgvs_p": "V600E"}'

# Expected:
{
  "answer": "BRAF V600E is a somatic mutation...",  # REAL ANSWER, not "I couldn't find"
  "evidence_level": "Strong",                       # Not "Insufficient"
  "confidence_score": 0.85,                         # Not 0.0
  "supporting_papers": [                            # Not []
    {"pmid": "12345", "title": "BRAF mutations in melanoma"}
  ]
}
```

---

### **PHASE 2: Backend Integration Tests**

#### **Test 2.1: Q2C Intent Classification**
**Test**: Verify Q2C Router classifies intents correctly

```javascript
// In browser console or test file
import { Q2C_ROUTER } from './Q2CRouter';

// Test 1: Variant Impact
const intent1 = Q2C_ROUTER.classifyIntent("What is the functional impact of BRAF V600E?");
// Expected: intent = "variant_impact", endpoint = "/api/evidence/deep_analysis"

// Test 2: Drug Efficacy
const intent2 = Q2C_ROUTER.classifyIntent("Will Vemurafenib work for me?");
// Expected: intent = "drug_efficacy", endpoint = "/api/efficacy/predict"

// Test 3: RAG Fallback
const intent3 = Q2C_ROUTER.classifyIntent("Explain how this mutation affects cells");
// Expected: intent = null (fallback to RAG)
```

#### **Test 2.2: Treatment History Context**
**Test**: Verify treatment history is passed to backend

```bash
curl -X POST http://127.0.0.1:8000/api/efficacy/predict \
  -d '{
    "mutations": [{"gene": "BRAF", "hgvs_p": "V600E"}],
    "disease": "melanoma",
    "treatment_history": [
      {"line": 1, "drugs": ["Pembrolizumab"], "outcome": "progression"}
    ]
  }'

# Expected:
{
  "drugs": [...],
  "sae_features": {
    "line_appropriateness": {...},  # Uses treatment history
    "cross_resistance": {...}       # Detects prior therapy resistance
  }
}
```

---

### **PHASE 3: Frontend Co-Pilot UI Tests**

#### **Test 3.1: Start Frontend**
```bash
cd oncology-coPilot/oncology-frontend
npm run dev

# Expected:
# ✅ Frontend starts on http://localhost:5173
# ✅ No console errors
```

#### **Test 3.2: Open Co-Pilot Drawer**
**Steps**:
1. Navigate to any page
2. Look for Co-Pilot icon/button (usually right side)
3. Click to open drawer

**Expected**:
- ✅ Drawer opens
- ✅ Welcome message appears
- ✅ Input field visible
- ✅ Suggested prompts shown

#### **Test 3.3: Ask Variant Impact Question**
**Steps**:
1. Set variant context: BRAF V600E, melanoma
2. Type: "What is the functional impact of this mutation?"
3. Send message

**Expected**:
- ✅ Q2C Router classifies as "variant_impact"
- ✅ Routes to `/api/evidence/deep_analysis`
- ✅ ClinVar classification shown
- ✅ Response displays in <5s

#### **Test 3.4: Ask RAG Question (Fallback)**
**Steps**:
1. Set variant context: BRAF V600E, melanoma
2. Type: "Explain how this mutation leads to cancer"
3. Send message

**Expected**:
- ✅ Q2C Router: no intent match → fallback to RAG
- ✅ Routes to `/api/evidence/rag-query`
- ✅ LLM answer with real content (not "I couldn't find")
- ✅ Supporting papers shown with PMIDs
- ✅ Evidence level: "Strong"
- ✅ Confidence: >0.8

#### **Test 3.5: Ask Drug Efficacy Question**
**Steps**:
1. Set variant: BRAF V600E, melanoma
2. Set treatment history: L1 Pembrolizumab (progression)
3. Type: "Will Vemurafenib work for me?"
4. Send message

**Expected**:
- ✅ Q2C Router classifies as "drug_efficacy"
- ✅ Routes to `/api/efficacy/predict`
- ✅ Treatment history passed in payload
- ✅ SAE features computed (line_appropriateness, cross_resistance)
- ✅ Drug ranking shown with confidence
- ✅ Response displays in <10s

#### **Test 3.6: Multi-Turn Conversation**
**Steps**:
1. Ask: "What is BRAF V600E?"
2. Wait for response
3. Ask: "What treatments target this mutation?"
4. Wait for response

**Expected**:
- ✅ First query: RAG response about BRAF
- ✅ Second query: Uses context from first query
- ✅ Conversation history maintained
- ✅ Proactive suggestions appear

---

### **PHASE 4: Context Awareness Tests**

#### **Test 4.1: Variant Context Propagation**
**Steps**:
1. Navigate to Myeloma Digital Twin page
2. Select variant: BRAF V600E
3. Open Co-Pilot
4. Ask: "What's the impact?"

**Expected**:
- ✅ Co-Pilot knows current variant (BRAF V600E)
- ✅ Co-Pilot knows current disease (melanoma)
- ✅ Query uses this context automatically

#### **Test 4.2: Treatment History Context**
**Steps**:
1. Set treatment history: L1 Carboplatin+Paclitaxel, L2 Olaparib
2. Open Co-Pilot
3. Ask: "What should I take next?"

**Expected**:
- ✅ Co-Pilot knows treatment history
- ✅ Recommendations reflect prior therapies
- ✅ Cross-resistance warnings shown

---

### **PHASE 5: Error Handling & Edge Cases**

#### **Test 5.1: Empty Variant Context**
**Steps**:
1. Clear all variant context
2. Ask: "What is BRAF V600E?"

**Expected**:
- ✅ Co-Pilot still responds
- ✅ Extracts gene from query
- ✅ Graceful fallback

#### **Test 5.2: Backend Unavailable**
**Steps**:
1. Stop backend server
2. Ask question in Co-Pilot

**Expected**:
- ✅ Error message shown (not crash)
- ✅ Retry button appears
- ✅ Conversation history preserved

#### **Test 5.3: Empty Knowledge Base**
**Steps**:
1. Clear knowledge base
2. Ask: "What is NOVEL_GENE ABC123?"

**Expected**:
- ✅ LLM responds with "couldn't find literature"
- ✅ Suggests alternative actions
- ✅ Confidence: 0.0
- ✅ No crash

---

## 🎯 **SUCCESS CRITERIA (REAL)**

| Test | Current Status | Target |
|------|----------------|--------|
| KB populated | ❌ 0 papers | ✅ 100+ papers |
| RAG answers with content | ❌ "I couldn't find" | ✅ Real answers |
| Evidence level | ❌ "Insufficient" | ✅ "Strong" |
| Confidence score | ❌ 0.0 | ✅ >0.8 |
| Supporting papers | ❌ [] | ✅ 3-5 papers |
| Frontend UI working | ❓ Not tested | ✅ Tested & verified |
| Q2C Router tested | ❓ Not tested | ✅ All intents verified |
| Treatment history | ❓ Not tested | ✅ Context propagates |
| Multi-turn conversation | ❓ Not tested | ✅ History maintained |

---

## 🔥 **WHAT I NEED TO DO NOW**

### **Immediate (P0)**
1. ✅ Admit previous tests were inadequate
2. ⏳ Populate knowledge base manually (100+ papers)
3. ⏳ Retest RAG endpoint with populated KB
4. ⏳ Start frontend and test UI
5. ⏳ Test Q2C Router classifications
6. ⏳ Test multi-turn conversations

### **Short-term (P1)**
1. ⏳ Fix async event loop bug
2. ⏳ Test treatment history context
3. ⏳ Test error handling
4. ⏳ Document real test results

### **Medium-term (P2)**
1. ⏳ Automated test suite
2. ⏳ Performance benchmarks
3. ⏳ Load testing

---

## 💥 **COMMANDER'S VERDICT**

**Previous Tests**: ❌ **INADEQUATE**  
**Reason**: Only tested that endpoints exist, not that Co-Pilot works  
**Required**: Real UI testing with populated knowledge base  

**Status**: 🔴 **NOT ACTUALLY TESTED** 🔴

---

**I need to stop pretending curl tests = real Co-Pilot tests and do the actual work!** ⚔️

