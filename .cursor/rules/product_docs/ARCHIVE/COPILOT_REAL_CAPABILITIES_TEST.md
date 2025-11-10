# ⚔️ CO-PILOT REAL CAPABILITIES TEST ⚔️

**Date**: November 4, 2025  
**Status**: 🔴 **TESTING THE WRONG THING - CORRECTED**  
**Commander's Feedback**: "this copilot is meant to work for everything we worked on not the rag"

---

## 💥 **WHAT I WAS TESTING (WRONG)**

❌ RAG system with PubMed  
❌ Knowledge base population  
❌ Generic literature search  

**Problem**: That's just ONE feature. I ignored 90% of what we built!

---

## 🎯 **WHAT THE CO-PILOT SHOULD ACTUALLY DO**

Based on `@ayesha_plan.mdc` and our actual work, the Co-Pilot should orchestrate **ALL** these capabilities:

### **1. DRUG EFFICACY (WIWFM)** 💊
**What we built**: 
- Endpoint: `POST /api/efficacy/predict`
- S/P/E + SAE features
- Treatment line context
- Biomarker gates (HRD+, TMB, TP53)

**Co-Pilot should**:
- User: "Will Olaparib work for me?"
- Co-Pilot: Routes to `/api/efficacy/predict` with variant + biomarkers
- Response: Efficacy score, confidence, SAE features, evidence tier

**Test Cases**:
```
1. "Will Olaparib work for BRCA1 mutation?"
   → Efficacy prediction with HRD context

2. "Should I use platinum after progression?"
   → Cross-resistance detection via SAE

3. "What's my best option for L3?"
   → Treatment line context (L1→L2→L3)
```

---

### **2. FOOD/SUPPLEMENT VALIDATOR** 🥗
**What we built**:
- Endpoint: `POST /api/hypothesis/validate_food_dynamic`
- Works for ANY compound
- Evidence-backed dosage extraction
- Treatment line appropriateness
- Drug interaction checking

**Co-Pilot should**:
- User: "Should I take Vitamin D?"
- Co-Pilot: Routes to food validator with biomarkers + treatment history
- Response: Verdict, evidence grade, dosage, timing, interactions

**Test Cases**:
```
1. "Should I take Vitamin D after platinum?"
   → Post-platinum oxidative stress → NAC recommendation

2. "Can I take curcumin with my medications?"
   → Drug interaction check (warfarin warning)

3. "What foods help with HRD?"
   → Biomarker-targeted food recommendations
```

---

### **3. COMPLETE CARE PLAN (HOLISTIC)** 🎯
**What we built**:
- Endpoint: `POST /api/ayesha/complete_care_plan`
- Unified drug + food recommendations
- Shared biomarkers (HRD+, TMB, TP53)
- Integrated confidence scoring

**Co-Pilot should**:
- User: "What should I do?"
- Co-Pilot: Routes to unified care plan
- Response: Drugs ranked + Foods ranked + Integrated confidence

**Test Cases**:
```
1. "Give me a complete care plan"
   → Drug recommendations + Food recommendations

2. "What's my holistic treatment strategy?"
   → L3 context → Olaparib + NAC + Vitamin D

3. "What works best for my biomarkers?"
   → HRD+ → PARP inhibitor + DNA repair support foods
```

---

### **4. CLINICAL TRIALS MATCHING** 🏥
**What we built**:
- Endpoint: `POST /api/trials/agent/search`
- Autonomous agent extracts from patient data
- Hybrid AstraDB + Neo4j search
- Graph-optimized ranking (PI reputation, sites)

**Co-Pilot should**:
- User: "Find trials for me"
- Co-Pilot: Routes to trials agent with patient context
- Response: Recruiting trials ranked by relevance

**Test Cases**:
```
1. "Find ovarian cancer trials"
   → Autonomous extraction of criteria → Ranked trials

2. "Any trials for HRD+ patients?"
   → Biomarker-filtered trial search

3. "Trials near me in New York?"
   → Location + graph optimization
```

---

### **5. TOXICITY RISK (PGx)** ⚠️
**What we built**:
- Endpoint: `POST /api/safety/toxicity_risk`
- Pharmacogene detection (DPYD, TPMT, DPYD)
- Pathway overlap (drug MoA ∩ variant pathways)
- Cross-resistance warnings

**Co-Pilot should**:
- User: "Will this be toxic?"
- Co-Pilot: Routes to toxicity risk with germline + drug
- Response: Risk score, factors, warnings, alternatives

**Test Cases**:
```
1. "Is platinum safe for me?"
   → Germline DPYD variant → Warning

2. "Can I tolerate Olaparib after platinum?"
   → Cross-resistance check → SAE features

3. "What are my drug interaction risks?"
   → PGx + current medications → Risk profile
```

---

### **6. RADIATION GUIDANCE (PrecisionRad)** ☢️
**What we built**:
- Endpoint: `POST /api/guidance/radonc`
- Radiosensitivity prediction
- Chromatin accessibility context
- Evidence tier classification

**Co-Pilot should**:
- User: "Should I get radiation?"
- Co-Pilot: Routes to radonc guidance
- Response: Tier, radiosensitivity score, confidence, chromatin context

**Test Cases**:
```
1. "Is radiation right for me?"
   → Radiosensitivity prediction + chromatin accessibility

2. "Will radiation work for TP53 mutation?"
   → Functionality + chromatin → Guidance tier

3. "What's my radiation risk?"
   → Toxicity + radiosensitivity → Risk profile
```

---

### **7. CHEMO GUIDANCE (Tiered)** 💉
**What we built**:
- Endpoint: `POST /api/guidance/chemo`
- On-label rules + evidence gates
- Tier classification (I/II/III/Research)
- MoA-aware highlighting

**Co-Pilot should**:
- User: "Should I use platinum?"
- Co-Pilot: Routes to chemo guidance with drug class
- Response: Tier, on-label status, strength, evidence, citations

**Test Cases**:
```
1. "Should I use platinum for ovarian cancer?"
   → HRD+ → Tier I (Yes GO) with guideline support

2. "Is PARP inhibitor on-label?"
   → On-label check + evidence tier

3. "What chemo class fits my mutations?"
   → Pathway alignment → Therapy mapping
```

---

### **8. SYNTHETIC LETHALITY** 🎯
**What we built**:
- Endpoint: `POST /api/guidance/synthetic_lethality`
- Damage + Dependency → Therapy mapping
- HRD → platinum/PARP logic
- A-B dependency technique

**Co-Pilot should**:
- User: "What's my vulnerability?"
- Co-Pilot: Routes to synthetic lethality
- Response: Damage report, essentiality report, suggested therapy

**Test Cases**:
```
1. "What therapies target my BRCA1 mutation?"
   → HRD detection → Platinum/PARP recommendation

2. "Exploit my cancer's weakness"
   → Synthetic lethality mapping → Therapy class

3. "What's the A-B dependency here?"
   → Damage (BRCA1) + Dependency (PARP1) → Olaparib
```

---

### **9. VARIANT IMPACT (Deep Analysis)** 🔬
**What we built**:
- Endpoint: `POST /api/evidence/deep_analysis`
- ClinVar classification
- Literature synthesis
- Functional impact prediction

**Co-Pilot should**:
- User: "What is BRAF V600E?"
- Co-Pilot: Routes to deep analysis
- Response: ClinVar classification, literature summary, functional impact

**Test Cases**:
```
1. "What is BRAF V600E?"
   → ClinVar + literature + functional impact

2. "Is this mutation pathogenic?"
   → ClinVar review status + evidence level

3. "What does this variant do?"
   → Functional impact + affected pathways
```

---

### **10. TREATMENT LINE CONTEXT (SAE)** 📊
**What we built**:
- SAE features in food + drug
- Line appropriateness scoring
- Cross-resistance detection
- Sequencing fitness

**Co-Pilot should**:
- User: Automatically knows treatment history (L1, L2, L3)
- Co-Pilot: Adjusts ALL recommendations based on line
- Response: Line-appropriate drugs + foods

**Test Cases**:
```
1. Set context: L3, post-platinum
   → Drug recommendations flag cross-resistance
   → Food recommendations boost NAC (line_appropriateness: 1.0)

2. Set context: L1, treatment-naive
   → Platinum-based chemo recommended
   → Basic supportive foods

3. Set context: L2, post-PARP inhibitor
   → Alternative DDR drugs
   → Cross-resistance warnings
```

---

## 🎯 **Q2C ROUTER INTENTS TO TEST**

Based on `intents.js`, verify these routes work:

| Intent | Pattern | Endpoint | Test Query |
|--------|---------|----------|------------|
| `variant_impact` | "what impact" | `/api/evidence/deep_analysis` | "What is the impact of BRAF V600E?" |
| `drug_efficacy` | "will drug work" | `/api/efficacy/predict` | "Will Olaparib work for me?" |
| `radonc_guidance` | "radiation" | `/api/guidance/radonc` | "Should I get radiation?" |
| `chemo_guidance` | "chemo" | `/api/guidance/chemo` | "Should I use platinum?" |
| `literature_retrieval` | "find papers" | `/api/evidence/literature` | "Find papers on BRCA1" |
| `clinvar_context` | "clinvar" | `/api/evidence/deep_analysis` | "What does ClinVar say?" |
| `design_request` | "design guide" | `/api/design/guide_rna` | "Design a guide RNA" |
| **NEW: food_validator** | "should I take" | `/api/hypothesis/validate_food_dynamic` | "Should I take Vitamin D?" |
| **NEW: trials** | "find trials" | `/api/trials/agent/search` | "Find trials for me" |
| **NEW: complete_care** | "complete plan" | `/api/ayesha/complete_care_plan` | "Give me a complete care plan" |
| **NEW: synthetic_lethality** | "exploit weakness" | `/api/guidance/synthetic_lethality` | "What's my cancer's weakness?" |

---

## 🧪 **REAL TEST PLAN (COMPREHENSIVE)**

### **PHASE 1: Q2C Router Intent Classification**

**Test**: Verify Q2C Router correctly classifies all intent types

```javascript
// Test file: .cursor/tests/copilot/test_q2c_router.js

import { Q2C_ROUTER } from 'Q2CRouter/intents.js';

// Test 1: Drug Efficacy
const intent1 = Q2C_ROUTER.classifyIntent("Will Olaparib work for me?");
assert(intent1.intent === "drug_efficacy");
assert(intent1.endpoint === "/api/efficacy/predict");

// Test 2: Food Validator
const intent2 = Q2C_ROUTER.classifyIntent("Should I take Vitamin D?");
assert(intent2.intent === "food_validator");
assert(intent2.endpoint === "/api/hypothesis/validate_food_dynamic");

// Test 3: Complete Care
const intent3 = Q2C_ROUTER.classifyIntent("Give me a complete care plan");
assert(intent3.intent === "complete_care");
assert(intent3.endpoint === "/api/ayesha/complete_care_plan");

// Test 4: Clinical Trials
const intent4 = Q2C_ROUTER.classifyIntent("Find trials for me");
assert(intent4.intent === "trials");
assert(intent4.endpoint === "/api/trials/agent/search");

// Test 5: RAG Fallback
const intent5 = Q2C_ROUTER.classifyIntent("Explain how cancer works");
assert(intent5.intent === null); // Falls back to RAG
```

---

### **PHASE 2: Backend Endpoint Integration**

**Test**: Verify all endpoints work with realistic payloads

#### **Test 2.1: Drug Efficacy**
```bash
curl -X POST http://127.0.0.1:8000/api/efficacy/predict \
  -d '{
    "mutations": [{"gene":"BRCA1","hgvs_p":"p.Cys61Gly"}],
    "disease": "ovarian_cancer",
    "treatment_history": [{"line": 2, "drugs": ["carboplatin", "paclitaxel"], "outcome": "progression"}],
    "biomarkers": {"HRD": "POSITIVE", "TMB": 8.2},
    "include_sae_features": true
  }'

# Expected:
{
  "drugs": [
    {
      "name": "Olaparib",
      "efficacy_score": 0.85,
      "confidence": 0.78,
      "evidence_tier": "SUPPORTED",
      "sae_features": {
        "line_appropriateness": 0.9,
        "cross_resistance": 0.2,
        "sequencing_fitness": 0.85
      }
    }
  ]
}
```

#### **Test 2.2: Food Validator**
```bash
curl -X POST http://127.0.0.1:8000/api/hypothesis/validate_food_dynamic \
  -d '{
    "compound": "Vitamin D",
    "disease_context": {"disease": "ovarian_cancer_hgs", "biomarkers": {"HRD": "POSITIVE"}},
    "treatment_history": {"current_line": "L3", "prior_therapies": ["carboplatin", "paclitaxel"]},
    "patient_medications": ["warfarin"],
    "use_llm": true
  }'

# Expected:
{
  "verdict": "SUPPORTED",
  "overall_score": 0.82,
  "evidence": {"grade": "STRONG", "papers": [...]},
  "recommendations": {
    "dosage": "2000-4000 IU daily",
    "timing": "Morning with breakfast",
    "interactions": ["Monitor INR with warfarin"]
  }
}
```

#### **Test 2.3: Complete Care Plan**
```bash
curl -X POST http://127.0.0.1:8000/api/ayesha/complete_care_plan \
  -d '{
    "patient_context": {
      "disease": "ovarian_cancer",
      "mutations": [{"gene":"BRCA1","hgvs_p":"p.Cys61Gly"}],
      "biomarkers": {"HRD": "POSITIVE", "TMB": 8.2},
      "treatment_history": {"line": 3, "prior_therapies": ["carboplatin", "olaparib"]}
    }
  }'

# Expected:
{
  "drug_recommendations": [...],
  "food_recommendations": [...],
  "integrated_confidence": 0.75,
  "provenance": {"run_id": "...", "drug_weight": 0.7, "food_weight": 0.3}
}
```

#### **Test 2.4: Clinical Trials**
```bash
curl -X POST http://127.0.0.1:8000/api/trials/agent/search \
  -d '{
    "patient_summary": "55yo female, ovarian cancer, HRD+, post-platinum progression"
  }'

# Expected:
{
  "trials": [
    {"nct_id": "NCT...", "title": "...", "relevance": 0.9, "status": "RECRUITING"}
  ]
}
```

---

### **PHASE 3: Frontend Co-Pilot UI**

**Test**: Real user interactions with Co-Pilot interface

#### **Test 3.1: Open Co-Pilot & Context Detection**
**Steps**:
1. Start frontend: `npm run dev`
2. Navigate to Ayesha Complete Care page
3. Set patient context:
   - Disease: Ovarian Cancer
   - Mutation: BRCA1 p.Cys61Gly
   - Treatment History: L3, post-platinum
4. Open Co-Pilot drawer

**Expected**:
- ✅ Co-Pilot knows disease = ovarian cancer
- ✅ Co-Pilot knows mutation = BRCA1
- ✅ Co-Pilot knows treatment line = L3
- ✅ Suggested questions appear

#### **Test 3.2: Drug Efficacy Question**
**Steps**:
1. Type: "Will Olaparib work for me?"
2. Send message

**Expected**:
- ✅ Q2C Router classifies as "drug_efficacy"
- ✅ Routes to `/api/efficacy/predict`
- ✅ Passes treatment history + biomarkers
- ✅ Response shows:
  - Efficacy score
  - Confidence
  - SAE features (line_appropriateness, cross_resistance)
  - Evidence tier
  - Citations

#### **Test 3.3: Food Validator Question**
**Steps**:
1. Type: "Should I take Vitamin D?"
2. Send message

**Expected**:
- ✅ Q2C Router classifies as "food_validator"
- ✅ Routes to `/api/hypothesis/validate_food_dynamic`
- ✅ Passes treatment history + biomarkers
- ✅ Response shows:
  - Verdict (SUPPORTED/WEAK_SUPPORT/NOT_SUPPORTED)
  - Evidence grade
  - Dosage recommendation
  - Timing recommendation
  - Drug interactions warning (warfarin)

#### **Test 3.4: Complete Care Plan**
**Steps**:
1. Type: "Give me a complete care plan"
2. Send message

**Expected**:
- ✅ Q2C Router classifies as "complete_care"
- ✅ Routes to `/api/ayesha/complete_care_plan`
- ✅ Response shows:
  - Drug recommendations ranked
  - Food recommendations ranked
  - Integrated confidence
  - Provenance (drug_weight + food_weight)

#### **Test 3.5: Clinical Trials**
**Steps**:
1. Type: "Find trials for me"
2. Send message

**Expected**:
- ✅ Q2C Router classifies as "trials"
- ✅ Routes to `/api/trials/agent/search`
- ✅ Autonomous agent extracts criteria from context
- ✅ Response shows:
  - Recruiting trials
  - Relevance scores
  - Contact information
  - Location

#### **Test 3.6: Multi-Turn Conversation**
**Steps**:
1. Ask: "Will Olaparib work for me?"
2. Wait for response
3. Ask: "What foods support PARP inhibitor therapy?"
4. Wait for response

**Expected**:
- ✅ First query: Drug efficacy prediction
- ✅ Second query: Food recommendations
- ✅ Context maintained (BRCA1, L3, post-platinum)
- ✅ Food recommendations aware of PARP inhibitor context

---

### **PHASE 4: Treatment Line Context Awareness**

**Test**: Verify treatment history propagates correctly

#### **Test 4.1: L1 (Treatment-Naive)**
**Context**: L1, no prior therapies
**Expected**:
- Drug recommendations: Platinum-based chemo
- Food recommendations: Basic supportive care
- SAE line_appropriateness: High for first-line drugs

#### **Test 4.2: L2 (Post-Platinum)**
**Context**: L2, post-carboplatin progression
**Expected**:
- Drug recommendations: PARP inhibitors
- Food recommendations: NAC (oxidative stress), Vitamin D
- SAE cross_resistance: Flags platinum resistance

#### **Test 4.3: L3 (Post-PARP)**
**Context**: L3, post-olaparib progression
**Expected**:
- Drug recommendations: Alternative DDR drugs
- Food recommendations: NAC (line_appropriateness: 1.0)
- SAE cross_resistance: Flags PARP resistance

---

### **PHASE 5: Biomarker-Aware Recommendations**

**Test**: Verify biomarkers affect recommendations

#### **Test 5.1: HRD+ Context**
**Biomarkers**: HRD = POSITIVE
**Expected**:
- Drug: PARP inhibitors ranked high
- Food: Vitamin D (DNA repair support)
- Evidence: HRD-specific citations

#### **Test 5.2: TMB-High Context**
**Biomarkers**: TMB = 15
**Expected**:
- Drug: Immunotherapy ranked high
- Food: Inflammation-targeting foods
- Evidence: TMB-specific citations

#### **Test 5.3: TP53 Mutation Context**
**Mutations**: TP53 p.R175H
**Expected**:
- Drug: TP53-targeting therapies
- Food: Apoptosis-supporting compounds
- Evidence: TP53-specific citations

---

## 🎯 **SUCCESS CRITERIA (REAL)**

| Capability | Current Status | Target |
|------------|----------------|--------|
| Q2C Router classifies all intents | ❓ Not tested | ✅ 100% accuracy |
| Drug efficacy with SAE | ❓ Not tested | ✅ Tested & verified |
| Food validator with evidence | ❓ Not tested | ✅ Tested & verified |
| Complete care plan | ❓ Not tested | ✅ Tested & verified |
| Clinical trials matching | ❓ Not tested | ✅ Tested & verified |
| Treatment line context | ❓ Not tested | ✅ Propagates correctly |
| Biomarker awareness | ❓ Not tested | ✅ Affects recommendations |
| Multi-turn conversation | ❓ Not tested | ✅ Context maintained |
| Frontend UI working | ❓ Not tested | ✅ All flows verified |

---

## 💥 **WHAT I NEED TO DO NOW**

### **Immediate (P0)**
1. ✅ Admit I was testing the wrong thing (RAG only)
2. ⏳ Add missing Q2C intents: `food_validator`, `trials`, `complete_care`, `synthetic_lethality`
3. ⏳ Test Q2C Router with all 10+ intent types
4. ⏳ Start frontend and test real Co-Pilot UI
5. ⏳ Test treatment line context propagation
6. ⏳ Test biomarker-aware recommendations

### **Short-term (P1)**
1. ⏳ Test multi-turn conversations
2. ⏳ Test all 10 use-cases end-to-end
3. ⏳ Document real capabilities (not just RAG)
4. ⏳ Create comprehensive test suite

### **Medium-term (P2)**
1. ⏳ Automated testing for all intents
2. ⏳ Performance benchmarks
3. ⏳ Load testing

---

## 💥 **COMMANDER'S VERDICT**

**Previous Tests**: ❌ **ONLY TESTED RAG (1 of 10 capabilities)**  
**What I Missed**: Drug efficacy, Food validator, Clinical trials, Treatment line context, Toxicity risk, Complete care, Synthetic lethality, Guidance, SAE features  
**Required**: Test ALL 10 capabilities the Co-Pilot orchestrates  

**Status**: 🔴 **COMPLETELY INADEQUATE TESTING** 🔴

---

**The Co-Pilot isn't just a RAG chatbot. It's an orchestration layer for 10+ distinct clinical capabilities. I need to test ALL of them, not just PubMed search!** ⚔️

