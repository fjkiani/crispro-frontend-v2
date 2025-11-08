# ⚔️ DYNAMIC FOOD VALIDATOR FRAMEWORK - COMPLETE BUILD DOCTRINE

**Mission:** Build end-to-end framework that works for ANY food/supplement, not just 6 hardcoded compounds

**Target Users:** Patients, Dieticians, Oncology Care Teams

**Timeline:** Full implementation (backup + frontend + testing)

**Status:** ✅ **READY TO BUILD**

---

## **🎯 ARCHITECTURE OVERVIEW**

### **Complete Data Flow:**

```
User Query: "Can [ANY FOOD] help [CANCER TYPE]?"
    ↓
[1] FAST CACHE CHECK (<100ms)
    → Check validated_claims_cache.json
    → If found → Return instantly
    ↓
[2] DYNAMIC TARGET EXTRACTION (2-5s)
    → ChEMBL API: Extract molecular targets
    → PubChem: Get compound structure/metadata
    → LLM: Extract targets from literature if APIs fail
    ↓
[3] PATHWAY MAPPING (1-2s)
    → Map targets → Cancer pathways:
      • Angiogenesis (VEGF, EGFR, PDGFR)
      • DNA Repair (BRCA1/2, PARP, TP53)
      • Inflammation (NF-κB, COX-2, IL-6)
      • Cell Cycle (CDK, Cyclin D1, p21)
      • Apoptosis (Bcl-2, Caspase, p53)
      • Metabolism (mTOR, PI3K, AKT)
    ↓
[4] EVIDENCE MINING (5-10s)
    → PubMed: Search compound + disease + pathway terms
    → LLM Synthesis: Extract mechanisms, doses, safety
    → Evidence Grading: STRONG/MODERATE/WEAK/INSUFFICIENT
    ↓
[5] S/P/E + SAE SCORING (1-2s)
    → Sequence (S): Evo2 if enabled, else neutral 0.5
    → Pathway (P): Alignment score (compound pathways ∩ disease pathways)
    → Evidence (E): Literature grade → 0-1 score
    → SAE: Treatment line + biomarker gating
    → Overall: 0.4×S + 0.3×P + 0.3×E
    ↓
[6] DIETICIAN RECOMMENDATIONS (LLM, 3-5s)
    → Dosage: Extract from literature + safety thresholds
    → Timing: Best time of day, with/without food
    → Interactions: Drug-food interactions from literature
    → Meal Planning: Foods to combine/avoid
    → Monitoring: What labs to track
    ↓
[7] CACHE STORAGE (background)
    → Store result in validated_claims_cache.json
    → TTL: 7 days (revalidate weekly)
    → Next query gets instant response
    ↓
[8] RETURN COMPLETE RESULT
```

---

## **📋 IMPLEMENTATION PLAN**

### **Phase 1: Dynamic Target Extraction Service (2 hours)**

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/dynamic_food_extraction.py`

**Capabilities:**
1. ChEMBL API integration (primary)
2. PubChem API fallback (compound metadata)
3. LLM literature extraction (backup if APIs fail)
4. Target → Pathway mapping
5. Cancer mechanism classification

---

### **Phase 2: Cancer Pathway Intelligence (1 hour)**

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/cancer_pathway_mapper.py`

**Pathway Database:**
```json
{
  "cancer_mechanisms": {
    "angiogenesis": {
      "targets": ["VEGF", "VEGFR", "EGFR", "PDGFR", "FGF", "Angiopoietin"],
      "pathways": ["VEGF signaling", "Angiogenesis", "Blood vessel formation"],
      "diseases": ["ovarian_cancer", "breast_cancer", "colon_cancer"],
      "intervention_types": ["anti-angiogenic", "VEGFR inhibitor", "tyrosine kinase inhibitor"]
    },
    "dna_repair": {
      "targets": ["BRCA1", "BRCA2", "PARP1", "TP53", "ATM", "ATR", "RAD51"],
      "pathways": ["Homologous recombination", "DNA repair", "Cell cycle checkpoint"],
      "diseases": ["ovarian_cancer", "breast_cancer", "pancreatic_cancer"],
      "intervention_types": ["PARP inhibitor", "DNA repair enhancer", "HRD-targeting"]
    },
    "inflammation": {
      "targets": ["NF-κB", "COX-2", "IL-6", "TNF-α", "STAT3", "iNOS"],
      "pathways": ["Inflammatory response", "NF-κB signaling", "JAK-STAT"],
      "diseases": ["ovarian_cancer", "colon_cancer", "prostate_cancer"],
      "intervention_types": ["anti-inflammatory", "NF-κB inhibitor", "COX-2 inhibitor"]
    },
    "cell_cycle": {
      "targets": ["CDK4", "CDK6", "Cyclin D1", "p21", "p27", "RB1"],
      "pathways": ["Cell cycle progression", "G1/S checkpoint", "CDK signaling"],
      "diseases": ["breast_cancer", "ovarian_cancer", "lung_cancer"],
      "intervention_types": ["CDK inhibitor", "cell cycle modulator"]
    },
    "apoptosis": {
      "targets": ["Bcl-2", "Bax", "Caspase-3", "p53", "Survivin", "XIAP"],
      "pathways": ["Apoptosis", "Programmed cell death", "p53 signaling"],
      "diseases": ["ovarian_cancer", "breast_cancer", "leukemia"],
      "intervention_types": ["pro-apoptotic", "Bcl-2 inhibitor", "p53 activator"]
    },
    "metabolism": {
      "targets": ["mTOR", "PI3K", "AKT", "GLUT1", "HK2", "LDH"],
      "pathways": ["mTOR signaling", "PI3K/AKT", "Glycolysis", "Warburg effect"],
      "diseases": ["ovarian_cancer", "breast_cancer", "pancreatic_cancer"],
      "intervention_types": ["mTOR inhibitor", "metabolic modulator", "glycolysis inhibitor"]
    }
  }
}
```

---

### **Phase 3: Enhanced Evidence Service (2 hours)**

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/enhanced_evidence_service.py`

**Features:**
1. Multi-source literature search (PubMed, OpenAlex, Semantic Scholar)
2. LLM synthesis of mechanisms
3. Dosage extraction from clinical trials
4. Safety/contraindication detection
5. Drug-food interaction checking
6. Evidence grading (RCT → STRONG, Observational → MODERATE, etc.)

---

### **Phase 4: Dietician Recommendations Engine (2 hours)**

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/dietician_recommendations.py`

**Outputs:**
- Dosage recommendations (with source citations)
- Timing guidance (best time of day, with/without food)
- Meal planning suggestions (foods to combine/avoid)
- Drug interactions (check against patient medication list)
- Lab monitoring (what to track: serum levels, liver function, etc.)
- Safety alerts (contraindications, precautions)

---

### **Phase 5: Unified Endpoint (1 hour)**

**File:** `oncology-coPilot/oncology-backend-minimal/api/routers/hypothesis_validator.py`

**New Endpoint:** `POST /api/hypothesis/validate_food_dynamic`

**Features:**
- Works for ANY food/supplement
- Fast cache check first
- Dynamic extraction if not cached
- Complete S/P/E + SAE scoring
- Dietician recommendations
- Caching for future queries

---

### **Phase 6: Frontend Enhancement (2 hours)**

**File:** `oncology-coPilot/oncology-frontend/src/pages/FoodValidatorDynamic.jsx`

**Features:**
- Free-text compound input (not just dropdown)
- Real-time validation feedback
- Complete results display:
  - Verdict + confidence
  - S/P/E breakdown
  - Target → Pathway mapping
  - Evidence cards with citations
  - Dietician recommendations
  - Safety alerts
- Save to session / export to PDF

---

## **🔬 DETAILED IMPLEMENTATION**

Let me build each component now...

