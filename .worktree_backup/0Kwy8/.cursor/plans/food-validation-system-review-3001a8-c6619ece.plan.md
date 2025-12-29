---
name: Food Validation Product - Cancer-Specific Treatment Line Intelligence
overview: ""
todos:
  - id: b119b8f9-80d8-463e-aba9-ade7eb572889
    content: supplement_treatment_rules.json EXISTS (verified at `.cursor/ayesha/hypothesis_validator/data/`)
    status: completed
  - id: 83532c6d-19d1-4c66-a478-788e95945a1d
    content: "**CRITICAL**: Add SAE structure adapter in hypothesis_validator.py:798-858 (copy from validate_food_ab_enhanced line 435-446)"
    status: pending
  - id: 09a76373-af5f-4de7-b4d0-9b7d607e97e4
    content: Add treatment line format normalization in food_treatment_line_service.py (all formats → "L1"/"L2"/"L3")
    status: pending
  - id: abe7eeb2-6d9d-432f-8ee0-7cc1c0241b4c
    content: Create cancer_type_food_recommendations.json mapping cancer types (ovarian_cancer_hgs, breast_cancer, etc.) to pathways and recommended foods with treatment line filters
    status: pending
  - id: 9ccada90-2527-4dbe-8fa5-9ccd8224a3d2
    content: Treatment line affects confidence only (by design) - verified in code (food_spe_integration.py:233-237, 458-463)
    status: completed
  - id: 367dc7fc-bd61-41d1-9a7d-de082a00c17c
    content: SAE boost is properly applied in food_spe_integration.py _compute_confidence (line_app + seq_fit) × 0.05 - verified
    status: completed
  - id: 1b1ae713-6fec-4912-a7d2-aa4da9c134d9
    content: Update FoodRankingPanel.jsx to display treatment line context (L1, L2, L3), cancer type, biomarker matches (HRD+ match)
    status: pending
  - id: 00d005d8-5719-4455-8cde-f82ce8c574cf
    content: "Fix SAE structure adapter so treatment line intelligence accordion displays correctly (nested structure) - **CRITICAL**: This is the same as Priority 1 #1"
    status: pending
  - id: 2010879e-285f-4fa1-a1fd-f291f57af936
    content: Create biomarker_food_mapping.json at `.cursor/ayesha/hypothesis_validator/data/` mapping HRD_POSITIVE, TMB_HIGH, MSI_HIGH to target pathways and recommended foods with mechanisms
    status: pending
  - id: 083d19ff-695a-4601-9d7a-c8ea827f9551
    content: Update enhanced_evidence_service.py get_complete_evidence to filter evidence by treatment line context (first-line vs maintenance)
    status: pending
  - id: f019b639-c58b-46e3-a9cd-60f3ad82c1cc
    content: "Test first-line chemotherapy for ovarian cancer (HGS) scenario: verify foods target TP53/DNA repair, have high line_appropriateness, show HRD+ match in FoodRankingPanel"
    status: pending
  - id: 56981c8b-a6bd-4279-aaf3-89cfee83a6e2
    content: "Verify SAE structure adapter works: flat structure from service → nested structure in frontend (test with validate_food_dynamic endpoint)"
    status: pending
  - id: fc94dbf4-c108-4c68-9697-b19d3b606424
    content: "Verify treatment line format normalization: \"first-line\", \"L1\", \"frontline\" → all normalize to \"L1\" (test all format variations)"
    status: pending
  - id: 4dcf6e7f-f42e-4a14-8469-3160ae9cf419
    content: "Test batch recommendation via orchestrator: verify ayesha_orchestrator.py call_food_validator can recommend multiple foods"
    status: pending
---

# Food Validation Product - Cancer-Specific Treatment Line Intelligence

**Product Vision**: Build a precision nutrition system that recommends foods/supplements based on cancer type, treatment line, and patient biomarkers, validated through SPE framework and displayed to patients via FoodRankingPanel.

**Review Methodology**: Read code, trace execution, map integration, document state, identify gaps, validate understanding with code references.

---

## 1. PRODUCT CAPABILITIES (Not Just Code Review)

### 1.1 Cancer Type-Specific Food Recommendations

**Current State**:

- Universal pathway database exists: `api/resources/universal_disease_pathway_database.json`
- TCGA-weighted pathways per disease (e.g., ovarian_cancer_hgs, breast_cancer)
- Code: `food_spe_integration.py:67-90` (_get_disease_pathway_weights)

**Gap**: No cancer type-specific food library or recommendations

- Need: Cancer type → recommended foods mapping
- Example: Ovarian cancer (HGS) → foods targeting TP53, HRD/DDR pathways
- Example: Breast cancer → foods targeting HER2, ER/PR pathways

**Action Required**:

- Create cancer type-specific food recommendation database
- Map foods to cancer-specific pathways (using TCGA weights)
- Integrate with existing pathway alignment logic

### 1.2 Treatment Line Intelligence Integration

**Current State** (Code-Verified):

- Treatment line service exists: `food_treatment_line_service.py:33-141`
- Computes: line_appropriateness, cross_resistance, sequencing_fitness
- Data source: `supplement_treatment_rules.json` ✅ **EXISTS** at `.cursor/ayesha/hypothesis_validator/data/` (22+ compounds verified)
- Code: `food_spe_integration.py:239-250` (SAE features integration)
- Code: `food_spe_integration.py:458-463` (SAE boost affects confidence only, not ranking)

**Gaps Identified**:

1. **SAE Structure Mismatch (CRITICAL)**: 

   - `validate_food_dynamic` (line 787) returns flat structure: `{"line_appropriateness": 0.9, ...}`
   - `validate_food_ab_enhanced` (line 435-446) ✅ HAS adapter that transforms flat → nested
   - Frontend expects nested: `{line_fitness: {score, status, reason}, ...}`
   - **Impact**: Frontend cannot display treatment line intelligence correctly

2. **Treatment Line Format Inconsistency**:

   - `validate_food_dynamic` (line 654) expects `"current_line": "L3"` format
   - `ayesha_orchestrator.py` (line 302) uses `"L{treatment_line}"` format
   - No normalization layer exists
   - **Impact**: Format confusion across endpoints

3. **Treatment Line in Ranking (By Design)**:

   - Code: `food_spe_integration.py:233-237` - `overall_score = 0.4×S + 0.3×P + 0.3×E` (NO treatment line)
   - Code: `food_spe_integration.py:458-463` - SAE boost only affects confidence: `sae_boost = (line_app + seq_fit) × 0.05`
   - **Design Decision**: Treatment line appropriateness boosts confidence, not ranking order
   - **Action**: Document this design (no code change needed)

**Action Required**:

- ✅ File exists - no creation needed
- **CRITICAL**: Add SAE structure adapter in `validate_food_dynamic` (copy from `validate_food_ab_enhanced`)
- Add treatment line format normalization (all formats → "L1"/"L2"/"L3")
- Document design: Treatment line affects confidence boost, not primary ranking

### 1.3 Biomarker Targeting

**Current State**:

- Biomarker gates exist: `food_treatment_line_service.py:102-122`
- Checks: HRD, TMB, MSI, DNA repair pathways
- Code: `food_spe_integration.py:465-476` (biomarker boosts)

**Gap**: Biomarker-specific food recommendations not systematic

- Need: HRD+ → foods targeting DNA repair
- Need: TMB-high → foods targeting immune pathways
- Need: MSI-high → foods targeting mismatch repair

**Action Required**:

- Create biomarker → food mapping database
- Enhance biomarker gates in treatment line service
- Add biomarker-specific rationale to patient display

### 1.4 Patient-Facing Display (FoodRankingPanel)

**Current State**:

- Component exists: `oncology-coPilot/oncology-frontend/src/components/ayesha/FoodRankingPanel.jsx`
- Displays: compound, pathways, efficacy_score, confidence, sae_features, dosage, rationale, citations
- Code: Lines 47-189 (rendering logic)

**Gap**: Missing treatment line context and cancer type display

- Need: Show "Recommended for: First-line chemotherapy, Ovarian Cancer (HGS)"
- Need: Show biomarker match (e.g., "HRD+ match")
- Need: Show treatment line intelligence breakdown

**Action Required**:

- Enhance FoodRankingPanel to display treatment line context
- Add cancer type and biomarker match indicators
- Show treatment line intelligence (line_appropriateness, cross_resistance, sequencing_fitness) prominently

---

## 2. EXECUTION FLOW (Enhanced for Product)

### 2.1 Patient Input → Food Recommendations

**Step 1: Patient Context Collection**

- Cancer type (e.g., "ovarian_cancer_hgs")
- Treatment line (e.g., "first-line", "second-line", "maintenance")
- Biomarkers (HRD, TMB, MSI, germline status)
- Treatment history (prior therapies, current medications)

**Step 2: Cancer Type-Specific Food Library**

- Load cancer type → foods mapping
- Filter foods by treatment line appropriateness
- Filter foods by biomarker match

**Step 3: SPE Validation for Each Food**

- Sequence (S): 0.4 weight (currently 0.5 neutral, Evo2 disabled)
- Pathway (P): 0.3 weight (TCGA-weighted alignment)
- Evidence (E): 0.3 weight (literature strength)
- Formula: `0.4×S + 0.3×P + 0.3×E`

**Step 4: Treatment Line Intelligence**

- Compute line_appropriateness (treatment line match)
- Compute cross_resistance (prior therapy conflicts)
- Compute sequencing_fitness (optimal timing)

**Step 5: Ranking & Display**

- Rank foods by: `efficacy_score` (SPE: 0.4×S + 0.3×P + 0.3×E)
- Confidence: `(S+P+E)/3 + SAE_boost` where `SAE_boost = (line_app + seq_fit) × 0.05`
- **Design**: Treatment line features boost confidence, not ranking order
- Display top N foods in FoodRankingPanel (with SAE structure adapter: flat → nested)
- Show rationale, dosage, safety, citations

### 2.2 Example: First-Line Chemotherapy for Ovarian Cancer

**Input**:

```json
{
  "cancer_type": "ovarian_cancer_hgs",
  "treatment_line": "first-line",
  "biomarkers": {
    "HRD": "POSITIVE",
    "TMB": 8,
    "germline_BRCA": "NEGATIVE"
  },
  "treatment_history": {
    "current_line": "first-line",
    "prior_therapies": []
  }
}
```

**Processing**:

1. Load ovarian cancer pathway weights (TP53: 0.95, HRD/DDR: 0.112, etc.)
2. Filter foods targeting TP53, DNA repair pathways
3. Apply treatment line filter (first-line appropriate foods)
4. Apply biomarker gates (HRD+ → boost DNA repair foods)
5. SPE validation for each food
6. Rank by composite score

**Output** (FoodRankingPanel):

- Top food: "Vitamin D" (efficacy: 0.72, confidence: 0.68)
  - Rationale: "Targets DNA repair pathways (HRD+ match), first-line appropriate"
  - Dosage: "2000-4000 IU daily"
  - SAE features: line_appropriateness: 0.9, sequencing_fitness: 0.85

---

## 3. DATA STRUCTURES & DATABASES

### 3.1 Cancer Type-Specific Food Library

**New File**: `api/resources/cancer_type_food_recommendations.json`

**Structure**:

```json
{
  "ovarian_cancer_hgs": {
    "primary_pathways": ["tp53", "hrd_ddr", "cell_cycle"],
    "recommended_foods": [
      {
        "compound": "Vitamin D",
        "pathways": ["dna_repair", "immune_modulation"],
        "treatment_lines": ["first-line", "maintenance"],
        "biomarker_gates": {
          "HRD": ["POSITIVE"],
          "TMB": [">=10"]
        },
        "mechanisms": ["dna_repair_support", "immune_boost"]
      }
    ]
  }
}
```

**Code Integration**:

- Load in `hypothesis_validator.py:711-737` (pathway loading section)
- Use in `food_spe_integration.py:216-224` (pathway alignment)

### 3.2 Treatment Line Rules Database

**File**: `.cursor/ayesha/hypothesis_validator/data/supplement_treatment_rules.json` ✅ **EXISTS**

**Status**: File exists with 22+ compounds (verified in codebase)

**Structure** (actual file structure):

```json
{
  "supplement_rules": {
    "Vitamin D": {
      "mechanism": "dna_repair_support",
      "high_appropriateness_contexts": ["hrd_positive", "dna_repair_deficient"],
      "default_scores": {
        "line_appropriateness": 0.9,
        "cross_resistance": 0.0,
        "sequencing_fitness": 0.85
      },
      "biomarker_gates": {
        "HRD": "POSITIVE"
      }
    }
  }
}
```

**Code Integration**:

- Load in `food_treatment_line_service.py:17` (file path verified)
- Use in `food_treatment_line_service.py:94-122` (biomarker gates)
- Service returns flat structure: `{"line_appropriateness": 0.9, "cross_resistance": 0.0, "sequencing_fitness": 0.85}`
- Frontend expects nested: `{line_fitness: {score: 0.9, status: "appropriate", reason: "..."}, ...}`
- **Gap**: `validate_food_dynamic` needs adapter (same as `validate_food_ab_enhanced` has)

### 3.3 Biomarker → Food Mapping

**New File**: `api/resources/biomarker_food_mapping.json`

**Structure**:

```json
{
  "HRD_POSITIVE": {
    "target_pathways": ["dna_repair", "hrd_ddr"],
    "recommended_foods": ["Vitamin D", "NAC", "Folate"],
    "mechanisms": ["dna_repair_support", "homologous_recombination"]
  },
  "TMB_HIGH": {
    "target_pathways": ["immune_surveillance", "checkpoint_inhibition"],
    "recommended_foods": ["Omega-3", "Curcumin", "Green Tea"],
    "mechanisms": ["immune_modulation", "anti_inflammatory"]
  }
}
```

**Code Integration**:

- Use in `food_treatment_line_service.py:102-122` (biomarker gates)
- Use in `food_spe_integration.py:465-476` (biomarker boosts)

---

## 4. INTEGRATION POINTS (Enhanced)

### 4.1 Backend → Frontend Integration

**Current**: FoodRankingPanel receives food objects with basic fields

**Enhanced**: FoodRankingPanel receives:

- Treatment line context (normalized to "L1", "L2", "L3")
- Cancer type (ovarian_cancer_hgs, breast_cancer, etc.)
- Biomarker matches (HRD+ match, TMB-high match)
- Treatment line intelligence breakdown (nested SAE structure)

**Code Changes**:

- `hypothesis_validator.py:798-858` (response assembly) - **CRITICAL**: Add SAE structure adapter (flat → nested)
  - Transform: `{"line_appropriateness": 0.9}` → `{"line_fitness": {"score": 0.9, "status": "appropriate", "reason": "..."}}`
  - **Reference**: Copy adapter logic from `validate_food_ab_enhanced` (line 435-446)
- `food_treatment_line_service.py:33-141` - add treatment line format normalization ("first-line" → "L1")
- `FoodRankingPanel.jsx:47-189` - display treatment line and biomarker info

### 4.2 Treatment Line Service → Food Ranking

**Current**: Treatment line features computed and used for confidence boost only (by design)

**Design Decision** (Code-Verified): Treatment line features affect confidence, not primary ranking

- Ranking: `efficacy_score = 0.4×S + 0.3×P + 0.3×E` (SPE only)
- Confidence: `(S+P+E)/3 + SAE_boost` where `SAE_boost = (line_app + seq_fit) × 0.05`
- Code: `food_spe_integration.py:233-237` (overall_score) - does NOT include treatment line
- Code: `food_spe_integration.py:458-463` (SAE boost) - only affects confidence

**Code Changes**:

- Document this design decision (treatment line appropriateness boosts confidence, not ranking)
- Ensure SAE boost is properly applied in `_compute_confidence` (already implemented)

### 4.3 Cancer Type → Pathway → Food Chain

**Current**: Pathway alignment exists but not cancer type-specific

**Enhanced**:

1. Load cancer type → pathways mapping
2. Load pathways → foods mapping
3. Filter foods by cancer type pathways
4. Rank by pathway alignment score

**Code Changes**:

- `food_spe_integration.py:67-90` (_get_disease_pathway_weights) - enhance with food filtering
- `hypothesis_validator.py:711-737` (pathway loading) - load cancer type food library

---

## 5. REAL GAPS (Product-Focused)

### 5.1 Critical Gaps

**Gap 1: Missing Cancer Type-Specific Food Library**

- Status: No database mapping cancer types to recommended foods
- Impact: Cannot provide cancer-specific recommendations
- Action: Create `cancer_type_food_recommendations.json`

**Gap 2: SAE Structure Mismatch (CRITICAL - Blocks Frontend Display)**

- Status: Inconsistent implementation - one endpoint has adapter, other doesn't
- Code Evidence:
  - `validate_food_ab_enhanced` (line 435-446): ✅ HAS adapter that transforms flat → nested
  - `validate_food_dynamic` (line 787): ❌ NO adapter - passes flat structure through
  - Service: `food_treatment_line_service.py:95-99` returns `{"line_appropriateness": 0.9, "cross_resistance": 0.0, "sequencing_fitness": 0.85}`
  - Frontend: `FoodRankingPanel.jsx:119-125` expects `sae_features.line_fitness.score`
  - Frontend: `SAEFeatureCards.jsx:31, 172-188` expects nested structure with `status` and `reason`
- Impact: FoodRankingPanel cannot display treatment line intelligence correctly
- Action: Add adapter function in `hypothesis_validator.py:798-858` to transform flat → nested (copy from `validate_food_ab_enhanced`)

**Gap 3: Treatment Line Format Inconsistency**

- Status: Multiple formats used: "L3", "first-line", "first-line chemotherapy", "frontline"
- Code Evidence: 
  - `validate_food_dynamic` (line 654): Expects `"current_line": "L3"` format
  - `ayesha_orchestrator.py` (line 302): Uses `"L{treatment_line}"` format
  - `trial_intelligence_universal/config.py` (line 111): Uses `['frontline', 'first-line', 'first line', 'primary', '1l']`
- User Decision: Use numeric format ("L1", "L2", "L3")
- Impact: Need normalization layer to handle all formats
- Action: Add treatment line normalization in `food_treatment_line_service.py`

**Gap 4: FoodRankingPanel Missing Treatment Line Display**

- Status: Component doesn't show treatment line context or biomarker matches
- Impact: Patients don't see why foods are recommended for their treatment line
- Action: Enhance FoodRankingPanel to display treatment line and biomarker context

### 5.2 Medium Priority Gaps

**Gap 5: Biomarker → Food Mapping Not Systematic**

- Status: Biomarker gates exist but not organized into database
- Impact: Cannot easily add new biomarker → food mappings
- Action: Create `biomarker_food_mapping.json`

**Gap 6: Evo2 Sequence Scoring Disabled**

- Status: Sequence (S) component always 0.5 (neutral)
- Impact: Missing biological plausibility signal
- Action: Implement Evo2 food plausibility service (Phase 2)

**Gap 7: No Treatment Line-Specific Evidence**

- Status: Evidence mining doesn't consider treatment line
- Impact: Evidence may not be relevant to patient's treatment line
- Action: Enhance evidence service to filter by treatment line context

### 5.3 Low Priority Gaps

**Gap 8: No Batch Food Recommendations**

- Status: Can validate one food at a time
- Impact: Cannot recommend multiple foods for a patient scenario
- Action: Create batch recommendation endpoint

**Gap 9: No Food Interaction Checking**

- Status: Drug interactions checked, but not food-food interactions
- Impact: May recommend conflicting foods
- Action: Add food-food interaction checking

---

## 6. CODE REFERENCES (Updated)

### 6.1 Entry Points

- API Endpoint: `oncology-coPilot/oncology-backend-minimal/api/routers/hypothesis_validator.py:630-858`
- Frontend Component: `oncology-coPilot/oncology-frontend/src/components/ayesha/FoodRankingPanel.jsx:26-192`
- Treatment Line Service: `oncology-coPilot/oncology-backend-minimal/api/services/food_treatment_line_service.py:33-141`

### 6.2 Core Services

- SPE Integration: `oncology-coPilot/oncology-backend-minimal/api/services/food_spe_integration.py:32-499`
- Treatment Line: `oncology-coPilot/oncology-backend-minimal/api/services/food_treatment_line_service.py:33-141`
- Pathway Weights: `oncology-coPilot/oncology-backend-minimal/api/services/food_spe_integration.py:67-90`

### 6.3 Data Files (Code-Verified Locations)

**EXISTS (Verified)**:

- Universal Pathway DB: `oncology-coPilot/oncology-backend-minimal/api/resources/universal_disease_pathway_database.json` ✅
- Treatment Line Rules: `.cursor/ayesha/hypothesis_validator/data/supplement_treatment_rules.json` ✅ (22+ compounds)
- Cancer Pathways: `.cursor/ayesha/hypothesis_validator/data/cancer_pathways.json` ✅
- Drug Interactions: `.cursor/ayesha/hypothesis_validator/data/drug_interactions.json` ✅
- Safety Database: `.cursor/ayesha/hypothesis_validator/data/safety_database.json` ✅

**NEEDS CREATION**:

- Cancer Type Foods: `.cursor/ayesha/hypothesis_validator/data/cancer_type_food_recommendations.json` (new)
- Biomarker Mapping: `.cursor/ayesha/hypothesis_validator/data/biomarker_food_mapping.json` (new)

**Pattern**: Food validator data files in `.cursor/ayesha/hypothesis_validator/data/`, universal/system data in `api/resources/`

---

## 7. IMMEDIATE ACTION ITEMS (Product-Focused)

### Priority 1 (Critical - Product Launch)

1. **Fix SAE Structure Adapter in validate_food_dynamic** ⚠️ **CRITICAL**

   - Issue: `validate_food_dynamic` (line 787) returns flat structure, frontend expects nested
   - Code: `hypothesis_validator.py:798-858` (response assembly)
   - Action: Add adapter function (copy from `validate_food_ab_enhanced` line 435-446)
   - Transform: `{"line_appropriateness": 0.9}` → `{"line_fitness": {"score": 0.9, "status": "appropriate", "reason": "..."}}`
   - Test: Verify FoodRankingPanel displays treatment line intelligence correctly

2. **Add Treatment Line Format Normalization**

   - Issue: Multiple formats ("L3", "first-line", "frontline")
   - User Decision: Use numeric format ("L1", "L2", "L3")
   - Code: `food_treatment_line_service.py:33-141` (add normalization)
   - Action: Add function to normalize all formats to "L1"/"L2"/"L3"
   - Test: Handle "first-line", "L3", "frontline" → all normalize correctly

3. **Create Cancer Type-Specific Food Library**

   - File: `.cursor/ayesha/hypothesis_validator/data/cancer_type_food_recommendations.json`
   - Structure: Cancer types → pathways → recommended foods with treatment line filters
   - Integration: `hypothesis_validator.py:711-737` (pathway loading section)
   - Test: Load ovarian_cancer_hgs → get foods targeting TP53, HRD/DDR pathways

4. **Enhance FoodRankingPanel Display**

   - Code: `FoodRankingPanel.jsx:47-189`
   - Add: 
     - Treatment line context display ("Recommended for: L1 chemotherapy, Ovarian Cancer (HGS)")
     - Biomarker match indicators ("HRD+ match", "TMB-high match")
     - Treatment line intelligence breakdown (already in accordion, needs structure fix from #1)
   - Test: Display first-line ovarian cancer recommendations with all context

### Priority 2 (Important - Product Enhancement)

5. **Create Biomarker → Food Mapping**

   - File: `.cursor/ayesha/hypothesis_validator/data/biomarker_food_mapping.json`
   - Structure: HRD_POSITIVE, TMB_HIGH, MSI_HIGH → target pathways → recommended foods
   - Integration: `food_treatment_line_service.py:102-122` (biomarker gates)
   - Test: HRD+ patient → get DNA repair foods (Vitamin D, NAC, Folate)

6. **Enhance Evidence Service for Treatment Line**

   - Code: `enhanced_evidence_service.py:955-1016` (get_complete_evidence)
   - Filter: Evidence by treatment line context
   - Test: First-line vs maintenance evidence

7. **Add Batch Food Recommendations**

   - Endpoint: `POST /api/hypothesis/recommend_foods_batch`
   - Input: Cancer type, treatment line, biomarkers
   - Output: Ranked list of foods with SPE scores

### Priority 3 (Enhancement - Future)

8. **Implement Evo2 Food Plausibility Service**

   - Code: `food_spe_integration.py:196-214` (currently disabled)
   - Service: `api/services/evo2_food_plausibility.py` (needs creation)
   - Impact: Enable Sequence (S) component (currently 0.5 neutral)

9. **Add Food-Food Interaction Checking**

   - Service: `api/services/food_interaction_service.py` (needs creation)
   - Integration: `dietician_recommendations.py:181-243` (currently only drug interactions)

---

## 8. TESTING STRATEGY (Product Validation)

### 8.1 Test Scenarios

**Scenario 1: First-Line Chemotherapy for Ovarian Cancer (HGS)**

- Input: Cancer type = ovarian_cancer_hgs, Treatment line = first-line, HRD = POSITIVE
- Expected: Foods targeting TP53, DNA repair pathways, first-line appropriate
- Validation: Top foods have high line_appropriateness, HRD+ biomarker match

**Scenario 2: Maintenance Therapy for Ovarian Cancer**

- Input: Cancer type = ovarian_cancer_hgs, Treatment line = maintenance, HRD = POSITIVE
- Expected: Different foods than first-line (maintenance-appropriate)
- Validation: Foods have high sequencing_fitness for maintenance

**Scenario 3: Breast Cancer with HER2+**

- Input: Cancer type = breast_cancer, Treatment line = first-line, HER2 = POSITIVE
- Expected: Foods targeting HER2 pathway
- Validation: Pathway alignment includes HER2 signaling

### 8.2 Validation Checklist

- [ ] Treatment line rules database created and loaded
- [ ] Cancer type food library created and integrated
- [ ] Treatment line used in ranking algorithm
- [ ] FoodRankingPanel displays treatment line context
- [ ] Biomarker gates working (HRD+, TMB-high, etc.)
- [ ] SPE scores computed correctly (0.4×S + 0.3×P + 0.3×E)
- [ ] Patient-facing display shows rationale, dosage, safety

---

## 9. PRODUCT SUCCESS CRITERIA (Code Review Validated)

### 9.1 Functional Requirements (Verified in Code)

- ✅ Can recommend foods for specific cancer types (ovarian, breast, etc.)
  - Code: `food_spe_integration.py:67-90` (TCGA-weighted pathways per disease)
  - Gap: Need cancer type→food mapping (Gap 1)
- ⚠️ Can filter foods by treatment line (L1, L2, L3)
  - Code: `food_treatment_line_service.py:33-141` (computes line_appropriateness)
  - Gap: Treatment line format normalization needed (Gap 3)
- ✅ Can target foods to patient biomarkers (HRD+, TMB-high, MSI-high)
  - Code: `food_treatment_line_service.py:102-122` (biomarker gates)
  - Code: `food_spe_integration.py:465-476` (biomarker boosts)
- ✅ Can validate foods through SPE framework (S/P/E scores)
  - Code: `food_spe_integration.py:233-237` (formula: 0.4×S + 0.3×P + 0.3×E)
- ⚠️ Can display recommendations to patients via FoodRankingPanel
  - Code: `FoodRankingPanel.jsx:26-192` (component exists)
  - Gap: SAE structure mismatch needs adapter (Gap 2, Gap 9)

### 9.2 Quality Requirements (Code-Validated Thresholds)

- Treatment line appropriateness > 0.7 for top recommendations
  - Code: `food_treatment_line_service.py:95-99` (default 0.6, boosted by gates)
  - Validation: Check supplement_treatment_rules.json for compound-specific rules
- Pathway alignment > 0.6 for recommended foods
  - Code: `food_spe_integration.py:329-415` (TCGA-weighted alignment)
  - Validation: Check pathway matching logic and weights
- Evidence tier >= "CONSIDER" for top recommendations
  - Code: `food_spe_integration.py:417-425` (_convert_evidence_grade)
  - Code: `food_spe_integration.py:482-497` (_classify_verdict)
  - Validation: Check evidence_grade conversion and verdict classification
- Patient-facing display shows clear rationale and dosage
  - Code: `FoodRankingPanel.jsx:90-107` (displays dosage and rationale)
  - Validation: Verify response includes dietician_recommendations

### 9.3 User Experience Requirements (Frontend Verified)

- ⚠️ Patient can see why foods are recommended (treatment line, biomarker match)
  - Code: `FoodRankingPanel.jsx:47-189` (missing treatment line context display)
  - Gap: Need to add treatment line label and biomarker match indicators (Gap 5)
- ⚠️ Patient can see treatment line intelligence (line_appropriateness, sequencing_fitness)
  - Code: `FoodRankingPanel.jsx:109-152` (accordion exists but structure mismatch)
  - Gap: SAE structure adapter needed (Gap 9)
- ✅ Patient can see SPE breakdown (S/P/E scores)
  - Code: `FoodRankingPanel.jsx:70-77` (displays efficacy_score and confidence)
  - Note: Full SPE breakdown not displayed, only overall score
- ✅ Patient can see safety and interaction warnings
  - Code: `FoodRankingPanel.jsx:90-107` (displays dosage)
  - Code: Response includes `dietician_recommendations` (safety, interactions)
  - Note: Safety/interactions not displayed in FoodRankingPanel (may be in details view)

---

**DOCTRINE STATUS**: Product-Focused Plan - Ready for Implementation

**LAST UPDATED**: 2025-01-XX (Comprehensive Code Review Complete)

**FOCUS**: Cancer-specific food recommendations with treatment line intelligence

---

## 10. CODE REVIEW FINDINGS (Comprehensive)

### 10.1 Critical Findings (Code-Verified)

**Finding 1: SAE Structure Mismatch (CRITICAL)** ⚠️

- **Status**: Inconsistent implementation
- **Code Evidence**:
  - `validate_food_ab_enhanced` (line 435-446): ✅ HAS adapter (flat → nested)
  - `validate_food_dynamic` (line 787): ❌ NO adapter (passes flat through)
  - Service returns: `{"line_appropriateness": 0.9, "cross_resistance": 0.0, "sequencing_fitness": 0.85}`
  - Frontend expects: `{line_fitness: {score: 0.9, status: "appropriate", reason: "..."}, ...}`
- **Impact**: Frontend cannot display treatment line intelligence correctly
- **Action**: Add adapter in `validate_food_dynamic` response assembly (line 798-858)

**Finding 2: Treatment Line Rules File EXISTS** ✅

- **Location**: `.cursor/ayesha/hypothesis_validator/data/supplement_treatment_rules.json`
- **Status**: File exists with 22+ compounds (verified)
- **Structure**: Matches expected format with `default_scores`, `biomarker_gates`, `high_appropriateness_contexts`
- **Action**: No creation needed - file is already there

**Finding 3: Treatment Line Format Inconsistency**

- **Multiple formats**: "L3", "first-line", "frontline", "first line"
- **User Decision**: Standardize to "L1", "L2", "L3"
- **Action**: Add normalization function to handle all formats

**Finding 4: Treatment Line NOT in Ranking (By Design)**

- **Code**: `food_spe_integration.py:233-237` (overall_score) does NOT include treatment line
- **Code**: `food_spe_integration.py:458-463` (SAE boost) only affects confidence
- **Design Decision**: Treatment line appropriateness boosts confidence, not ranking
- **Action**: Document this design (no code change needed)

### 10.2 Data File Locations (Verified)

- **Food validator data**: `.cursor/ayesha/hypothesis_validator/data/` (consistent pattern)
- **Universal/system data**: `api/resources/` (different type of data)
- **All food validator data files use 5-level parent traversal pattern** (verified working)

### 10.3 Questions Answered

1. ✅ **supplement_treatment_rules.json EXISTS** (verified at `.cursor/ayesha/hypothesis_validator/data/`)
2. ✅ **SAE structure mismatch identified** (inconsistent implementation - one endpoint has adapter, other doesn't)
3. ✅ **Treatment line format inconsistency identified** (no normalization layer)
4. ✅ **Treatment line affects confidence only** (by design, verified in code)
5. ✅ **Cancer type food library missing** (confirmed - does not exist)
6. ✅ **Batch recommendation partial** (orchestrator has it, no dedicated endpoint)

### 10.4 Implementation Phases

**Phase 1: Critical Fixes (P0 - Blocks Product Launch)**

1. Fix SAE structure adapter in `validate_food_dynamic`
2. Add treatment line format normalization

**Phase 2: Core Capabilities (P1 - Product Features)**

3. Create cancer type food library
4. Enhance FoodRankingPanel display

**Phase 3: Enhanced Features (P2 - Product Enhancement)**

5. Create biomarker → food mapping
6. Add batch recommendation endpoint
7. Enhance evidence service for treatment line