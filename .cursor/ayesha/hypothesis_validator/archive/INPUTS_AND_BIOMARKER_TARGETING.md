# 📋 INPUTS & BIOMARKER TARGETING GUIDE

**What makes our Dynamic Food Validator unique vs. googling?**

---

## **🎯 COMPLETE INPUT SPECIFICATION**

### **Required Inputs:**

```json
{
  "compound": "Resveratrol",  // ANY food/supplement name
  "disease_context": {
    "disease": "ovarian_cancer_hgs",
    "mutations": [{"gene": "TP53", "hgvs_p": "R248Q"}],
    "biomarkers": {
      "HRD": "POSITIVE",      // ⭐ KEY: Homologous Recombination Deficiency
      "TMB": 8.2,             // Tumor Mutational Burden
      "MSI": "STABLE",        // Microsatellite Instability
      "PIK3CA": "MUTANT",     // Specific gene mutations
      "BRCA1": "NEGATIVE",    // Germline status
      "BRCA2": "NEGATIVE"
    },
    "pathways_disrupted": [
      "DNA repair",
      "Angiogenesis",
      "Inflammation",
      "Cell cycle"
    ]
  }
}
```

### **Optional Inputs (Recommended for Best Results):**

```json
{
  "treatment_history": {
    "current_line": "L3",                    // Treatment line (L1, L2, L3+)
    "prior_therapies": [                     // Prior treatments
      "carboplatin",
      "paclitaxel",
      "bevacizumab"
    ]
  },
  "patient_medications": [                   // Current medications (for interaction checking)
    "warfarin",
    "metformin"
  ],
  "use_evo2": false                          // Phase 1: disabled, Phase 2: experimental
}
```

---

## **🔬 HOW BIOMARKER TARGETING WORKS**

### **Step-by-Step Biomarker Integration:**

#### **1. Biomarker Extraction from Disease Context**
```python
biomarkers = disease_context.get("biomarkers", {})
# Extracts: HRD, TMB, MSI, PIK3CA, BRCA1/2, etc.
```

#### **2. Compound Target Extraction**
```python
# System dynamically extracts targets for compound
# Example: Vitamin D → ["VDR", "TP53 pathway", "BRCA1", "DNA repair"]
targets = extractor.extract_all("Vitamin D", disease)
```

#### **3. Pathway Mapping**
```python
# Maps compound targets to cancer pathways
# Example: VDR, BRCA1 → "DNA repair" pathway
pathways = mapper.map_targets_to_pathways(targets)
# Result: ["DNA repair", "Cell cycle regulation"]
```

#### **4. Biomarker Gating (SAE Treatment Line Service)**

**File:** `supplement_treatment_rules.json`

```json
{
  "Vitamin D": {
    "biomarker_gates": {
      "HRD": "POSITIVE"        // ⭐ Gate: Only boost if HRD+
    },
    "high_appropriateness_contexts": [
      "hrd_positive",           // Context: HRD+ patient
      "dna_repair_deficient"    // Context: DNA repair pathway disrupted
    ]
  }
}
```

**Logic:**
```python
# Check if biomarker matches gate
if biomarkers.get("HRD") == "POSITIVE":
    # Vitamin D is especially appropriate for HRD+ patients
    line_appropriateness += 0.1  # Boost from 0.9 → 1.0
```

#### **5. Pathway Alignment Scoring**

**Compound Pathways:** `["DNA repair", "Cell cycle regulation"]`  
**Disease Pathways:** `["DNA repair", "Angiogenesis"]`

**Calculation:**
```python
# Intersection: {"DNA repair"}
aligned_count = 1
alignment_ratio = 1 / 2 = 0.5

# Score: aligned=1.0, misaligned=0.2
pathway_score = 0.5 * 1.0 + 0.5 * 0.2 = 0.6
```

**Biomarker Impact:**
- If `HRD="POSITIVE"` AND compound targets DNA repair → Higher pathway score
- If `TMB >= 10` → Additional confidence boost (+0.03)

#### **6. Confidence Modulation**

**Formula:**
```python
base_confidence = (S + P + E) / 3.0

# Biomarker boosts
biomarker_boost = 0.0
if biomarkers.get("HRD") == "POSITIVE":
    if "DNA repair" in pathways_disrupted:
        if compound_targets_dna_repair:
            biomarker_boost += 0.05  # HRD+ + DNA repair compound

if biomarkers.get("TMB", 0) >= 10:
    biomarker_boost += 0.03  # High TMB

final_confidence = min(base + biomarker_boost, 0.95)
```

---

## **⚡ WHAT MAKES THIS UNIQUE VS. GOOGLING**

### **❌ Google Search Limitations:**

1. **Generic Information Only**
   - Search: "Vitamin D ovarian cancer"
   - Result: General information, not personalized

2. **No Biomarker Context**
   - Doesn't know if you're HRD+ or HRD-
   - Doesn't adjust recommendations based on biomarkers

3. **No Treatment Line Intelligence**
   - Doesn't consider where you are in treatment journey
   - Doesn't account for prior therapies

4. **No Drug Interaction Checking**
   - Doesn't check YOUR medication list
   - Generic warnings only

5. **No Integrated Scoring**
   - Just papers, no unified score
   - No confidence assessment
   - No verdict (SUPPORTED/WEAK_SUPPORT/NOT_SUPPORTED)

6. **No Pathway Analysis**
   - Doesn't map compound targets to YOUR cancer pathways
   - Doesn't show pathway alignment score

7. **No Dietician-Grade Recommendations**
   - No dosage with citations
   - No timing guidance
   - No meal planning
   - No lab monitoring recommendations

### **✅ Our System Provides:**

#### **1. PERSONALIZED TO YOUR PROFILE**
```
Input: HRD+, L3 post-platinum, on warfarin
Output: Vitamin D recommendation adjusted for:
  - HRD+ biomarker (boosted appropriateness)
  - Post-platinum context (oxidative stress support)
  - Warfarin interaction (Vitamin K warning)
```

#### **2. BIOMARKER-SPECIFIC TARGETING**
```
Example: Vitamin D for HRD+ vs HRD- patient

HRD+ Patient:
  → Pathway alignment: HIGH (DNA repair pathway active)
  → Line appropriateness: 1.0 (biomarker gate matches)
  → Confidence boost: +0.05 (HRD+ + DNA repair compound)
  → Verdict: WEAK_SUPPORT → SUPPORTED (threshold crossed)

HRD- Patient:
  → Pathway alignment: MODERATE (no DNA repair pathway)
  → Line appropriateness: 0.9 (no biomarker boost)
  → Confidence: Lower (no biomarker boost)
  → Verdict: WEAK_SUPPORT
```

#### **3. INTEGRATED S/P/E + SAE SCORING**
```
Sequence (S): Evo2 plausibility (0.4 weight)
Pathway (P): Alignment to YOUR cancer pathways (0.3 weight)
Evidence (E): Literature grade → score (0.3 weight)
SAE: Treatment line + biomarker intelligence

Final Score = 0.4×S + 0.3×P + 0.3×E
Confidence = Base + SAE boost + Biomarker boost
Verdict = Classify(score, confidence)
```

#### **4. TREATMENT LINE INTELLIGENCE**
```
SAE Features:
  - Line Appropriateness: How appropriate for L3?
  - Cross-Resistance Risk: Risk of interfering with prior therapies
  - Sequencing Fitness: Safe to sequence with other treatments

Biomarker Impact:
  - HRD+ + DNA repair compound → Higher appropriateness
  - Post-platinum + NAC → Higher appropriateness (oxidative stress recovery)
```

#### **5. DRUG INTERACTION CHECKING**
```
Input: patient_medications = ["warfarin"]
System checks: drug_interactions.json
Result: 
  ⚠️ WARNING: Vitamin D + warfarin
    - Severity: Moderate
    - Action: Monitor INR closely
```

#### **6. EVIDENCE GRADING (Not Just Papers)**
```
Google: Returns 50 papers, user has to figure out quality
Our System:
  - Grades evidence: STRONG/MODERATE/WEAK/INSUFFICIENT
  - Counts RCTs vs. observational studies
  - Synthesizes mechanisms from papers
  - Provides unified evidence score (0-1)
```

#### **7. DIETICIAN-GRADE OUTPUT**
```
Complete Package:
  ✅ Dosage: "2000-4000 IU daily (target serum 25(OH)D: 40-60 ng/mL)"
  ✅ Timing: "Morning with breakfast (fat-soluble, needs dietary fat)"
  ✅ Meal Planning: "Combine with: Eggs, Avocado, Nuts"
  ✅ Drug Interactions: "Avoid or monitor with warfarin"
  ✅ Lab Monitoring: "Serum 25(OH)D q3-6 months, Calcium levels"
  ✅ Safety: "Max dose: 10000 IU/day, Contraindications: Hypercalcemia"
  ✅ Patient Instructions: Plain-language summary
```

---

## **📊 EXAMPLE: COMPARISON**

### **Scenario: "Can Vitamin D help my ovarian cancer?"**

#### **Google Search Result:**
```
- Generic articles about Vitamin D and cancer
- Not personalized to ovarian cancer
- Doesn't know you're HRD+
- Doesn't know you're L3 post-platinum
- Doesn't know you're on warfarin
- You have to read 20 papers yourself
- No verdict or recommendation
```

#### **Our System Result:**
```
Compound: Vitamin D
Verdict: WEAK_SUPPORT
Score: 0.635 (63.5%)
Confidence: 0.82 (82%)

S/P/E Breakdown:
  - Sequence (S): 0.500 (neutral, Phase 1)
  - Pathway (P): 0.850 (HIGH - DNA repair alignment)
  - Evidence (E): 0.600 (MODERATE grade)

Targets: VDR, TP53 pathway, BRCA1, DNA repair
Pathways: DNA repair, Cell cycle regulation

SAE Features:
  - Line Appropriateness: 1.0 (100% - HRD+ biomarker match)
  - Cross-Resistance: 0.0 (safe with prior therapies)
  - Sequencing Fitness: 0.85 (high - good for sequencing)

Evidence:
  - Grade: MODERATE
  - Papers: 15 total, 3 RCTs
  - Key finding: HR 0.77 for mortality with serum >30 ng/mL

Dietician Recommendations:
  - Dosage: 2000-4000 IU daily
  - Target Level: 40-60 ng/mL (serum 25(OH)D)
  - Timing: Morning with breakfast (fat-soluble)
  - ⚠️ Interaction: Monitor INR if on warfarin
  - Monitoring: Serum 25(OH)D q3-6 months, Calcium levels
  - Safety: Max 10000 IU/day, avoid if hypercalcemia

Biomarker Targeting:
  ✅ HRD+ detected → Line appropriateness boosted
  ✅ DNA repair pathway active → Pathway alignment HIGH
  ✅ Biomarker confidence boost applied (+0.05)
```

---

## **🎯 KEY DIFFERENTIATORS SUMMARY**

| Feature | Google Search | Our System |
|---------|---------------|------------|
| **Personalization** | ❌ Generic | ✅ Your cancer type, biomarkers, treatment history |
| **Biomarker Targeting** | ❌ None | ✅ HRD+, TMB, MSI-specific recommendations |
| **Treatment Line Intelligence** | ❌ None | ✅ SAE features (appropriateness, cross-resistance, sequencing) |
| **Drug Interactions** | ❌ Generic warnings | ✅ Checks YOUR medication list |
| **Integrated Scoring** | ❌ Just papers | ✅ S/P/E + SAE unified score |
| **Evidence Grading** | ❌ User must evaluate | ✅ STRONG/MODERATE/WEAK classification |
| **Pathway Analysis** | ❌ None | ✅ Target → Pathway mapping + alignment scores |
| **Dietician Guidance** | ❌ None | ✅ Complete dosage/timing/interactions/safety package |
| **Verdict** | ❌ User decides | ✅ SUPPORTED/WEAK_SUPPORT/NOT_SUPPORTED |
| **Confidence** | ❌ Unknown | ✅ 0-100% with biomarker modulation |

---

**⚔️ BOTTOM LINE: Google gives you information. Our system gives you PERSONALIZED, EVIDENCE-GRADED, TREATMENT-LINE-AWARE, BIOMARKER-TARGETED RECOMMENDATIONS.**

