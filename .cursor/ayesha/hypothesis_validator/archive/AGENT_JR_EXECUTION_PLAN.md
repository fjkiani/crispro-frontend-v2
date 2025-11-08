# ⚔️ AGENT JR EXECUTION PLAN - FOOD VALIDATOR FOR AYESHA

**Status:** 🚨 **NEEDS FIXES** - Dynamic framework built but has critical gaps

**Mission:** Fix stubbed components to make Food Validator truly dynamic for ANY compound

**Current State:** 60% dynamic, 40% hardcoded/stubbed (see CODE_AUDIT_REPORT.md)

**Timeline:** 2-3 hours to fix critical gaps, then demo-ready

---

## **📋 WHAT AGENT JR BUILT (COMPLETED)**

✅ **Dynamic Target Extraction** - ChEMBL/PubChem APIs working  
✅ **Pathway Mapping** - Dynamic mapping to cancer mechanisms  
✅ **S/P/E Integration** - Formula and scoring working  
✅ **Biomarker Targeting** - HRD+, TMB boosts functional  
✅ **SAE Treatment Line** - Biomarker gates working  
✅ **Drug Interactions** - Patient medication checking working  

## **🚨 CRITICAL GAPS FOUND (NEED FIXES)**

❌ **Evidence Synthesis is STUB** - Always returns "MODERATE" (hardcoded)  
❌ **Dosage Extraction is STUB** - Returns empty strings  
❌ **Timing/Meal Planning LIMITED** - Only ~10 compounds (hardcoded patterns)  
❌ **SAE Coverage LIMITED** - Requires supplement_treatment_rules.json entry  

**See:** `.cursor/ayesha/hypothesis_validator/CODE_AUDIT_REPORT.md` for complete audit

---

## **🎯 COMMANDER'S STRATEGIC DECISIONS**

Based on manager's assessment and Zo's tactical analysis:

### **Critical Decision: Build P/E/SAE FIRST (Phased Approach)**

**Rationale:**
- P/E/SAE alone: 90% confidence ("definitely work")
- Evo2 variant proxy: 60% confidence ("scientifically questionable but novel")
- **Value**: P/E/SAE + treatment line intelligence > PubMed manual search
- **Differentiation**: Integrated scoring, biomarker gating, SAE confidence

**Path:** Build guaranteed value first, add Evo2 as experimental enhancement later

---

## **🔧 PRIORITY FIXES (2-3 HOURS) - DO THIS FIRST**

### **Fix 1: Implement Real Evidence Synthesis (1 hour)** 🚨

**Problem:** `enhanced_evidence_service.py` lines 148-197 returns hardcoded "MODERATE"

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/enhanced_evidence_service.py`

**Current Code (STUB):**
```python
# Line 180: HARDCODED
synthesis = {
    "evidence_grade": "MODERATE",  # ← ALWAYS THE SAME
    "mechanisms": [],
    "dosage": "Based on literature review",
    "safety": "Generally safe with precautions"
}
```

**Required Fix:**
```python
async def synthesize_evidence_llm(self, compound: str, disease: str, papers: List[Dict]) -> Dict[str, Any]:
    """Use LLM to synthesize evidence from papers."""
    if not LLM_AVAILABLE or not papers:
        return {
            "evidence_grade": "INSUFFICIENT",
            "mechanisms": [],
            "dosage": "",
            "safety": "",
            "outcomes": []
        }
    
    try:
        llm_service = get_llm_service()
        
        # Build context from papers
        abstracts_text = "\n\n".join([
            f"PMID: {p['pmid']}\nTitle: {p['title']}\nAbstract: {p['abstract'][:800]}"
            for p in papers[:10]
        ])
        
        # LLM prompt for synthesis
        prompt = f"""Analyze these research papers about {compound} for {disease}.

Papers:
{abstracts_text}

Extract:
1. Evidence grade (STRONG if ≥3 RCTs, MODERATE if ≥5 papers, WEAK if ≥2, INSUFFICIENT otherwise)
2. Mechanisms of action
3. Typical dosages used in studies
4. Safety concerns mentioned

Return JSON format:
{{
  "evidence_grade": "STRONG|MODERATE|WEAK|INSUFFICIENT",
  "mechanisms": ["mechanism1", "mechanism2"],
  "dosage": "typical dose range",
  "safety": "key safety concerns"
}}
"""
        
        # Call LLM
        response = await llm_service.query(prompt)
        synthesis = json.loads(response)
        
        return synthesis
        
    except Exception as e:
        print(f"⚠️ LLM synthesis error: {e}")
        # Fallback to heuristic grading
        return {
            "evidence_grade": self._heuristic_grade(papers),
            "mechanisms": [],
            "dosage": "",
            "safety": "",
            "outcomes": []
        }

def _heuristic_grade(self, papers: List[Dict]) -> str:
    """Fallback heuristic grading when LLM unavailable."""
    paper_count = len(papers)
    has_rct = any(
        "randomized" in p.get("title", "").lower() or 
        "randomized" in p.get("abstract", "").lower() or
        "RCT" in p.get("title", "").upper()
        for p in papers
    )
    
    if has_rct and paper_count >= 3:
        return "STRONG"
    elif paper_count >= 5:
        return "MODERATE"
    elif paper_count >= 2:
        return "WEAK"
    else:
        return "INSUFFICIENT"
```

**Acceptance:**
- [ ] Evidence grade changes based on paper analysis (not always "MODERATE")
- [ ] Mechanisms extracted from abstracts
- [ ] Heuristic fallback works when LLM unavailable

---

### **Fix 2: Implement Dosage Extraction (45 min)** 🚨

**Problem:** `dietician_recommendations.py` lines 75-110 returns empty dosage

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/dietician_recommendations.py`

**Required Fix:**
```python
def extract_dosage_from_evidence(self, papers: List[Dict], compound: str) -> Dict[str, Any]:
    """Extract dosage from paper abstracts using regex + LLM."""
    
    dosage_info = {
        "recommended_dose": "",
        "dose_range": {},
        "frequency": "daily",
        "duration": "ongoing",
        "citations": [],
        "target_level": ""
    }
    
    # Try LLM extraction first
    if LLM_AVAILABLE and len(papers) > 0:
        try:
            llm_service = get_llm_service()
            
            abstracts_text = "\n\n".join([
                f"PMID: {p['pmid']}\n{p.get('abstract', '')[:500]}"
                for p in papers[:5]
            ])
            
            prompt = f"""Extract dosage information for {compound} from these abstracts:

{abstracts_text}

Return JSON:
{{
  "recommended_dose": "dose range with units",
  "frequency": "daily|weekly|etc",
  "citations": ["PMID1", "PMID2"]
}}

If no dosage found, return empty strings."""
            
            response = await llm_service.query(prompt)
            llm_result = json.loads(response)
            
            if llm_result.get("recommended_dose"):
                dosage_info.update(llm_result)
                return dosage_info
        except Exception as e:
            print(f"⚠️ LLM dosage extraction failed: {e}")
    
    # Fallback: Regex extraction
    import re
    
    for paper in papers[:5]:
        abstract = paper.get("abstract", "").lower()
        
        # Common dosage patterns
        patterns = [
            r'(\d+[-–]\d+)\s*(mg|iu|g|mcg)',  # Range: 2000-4000 IU
            r'(\d+)\s*(mg|iu|g|mcg)',          # Single: 2000 IU
            r'(\d+\.?\d*)\s*(mg|iu|g)',        # Decimal: 2.5 mg
        ]
        
        for pattern in patterns:
            matches = re.findall(pattern, abstract)
            if matches:
                dose, unit = matches[0]
                dosage_info["recommended_dose"] = f"{dose} {unit}"
                dosage_info["citations"].append(paper.get("pmid", ""))
                break
        
        if dosage_info["recommended_dose"]:
            break
    
    return dosage_info
```

**Acceptance:**
- [ ] Dosage extracted from papers when available
- [ ] Regex fallback works
- [ ] Returns empty gracefully when no dose found

---

### **Fix 3: Expand SAE Coverage (30 min)** 🚨

**Problem:** Only 4 compounds have SAE rules, others get neutral defaults

**File:** `.cursor/ayesha/hypothesis_validator/data/supplement_treatment_rules.json`

**Required Fix:** Add 16 more compounds

```json
{
  "supplement_rules": {
    // ... existing 6 compounds ...
    
    "Resveratrol": {
      "mechanism": "sirtuin_activation",
      "high_appropriateness_contexts": ["metabolic_stress", "oxidative_stress"],
      "default_scores": {"line_appropriateness": 0.7, "cross_resistance": 0.0, "sequencing_fitness": 0.75}
    },
    "Quercetin": {
      "mechanism": "antioxidant_anti_inflammatory",
      "high_appropriateness_contexts": ["chronic_inflammation", "oxidative_stress"],
      "default_scores": {"line_appropriateness": 0.7, "cross_resistance": 0.0, "sequencing_fitness": 0.70}
    },
    "Sulforaphane": {
      "mechanism": "nrf2_activation",
      "high_appropriateness_contexts": ["oxidative_stress_high", "post_chemotherapy"],
      "default_scores": {"line_appropriateness": 0.75, "cross_resistance": 0.0, "sequencing_fitness": 0.80}
    },
    "Genistein": {
      "mechanism": "tyrosine_kinase_inhibition",
      "high_appropriateness_contexts": ["hormone_receptor_positive"],
      "default_scores": {"line_appropriateness": 0.65, "cross_resistance": 0.0, "sequencing_fitness": 0.70}
    },
    "Lycopene": {
      "mechanism": "antioxidant",
      "high_appropriateness_contexts": ["oxidative_stress"],
      "default_scores": {"line_appropriateness": 0.6, "cross_resistance": 0.0, "sequencing_fitness": 0.65}
    },
    "Beta-glucan": {
      "mechanism": "immune_modulation",
      "high_appropriateness_contexts": ["post_immunotherapy", "immune_suppressed"],
      "default_scores": {"line_appropriateness": 0.7, "cross_resistance": 0.0, "sequencing_fitness": 0.75}
    },
    "Selenium": {
      "mechanism": "antioxidant_selenoprotein",
      "high_appropriateness_contexts": ["oxidative_stress"],
      "default_scores": {"line_appropriateness": 0.65, "cross_resistance": 0.0, "sequencing_fitness": 0.70}
    },
    "Zinc": {
      "mechanism": "immune_support",
      "high_appropriateness_contexts": ["immune_suppressed"],
      "default_scores": {"line_appropriateness": 0.7, "cross_resistance": 0.0, "sequencing_fitness": 0.75}
    },
    "Melatonin": {
      "mechanism": "circadian_antioxidant",
      "high_appropriateness_contexts": ["sleep_disruption", "oxidative_stress"],
      "default_scores": {"line_appropriateness": 0.65, "cross_resistance": 0.0, "sequencing_fitness": 0.70}
    },
    "CoQ10": {
      "mechanism": "mitochondrial_support",
      "high_appropriateness_contexts": ["post_chemotherapy", "cardiac_protection"],
      "default_scores": {"line_appropriateness": 0.75, "cross_resistance": 0.0, "sequencing_fitness": 0.80}
    },
    "Probiotics": {
      "mechanism": "microbiome_modulation",
      "high_appropriateness_contexts": ["post_chemotherapy", "gi_toxicity"],
      "default_scores": {"line_appropriateness": 0.8, "cross_resistance": 0.0, "sequencing_fitness": 0.85}
    },
    "Magnesium": {
      "mechanism": "electrolyte_muscle_support",
      "high_appropriateness_contexts": ["post_platinum", "neuropathy"],
      "default_scores": {"line_appropriateness": 0.75, "cross_resistance": 0.0, "sequencing_fitness": 0.80}
    },
    "Vitamin E": {
      "mechanism": "antioxidant",
      "high_appropriateness_contexts": ["oxidative_stress"],
      "default_scores": {"line_appropriateness": 0.6, "cross_resistance": 0.0, "sequencing_fitness": 0.65}
    },
    "Vitamin C": {
      "mechanism": "antioxidant_immune",
      "high_appropriateness_contexts": ["oxidative_stress", "immune_support"],
      "default_scores": {"line_appropriateness": 0.7, "cross_resistance": 0.0, "sequencing_fitness": 0.75}
    },
    "Alpha-lipoic acid": {
      "mechanism": "antioxidant_neuropathy",
      "high_appropriateness_contexts": ["neuropathy", "post_platinum"],
      "default_scores": {"line_appropriateness": 0.8, "cross_resistance": 0.0, "sequencing_fitness": 0.85}
    },
    "L-glutamine": {
      "mechanism": "mucosal_protection",
      "high_appropriateness_contexts": ["gi_toxicity", "post_chemotherapy"],
      "default_scores": {"line_appropriateness": 0.8, "cross_resistance": 0.0, "sequencing_fitness": 0.85}
    }
  },
  "default_supplement": {
    "line_appropriateness": 0.6,
    "cross_resistance": 0.0,
    "sequencing_fitness": 0.6
  }
}
```

**Acceptance:**
- [ ] 20+ compounds have SAE rules
- [ ] Coverage for common supplements (Vitamin C/E, Selenium, Zinc, CoQ10, etc.)

---

### **Fix 4: Add Dynamic Recommendations Fallback (30 min)** 🚨

**Problem:** Timing/meal planning only works for ~10 compounds (hardcoded patterns)

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/dietician_recommendations.py`

**Required Fix:** Add LLM fallback for unknown compounds

```python
def generate_timing_recommendations(self, compound: str, evidence: Dict[str, Any]) -> Dict[str, Any]:
    """Generate timing recommendations with LLM fallback."""
    
    timing = {
        "best_time": "As directed",
        "with_food": True,
        "timing_rationale": "",
        "meal_suggestions": []
    }
    
    compound_lower = compound.lower()
    
    # Hardcoded patterns for common compounds (existing code)
    if any(x in compound_lower for x in ["vitamin d", "vitamin a", "vitamin e", "vitamin k"]):
        timing.update({
            "best_time": "Morning with breakfast",
            "with_food": True,
            "timing_rationale": "Fat-soluble vitamins require dietary fat for optimal absorption",
            "meal_suggestions": ["Eggs", "Avocado", "Nuts", "Oily fish"]
        })
        return timing
    
    # ... (existing patterns) ...
    
    # LLM FALLBACK for unknown compounds
    if timing["best_time"] == "As directed" and LLM_AVAILABLE:
        try:
            llm_service = get_llm_service()
            
            prompt = f"""Generate timing recommendations for {compound}:

Return JSON:
{{
  "best_time": "morning|evening|with meals",
  "with_food": true|false,
  "timing_rationale": "brief explanation",
  "meal_suggestions": ["food1", "food2"]
}}

If unknown, return generic defaults."""
            
            response = await llm_service.query(prompt)
            llm_timing = json.loads(response)
            
            if llm_timing.get("best_time"):
                timing.update(llm_timing)
        except Exception as e:
            print(f"⚠️ LLM timing generation failed: {e}")
    
    return timing
```

**Acceptance:**
- [ ] LLM generates recommendations for unknown compounds
- [ ] Fallback to generic "As directed" if LLM fails
- [ ] Hardcoded patterns still work for common compounds

---

## **✅ PHASE 1: P/E/SAE MVP (ALREADY BUILT)**

### **What Agent Jr Built (Status: ✅ COMPLETE):**

#### **1. Data Files (30 min)**
Create these files with complete data:

**File:** `.cursor/ayesha/hypothesis_validator/data/food_targets.json`
```json
{
  "compounds": [
    {
      "compound": "Vitamin D",
      "targets": ["VDR", "TP53", "BRCA1"],
      "pathways": ["DNA repair", "Cell cycle regulation"],
      "mechanism": "VDR-mediated gene regulation + DNA repair support"
    },
    {
      "compound": "NAC",
      "targets": ["GSH synthesis", "NRF2"],
      "pathways": ["Oxidative stress response", "Glutathione metabolism"],
      "mechanism": "Precursor to glutathione, antioxidant support post-chemo"
    },
    {
      "compound": "Curcumin",
      "targets": ["NFKB", "COX2", "STAT3"],
      "pathways": ["Inflammation", "NFkB signaling"],
      "mechanism": "NFkB inhibition, anti-inflammatory"
    },
    {
      "compound": "Omega-3",
      "targets": ["COX2", "PPARG"],
      "pathways": ["Inflammation", "Lipid metabolism"],
      "mechanism": "Anti-inflammatory via prostaglandin modulation"
    },
    {
      "compound": "Green Tea (EGCG)",
      "targets": ["EGFR", "VEGF", "MMP"],
      "pathways": ["Angiogenesis", "Cell proliferation"],
      "mechanism": "EGFR inhibition, anti-angiogenic"
    },
    {
      "compound": "Folate",
      "targets": ["DHFR", "TYMS"],
      "pathways": ["One-carbon metabolism", "DNA synthesis"],
      "mechanism": "DNA synthesis and repair support"
    }
  ]
}
```

**File:** `.cursor/ayesha/hypothesis_validator/data/supplement_treatment_rules.json`
```json
{
  "supplement_rules": {
    "NAC": {
      "mechanism": "oxidative_stress_recovery",
      "high_appropriateness_contexts": [
        "post_platinum",
        "post_chemotherapy",
        "oxidative_stress_high"
      ],
      "default_scores": {
        "line_appropriateness": 1.0,
        "cross_resistance": 0.0,
        "sequencing_fitness": 0.95
      },
      "biomarker_gates": {
        "chemotherapy_history": "platinum-based"
      }
    },
    "Vitamin D": {
      "mechanism": "dna_repair_support",
      "high_appropriateness_contexts": [
        "hrd_positive",
        "dna_repair_deficient",
        "brca1_brca2_mutated"
      ],
      "default_scores": {
        "line_appropriateness": 0.9,
        "cross_resistance": 0.0,
        "sequencing_fitness": 0.85
      },
      "biomarker_gates": {
        "HRD": "POSITIVE"
      }
    },
    "Omega-3": {
      "mechanism": "anti_inflammatory",
      "high_appropriateness_contexts": [
        "post_chemotherapy",
        "chronic_inflammation",
        "post_immunotherapy"
      ],
      "default_scores": {
        "line_appropriateness": 0.85,
        "cross_resistance": 0.0,
        "sequencing_fitness": 0.80
      }
    },
    "Curcumin": {
      "mechanism": "nfkb_inhibition",
      "high_appropriateness_contexts": [
        "chronic_inflammation",
        "nfkb_active"
      ],
      "default_scores": {
        "line_appropriateness": 0.7,
        "cross_resistance": 0.0,
        "sequencing_fitness": 0.75
      }
    },
    "Green Tea (EGCG)": {
      "mechanism": "egfr_inhibition",
      "high_appropriateness_contexts": [
        "egfr_active",
        "angiogenesis_high"
      ],
      "default_scores": {
        "line_appropriateness": 0.65,
        "cross_resistance": 0.0,
        "sequencing_fitness": 0.70
      }
    },
    "Folate": {
      "mechanism": "dna_synthesis_support",
      "high_appropriateness_contexts": [
        "dna_repair_deficient",
        "folate_deficiency"
      ],
      "default_scores": {
        "line_appropriateness": 0.80,
        "cross_resistance": 0.0,
        "sequencing_fitness": 0.75
      }
    }
  },
  "default_supplement": {
    "line_appropriateness": 0.6,
    "cross_resistance": 0.0,
    "sequencing_fitness": 0.6
  }
}
```

#### **2. Backend Services (3 hours)**

**Service 1:** `oncology-coPilot/oncology-backend-minimal/api/services/food_treatment_line_service.py`
- Copy implementation from `EXECUTION_DECISIONS.md` Q3/Q6
- Load supplement rules JSON
- Apply biomarker + treatment history gating
- Return SAE features: `line_appropriateness`, `cross_resistance`, `sequencing_fitness`

**Service 2:** `oncology-coPilot/oncology-backend-minimal/api/services/compound_target_extraction.py`
```python
import json
from typing import Dict, Any

def extract_targets_phase1(compound: str) -> Dict[str, Any]:
    """
    Phase 1: Use food_targets.json only.
    Fast, deterministic, works for 6 compounds.
    """
    with open('.cursor/ayesha/hypothesis_validator/data/food_targets.json') as f:
        food_targets = json.load(f)
    
    for item in food_targets['compounds']:
        if compound.lower() in item['compound'].lower():
            return {
                "compound": item['compound'],
                "targets": item['targets'],
                "pathways": item.get('pathways', []),
                "mechanism": item.get('mechanism', ''),
                "source": "knowledge_base",
                "confidence": 0.95
            }
    
    return {
        "compound": compound,
        "targets": [],
        "source": "not_found",
        "confidence": 0.0,
        "error": f"Compound '{compound}' not in knowledge base"
    }
```

**Service 3:** `oncology-coPilot/oncology-backend-minimal/api/services/food_spe_integration.py`
- Copy full implementation from `EXECUTION_DECISIONS.md` Q4
- Class: `FoodSPEIntegrationService`
- Formula: `0.4×S + 0.3×P + 0.3×E` (S=0.5 neutral for Phase 1)
- Methods:
  - `compute_spe_score()` - main orchestrator
  - `_compute_pathway_alignment()` - keyword matching
  - `_convert_evidence_grade()` - STRONG/MODERATE/WEAK → 0-1
  - `_compute_confidence()` - multi-stage boosts
  - `_classify_verdict()` - SUPPORTED/WEAK_SUPPORT/NOT_SUPPORTED

#### **3. Endpoint (1 hour)**

**File:** `oncology-coPilot/oncology-backend-minimal/api/routers/hypothesis_validator.py`

Add this endpoint:
```python
@router.post("/validate_food_complete")
async def validate_food_complete(request: Dict[str, Any]):
    """
    Complete food/supplement validation with P/E/SAE.
    
    Request:
      {
        "compound": "Vitamin D",
        "disease_context": {
          "disease": "ovarian_cancer_hgs",
          "mutations": [{"gene": "TP53", "hgvs_p": "R248Q"}],
          "biomarkers": {"HRD": "POSITIVE", "TMB": 8.2},
          "pathways_disrupted": ["DNA repair", "Cell cycle"]
        },
        "treatment_history": {
          "current_line": "L3",
          "prior_therapies": ["carboplatin", "paclitaxel"]
        },
        "use_evo2": false
      }
    
    Response:
      {
        "overall_score": 0.615,
        "confidence": 0.82,
        "verdict": "WEAK_SUPPORT",
        "spe_breakdown": {"sequence": 0.5, "pathway": 0.85, "evidence": 0.60},
        "sae_features": {...},
        "evo2_analysis": {"enabled": false},
        "provenance": {...}
      }
    """
    
    # 1. Extract targets
    compound = request['compound']
    target_data = extract_targets_phase1(compound)
    
    if not target_data['targets']:
        return {
            "status": "ERROR",
            "error": f"Compound '{compound}' not in knowledge base"
        }
    
    # 2. Get evidence grade from LLM literature service
    # (fallback to "INSUFFICIENT" if unavailable)
    try:
        evidence_grade = await get_evidence_grade(compound, request['disease_context'])
    except Exception:
        evidence_grade = "INSUFFICIENT"
    
    # 3. Compute SAE features
    sae_features = None
    if request.get('treatment_history'):
        sae_features = compute_food_treatment_line_features(
            compound=compound,
            disease_context=request['disease_context'],
            treatment_history=request['treatment_history']
        )
    
    # 4. Compute S/P/E score
    spe_service = FoodSPEIntegrationService()
    result = await spe_service.compute_spe_score(
        compound=compound,
        targets=target_data['targets'],
        pathways=target_data['pathways'],
        disease_context=request['disease_context'],
        evidence_grade=evidence_grade,
        treatment_history=request.get('treatment_history'),
        evo2_enabled=request.get('use_evo2', False)  # Phase 1: always False
    )
    
    # 5. Add provenance
    result['provenance'] = {
        "run_id": str(uuid.uuid4()),
        "profile": {"use_evo2": request.get('use_evo2', False)},
        "sources": ["food_targets.json", "llm_literature_service"],
        "timestamp": datetime.utcnow().isoformat()
    }
    
    return result
```

#### **4. Testing (1 hour)**

**Test File:** `.cursor/ayesha/hypothesis_validator/evo2_food_validator_doctrine/testing/test_phase1_mvp.py`
```python
import requests
import json

BASE_URL = "http://127.0.0.1:8000"

def test_vitamin_d_ayesha_case():
    """Test Vitamin D for Ayesha's ovarian cancer case."""
    
    payload = {
        "compound": "Vitamin D",
        "disease_context": {
            "disease": "ovarian_cancer_hgs",
            "mutations": [{"gene": "TP53", "hgvs_p": "R248Q"}],
            "biomarkers": {"HRD": "POSITIVE", "TMB": 8.2},
            "pathways_disrupted": ["DNA repair", "Cell cycle"]
        },
        "treatment_history": {
            "current_line": "L3",
            "prior_therapies": ["carboplatin", "paclitaxel"]
        },
        "use_evo2": False
    }
    
    response = requests.post(
        f"{BASE_URL}/api/hypothesis/validate_food_complete",
        json=payload
    )
    
    assert response.status_code == 200
    result = response.json()
    
    # Verify structure
    assert "overall_score" in result
    assert "confidence" in result
    assert "verdict" in result
    assert "spe_breakdown" in result
    assert "sae_features" in result
    
    # Verify values make sense
    assert 0 <= result['overall_score'] <= 1
    assert 0 <= result['confidence'] <= 1
    assert result['verdict'] in ["SUPPORTED", "WEAK_SUPPORT", "NOT_SUPPORTED"]
    
    # Verify SAE features present for Ayesha's case
    assert result['sae_features']['line_appropriateness'] >= 0.8  # Post-platinum, HRD+
    
    print(f"✅ Vitamin D Test PASSED")
    print(f"   Score: {result['overall_score']:.3f}, Confidence: {result['confidence']:.3f}")
    print(f"   Verdict: {result['verdict']}")
    
    return result

def test_nac_post_chemo():
    """Test NAC for post-chemotherapy support."""
    
    payload = {
        "compound": "NAC",
        "disease_context": {
            "disease": "ovarian_cancer_hgs",
            "biomarkers": {},
            "pathways_disrupted": ["Oxidative stress"]
        },
        "treatment_history": {
            "current_line": "L3",
            "prior_therapies": ["carboplatin"]
        },
        "use_evo2": False
    }
    
    response = requests.post(
        f"{BASE_URL}/api/hypothesis/validate_food_complete",
        json=payload
    )
    
    assert response.status_code == 200
    result = response.json()
    
    # NAC should score HIGH for post-platinum
    assert result['sae_features']['line_appropriateness'] >= 0.9
    
    print(f"✅ NAC Test PASSED")
    print(f"   Score: {result['overall_score']:.3f}, Confidence: {result['confidence']:.3f}")
    
    return result

if __name__ == "__main__":
    print("🧪 PHASE 1 MVP TESTING")
    print("=" * 50)
    
    try:
        test_vitamin_d_ayesha_case()
        print()
        test_nac_post_chemo()
        
        print("\n✅ ALL PHASE 1 TESTS PASSED!")
        print("\n📋 READY FOR FRONTEND INTEGRATION")
        
    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
```

#### **5. Frontend Update (1 hour)**

**File:** `oncology-coPilot/oncology-frontend/src/pages/FoodValidatorAB.jsx`

Update to use new endpoint:
```jsx
const handleValidate = async () => {
  setLoading(true);
  
  const payload = {
    compound: selectedCompound,
    disease_context: {
      disease: "ovarian_cancer_hgs",
      mutations: [{ gene: "TP53", hgvs_p: "R248Q" }],
      biomarkers: { HRD: "POSITIVE", TMB: 8.2 },
      pathways_disrupted: ["DNA repair", "Cell cycle"]
    },
    treatment_history: {
      current_line: "L3",
      prior_therapies: ["carboplatin", "paclitaxel"]
    },
    use_evo2: false  // Phase 1: Evo2 disabled
  };
  
  try {
    const response = await fetch(
      `${API_BASE}/api/hypothesis/validate_food_complete`,
      {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload)
      }
    );
    
    const result = await response.json();
    setResults(result);
    
  } catch (error) {
    console.error("Validation error:", error);
  } finally {
    setLoading(false);
  }
};
```

Display results:
```jsx
{results && (
  <div className="results-section">
    <div className="verdict-card">
      <h3>{results.verdict}</h3>
      <p>Score: {(results.overall_score * 100).toFixed(1)}%</p>
      <p>Confidence: {(results.confidence * 100).toFixed(1)}%</p>
    </div>
    
    <div className="spe-breakdown">
      <h4>S/P/E Breakdown</h4>
      <div className="score-bar">
        <span>Sequence (S):</span>
        <span>{(results.spe_breakdown.sequence * 100).toFixed(0)}%</span>
      </div>
      <div className="score-bar">
        <span>Pathway (P):</span>
        <span>{(results.spe_breakdown.pathway * 100).toFixed(0)}%</span>
      </div>
      <div className="score-bar">
        <span>Evidence (E):</span>
        <span>{(results.spe_breakdown.evidence * 100).toFixed(0)}%</span>
      </div>
    </div>
    
    <div className="sae-features">
      <h4>Treatment Line Intelligence</h4>
      <div className="sae-chip">
        Line Appropriateness: {(results.sae_features.line_appropriateness * 100).toFixed(0)}%
      </div>
      <div className="sae-chip">
        Sequencing Fitness: {(results.sae_features.sequencing_fitness * 100).toFixed(0)}%
      </div>
    </div>
  </div>
)}
```

---

## **🧬 PHASE 2: EVO2 TOGGLE (2-3 HOURS) - BUILD AFTER PHASE 1 VALIDATES**

### **Only build this if Phase 1 works!**

Add Evo2 service:

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/evo2_food_plausibility.py`
```python
async def score_compound_plausibility(compound, target_gene, disease_context):
    """
    Score compound effect using promoter variant proxy.
    
    EXPERIMENTAL: 60% confidence this approach works.
    """
    
    # Step 1: Get gene TSS (transcription start site)
    gene_info = await fetch_ensembl_gene_info(target_gene)
    tss = gene_info['start']  # GRCh38 coordinate
    chrom = gene_info['seq_region_name']
    
    # Step 2: Create synthetic variant at promoter (TSS - 500bp)
    promoter_pos = tss - 500
    
    # Step 3: Get reference allele at that position
    ref_allele = await fetch_reference_allele(chrom, promoter_pos)
    alt_allele = "G" if ref_allele != "G" else "A"
    
    # Step 4: Score variant using Evo2
    response = await http_client.post(
        f"{API_BASE}/api/evo/score_variant_multi",
        json={
            "assembly": "GRCh38",
            "chrom": chrom,
            "pos": promoter_pos,
            "ref": ref_allele,
            "alt": alt_allele,
            "model_id": "evo2_1b"
        }
    )
    
    result = response.json()
    delta = abs(result.get('min_delta', 0))
    
    # Classify plausibility
    if delta > 0.5:
        plausibility = "HIGH"
    elif delta > 0.2:
        plausibility = "MODERATE"
    else:
        plausibility = "LOW"
    
    return {
        "gene": target_gene,
        "delta": delta,
        "plausibility": plausibility,
        "proxy_method": "promoter_variant",
        "coordinates": {"chrom": chrom, "pos": promoter_pos}
    }
```

Update `food_spe_integration.py` to use Evo2 when `evo2_enabled=True`.

Update frontend to add Evo2 toggle (default OFF).

---

## **🎯 AYESHA-SPECIFIC RECOMMENDATIONS**

### **For Ayesha's Case (Ovarian HGS, TP53 mutant, HRD+, L3 post-platinum):**

**Top Recommendations (expected from MVP):**

1. **Vitamin D** - WEAK_SUPPORT
   - Score: ~0.60-0.65
   - Confidence: ~0.80-0.85
   - Rationale: HRD+ → DNA repair support; line appropriateness HIGH (0.9)
   - Dosage: 2000-4000 IU daily
   - Timing: During/after L3 therapy

2. **NAC** - SUPPORTED
   - Score: ~0.70-0.75
   - Confidence: ~0.85-0.90
   - Rationale: Post-platinum oxidative stress; line appropriateness VERY HIGH (1.0)
   - Dosage: 600-1200mg daily
   - Timing: Post-chemotherapy

3. **Omega-3** - WEAK_SUPPORT
   - Score: ~0.55-0.60
   - Confidence: ~0.70-0.75
   - Rationale: Anti-inflammatory support; moderate line appropriateness (0.85)
   - Dosage: 2-3g EPA+DHA daily
   - Timing: Ongoing

**Not Recommended:**

4. **Curcumin** - NOT_SUPPORTED
   - Score: ~0.40-0.45
   - Confidence: ~0.50-0.60
   - Rationale: NFkB not primary pathway for ovarian HGS; lower appropriateness (0.7)

---

## **✅ ACCEPTANCE CRITERIA**

### **Phase 1 MVP (Must Have):**
- [ ] All 6 compounds return valid scores (0-1 range)
- [ ] Vitamin D scores HIGH for Ayesha's case (HRD+, post-platinum)
- [ ] NAC scores VERY HIGH for post-chemo context
- [ ] SAE features boost confidence for appropriate contexts
- [ ] P/E/SAE breakdown displayed correctly
- [ ] Frontend shows verdict + scores + SAE chips
- [ ] Response time < 2 seconds

### **Phase 2 Evo2 (Nice to Have):**
- [ ] Evo2 toggle works without breaking Phase 1
- [ ] Promoter variant scoring returns non-zero deltas
- [ ] Evo2 S component changes overall score
- [ ] Provenance tracks Evo2 usage

---

## **🚨 CRITICAL NOTES FOR AGENT JR**

### **DO:**
- ✅ Build Phase 1 first (P/E/SAE MVP)
- ✅ Test with Ayesha's case (HRD+, L3, post-platinum)
- ✅ Use exact implementations from `EXECUTION_DECISIONS.md`
- ✅ Keep `use_evo2=false` default
- ✅ Add proper error handling and fallbacks
- ✅ Track provenance (run_id, profile, sources)

### **DON'T:**
- ❌ Build Evo2 before Phase 1 validates
- ❌ Overcomplicate the pathway alignment (keyword matching is fine)
- ❌ Reuse drug efficacy orchestrator (too complex)
- ❌ Skip SAE features (critical differentiation)
- ❌ Hardcode Ayesha's data in backend (parameterize)

---

## **📊 EXPECTED TIMELINE**

**Phase 1:**
- Data files: 30 min
- Services: 3 hours
- Endpoint: 1 hour
- Testing: 1 hour
- Frontend: 1 hour
**Total: 4-5 hours**

**Phase 2 (if Phase 1 passes):**
- Evo2 service: 1.5 hours
- Integration: 1 hour
- Testing: 30 min
**Total: 2-3 hours**

---

## **🎯 SUCCESS = AYESHA GETS ACTIONABLE RECOMMENDATIONS**

**The prize:** Ayesha receives:
1. Clear verdict (SUPPORTED/WEAK_SUPPORT/NOT_SUPPORTED)
2. Confidence score (how sure are we?)
3. S/P/E breakdown (why this score?)
4. Treatment line intelligence (why now?)
5. Dosage + timing guidance

**Differentiation from PubMed:**
- Integrated P/E/SAE scoring ✅
- Biomarker gating ✅
- Treatment line intelligence ✅
- SAE confidence modulation ✅
- Automated verdict classification ✅

---

**AGENT JR: This is your build plan. Execute Phase 1 first. Report results. Then we decide on Phase 2.** ⚔️

**Commander Zo - Battle Plan Complete**

