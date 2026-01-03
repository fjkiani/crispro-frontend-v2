# Food Validation Query Format: Current State vs. Desired Format

## Your Query Format

```json
{
  "compounds": [
    "Vitamin D3 5000 IU",
    "Vitamin C 2000mg",
    "Alpha-lipoic acid 600mg",
    "Magnesium glycinate 400mg",
    "Fish oil 3g omega-3",
    "Probiotics 50 billion CFU",
    "Black seed oil (thymoquinone)"
  ],
  "disease": "ovarian_cancer_hgs",
  "patient_context": {
    "treatment_line": "first-line",
    "current_drugs": ["carboplatin_AUC5", "paclitaxel_175mg"],
    "biomarkers": {
      "TP53": "MUTANT",
      "MBD4": "HOMOZYGOUS_LOSS",
      "HRD": "UNKNOWN",
      "TMB": "PRESUMED_HIGH"
    },
    "cycle": 1,
    "day_post_chemo": 7
  }
}
```

## Current API Format (Single Compound)

**Endpoint**: `POST /api/hypothesis/validate_food_dynamic`

**Request** (per compound):
```json
{
  "compound": "Vitamin D3",  // Base name (dosage extracted separately)
  "disease_context": {
    "disease": "ovarian_cancer_hgs",
    "mutations": [
      {"gene": "TP53", "status": "MUTANT"},
      {"gene": "MBD4", "status": "HOMOZYGOUS_LOSS"}
    ],
    "biomarkers": {
      "HRD": "UNKNOWN",
      "TMB": "PRESUMED_HIGH"
    }
  },
  "treatment_history": {
    "current_line": "first-line",  // Will normalize to "L1"
    "prior_therapies": []
  },
  "patient_medications": ["carboplatin_AUC5", "paclitaxel_175mg"],
  "use_evo2": false
}
```

**Response** (per compound):
```json
{
  "status": "SUCCESS",
  "compound": "Vitamin D3",
  "alignment_score": 0.72,
  "overall_score": 0.72,  // Backward compatibility
  "confidence": 0.68,
  "verdict": "SUPPORTED",
  
  "spe_breakdown": {
    "sequence": 0.5,      // Neutral (Evo2 disabled)
    "pathway": 0.75,      // TCGA-weighted alignment
    "evidence": 0.70      // Literature strength
  },
  
  "sae_features": {
    "line_fitness": {
      "score": 0.9,
      "status": "appropriate",
      "reason": "Compatible with first-line platinum chemotherapy"
    },
    "cross_resistance": {
      "risk": "LOW",
      "score": 0.0,
      "reason": "No prior therapies to assess cross-resistance"
    },
    "sequencing_fitness": {
      "score": 0.85,
      "optimal": true,
      "reason": "Appropriate timing for first-line treatment phase"
    }
  },
  
  "targets": ["VDR", "DNA repair pathways"],
  "pathways": ["dna_repair", "immune_modulation"],
  "mechanisms": ["dna_repair_support", "immune_boost"],
  
  "evidence": {
    "papers": [...],
    "evidence_grade": "MODERATE",
    "total_papers": 15,
    "rct_count": 2
  },
  
  "dietician_recommendations": {
    "dosage": {
      "recommended": "2000-4000 IU daily",
      "user_specified": "5000 IU",  // Extracted from input
      "note": "User specified 5000 IU - monitor for hypercalcemia"
    },
    "timing": {
      "relative_to_chemo": "Day 7 post-infusion: Safe to resume",
      "optimal_timing": "With meals, preferably morning"
    },
    "interactions": {
      "carboplatin": "No known interaction",
      "paclitaxel": "No known interaction"
    },
    "monitoring": {
      "labs": ["25-OH vitamin D", "Calcium"],
      "frequency": "Every 3 months"
    }
  },
  
  "llm_enhancements": {
    "personalized_rationale": "...",
    "mechanism_synthesis": "...",
    "evidence_interpretation": "..."
  }
}
```

## Format Differences

| Your Format | Current API | Notes |
|-------------|-------------|-------|
| `compounds` (array) | `compound` (string) | **Need batch endpoint** |
| `current_drugs` | `patient_medications` | ✅ Same concept |
| `treatment_line` | `current_line` | ✅ Same, will normalize |
| `cycle` | ❌ Not used | Could add to timing logic |
| `day_post_chemo` | ❌ Not used | Could add to timing logic |
| Dosage in compound name | Extracted separately | Need to parse "5000 IU" from name |

## What We'd Need to Add

### 1. Batch Endpoint (High Priority)

**New Endpoint**: `POST /api/hypothesis/validate_foods_batch`

```json
{
  "compounds": ["Vitamin D3 5000 IU", "Vitamin C 2000mg", ...],
  "disease": "ovarian_cancer_hgs",
  "patient_context": {...}
}
```

**Response**:
```json
{
  "status": "SUCCESS",
  "results": [
    {
      "compound": "Vitamin D3",
      "user_specified_dosage": "5000 IU",
      "alignment_score": 0.72,
      ...
    },
    {
      "compound": "Vitamin C",
      "user_specified_dosage": "2000mg",
      "alignment_score": 0.65,
      ...
    }
  ],
  "summary": {
    "total_compounds": 7,
    "high_alignment": 3,
    "interactions_detected": 0,
    "warnings": []
  }
}
```

### 2. Dosage Parsing

Extract dosage from compound name:
- "Vitamin D3 5000 IU" → compound: "Vitamin D3", dosage: "5000 IU"
- "Alpha-lipoic acid 600mg" → compound: "Alpha-lipoic acid", dosage: "600mg"

### 3. Cycle/Day Post-Chemo Logic

Use `cycle` and `day_post_chemo` for timing recommendations:
- Day 0-3: Avoid certain supplements (chemo active)
- Day 4-7: Resume with caution
- Day 8+: Full recommendations

## Current Workaround

**For now**, you can:

1. **Loop through compounds**:
```python
results = []
for compound in query["compounds"]:
    api_request = convert_to_api_format(query, compound)
    response = await validate_food_dynamic(api_request)
    results.append(response)
```

2. **Extract dosage manually**:
```python
import re
def parse_compound(compound_str):
    # "Vitamin D3 5000 IU" → ("Vitamin D3", "5000 IU")
    match = re.match(r"(.+?)\s+(\d+\s*(?:mg|IU|g|mcg|CFU))", compound_str)
    if match:
        return match.group(1), match.group(2)
    return compound_str, None
```

## Recommendation

**Ship V1**: Current single-compound endpoint works, just needs client-side batching

**Add V2**: Batch endpoint + dosage parsing + cycle/day logic

This would make it a true "patient supplement stack validator" rather than "single compound checker."

