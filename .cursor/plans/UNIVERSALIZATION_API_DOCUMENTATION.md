# Universal Complete Care API Documentation

**Date:** January 28, 2025  
**Version:** 1.0  
**Status:** ✅ **ENDPOINT ENABLED**

---

## Overview

The Universal Complete Care API provides comprehensive care planning for **any patient** with **any cancer type**. This is the universalized version of the Ayesha orchestrator, making all capabilities available to any user.

**Base URL:** `http://localhost:8000`  
**API Prefix:** `/api/complete_care`

---

## Endpoints

### 1. Complete Care Plan v2

**Endpoint:** `POST /api/complete_care/v2`

**Description:**  
Get a complete care plan for any patient. Orchestrates all care plan components:
- Clinical trials (universal pipeline)
- SOC recommendation (disease-specific)
- Biomarker intelligence (disease-specific)
- Drug efficacy (WIWFM)
- Food validator (optional)
- Resistance playbook (optional)
- Resistance Prophet (optional)
- SAE services (Phase 1 & 2)

**Request Schema:**
```json
{
  "patient_profile": {
    // Simple or full format (see Profile Formats below)
  },
  "drug_query": "optional_drug_name",
  "food_query": "optional_food_name",
  "include_trials": true,
  "include_soc": true,
  "include_biomarker": true,
  "include_wiwfm": true,
  "include_food": false,
  "include_resistance": false,
  "include_resistance_prediction": false
}
```

**Response Schema:**
```json
{
  "summary": {
    "components_included": ["trials", "soc", "biomarker", "wiwfm"],
    "ngs_status": "available",
    "confidence_level": "high"
  },
  "trials": {
    // Clinical trials matching
  },
  "soc_recommendation": {
    "regimen": "Carboplatin + Paclitaxel + Bevacizumab",
    "nccn_category": "Category 1"
  },
  "biomarker_intelligence": {
    "burden_assessment": {},
    "response_forecast": {},
    "resistance_risk": {}
  },
  "wiwfm": {
    // Drug efficacy analysis
  },
  "food_validation": null,
  "resistance_playbook": null,
  "resistance_prediction": null,
  "sae_features": {},
  "provenance": {
    "orchestrator": "universal_complete_care_v2",
    "run_id": "uuid",
    "timestamp": "2025-01-28T..."
  }
}
```

**Example Request:**
```bash
curl -X POST http://localhost:8000/api/complete_care/v2 \
  -H "Content-Type: application/json" \
  -d '{
    "patient_profile": {
      "patient_id": "test_001",
      "name": "Test Patient",
      "disease": "ovarian_cancer_hgs",
      "stage": "IVB",
      "treatment_line": "first-line",
      "location": "New York",
      "zip_code": "10001",
      "biomarkers": {
        "ca125_value": 1500.0
      },
      "tumor_context": {
        "somatic_mutations": [
          {
            "gene": "BRCA1",
            "hgvs_p": "p.Arg1835Ter"
          }
        ]
      }
    },
    "include_trials": true,
    "include_soc": true,
    "include_biomarker": true,
    "include_wiwfm": true
  }'
```

---

### 2. Universal Endpoint (Alias)

**Endpoint:** `POST /api/complete_care/universal`

**Description:**  
Alias for `/api/complete_care/v2`. Same functionality.

---

### 3. Biomarker Analysis

**Endpoint:** `POST /api/biomarker/analyze`

**Description:**  
Analyze biomarker values for any cancer type. Supports CA-125 (ovarian), PSA (prostate), CEA (colorectal), etc.

**Request Schema:**
```json
{
  "disease_type": "ovarian_cancer_hgs",
  "biomarker_type": "ca125",
  "current_value": 1500.0,
  "baseline_value": 35.0,
  "cycle": 2,
  "treatment_ongoing": true
}
```

**Response Schema:**
```json
{
  "burden_assessment": {
    "level": "high",
    "value": 1500.0,
    "threshold": 35.0
  },
  "response_forecast": {
    "trend": "increasing",
    "risk_level": "high"
  },
  "resistance_risk": {
    "probability": 0.75,
    "factors": []
  }
}
```

---

## Profile Formats

### Simple Profile Format

```json
{
  "patient_id": "string",
  "name": "string (optional)",
  "disease": "string | object",
  "stage": "string (optional)",
  "treatment_line": "string",
  "location": "string",
  "zip_code": "string (optional)",
  "age": "number (optional)",
  "sex": "string (optional)",
  "biomarkers": {
    "ca125_value": "number (optional)",
    "psa_value": "number (optional)",
    // ... other biomarkers
  },
  "tumor_context": {
    "somatic_mutations": [
      {
        "gene": "string",
        "hgvs_p": "string (optional)",
        "consequence": "string (optional)"
      }
    ],
    "hrd_score": "number (optional)"
  }
}
```

### Full Profile Format

```json
{
  "patient_id": "string",
  "demographics": {
    "name": "string",
    "age": "number",
    "sex": "string"
  },
  "disease": {
    "type": "string",
    "stage": "string"
  },
  "treatment": {
    "line": "string"
  },
  "logistics": {
    "location": "string",
    "zip_code": "string"
  },
  "biomarkers": {
    // Disease-specific biomarkers
  },
  "tumor_context": {
    "somatic_mutations": [],
    "hrd_score": "number"
  }
}
```

**Note:** The API accepts both formats. Simple profiles are automatically converted to full format.

---

## Supported Disease Types

- `ovarian_cancer_hgs` - Ovarian cancer (high-grade serous)
- `melanoma` - Melanoma
- `multiple_myeloma` - Multiple myeloma
- `breast_cancer` - Breast cancer
- `colorectal_cancer` - Colorectal cancer
- `prostate_cancer` - Prostate cancer
- (More can be added via configuration)

---

## Supported Biomarkers

- `ca125` - CA-125 (ovarian cancer)
- `psa` - PSA (prostate cancer)
- `cea` - CEA (colorectal cancer)
- `afp` - AFP (liver cancer)
- `hcg` - hCG (germ cell tumors)
- (More can be added via configuration)

---

## Error Handling

**400 Bad Request:**
- Invalid request format
- Missing required fields

**422 Unprocessable Entity:**
- Validation errors
- Invalid disease type
- Invalid biomarker type

**500 Internal Server Error:**
- Service unavailable
- Unexpected errors

**Example Error Response:**
```json
{
  "detail": "Invalid disease type: unknown_cancer"
}
```

---

## Rate Limiting

No rate limiting currently implemented. Consider implementing for production.

---

## Authentication

No authentication currently required. Consider adding for production.

---

## Examples

### Example 1: Simple Profile with Ovarian Cancer

```python
import httpx
import asyncio

async def get_care_plan():
    async with httpx.AsyncClient() as client:
        response = await client.post(
            "http://localhost:8000/api/complete_care/v2",
            json={
                "patient_profile": {
                    "patient_id": "patient_001",
                    "name": "Jane Doe",
                    "disease": "ovarian_cancer_hgs",
                    "stage": "IVB",
                    "treatment_line": "first-line",
                    "location": "New York",
                    "zip_code": "10001",
                    "biomarkers": {
                        "ca125_value": 1500.0
                    },
                    "tumor_context": {
                        "somatic_mutations": [
                            {"gene": "BRCA1", "hgvs_p": "p.Arg1835Ter"}
                        ]
                    }
                },
                "include_trials": True,
                "include_soc": True,
                "include_biomarker": True,
                "include_wiwfm": True
            }
        )
        return response.json()

result = asyncio.run(get_care_plan())
print(result)
```

### Example 2: Full Profile with Melanoma

```python
async def get_melanoma_care_plan():
    async with httpx.AsyncClient() as client:
        response = await client.post(
            "http://localhost:8000/api/complete_care/v2",
            json={
                "patient_profile": {
                    "patient_id": "patient_002",
                    "demographics": {
                        "name": "John Smith",
                        "age": 55,
                        "sex": "male"
                    },
                    "disease": {
                        "type": "melanoma",
                        "stage": "III"
                    },
                    "treatment": {
                        "line": "first-line"
                    },
                    "logistics": {
                        "location": "Los Angeles",
                        "zip_code": "90001"
                    }
                }
            }
        )
        return response.json()
```

---

## Testing

See `tests/test_complete_care_universal_integration.py` for integration test examples.

---

**Created:** January 28, 2025  
**Status:** ✅ **ENDPOINT ENABLED**  
**Last Updated:** January 28, 2025


