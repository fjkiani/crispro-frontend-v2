# ⚔️ API SPECIFICATIONS - BACKEND ENDPOINTS ⚔️

**Purpose**: Define all API endpoints for dossier generation and review

---

## 📋 **DOSSIER GENERATION API**

### **POST /api/dossiers/generate**

**Request**:
```json
{
  "nct_id": "NCT06819007",
  "patient_id": "ayesha_001",
  "scrape_full": true
}
```

**Response**:
```json
{
  "dossier_id": "uuid",
  "nct_id": "NCT06819007",
  "patient_id": "ayesha_001",
  "generated_at": "2025-01-13T22:00:00Z",
  "sections": {
    "intelligence": {...},
    "why_matters": {...},
    "mechanism": {...},
    "eligibility_table": [...],
    "decision_tree": "...",
    "strategic_implications": {...},
    "tactical_recommendations": [...],
    "clinical_evidence": "...",
    "competitive_positioning": [...],
    "final_recommendation": {...}
  },
  "markdown": "# Dossier content...",
  "pdf_ready": true,
  "confidence_score": 0.92,
  "manual_review_required": false
}
```

**Implementation**:
```python
# File: oncology-backend-minimal/api/routers/dossiers.py
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from api.services.client_dossier_generator import generate_client_dossier

router = APIRouter(prefix="/api/dossiers", tags=["dossiers"])

class DossierGenerateRequest(BaseModel):
    nct_id: str
    patient_id: str
    scrape_full: bool = True

@router.post("/generate")
async def generate_dossier(request: DossierGenerateRequest):
    """Generate trial dossier for patient."""
    # Get patient profile from database
    patient_profile = get_patient_profile(request.patient_id)
    
    # Get trial from AstraDB
    astradb_record = get_trial_from_astradb(request.nct_id)
    
    # Generate dossier
    dossier = await generate_client_dossier(
        nct_id=request.nct_id,
        patient_profile=patient_profile,
        astradb_record=astradb_record,
        scrape_full_page=request.scrape_full
    )
    
    # Store in database (AstraDB collection: clinical_dossiers)
    dossier_id = store_dossier(dossier)
    
    return {
        "dossier_id": dossier_id,
        "nct_id": request.nct_id,
        "patient_id": request.patient_id,
        **dossier
    }
```

---

## 🔍 **DOSSIER REVIEW API (FOR ZO)**

### **GET /api/dossiers/{dossier_id}**

**Response**: Full dossier dict (same as generate response)

### **POST /api/dossiers/{dossier_id}/approve**

**Request**:
```json
{
  "approved": true,
  "edits": [],
  "notes": "Looks good, ready for oncologist"
}
```

**Response**:
```json
{
  "status": "approved",
  "ready_for_oncologist": true,
  "reviewed_at": "2025-01-13T23:00:00Z",
  "reviewed_by": "zo"
}
```

---

## 🔄 **BATCH FILTER API**

### **POST /api/trials/filter-batch**

**Request**:
```json
{
  "trials": [...],  # 50 candidates from Zo
  "patient_profile": {
    "stage": "IVB",
    "treatment_line": "first-line",
    "location": "New York, NY"
  },
  "filters": {
    "stage": "IV",
    "treatment_line": "first-line",
    "recruiting_only": true,
    "geography": "USA"
  }
}
```

**Response**:
```json
{
  "top_tier": [...],      # 5-10 trials (pass ALL filters)
  "good_tier": [...],     # 10-15 trials (pass MOST filters)
  "ok_tier": [...],       # 15-20 trials (interesting but not immediate)
  "rejected": [...]       # 10-20 trials (not eligible)
}
```

---

## 📁 **DATABASE STORAGE**

**Collection**: `clinical_dossiers` (AstraDB)

**Schema**:
```json
{
  "dossier_id": "uuid",
  "nct_id": "NCT06819007",
  "patient_id": "ayesha_001",
  "generated_at": "2025-01-13T22:00:00Z",
  "sections": {...},
  "markdown": "...",
  "approved_by_zo": false,
  "sent_to_oncologist": false,
  "oncologist_response": null,
  "version": "1.0"
}
```

---

**Next**: See [10_FRONTEND_REQUIREMENTS.md](./10_FRONTEND_REQUIREMENTS.md) for UI components

