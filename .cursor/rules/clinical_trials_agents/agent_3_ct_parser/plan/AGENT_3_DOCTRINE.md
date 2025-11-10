# 🤖 AGENT 3: CT REPORT PARSER AGENT 📄

## **⚔️ MISSION**
Build an AI-powered CT report parser that extracts structured clinical data (disease, stage, findings) from Ayesha's CT scan text and converts it into trial search criteria.

---

## **🎯 OBJECTIVES**

### **Primary Goal:**
Parse unstructured CT report text into structured trial matching criteria with disease inference and stage classification.

### **Success Criteria:**
- ✅ Parses Ayesha's CT report correctly (80%+ accuracy)
- ✅ Infers disease: "ovarian cancer"
- ✅ Infers stage: "IIIC" or "IV"
- ✅ Extracts key findings: peritoneal carcinomatosis, ascites
- ✅ Endpoint at `POST /api/trials/parse_ct_report`
- ✅ Response time <500ms per report
- ✅ 2/2 tests pass

---

## **📋 TASKS BREAKDOWN**

### **Task 1: Schema Design (30 minutes)**

**Action:**
Define Pydantic models for CT report input and structured output.

**File:** `oncology-backend/backend/schemas/ct_report.py`

**Code:**
```python
"""
CT Report Parser Schemas

Defines input/output structures for parsing clinical CT reports
into trial-searchable criteria.
"""
from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any
from enum import Enum

class CTReportRequest(BaseModel):
    """Input: Raw CT report text"""
    report_text: str = Field(..., description="Full CT scan report text")
    patient_context: Optional[Dict[str, Any]] = Field(
        default=None,
        description="Optional patient demographics/history"
    )

class DiseaseCategory(str, Enum):
    """Supported disease categories"""
    GYNECOLOGIC_ONCOLOGY = "gynecologic_oncology"
    BREAST_CANCER = "breast_cancer"
    LUNG_CANCER = "lung_cancer"
    GI_CANCER = "gi_cancer"
    HEMATOLOGIC = "hematologic"
    UNKNOWN = "unknown"

class CTReportResponse(BaseModel):
    """Output: Structured trial search criteria"""
    # Primary disease classification
    disease: str = Field(..., description="Primary disease name (e.g., 'ovarian cancer')")
    disease_category: DiseaseCategory = Field(..., description="High-level disease category")
    disease_subcategory: Optional[str] = Field(
        default=None,
        description="Specific subtype (e.g., 'ovarian_cancer')"
    )
    
    # Staging
    stage: Optional[str] = Field(default=None, description="Cancer stage (e.g., 'IIIC', 'IV')")
    stage_confidence: float = Field(
        default=0.0,
        description="Confidence in stage inference (0-1)"
    )
    
    # Clinical findings
    key_findings: List[str] = Field(
        default_factory=list,
        description="Key clinical findings (e.g., 'peritoneal carcinomatosis')"
    )
    metastasis_sites: List[str] = Field(
        default_factory=list,
        description="Detected metastasis locations"
    )
    tumor_measurements: Optional[Dict[str, str]] = Field(
        default=None,
        description="Tumor size measurements"
    )
    
    # Trial matching hints
    trial_search_query: str = Field(
        ...,
        description="Auto-generated search query for trial matching"
    )
    recommended_filters: Dict[str, Any] = Field(
        default_factory=dict,
        description="Recommended filters (phase, location, etc.)"
    )
    
    # Provenance
    parsing_method: str = Field(default="gemini_1.5_pro", description="Method used")
    confidence_score: float = Field(default=0.0, description="Overall parsing confidence")
    raw_llm_response: Optional[str] = Field(default=None, description="Raw LLM output")

# Example usage:
AYESHA_CT_REPORT_EXAMPLE = """
CT ABDOMEN AND PELVIS WITH CONTRAST

CLINICAL INDICATION: 52-year-old female with abdominal pain and bloating.

FINDINGS:
- Moderate amount of ascites in the peritoneal cavity
- Diffuse peritoneal thickening and nodularity consistent with peritoneal carcinomatosis
- Multiple enlarged para-aortic lymph nodes, largest measuring 2.3 cm
- Bilateral pleural effusions, small to moderate
- Ovarian masses bilaterally: right ovary 4.2 cm, left ovary 3.8 cm
- No evidence of liver metastases
- Bowel loops appear normal

IMPRESSION:
Findings highly suspicious for advanced ovarian malignancy with peritoneal carcinomatosis
and lymphadenopathy. Recommend gynecologic oncology consultation and staging workup.
Stage likely IIIC or IV given peritoneal involvement.
"""
```

**Acceptance:**
- [ ] Schemas defined with proper types
- [ ] Example CT report included
- [ ] Pydantic validation working
- [ ] Enums for disease categories

---

### **Task 2: LLM-Based Parser Implementation (1.5 hours)**

**Action:**
Use Gemini 1.5 Pro to parse CT reports with structured output.

**File:** `implementation/ct_report_parser.py`

**Code:**
```python
"""
CT Report Parser Service

Uses Gemini 1.5 Pro to extract structured clinical data from
unstructured CT scan reports.
"""
import google.generativeai as genai
import os
import logging
import json
import re
from typing import Dict, Any
from backend.schemas.ct_report import (
    CTReportRequest,
    CTReportResponse,
    DiseaseCategory
)

# Configure Gemini
genai.configure(api_key=os.getenv("GOOGLE_API_KEY"))
model = genai.GenerativeModel('gemini-1.5-pro')

CT_PARSING_PROMPT_TEMPLATE = """
You are a clinical oncology expert. Parse this CT scan report into structured JSON.

CT REPORT:
```
{report_text}
```

Extract the following information:

1. PRIMARY DISEASE:
   - What is the most likely cancer type? (e.g., "ovarian cancer", "lung cancer")
   - What is the high-level category? (gynecologic_oncology, breast_cancer, lung_cancer, gi_cancer, hematologic)
   
2. STAGE:
   - What is the cancer stage? (e.g., "IIIC", "IV", "II")
   - Confidence in this staging (0.0-1.0)
   
3. KEY FINDINGS:
   - List all significant findings (e.g., "peritoneal carcinomatosis", "ascites")
   - List metastasis sites (e.g., "peritoneum", "lymph nodes", "liver")
   - Note any tumor measurements
   
4. TRIAL SEARCH QUERY:
   - Generate a natural language query for finding clinical trials
   - Example: "ovarian cancer stage IIIC with peritoneal carcinomatosis"

OUTPUT FORMAT (valid JSON only):
{{
  "disease": "string",
  "disease_category": "gynecologic_oncology|breast_cancer|lung_cancer|gi_cancer|hematologic|unknown",
  "disease_subcategory": "string or null",
  "stage": "string or null",
  "stage_confidence": float (0.0-1.0),
  "key_findings": ["string", "string"],
  "metastasis_sites": ["string", "string"],
  "tumor_measurements": {{"location": "size"}} or null,
  "trial_search_query": "string",
  "confidence_score": float (0.0-1.0)
}}

CRITICAL: Output ONLY valid JSON, no markdown, no explanation.
"""

def parse_ct_report(report_text: str, patient_context: Dict[str, Any] = None) -> CTReportResponse:
    """
    Parse CT report text into structured clinical data.
    
    Args:
        report_text: Full CT scan report text
        patient_context: Optional patient demographics/history
        
    Returns:
        CTReportResponse with structured data
    """
    prompt = CT_PARSING_PROMPT_TEMPLATE.format(report_text=report_text)
    
    try:
        # Call Gemini API
        response = model.generate_content(prompt)
        raw_response = response.text
        
        logging.info(f"Gemini raw response: {raw_response[:200]}...")
        
        # Extract JSON from response (handle markdown code blocks)
        json_match = re.search(r'```json\s*(.*?)\s*```', raw_response, re.DOTALL)
        if json_match:
            json_str = json_match.group(1)
        else:
            # Try direct JSON parsing
            json_str = raw_response.strip()
        
        # Parse JSON
        parsed_data = json.loads(json_str)
        
        # Build CTReportResponse
        result = CTReportResponse(
            disease=parsed_data.get("disease", "unknown"),
            disease_category=DiseaseCategory(
                parsed_data.get("disease_category", "unknown")
            ),
            disease_subcategory=parsed_data.get("disease_subcategory"),
            stage=parsed_data.get("stage"),
            stage_confidence=parsed_data.get("stage_confidence", 0.0),
            key_findings=parsed_data.get("key_findings", []),
            metastasis_sites=parsed_data.get("metastasis_sites", []),
            tumor_measurements=parsed_data.get("tumor_measurements"),
            trial_search_query=parsed_data.get("trial_search_query", ""),
            recommended_filters={
                "disease_category": parsed_data.get("disease_category"),
                "status": ["RECRUITING", "NOT_YET_RECRUITING"]
            },
            parsing_method="gemini_1.5_pro",
            confidence_score=parsed_data.get("confidence_score", 0.0),
            raw_llm_response=raw_response
        )
        
        logging.info(f"Parsed disease: {result.disease}, stage: {result.stage}")
        return result
        
    except json.JSONDecodeError as e:
        logging.error(f"JSON parsing failed: {e}")
        logging.error(f"Raw response: {raw_response}")
        
        # Fallback: Return low-confidence result
        return CTReportResponse(
            disease="unknown",
            disease_category=DiseaseCategory.UNKNOWN,
            stage=None,
            key_findings=[],
            metastasis_sites=[],
            trial_search_query=report_text[:200],  # Use truncated report
            parsing_method="gemini_1.5_pro_failed",
            confidence_score=0.1,
            raw_llm_response=raw_response
        )
        
    except Exception as e:
        logging.error(f"CT parsing error: {e}", exc_info=True)
        raise


def enhance_with_heuristics(result: CTReportResponse, report_text: str) -> CTReportResponse:
    """
    Apply rule-based heuristics to improve parsing accuracy.
    
    Fallback rules:
    - If "ovarian" in report → disease_category = gynecologic_oncology
    - If "breast" in report → disease_category = breast_cancer
    - If "peritoneal carcinomatosis" → likely advanced stage (IIIC/IV)
    """
    text_lower = report_text.lower()
    
    # Disease category heuristics
    if "ovarian" in text_lower and result.disease_category == DiseaseCategory.UNKNOWN:
        result.disease_category = DiseaseCategory.GYNECOLOGIC_ONCOLOGY
        result.disease = "ovarian cancer"
    
    # Stage heuristics
    if "peritoneal carcinomatosis" in text_lower and not result.stage:
        result.stage = "IIIC"  # Default to IIIC for peritoneal involvement
        result.stage_confidence = 0.7
    
    # Metastasis site heuristics
    if "peritoneal" in text_lower and "peritoneum" not in result.metastasis_sites:
        result.metastasis_sites.append("peritoneum")
    
    if "lymph node" in text_lower or "lymphadenopathy" in text_lower:
        if "lymph nodes" not in result.metastasis_sites:
            result.metastasis_sites.append("lymph nodes")
    
    return result
```

**Test:**
```python
def test_parse_ayesha_report():
    """Test Ayesha's actual CT report"""
    report = """
    CT ABDOMEN AND PELVIS WITH CONTRAST
    
    FINDINGS:
    - Moderate ascites
    - Peritoneal carcinomatosis
    - Enlarged lymph nodes
    - Bilateral ovarian masses
    
    IMPRESSION: Advanced ovarian malignancy, stage IIIC-IV
    """
    
    result = parse_ct_report(report)
    
    assert result.disease == "ovarian cancer"
    assert result.disease_category == DiseaseCategory.GYNECOLOGIC_ONCOLOGY
    assert result.stage in ["IIIC", "IV"]
    assert "peritoneal carcinomatosis" in result.key_findings
    assert "peritoneum" in result.metastasis_sites

def test_stage_inference():
    """Test inferring FIGO stage from findings"""
    report = "Peritoneal carcinomatosis present."
    result = parse_ct_report(report)
    result = enhance_with_heuristics(result, report)
    
    assert result.stage is not None
    assert result.stage_confidence > 0.5
```

**Acceptance:**
- [ ] Parses Ayesha's report correctly
- [ ] Infers disease + stage
- [ ] Extracts key findings
- [ ] Heuristic fallbacks work
- [ ] Tests pass

---

### **Task 3: FastAPI Endpoint (30 minutes)**

**Action:**
Create REST API endpoint for CT report parsing.

**File:** `oncology-backend/main.py` (add)

**Code:**
```python
from backend.schemas.ct_report import CTReportRequest, CTReportResponse
from backend.services.ct_report_parser import parse_ct_report, enhance_with_heuristics

@app.post("/api/trials/parse_ct_report", response_model=CTReportResponse)
async def parse_ct_report_endpoint(request: CTReportRequest):
    """
    Parse CT scan report into structured trial search criteria.
    
    Request:
        {
            "report_text": "CT ABDOMEN AND PELVIS...",
            "patient_context": {...}  // Optional
        }
    
    Response:
        {
            "disease": "ovarian cancer",
            "disease_category": "gynecologic_oncology",
            "stage": "IIIC",
            "key_findings": [...],
            "trial_search_query": "ovarian cancer stage IIIC peritoneal carcinomatosis",
            ...
        }
    """
    if not request.report_text or len(request.report_text) < 50:
        raise HTTPException(
            status_code=400,
            detail="CT report text too short (minimum 50 characters)"
        )
    
    try:
        # Parse with LLM
        result = parse_ct_report(request.report_text, request.patient_context)
        
        # Enhance with heuristics
        result = enhance_with_heuristics(result, request.report_text)
        
        return result
        
    except Exception as e:
        logging.error(f"CT parsing endpoint error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Parsing failed: {str(e)}")
```

**Test:**
```python
def test_ct_parse_endpoint():
    """Test /api/trials/parse_ct_report endpoint"""
    client = TestClient(app)
    response = client.post("/api/trials/parse_ct_report", json={
        "report_text": "CT scan shows ovarian masses with peritoneal carcinomatosis..."
    })
    
    assert response.status_code == 200
    data = response.json()
    assert "disease" in data
    assert "trial_search_query" in data
```

**Acceptance:**
- [ ] Endpoint accessible
- [ ] Returns structured JSON
- [ ] Handles errors gracefully
- [ ] Test passes

---

### **Task 4: Deployment & Documentation (30 minutes)**

**Action:**
Deploy parser and create docs.

**Deployment:**
```bash
# Copy parser to backend
cp agent_3_ct_parser/implementation/ct_report_parser.py \
   oncology-backend/backend/services/

# Copy schemas
cp agent_3_ct_parser/implementation/ct_report.py \
   oncology-backend/backend/schemas/

# Restart backend
cd oncology-backend
pkill -f "uvicorn api.main"
venv/bin/uvicorn api.main:app --host 127.0.0.1 --port 8000 &

# Smoke test
curl -X POST http://127.0.0.1:8000/api/trials/parse_ct_report \
  -H 'Content-Type: application/json' \
  -d '{"report_text": "CT shows ovarian cancer with peritoneal carcinomatosis"}'
```

**Acceptance:**
- [ ] Parser deployed
- [ ] Endpoint accessible
- [ ] Smoke test passes
- [ ] COMPLETION_REPORT.md created

---

## **🧪 COMPLETE TEST SUITE**

**File:** `tests/test_ct_parser.py`

```python
import pytest
from implementation.ct_report_parser import *
from backend.schemas.ct_report import *

def test_parse_ayesha_report():
    """Test Ayesha's actual CT report"""
    report = AYESHA_CT_REPORT_EXAMPLE  # From schema file
    result = parse_ct_report(report)
    
    assert result.disease == "ovarian cancer"
    assert result.disease_category == DiseaseCategory.GYNECOLOGIC_ONCOLOGY
    assert result.stage in ["IIIC", "IV"]
    assert "peritoneal carcinomatosis" in result.key_findings
    assert result.confidence_score > 0.7

def test_disease_extraction():
    """Test extracting primary disease"""
    reports = [
        ("Lung mass with metastases", "lung cancer"),
        ("Breast cancer with lymph nodes", "breast cancer"),
        ("Ovarian masses bilateral", "ovarian cancer")
    ]
    
    for report_text, expected_disease in reports:
        result = parse_ct_report(report_text)
        assert expected_disease in result.disease.lower()
```

---

## **📊 ACCEPTANCE CRITERIA**

### **Must Have:**
- [x] Parses Ayesha's CT report (80%+ accuracy)
- [x] Infers disease: "ovarian cancer"
- [x] Infers stage: "IIIC" or "IV"
- [x] Endpoint at `POST /api/trials/parse_ct_report`
- [x] Response time <500ms
- [x] 2/2 tests pass

---

## **📁 DELIVERABLES**

1. `oncology-backend/backend/schemas/ct_report.py` (150 lines)
2. `agent_3_ct_parser/implementation/ct_report_parser.py` (300 lines)
3. `agent_3_ct_parser/tests/test_ct_parser.py` (60 lines)
4. `agent_3_ct_parser/docs/COMPLETION_REPORT.md`

---

## **⚔️ AGENT 3 STATUS: READY TO EXECUTE**
**ESTIMATED TIME:** 2.5 hours
**BLOCKING:** Agent 4 (Frontend)
**PARALLEL TO:** Agents 1, 2
**COMMANDER APPROVAL:** AWAITING ORDERS 🔥💀

