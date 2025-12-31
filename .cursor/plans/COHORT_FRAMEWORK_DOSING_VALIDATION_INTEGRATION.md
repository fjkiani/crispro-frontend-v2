# Cohort Context Framework → Dosing Guidance Validation Integration

**Created:** January 2025  
**Purpose:** Articulate how the existing cohort context framework supports the dosing guidance validation plan

---

## 🎯 Executive Summary

The **Cohort Context Framework** (`cohort_context_concept.mdc`) provides a complete data acquisition infrastructure that directly supports **Phase 1** (Literature Cases) and **Phase 2** (TCGA Cases) of the dosing guidance validation plan. Instead of building new extraction scripts from scratch, we can **reuse and adapt** existing clients and patterns.

**Key Insight:** The framework's integration pattern (Discovery → Extraction → Transformation → Validation → Integration) maps directly to the validation plan's data acquisition phases.

---

## 📊 Framework Capabilities → Validation Needs Mapping

### Phase 1: Lite20-30)

| Validation Need | Framework Capability | Integration Action |
|----------------|---------------------|-------------------|
| PubMed case reports | `EnhancedPubMedPortal` (already built) | Adapt `search_kelim_researchers()` pattern for pharmacogenomics queries |
| PharmGKB annotations | PubMed search (via framework) | Use framework's PubMed client for PharmGKB literature search |
| CPIC case studies | Manual extraction (not automated) | Framework provides extraction pattern, manual curation required |
| Case parsing | Framework's extraction pattern | Reuse transformation logic from cohort framework |

**Current State:**
- ✅ `EnhancedPubMedPortal` exists at `api/services/research_intelligence/portals/pubmed_enhanced.py`
- ✅ PubMed search infrastructure ready
- ⚠️ Need to adapt queries for pharmacogenomics (DPYD, UGT1A1, TPMT)

**Action Required:**
```python
# Instead of building new extract_literature_cases.py from scratch:
# REUSE: EnhancedPubMedPortal.search() method

from api.services.restelligence.portals.pubmed_enhanced import EnhancedPubMedPortal

portal = EnhancedPubMedPortal()

# Adapt existing search pattern for pharmacogenomics
dpyd_cases = portal.search(
    query='"DPYD deficiency" AND "fluoropyrimidine" AND "case report"',
    max_results=50
)
```

---

### Phase 2: TCGA Cases (N=30-50)

| Validation Need | Framework Capability | Integration Action |
|----------------|---------------------|-------------------|
| TCGA germline variants | GDC API (can be added to framework) | Extend cohort framework with GDC client pattern |
| cBioPortal pharmacogenomics | `CBioportalClient` (already built) | Use existing client for MSK-IMPACT, Foundation cohorts |
| Treatment history | cBioPortal clinical data | Framework's `get_clinical_data()` method |
| Project Data Sphere | `ProjectDataSphereClient` (already connected) | Explore PDS for colorectal/breast trials with DPYD/UGT1A1 data |

**Current State:**
- ✅ `CBioportalClient` exists at `scripts/data_acquisition/utils/cbioportal_client.py`
- `ProjectDataSphereClient` exists and connected
- ⚠️ GDC API client not yet in framework (can be added)

**Action Required:**
```python
# Instead of building new extract_tcga_pharmacogenes.py from scratch:
# REUSE: CBioportalClient for studies with pharmacogenomics data

from scripts.data_acquisition.utils.cbioportal_client import CBioportalClient

client = CBioportalClient()

# Use existing client for MSK-IMPACT (has germline data)
studies = client.list_studies()
msk_studies = [s for s in studies if 'msk_impact' in s['study_id'].lower()]

# Extract clinical data (includes pharmacogenomics if available)
for study in msk_studies:
    clinical = client.get_clinical_data(study['study_id'], entity_type="PATIENT")
    # Filter for DPYD/UGT1A1/TPMT variants
```

---

### Phase 3: Validation Metrics

| Validation Need | Framework Capability | Integration Action |
|----------------|---------------------|-------------------|
| Standardized schema | Cohort JSON schema pattern | Extend cohort schema to include phargenomics fields |
| Quality validation | Framework's validation pattern | Reuse data quality checks |
| Integration harness | Framework's integration pattern | Feed validation cases into metrics calculator |

**Current State:**
- ✅ Framework defines JSON schema pattern
- ✅ Framework provides validation pattern
- ⚠️ Need to extend schema for pharmacogenomics fields

**Action Required:**
```json
// Extend cohort schema to include pharmacogenomics:
{
  "cohort": {
    "source": "cbioportal",
    "study_id": "msk_impact_2017",
    "patients": [
      {
        "patient_id": "P001",
        "pharmacogenomics": {
          "gene": "DPYD",
          "variant": "c.1905+1G>A (*2A)",
          "zygosity": "heterozygous",
          "predicted_phenotype": "Intermediate Metabolizer"
        },
        "treatment": {
          "drug": "5-fluorouracil",
          "standard_dose": "400 mg/m²",
          "actual_dose_given": "400 mg/m²"
        },
        "outcome": {
          "toxicity_occurred": true,
         _grade": 4,
          "toxicity_type": "neutropenia"
        }
      }
    ],
    "metadata": {
      "extraction_date": "2025-01-28",
      "data_quality": "high"
    }
  }
}
```

---

## 🔄 Integration Pattern Application

The framework's **5-step integration pattern** maps directly to validation phases:

### Step 1: Discovery
**Framework Pattern:**
```python
# List/search available data
studies = client.list_studies()
```

**Validation Application:**
- **Phase 1:** Search PubMed for case reports (`portal.search()`)
- **Phase 2:** List TCGA projects with germline data (`query_gdc_projects()`)
- **Phase 2:** List cBioPortal studies with pharmacogenomics (`client.list_studies()`)

### Step 2: Extraction
**Framework Pattern:**
```python
# Load and parse data
clinical = client.get_clinical_data(study_id, entity_type="PATIENT")
```

**Validation Application:**
- **Phase 1:** Extract case details from PubMed abstracts (`parse_case_from_abstract()`)
- **Phase 2:** Extract germline variants from TCGA VCF files query_gdc_variants()`)
- **Phase 2:** Extract treatment history from cBioPortal (`get_clinical_data()`)

### Step 3: Transformation
**Framework Pattern:**
```python
# Map to validation schema
cohort_data = transform_to_schema(raw_data)
```

**Validation Application:**
- **Phase 1:** Map PubMed cases to validation case schema (gene, variant, dose, outcome)
- **Phase 2:** Map TCGA data to validation schema (germline variant + treatment + outcome)
- **Phase 2:** Map cBioPortal data to validation schema

### Step 4: Validation
**Framework Pattern:**
```python
# Quality checks
validate_data_quality(cohort_data)
```

**Validation Application:**
- **Phase 1:** Validate case completeness (all required fields present)
- **Phase 2:** Validate variant nomenclature (PharmVar alignment)
- **Phase 2:** Validate treatment-outcome linkage (same patient, same time period)

### Step 5: Integration
**Framework Pattern:**
```python
# Feed into validation harness
validation_harness.add_cohort(cohort_data)
```

**Validation Application:**
- **Phase 3:** Feed validation cases into metrics calculator
- **Phase 3:** Run dosing guidance API on each case
- **Phase 3:** Compare predictions to actual outcomes

---

## 🛠️ Specific Code Reuse Opportunities

### Opportunity 1: Reuse PubMed Infrastructure

**Current Plan (Dosing Validation):**
```python
# scripts/validation/extract_literature_cases.py (NEW)
import requests
from ncbi_api import EUtils

def extract_dpyd_cases():
    query = '"DPYD deficiency" AND "fluoropyrimidine" AND "case report"'
    pmids = search_pubmed(query, max_results=50)
    # ... manual parsing ...
```

**Framework Has (Cohort Context):**
```python
# api/services/research_intelligence/portals/pubmed_enhanced.py (EXISTS)
from api.services.research_intelligence.portals.pubmed_enhanced import EnhancedPubMedPortal

portal = EnhancedPubMedPortal()
results = portal.search(query, max_results=50)
# Already handles: rate limiting, retries, error handling
```

**Action:**
- ✅ **REUSE** `EnhancedPubMedPortal` instead oding new PubMed client
- ✅ **ADAPT** search queries for pharmacogenomics (DPYD, UGT1A1, TPMT)
- ✅ **EXTEND** parsing logic for case report extraction (new, but uses existing infrastructure)

---

### Opportunity 2: Reuse cBioPortal Client

**Current Plan (Dosing Validation):**
```python
# scripts/validation/extract_tcga_pharmacogenes.py (NEW)
def extract_pharmacogene_cases():
    # New GDC client code
    dpyd_variants = query_gdc_variants(gene="DPYD")
    # ...
```

**Framework Has (Cohort Context):**
```python
# scripts/data_acquisition/utils/cbioportal_client.py (EXISTS)
from scripts.data_acquisition.utils.cbioportal_client import CBioportalClient

client = CBioportalClient()
studies = client.list_studies()
clinical = client.get_clinical_data(study_id, entity_type="PATIENT")
```

**Action:**
- ✅ **REUSE** `CBioportalClient` for MSK-IMPACT, Foundation cohorts (have pharmacogenomics data)
- ✅ **ADD** GDC client to framework (new, but follows same pattern)
- ✅ **EXTEND** client methods for pharmacering (new logic, existing infrastructure)

---

### Opportunity 3: Extend JSON Schema

**Current Plan (Dosing Validation):**
```json
// data/validation/dosing_guidance/sample_validation_cases.json (CUSTOM)
{
  "case_id": "LIT-001",
  "pharmacogenomics": {...},
  "treatment": {...},
  "outcome": {...}
}
```

**Framework Has (Cohort Context):**
```json
// Framework's cohort schema pattern
{
  "cohort": {
    "source": "cbioportal",
    "patients": [...],
    "metadata": {...}
  }
}
```

**Action:**
- ✅ **EXTEND** cohort schema to include pharmacogenomics fields
- ✅ **REUSE** metadata structure (extraction_date, data_quality, completeness)
- ✅ **ALIGN** validation case format with cohort schema (enables framework integration)

---

### Opportunity 4: Add Project Data Sphere as Data Source

**Current Plan (Dosing Validation):**
- ❌ **NOT MENTIONED** in dosing validation plan

**Framework Has (Cohort Context):**
```python
# scripts/data_acquisition/utils/project_data_sphere_client.py (EXISTS & CONNECTEDcripts.data_acquisition.utils.project_data_sphere_client import ProjectDataSphereClient

client = ProjectDataSphereClient(...)
if client.connect(...):
    caslibs = client.list_caslibs()
    # 102 caslibs available
    # Colorectal (16), Breast (17) - may have DPYD/UGT1A1 data
```

**Action:**
- ✅ **EXPLORE** PDS caslibs for pharmacogenomics datasets
- ✅ **PRIORITIZE** Colorectal and Breast caslibs (common 5-FU/irinotecan use)
- ✅ **EXTRACT** patient-level data with germline variants + treatment + outcomes

---

## 📋 Implementation Roadmap

### Week 1: Adapt Existing Infrastructure (Zo)

**Task 1.1: Extend PubMed Portal for Pharmacogenomics**
- **File:** `api/services/research_intelligence/portals/pubmed_enhanced.py`
- **Action:** Add `search_pharmacogenomics_cases(gene, drug, max_results)` method
- **Reuses:** Existing `portal.search()` infrastructure
- **New:** Pharmacogenomics-specific query building

**Task 1.2: Extend cBioPortal Client for Pharmacogenes**
- **File:** `scripts/data_acquisition/ioportal_client.py`
- **Action:** Add `filter_pharmacogenes(study_id, genes=['DPYD', 'UGT1A1', 'TPMT'])` method
- **Reuses:** Existing `get_clinical_data()` method
- **New:** Pharmacogene filtering logic

**Task 1.3: Extend Cohort Schema for Pharmacogenomics**
- **File:** `.cursor/rules/research/cohort_context_concept.mdc`
- **Action:** Add pharmacogenomics fields to JSON schema
- **Reuses:** Existing schema structure
- **New:** `pharmacogenomics`, `treatment`, `outcome` fields

---

### Week 2: Build Extraction Scripts (Agent Jr)

**Task 2.1: Literature Case Extraction (REUSES Framework)**
```python
# scripts/validation/extract_literature_cases.py
from api.services.research_intelligence.portals.pubmed_enhanced import EnhancedPubMedPortal

portal = EnhancedPubMedPortal()

# REUSE framework's PubMed infrastructure
dpyd_cases = portal.search_pharmacogenomics_cases(
    gene="DPYD",
    drug="fluoropyrimidine",
    max_results=50
)

# NEW: Case parsing logic (but uses framework's extraction pattern)
for case in dpyd_cases:
    parsed = parse_case_from_abstract(case)
    validation_cases.append(parsed)
```

**Task 2.2: TCGA Extraction (EXTENDS Framework)**
```python
# scripts/validation/extract_tcga_pharmacogenes.py
from scripts.data_acquisition.utils.cbioportal_client import CBioportalClient

client = CBioportalClient()

# REUSE framework's cBioPortal client
msk_studies = [s for s in client.list_studies() if 'msk_impact' in s['study_id'].lower()]

# NEW: GDC client (adds to framework)
from scripts.data_acquisition.utils.gdc_client import GDCClient  # NEW, but follows pattern
gdc = GDCClient()
tcga_variants = gdc.query_variants(gene="DPYD", project="TCGA-COAD")
```

**Task 2.3: Project Data Sphere Exploration (NEW DATA SOURCE)**
```python
# scripts/validation/explore_pds_pharmacogenes.py
from scripts.data_acquisition.utils.project_data_sphere_client import ProjectDataSphereClient

client = ProjectDataSphereClient(...)
if client.connect(...):
    # REUSE framework's PDS client
    colorectal_caslibs = [c for c in client.list_caslibs() if 'colorectal' in c.lower()]
    
    # NEW: Extract pharmacogenomics data
    for caslib in colorectal_caslibs:
        files = client.list_files_in_caslib(caslib)
        # Search for DPYD/UGT1A1 data
```

---

### Week 3: Integration & Validation (Zo)

**Task 3.1: Schema Alignment**
- Map all extracted cases to extended cohort schema
- Ensure consistency across data sources

**Task 3.2: Validation Harness Integration**
- Feed validation cases into metrics calculator
- Use framework's integration pattern

**Task 3.3: Metrics Calculation**
- Run dosing guidance API on each case
- Compare predictions to actual outcomes
- Generate validation report

---

## 🎯 Key Benefits of Framework Integration

### 1. Code Reuse (Time Savings)
- **Without Framework:** Build PubMed client, cBioPortal client, PDS client from scratch (~20 hours)
- **With Framework:** Adapt existing clients (~5 hours)
- **Savings:** 15 hours

### 2. Consistency (Quality)
- **Without Framework:** Each script has diffent error handling, rate limiting, retry logic
- **With Framework:** Standardized patterns across all data sources
- **Benefit:** More reliable, maintainable code

### 3. Extensibility (Future)
- **Without Framework:** New data sources require new infrastructure
- **With Framework:** Add new sources following established pattern
- **Benefit:** Framework grows with validation needs

### 4. Documentation (Clarity)
- **Without Framework:** Validation scripts are isolated, undocumented
- **With Framework:** All data sources documented in one place (`cohort_context_concept.mdc`)
- **Benefit:** Easier onboarding, maintenance

---

## 📝 Next Steps

### Immediate (Zo):
1. ✅ Review framework capabilities (this document)
2. ⏳ Extend `EnhancedPubMedPortal` with pharmacogenomics search method
3. ⏳ Extend `CBioportalClient` with pharmacogene filtering
4. ⏳ Update cohort schema to include pharmacogenomics fields

### Next Sprint (Agent Jr):
1. Execute literature extraction using extended PubMed portal
2. Executraction using cBioPortal client + new GDC client
3. Explore Project Data Sphere for pharmacogenomics data
4. Report case yields and data quality

### Manager Review:
- Approve framework integration approach
- Confirm data source prioritization
- Validate extended schema design

---

**Status:** ANALYSIS COMPLETE - READY FOR IMPLEMENTATION  
**Last Updated:** January 2025  
**Author:** Zo (Agent)
