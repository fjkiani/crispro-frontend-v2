# 🔬 RESEARCH INTELLIGENCE FRAMEWORK - FULL LLM-BASED RESEARCH

**Date**: January 15, 2025  
**Status**: 🚀 **ARCHITECTURE DESIGN**  
**Alpha's Vision**: Full machine learning/LLM-based research queries across multiple portals, then run mechanisms/analysis

---

## 🎯 WHAT ALPHA WANTS

**NOT**: Keyword matching, simple queries  
**YES**: Deep LLM-based research intelligence that:
1. **Formulates research questions** using LLM (like a researcher would)
2. **Queries multiple portals** (PubMed, ClinVar, cBioPortal, etc.)
3. **Synthesizes results** with LLM comprehension
4. **Runs mechanism analysis** on synthesized findings

**Example**: "How do purple potatoes help with ovarian cancer?"
- LLM formulates: "What are the active compounds in purple potatoes? What mechanisms do they target? What evidence exists for ovarian cancer?"
- Queries PubMed, cBioPortal, ClinVar
- LLM synthesizes: "Anthocyanins → VEGF inhibition → angiogenesis suppression"
- Runs MOAT mechanism analysis

---

## ✅ WHAT WE ALREADY HAVE

### **1. pyBioPortal** ✅
**Location**: `oncology-coPilot/oncology-backend/tests/pyBioPortal-master/`

**Capabilities**:
- Query cBioPortal API for cancer genomics data
- Get mutations, copy number alterations, clinical data
- Transform to pandas DataFrames

**Modules Available**:
- `mutations.py` - Mutation data
- `clinical_data.py` - Clinical attributes
- `molecular_profiles.py` - Molecular data
- `studies.py` - Study metadata
- `genes.py` - Gene information

**Status**: ✅ **READY TO USE**

---

### **2. PubMed-LLM-Agent** ✅
**Location**: `oncology-coPilot/oncology-backend-minimal/Pubmed-LLM-Agent-main/`

**Capabilities**:
- Natural language → PubMed query conversion
- LLM-based paper ranking (relevance scoring)
- Full-text PMC article retrieval
- Knowledge base caching

**Status**: ✅ **INTEGRATED** (used in `llm_literature_service.py`)

---

### **3. ClinVar Integration** ✅
**Location**: Multiple files

**Capabilities**:
- Variant pathogenicity lookup
- Clinical significance extraction
- E-utils API integration

**Files**:
- `api/routers/evidence2.py:289-476` (deep_analysis)
- `api/routers/evidence/clinvar.py`
- `src/tools/mutation_validator.py` (query_clinvar)

**Status**: ✅ **WORKING**

---

### **4. Evidence Endpoints** ✅
**Location**: `api/routers/evidence2.py`

**Endpoints**:
- `/api/evidence/literature` - PubMed search + LLM synthesis
- `/api/evidence/deep_analysis` - ClinVar + PubMed integration
- `/api/evidence/clinvar` - ClinVar lookup

**Status**: ✅ **WORKING**

---

## 🚀 WHAT WE NEED TO BUILD

### **ARCHITECTURE: Research Intelligence Orchestrator**

```
User Question: "How do purple potatoes help with ovarian cancer?"
    ↓
[1] LLM RESEARCH QUESTION FORMULATOR
    - Breaks down question into sub-questions
    - Identifies key entities (compound, disease, mechanisms)
    - Generates portal-specific queries
    ↓
[2] MULTI-PORTAL QUERY EXECUTOR
    - PubMed: "anthocyanins AND ovarian cancer AND angiogenesis"
    - cBioPortal: "ovarian cancer AND VEGF pathway"
    - ClinVar: "BRCA1 variants AND ovarian cancer"
    ↓
[3] LLM SYNTHESIS ENGINE
    - Reads all results from portals
    - Extracts mechanisms, evidence, confidence
    - Identifies knowledge gaps
    ↓
[4] MOAT MECHANISM ANALYSIS
    - Maps synthesized findings to pathways
    - Calculates pathway alignment
    - Generates personalized recommendations
```

---

## 📋 EXACT REQUIREMENTS

### **1. Research Question Formulator** (NEW)

**Input**: Natural language question  
**Output**: Structured research plan

**Needs**:
- LLM service (OpenAI/Gemini) ✅ **HAVE**
- Prompt engineering for question decomposition
- Entity extraction (compound, disease, mechanism, etc.)

**File to Create**: `api/services/research_intelligence/question_formulator.py`

---

### **2. Multi-Portal Query Executor** (NEW)

**Portals to Support**:
- ✅ **PubMed** - Already integrated (PubMed-LLM-Agent)
- ✅ **ClinVar** - Already integrated (evidence2.py)
- ✅ **cBioPortal** - Have pyBioPortal, need integration
- ❌ **COSMIC** - Need to add
- ❌ **UniProt** - Need to add (protein data)
- ❌ **ChEMBL** - Need to add (compound data)

**File to Create**: `api/services/research_intelligence/portal_executor.py`

---

### **3. LLM Synthesis Engine** (ENHANCE EXISTING)

**Current**: Basic LLM extraction in `enhanced_evidence_service.py`  
**Enhancement**: Deep comprehension across multiple portals

**Needs**:
- Multi-source synthesis (combine PubMed + cBioPortal + ClinVar)
- Mechanism extraction with confidence scoring
- Evidence quality assessment
- Knowledge gap identification

**File to Enhance**: `api/services/research_intelligence/synthesis_engine.py`

---

### **4. MOAT Integration** (EXISTS, NEEDS WIRING)

**Current**: MOAT framework exists (treatment line, biomarkers, pathways)  
**Enhancement**: Connect research findings → MOAT analysis

**File**: `api/services/research_intelligence/moat_integrator.py`

---

## 🔧 IMPLEMENTATION PLAN

### **Phase 1: Research Question Formulator** (1-2 days)

**Create**: `api/services/research_intelligence/question_formulator.py`

```python
class ResearchQuestionFormulator:
    """
    Uses LLM to decompose natural language questions into structured research plan.
    """
    
    async def formulate_research_plan(
        self, 
        question: str,
        context: Dict[str, Any]  # disease, treatment_line, biomarkers
    ) -> ResearchPlan:
        """
        Input: "How do purple potatoes help with ovarian cancer?"
        
        Output:
        {
            "primary_question": "How do purple potatoes help with ovarian cancer?",
            "entities": {
                "compound": "purple potatoes",
                "active_compounds": ["anthocyanins", "cyanidin"],
                "disease": "ovarian cancer",
                "mechanisms_of_interest": ["angiogenesis", "inflammation", "DNA repair"]
            },
            "sub_questions": [
                "What active compounds are in purple potatoes?",
                "What mechanisms do anthocyanins target?",
                "What evidence exists for anthocyanins in ovarian cancer?",
                "What pathways are disrupted in ovarian cancer?"
            ],
            "portal_queries": {
                "pubmed": [
                    "anthocyanins AND ovarian cancer AND angiogenesis",
                    "purple potatoes AND cancer AND mechanisms"
                ],
                "cbioportal": [
                    "ovarian cancer AND VEGF pathway",
                    "ovarian cancer AND angiogenesis genes"
                ],
                "clinvar": [
                    "BRCA1 AND ovarian cancer",
                    "TP53 AND ovarian cancer"
                ]
            }
        }
        """
```

**Dependencies**:
- ✅ LLM service (OpenAI/Gemini) - HAVE
- ✅ Entity extraction - Can use LLM

---

### **Phase 2: Portal Executor** (2-3 days)

**Create**: `api/services/research_intelligence/portal_executor.py`

```python
class MultiPortalExecutor:
    """
    Executes queries across multiple research portals.
    """
    
    def __init__(self):
        # Initialize portal clients
        self.pubmed_agent = PubMedLLMAgent()  # ✅ HAVE
        self.clinvar_client = ClinVarClient()  # ✅ HAVE
        self.cbioportal = pybioportal.Client()  # ✅ HAVE (pyBioPortal)
        self.cosmic_client = None  # ❌ NEED TO ADD
        self.chembl_client = None  # ❌ NEED TO ADD
    
    async def execute_research_plan(
        self, 
        plan: ResearchPlan
    ) -> PortalResults:
        """
        Execute queries across all portals in parallel.
        
        Returns:
        {
            "pubmed": [...papers...],
            "cbioportal": {...genomics_data...},
            "clinvar": [...variants...],
            "cosmic": [...mutations...],
            "chembl": [...compounds...]
        }
        """
```

**Dependencies**:
- ✅ PubMed-LLM-Agent - HAVE
- ✅ ClinVar client - HAVE
- ✅ pyBioPortal - HAVE (need to integrate)
- ❌ COSMIC API client - NEED
- ❌ ChEMBL API client - NEED

---

### **Phase 3: LLM Synthesis Engine** (2-3 days)

**Create**: `api/services/research_intelligence/synthesis_engine.py`

```python
class ResearchSynthesisEngine:
    """
    Uses LLM to synthesize findings from multiple portals.
    """
    
    async def synthesize_findings(
        self,
        portal_results: PortalResults,
        research_plan: ResearchPlan
    ) -> SynthesizedFindings:
        """
        LLM reads all portal results and synthesizes:
        - Mechanisms of action
        - Evidence strength
        - Confidence scores
        - Knowledge gaps
        
        Returns:
        {
            "mechanisms": [
                {
                    "mechanism": "angiogenesis_inhibition",
                    "target": "VEGF",
                    "evidence": "Strong (3 RCTs, 15 in vitro studies)",
                    "confidence": 0.85,
                    "sources": ["pubmed", "cbioportal"]
                }
            ],
            "evidence_summary": "...",
            "knowledge_gaps": [...],
            "confidence": 0.78
        }
        """
```

**Dependencies**:
- ✅ LLM service - HAVE
- ✅ Portal results - Phase 2

---

### **Phase 4: MOAT Integration** (1-2 days)

**Create**: `api/services/research_intelligence/moat_integrator.py`

```python
class MOATIntegrator:
    """
    Connects research findings to MOAT framework.
    """
    
    async def integrate_with_moat(
        self,
        synthesized_findings: SynthesizedFindings,
        patient_context: Dict[str, Any]  # disease, treatment_line, biomarkers
    ) -> MOATAnalysis:
        """
        Maps research findings to MOAT capabilities:
        - Pathway alignment
        - Treatment line appropriateness
        - Biomarker matching
        - Toxicity mitigation
        
        Returns MOAT-formatted analysis ready for frontend.
        """
```

**Dependencies**:
- ✅ MOAT framework - HAVE
- ✅ Synthesized findings - Phase 3

---

## 📦 EXTERNAL DEPENDENCIES NEEDED

### **1. COSMIC API Client** ❌

**What**: Cancer mutation database  
**API**: https://cancer.sanger.ac.uk/cosmic/api  
**Status**: Need to create client

**File**: `api/services/research_intelligence/portals/cosmic_client.py`

---

### **2. ChEMBL API Client** ❌

**What**: Compound/target database  
**API**: https://www.ebi.ac.uk/chembl/api  
**Status**: Need to create client (or use existing if we have it)

**File**: `api/services/research_intelligence/portals/chembl_client.py`

---

### **3. UniProt API Client** ❌ (Optional)

**What**: Protein function/domain data  
**API**: https://www.uniprot.org/help/api  
**Status**: Nice to have

---

## 🎯 INTEGRATION POINTS

### **1. pyBioPortal Integration**

**Current**: Exists in `oncology-coPilot/oncology-backend/tests/pyBioPortal-master/`

**Action**: 
- Move to `api/services/research_intelligence/portals/`
- Create wrapper: `cbioportal_client.py`
- Integrate with portal executor

**Example Usage**:
```python
from pybioportal import Client

client = Client()
studies = client.studies.get_studies()
mutations = client.mutations.get_mutations(
    study_id="ov_tcga",
    molecular_profile_id="ov_tcga_mutations",
    entrez_gene_id=672  # BRCA1
)
```

---

### **2. PubMed-LLM-Agent Integration**

**Current**: Used in `llm_literature_service.py`

**Action**: 
- Enhance to support multi-query execution
- Add mechanism extraction prompts
- Integrate with synthesis engine

---

### **3. Frontend Integration**

**Current**: `MyelomaDigitalTwin.jsx` has "Synthesize Literature" toggle

**Action**:
- Add "Deep Research" button
- Call new `/api/research/intelligence` endpoint
- Display synthesized findings + MOAT analysis

---

## 📝 EXACT FILES TO CREATE

### **New Files**:

1. `api/services/research_intelligence/__init__.py`
2. `api/services/research_intelligence/question_formulator.py`
3. `api/services/research_intelligence/portal_executor.py`
4. `api/services/research_intelligence/synthesis_engine.py`
5. `api/services/research_intelligence/moat_integrator.py`
6. `api/services/research_intelligence/portals/cbioportal_client.py`
7. `api/services/research_intelligence/portals/cosmic_client.py`
8. `api/services/research_intelligence/portals/chembl_client.py`
9. `api/routers/research_intelligence.py` (API endpoint)

### **Files to Enhance**:

1. `api/services/llm_literature_service.py` - Add multi-query support
2. `api/services/enhanced_evidence_service.py` - Integrate with synthesis engine
3. `oncology-coPilot/oncology-frontend/src/pages/MyelomaDigitalTwin.jsx` - Add UI

---

## ✅ WHAT I NEED FROM YOU, ALPHA

### **1. Confirm Portal Priorities**

**Must Have**:
- ✅ PubMed (have)
- ✅ ClinVar (have)
- ✅ cBioPortal (have pyBioPortal)

**Should Have**:
- ❓ COSMIC (need API key?)
- ❓ ChEMBL (public API, no key needed)

**Nice to Have**:
- UniProt
- DrugBank
- PharmGKB (we might already have this)

---

### **2. LLM Provider Preference**

**Current**: OpenAI (you provided key)  
**Options**:
- OpenAI GPT-4o (current)
- Gemini (quota-exhausted)
- Anthropic Claude (if you have key)

**Recommendation**: Use OpenAI for now, add fallback to others.

---

### **3. pyBioPortal Location**

**Current**: `oncology-coPilot/oncology-backend/tests/pyBioPortal-master/`

**Question**: Should I:
- Move it to `api/services/research_intelligence/portals/`?
- Or keep it where it is and import from there?

---

### **4. API Endpoint Design**

**Proposed**: `POST /api/research/intelligence`

**Request**:
```json
{
  "question": "How do purple potatoes help with ovarian cancer?",
  "context": {
    "disease": "ovarian_cancer_hgs",
    "treatment_line": "L2",
    "biomarkers": {"HRD": "POSITIVE"}
  },
  "portals": ["pubmed", "cbioportal", "clinvar"],
  "synthesize": true,
  "run_moat_analysis": true
}
```

**Response**:
```json
{
  "research_plan": {...},
  "portal_results": {...},
  "synthesized_findings": {...},
  "moat_analysis": {...}
}
```

**Does this work for you, Alpha?**

---

## 🚀 NEXT STEPS

1. **You confirm**: Portal priorities, LLM provider, pyBioPortal location
2. **I build**: Research Intelligence Framework (Phases 1-4)
3. **We test**: "How do purple potatoes help with ovarian cancer?"
4. **We integrate**: Frontend UI in MyelomaDigitalTwin

**Ready to build when you confirm, Alpha!** 🔥





