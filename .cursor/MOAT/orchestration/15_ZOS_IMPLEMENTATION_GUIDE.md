# Zo's Implementation Guide: Core Orchestrator Integration

**Date:** January 28, 2025  
**Purpose:** Step-by-step implementation guide for Zo's 9 core deliverables  
**Status:** ✅ **READY TO START**

---

## 🎯 ZO'S MISSION

**Focus:** Core orchestrator agent integration and infrastructure  
**Time:** ~2-3 weeks  
**Outcome:** Fully functional orchestrator with all blocking agents integrated

---

## 📋 DELIVERABLE 1: DATA EXTRACTION AGENT

### **Objective:**
Parse uploaded files (VCF, PDF, MAF) and extract structured patient data into PatientProfile.

### **Current State:**
- ⏳ Skeleton exists: `_run_data_extraction_agent()` in orchestrator.py (placeholder)
- ⏳ PatientProfile dataclass exists in `state.py`
- ❌ No actual parsing logic implemented

### **Implementation Steps:**

#### **Step 1: Create Data Extraction Agent Module** (1-2 hours)
**File:** `api/services/orchestrator/agents/data_extraction_agent.py` (NEW)

```python
from typing import Dict, List, Optional
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)

@dataclass
class ExtractedMutation:
    """Extracted mutation data"""
    gene: str
    chromosome: str
    position: int
    ref_allele: str
    alt_allele: str
    variant_type: str  # "SNV", "INDEL", "CNV", etc.
    vaf: Optional[float] = None
    depth: Optional[int] = None
    zygosity: Optional[str] = None  # "HET", "HOM", "UNKNOWN"

@dataclass
class ExtractedClinicalData:
    """Extracted clinical data"""
    cancer_type: Optional[str] = None
    stage: Optional[str] = None
    histology: Optional[str] = None
    age: Optional[int] = None
    sex: Optional[str] = None
    ecog: Optional[int] = None
    biomarkers: Dict[str, any] = None  # TMB, MSI, HRD, etc.

class DataExtractionAgent:
    """Extracts structured data from uploaded files"""
    
    def __init__(self):
        self.vcf_parser = None  # Will use PyVCF or cyvcf2
        self.pdf_parser = None  # Will use LLM-based extraction
        self.maf_parser = None  # Tab-delimited parser
    
    async def extract_from_vcf(self, vcf_path: str) -> List[ExtractedMutation]:
        """Parse VCF file and extract mutations"""
        # TODO: Implement VCF parsing
        # Use PyVCF or cyvcf2
        # Extract: CHROM, POS, REF, ALT, INFO (VAF, DP, etc.)
        pass
    
    async def extract_from_pdf(self, pdf_path: str) -> tuple[List[ExtractedMutation], ExtractedClinicalData]:
        """Parse PDF report using LLM extraction"""
        # TODO: Implement PDF parsing
        # Use LLM service to extract:
        # - Mutations from tables/text
        # - Clinical data (stage, histology, biomarkers)
        # - Demographics (age, sex, ECOG)
        pass
    
    async def extract_from_maf(self, maf_path: str) -> List[ExtractedMutation]:
        """Parse MAF file (tab-delimited)"""
        # TODO: Implement MAF parsing
        # Standard MAF format columns
        pass
    
    async def extract(self, file_path: str, file_type: str) -> Dict:
        """Main extraction method - routes to appropriate parser"""
        if file_type.lower() == 'vcf':
            mutations = await self.extract_from_vcf(file_path)
            return {'mutations': mutations, 'clinical_data': None}
        elif file_type.lower() == 'pdf':
            mutations, clinical_data = await self.extract_from_pdf(file_path)
            return {'mutations': mutations, 'clinical_data': clinical_data}
        elif file_type.lower() == 'maf':
            mutations = await self.extract_from_maf(file_path)
            return {'mutations': mutations, 'clinical_data': None}
        else:
            raise ValueError(f"Unsupported file type: {file_type}")
```

#### **Step 2: Wire to Orchestrator** (1 hour)
**File:** `api/services/orchestrator/orchestrator.py`

**Find:** `_run_data_extraction_agent()` method (around line 200-300)

**Replace with:**
```python
async def _run_data_extraction_agent(self, state: PatientState) -> Dict:
    """Run the data extraction agent."""
    execution = state.start_agent('data_extraction')
    
    try:
        # Import agent directly (not HTTP)
        from .agents.data_extraction_agent import DataExtractionAgent
        
        # Get file info from state (uploaded file path/type)
        # TODO: How are files passed? Check state structure
        file_path = state.uploaded_file_path  # Or however files are stored
        file_type = state.uploaded_file_type  # 'vcf', 'pdf', 'maf'
        
        # Create agent and extract
        agent = DataExtractionAgent()
        extracted_data = await agent.extract(file_path, file_type)
        
        # Convert to PatientProfile format
        # TODO: Map ExtractedMutation to PatientProfile.mutations format
        # TODO: Map ExtractedClinicalData to PatientProfile fields
        
        # Update state
        state.patient_profile = self._build_patient_profile(extracted_data)
        
        # Complete execution
        execution.complete({
            'mutations_count': len(extracted_data['mutations']),
            'clinical_data': extracted_data['clinical_data']
        })
        
        return extracted_data
    except Exception as e:
        logger.error(f"Data extraction failed: {e}")
        execution.fail(str(e))
        raise
```

#### **Step 3: Implement VCF Parser** (1-2 hours)
**Dependencies:** Install `pyvcf` or `cyvcf2`

```python
# In data_extraction_agent.py
import vcf  # or cyvcf2

async def extract_from_vcf(self, vcf_path: str) -> List[ExtractedMutation]:
    """Parse VCF file and extract mutations"""
    mutations = []
    
    with open(vcf_path, 'r') as f:
        vcf_reader = vcf.Reader(f)
        for record in vcf_reader:
            # Extract VAF from INFO or FORMAT
            vaf = None
            if 'VAF' in record.INFO:
                vaf = record.INFO['VAF']
            elif 'AF' in record.INFO:
                vaf = record.INFO['AF']
            
            # Extract depth
            depth = None
            if 'DP' in record.INFO:
                depth = record.INFO['DP']
            
            # Extract zygosity
            zygosity = None
            if record.genotype:
                gt = record.genotype(record.samples[0].sample) if record.samples else None
                if gt:
                    if gt.is_het:
                        zygosity = "HET"
                    elif gt.is_hom_alt:
                        zygosity = "HOM"
                    else:
                        zygosity = "UNKNOWN"
            
            mutation = ExtractedMutation(
                gene=record.CHROM,  # Or extract from INFO['GENE'] if available
                chromosome=record.CHROM,
                position=record.POS,
                ref_allele=record.REF,
                alt_allele=str(record.ALT[0]) if record.ALT else None,
                variant_type=self._classify_variant(record),
                vaf=vaf,
                depth=depth,
                zygosity=zygosity
            )
            mutations.append(mutation)
    
    return mutations
```

#### **Step 4: Implement PDF Parser (LLM-based)** (2-3 hours)
**Dependencies:** Use existing LLM service

```python
# In data_extraction_agent.py
from ..gpt_service import get_gpt_service  # Or appropriate LLM service

async def extract_from_pdf(self, pdf_path: str) -> tuple[List[ExtractedMutation], ExtractedClinicalData]:
    """Parse PDF report using LLM extraction"""
    
    # Step 1: Extract text from PDF
    # Use existing PDF extraction service or PyPDF2/pdfplumber
    
    # Step 2: Use LLM to extract structured data
    llm_service = get_gpt_service()
    
    prompt = f"""
    Extract the following from this NGS report:
    
    1. Mutations: List all mutations with gene, chromosome, position, ref, alt, VAF, depth
    2. Clinical Data:
       - Cancer type
       - Stage
       - Histology
       - Age
       - Sex
       - ECOG score
    3. Biomarkers: TMB, MSI status, HRD status, PD-L1, etc.
    
    Report text:
    {pdf_text}
    
    Return as JSON.
    """
    
    response = await llm_service.generate(prompt)
    # Parse JSON response
    # Map to ExtractedMutation and ExtractedClinicalData
    
    return mutations, clinical_data
```

#### **Step 5: Implement MAF Parser** (1 hour)
```python
async def extract_from_maf(self, maf_path: str) -> List[ExtractedMutation]:
    """Parse MAF file (tab-delimited)"""
    import csv
    
    mutations = []
    with open(maf_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            mutation = ExtractedMutation(
                gene=row.get('Hugo_Symbol', ''),
                chromosome=row.get('Chromosome', ''),
                position=int(row.get('Start_Position', 0)),
                ref_allele=row.get('Reference_Allele', ''),
                alt_allele=row.get('Tumor_Seq_Allele2', ''),
                variant_type=row.get('Variant_Classification', ''),
                vaf=float(row.get('tumor_f', 0)) if row.get('tumor_f') else None,
                depth=int(row.get('t_depth', 0)) if row.get('t_depth') else None
            )
            mutations.append(mutation)
    
    return mutations
```

#### **Step 6: Data Quality Validation** (1 hour)
```python
def validate_extraction(self, extracted_data: Dict) -> Dict[str, any]:
    """Validate extracted data quality"""
    warnings = []
    errors = []
    
    # Check mutation count
    if len(extracted_data['mutations']) == 0:
        warnings.append("No mutations found in file")
    
    # Check coverage thresholds
    for mut in extracted_data['mutations']:
        if mut.depth and mut.depth < 100:
            warnings.append(f"Low depth for {mut.gene}:{mut.position} ({mut.depth}x)")
        if mut.vaf and mut.vaf < 0.05:
            warnings.append(f"Low VAF for {mut.gene}:{mut.position} ({mut.vaf})")
    
    # Check required clinical fields
    if extracted_data['clinical_data']:
        if not extracted_data['clinical_data'].cancer_type:
            warnings.append("Cancer type not found")
    
    return {
        'valid': len(errors) == 0,
        'warnings': warnings,
        'errors': errors
    }
```

### **Acceptance Criteria:**
- ✅ Can parse VCF files (test with sample VCF)
- ✅ Can parse PDF reports (test with sample PDF)
- ✅ Can parse MAF files (test with sample MAF)
- ✅ Outputs PatientProfile object
- ✅ Validates data quality (coverage, VAF thresholds)
- ✅ Flags ambiguities for human review
- ✅ Unit tests pass
- ✅ End-to-end test: Upload file → Extract → PatientProfile created

### **Test Cases:**
1. **VCF Test:** Upload sample VCF, verify mutations extracted
2. **PDF Test:** Upload sample NGS PDF, verify mutations + clinical data extracted
3. **MAF Test:** Upload sample MAF, verify mutations extracted
4. **Validation Test:** Upload low-quality file, verify warnings generated
5. **Integration Test:** Full pipeline: Upload → Extract → Biomarker agent runs

---

## 📋 DELIVERABLE 2: DRUG EFFICACY INTEGRATION

### **Objective:**
Wire S/P/E framework to orchestrator for drug ranking.

### **Current State:**
- ✅ S/P/E framework exists: `api/services/efficacy_orchestrator/orchestrator.py`
- ✅ Framework validated and working
- ⏳ Needs orchestrator wiring

### **Implementation Steps:**

#### **Step 1: Understand S/P/E Framework Interface** (30 min)
**File:** `api/services/efficacy_orchestrator/orchestrator.py`

**Check:**
- What's the class name? (`DrugEfficacyOrchestrator`?)
- What's the method signature? (`predict()` or `rank_drugs()`?)
- What inputs does it need? (mutations, biomarker_profile, resistance_prediction?)
- What does it return? (Drug ranking list? Efficacy scores?)

#### **Step 2: Wire to Orchestrator** (2-3 hours)
**File:** `api/services/orchestrator/orchestrator.py`

**Find:** `_run_drug_efficacy_agent()` method

**Replace with:**
```python
async def _run_drug_efficacy_agent(self, state: PatientState) -> Dict:
    """Run the drug efficacy agent (S/P/E framework)."""
    execution = state.start_agent('drug_efficacy')
    
    try:
        # Import service DIRECTLY (not HTTP)
        from ..efficacy_orchestrator.orchestrator import DrugEfficacyOrchestrator
        
        # Build request from PatientState
        request = {
            'mutations': [m.to_dict() for m in state.patient_profile.mutations],
            'biomarker_profile': state.biomarker_profile,
            'resistance_prediction': state.resistance_prediction
        }
        
        # Call service method
        orchestrator = DrugEfficacyOrchestrator()
        result = await orchestrator.predict(request)  # Or appropriate method name
        
        # Convert to dict if dataclass
        if hasattr(result, 'to_dict'):
            result_dict = result.to_dict()
        elif hasattr(result, '__dict__'):
            result_dict = result.__dict__
        else:
            result_dict = result
        
        # Update state
        state.drug_ranking = result_dict.get('drugs', [])
        state.efficacy_scores = result_dict.get('scores', {})
        
        # Complete execution
        execution.complete(result_dict)
        
        return result_dict
    except Exception as e:
        logger.error(f"Drug efficacy prediction failed: {e}")
        execution.fail(str(e))
        raise
```

#### **Step 3: Test Integration** (1-2 hours)
- Test with sample patient data
- Verify drug ranking output
- Verify performance (<2 seconds)
- Test error handling

### **Acceptance Criteria:**
- ✅ S/P/E framework integrated into orchestrator
- ✅ Direct service imports (no HTTP calls)
- ✅ Drug ranking output in PatientState
- ✅ End-to-end test passes
- ✅ Performance: <2 seconds for drug ranking

---

## 📋 DELIVERABLE 3: NUTRITION INTEGRATION

### **Objective:**
Wire nutrition services to orchestrator for toxicity-aware nutrition recommendations.

### **Current State:**
- ✅ Nutrition services exist (check `api/services/nutrition/` or `api/services/food_*`)
- ⏳ Needs orchestrator wiring

### **Implementation Steps:**

#### **Step 1: Find Nutrition Services** (30 min)
**Check:**
- `api/services/nutrition/` directory
- `api/services/food_*` services
- `api/services/toxicity_pathway_mappings.py`
- `api/services/dietician_recommendations.py`

#### **Step 2: Wire to Orchestrator** (2-3 hours)
**File:** `api/services/orchestrator/orchestrator.py`

**Find:** `_run_nutrition_agent()` method

**Replace with:**
```python
async def _run_nutrition_agent(self, state: PatientState) -> Dict:
    """Run the nutrition agent (toxicity-aware recommendations)."""
    execution = state.start_agent('nutrition')
    
    try:
        # Import service DIRECTLY (not HTTP)
        from ..nutrition import NutritionAgent  # Or appropriate service
        # OR
        from ..toxicity_pathway_mappings import get_nutrition_recommendations
        
        # Build request from PatientState
        request = {
            'mutations': [m.to_dict() for m in state.patient_profile.mutations],
            'drug_ranking': state.drug_ranking,  # For toxicity-aware recommendations
            'biomarker_profile': state.biomarker_profile
        }
        
        # Call service method
        agent = NutritionAgent()  # Or appropriate service class
        result = await agent.recommend(request)  # Or appropriate method
        
        # Update state
        state.nutrition_plan = result.get('recommendations', [])
        state.food_restrictions = result.get('restrictions', [])
        
        # Complete execution
        execution.complete(result)
        
        return result
    except Exception as e:
        logger.error(f"Nutrition recommendation failed: {e}")
        execution.fail(str(e))
        raise
```

### **Acceptance Criteria:**
- ✅ Nutrition services integrated into orchestrator
- ✅ Direct service imports (no HTTP calls)
- ✅ Nutrition plan output in PatientState
- ✅ End-to-end test passes

---

## 📋 DELIVERABLE 5: TRIGGER SYSTEM

### **Objective:**
Implement event detection and automated actions.

### **Current State:**
- ✅ TriggerEngine exists (check `api/services/triggers/`)
- ⏳ Needs orchestrator integration

### **Implementation Steps:**

#### **Step 1: Check Existing Trigger System** (30 min)
**Check:** `api/services/triggers/` directory

#### **Step 2: Integrate with Orchestrator** (2-3 hours)
**File:** `api/services/orchestrator/orchestrator.py`

**Add after pipeline completes:**
```python
async def _run_trigger_system(self, state: PatientState) -> Dict:
    """Run trigger system for event detection and automated actions."""
    execution = state.start_agent('trigger_system')
    
    try:
        from ..triggers.trigger_engine import TriggerEngine
        
        engine = TriggerEngine()
        
        # Detect events from state
        events = engine.detect_events(state)
        
        # Execute automated actions
        actions = []
        for event in events:
            action = engine.execute_action(event, state)
            actions.append(action)
        
        # Update state
        state.triggered_events = events
        state.automated_actions = actions
        
        execution.complete({'events': events, 'actions': actions})
        return {'events': events, 'actions': actions}
    except Exception as e:
        logger.error(f"Trigger system failed: {e}")
        execution.fail(str(e))
        raise
```

### **Acceptance Criteria:**
- ✅ Event detection working
- ✅ Automated actions triggered
- ✅ Escalation protocols implemented
- ✅ Integration with monitoring agent working

---

## 📋 DELIVERABLE 10: ERROR HANDLING & RECOVERY

### **Objective:**
Implement robust error handling for agent failures.

### **Implementation Steps:**

#### **Step 1: Agent Failure Recovery** (2-3 hours)
**File:** `api/services/orchestrator/orchestrator.py`

**Add error handling to each agent:**
```python
async def _run_pipeline(self, state: PatientState):
    """Run complete pipeline with error handling"""
    
    agents = [
        ('data_extraction', self._run_data_extraction_agent),
        ('biomarker', self._run_biomarker_agent),
        ('resistance', self._run_resistance_agent),
        ('drug_efficacy', self._run_drug_efficacy_agent),
        ('nutrition', self._run_nutrition_agent),
        ('trial_matching', self._run_trial_matching_agent),
        ('care_plan', self._run_care_plan_agent)
    ]
    
    failed_agents = []
    partial_results = {}
    
    for agent_name, agent_func in agents:
        try:
            result = await agent_func(state)
            partial_results[agent_name] = result
        except Exception as e:
            logger.error(f"Agent {agent_name} failed: {e}")
            failed_agents.append(agent_name)
            
            # Decide: continue or stop?
            if agent_name in ['data_extraction']:  # Blocking agents
                raise  # Stop pipeline
            else:
                # Continue with partial results
                partial_results[agent_name] = {'error': str(e)}
    
    return {
        'success': len(failed_agents) == 0,
        'failed_agents': failed_agents,
        'partial_results': partial_results
    }
```

#### **Step 2: Retry Logic** (1-2 hours)
```python
async def _run_with_retry(self, agent_func, state: PatientState, max_retries: int = 3):
    """Run agent with retry logic"""
    for attempt in range(max_retries):
        try:
            return await agent_func(state)
        except Exception as e:
            if attempt == max_retries - 1:
                raise
            await asyncio.sleep(2 ** attempt)  # Exponential backoff
```

#### **Step 3: Circuit Breaker** (1-2 hours)
```python
class CircuitBreaker:
    """Circuit breaker pattern for agent calls"""
    
    def __init__(self, failure_threshold: int = 5, timeout: int = 60):
        self.failure_count = 0
        self.failure_threshold = failure_threshold
        self.timeout = timeout
        self.last_failure_time = None
        self.state = 'CLOSED'  # CLOSED, OPEN, HALF_OPEN
    
    def call(self, func, *args, **kwargs):
        if self.state == 'OPEN':
            if time.time() - self.last_failure_time > self.timeout:
                self.state = 'HALF_OPEN'
            else:
                raise CircuitBreakerOpenError("Circuit breaker is OPEN")
        
        try:
            result = func(*args, **kwargs)
            if self.state == 'HALF_OPEN':
                self.state = 'CLOSED'
                self.failure_count = 0
            return result
        except Exception as e:
            self.failure_count += 1
            self.last_failure_time = time.time()
            if self.failure_count >= self.failure_threshold:
                self.state = 'OPEN'
            raise
```

### **Acceptance Criteria:**
- ✅ Agent failures handled gracefully
- ✅ Partial failures don't break pipeline
- ✅ Retry logic prevents transient failures
- ✅ Circuit breakers prevent cascade failures

---

## 📋 DELIVERABLE 13: STATE PERSISTENCE & RECOVERY

### **Objective:**
Persist orchestrator state to database and enable recovery.

### **Current State:**
- ✅ StateStore exists: `api/services/orchestrator/state_store.py`
- ⏳ Needs database persistence

### **Implementation Steps:**

#### **Step 1: Database Persistence** (2-3 hours)
**File:** `api/services/orchestrator/state_store.py` (enhance)

**Add:**
```python
async def persist_state(self, state: PatientState) -> str:
    """Persist state to database"""
    # Serialize state to JSON
    state_json = state.to_json()
    
    # Store in database (SQLite or PostgreSQL)
    # Use existing database connection
    # Store: patient_id, state_json, version, timestamp
    
    return state_id

async def recover_state(self, patient_id: str, state_id: Optional[str] = None) -> PatientState:
    """Recover state from database"""
    # Load from database
    # Deserialize JSON to PatientState
    # Return recovered state
    pass
```

#### **Step 2: State Versioning** (1-2 hours)
```python
def version_state(self, state: PatientState) -> str:
    """Create version identifier for state"""
    # Hash state content
    # Return version string
    pass
```

### **Acceptance Criteria:**
- ✅ State persisted to database
- ✅ Recovery works after crash
- ✅ State versioning implemented
- ✅ Cleanup/archival working

---

## 📋 DELIVERABLE 15: DATA VALIDATION & QUALITY

### **Objective:**
Validate extracted data quality and coverage.

### **Implementation Steps:**

#### **Step 1: Input Validation** (2-3 hours)
**File:** `api/services/orchestrator/agents/data_extraction_agent.py` (add validation)

```python
def validate_mutation(self, mutation: ExtractedMutation) -> Dict[str, any]:
    """Validate single mutation"""
    errors = []
    warnings = []
    
    # Required fields
    if not mutation.gene:
        errors.append("Gene name required")
    if not mutation.chromosome:
        errors.append("Chromosome required")
    if mutation.position <= 0:
        errors.append("Position must be positive")
    
    # Quality checks
    if mutation.depth and mutation.depth < 100:
        warnings.append(f"Low depth: {mutation.depth}x (recommended: ≥100x)")
    if mutation.vaf and mutation.vaf < 0.05:
        warnings.append(f"Low VAF: {mutation.vaf} (recommended: ≥5%)")
    
    return {'valid': len(errors) == 0, 'errors': errors, 'warnings': warnings}

def validate_coverage(self, mutations: List[ExtractedMutation]) -> Dict[str, any]:
    """Validate overall coverage"""
    total_mutations = len(mutations)
    high_quality = sum(1 for m in mutations if m.depth and m.depth >= 100)
    
    coverage_ratio = high_quality / total_mutations if total_mutations > 0 else 0
    
    return {
        'total_mutations': total_mutations,
        'high_quality_count': high_quality,
        'coverage_ratio': coverage_ratio,
        'meets_threshold': coverage_ratio >= 0.8  # 80% high quality
    }
```

### **Acceptance Criteria:**
- ✅ Input validation rules implemented
- ✅ Data quality checks working
- ✅ Coverage thresholds enforced
- ✅ Quality scoring working
- ✅ Validation errors reported clearly

---

## 📋 DELIVERABLE 7 & 8: AGENT INTEGRATIONS (Shared)

### **Objective:**
Wire Access & Advocacy and Toxicity Risk agents to orchestrator.

### **Implementation Steps:**

#### **Step 1: Create Agent Interfaces** (1-2 hours each)
**Files:**
- `api/services/orchestrator/agents/access_advocacy_agent.py` (interface/skeleton)
- `api/services/orchestrator/agents/toxicity_risk_agent.py` (interface/skeleton)

#### **Step 2: Wire to Orchestrator** (1-2 hours each)
**File:** `api/services/orchestrator/orchestrator.py`

**Add methods:**
```python
async def _run_access_advocacy_agent(self, state: PatientState) -> Dict:
    """Run access & advocacy agent."""
    execution = state.start_agent('access_advocacy')
    
    try:
        from .agents.access_advocacy_agent import AccessAdvocacyAgent
        
        agent = AccessAdvocacyAgent()
        result = await agent.generate_packet(state)
        
        state.access_advocacy = result
        execution.complete(result)
        return result
    except Exception as e:
        execution.fail(str(e))
        raise

async def _run_toxicity_risk_agent(self, state: PatientState) -> Dict:
    """Run toxicity risk agent."""
    execution = state.start_agent('toxicity_risk')
    
    try:
        from .agents.toxicity_risk_agent import ToxicityRiskAgent
        
        agent = ToxicityRiskAgent()
        result = await agent.assess_risk(state)
        
        state.toxicity_risk = result
        execution.complete(result)
        return result
    except Exception as e:
        execution.fail(str(e))
        raise
```

### **Acceptance Criteria:**
- ✅ Agent wired to orchestrator
- ✅ Integration points defined
- ✅ Orchestrator integration tested
- ✅ Interface documented for plumber

---

## 🚀 IMPLEMENTATION ORDER

### **Week 1: Critical Blockers**
1. **Deliverable 1: Data Extraction** (4-6h)
   - Day 1: Create agent module, implement VCF parser
   - Day 2: Implement PDF parser, MAF parser, wire to orchestrator

2. **Deliverable 2: Drug Efficacy** (8-10h)
   - Day 2-3: Wire S/P/E framework, test integration

### **Week 2: Core Intelligence**
3. **Deliverable 3: Nutrition** (4-6h)
   - Day 1: Find services, wire to orchestrator

15. **Deliverable 15: Data Validation** (1d)
   - Day 2: Add validation to data extraction agent

### **Week 3: Core Infrastructure**
10. **Deliverable 10: Error Handling** (1-2d)
   - Day 1-2: Implement error handling, retry logic, circuit breakers

13. **Deliverable 13: State Persistence** (1-2d)
   - Day 2-3: Implement database persistence, recovery

### **Week 4: Automation**
5. **Deliverable 5: Trigger System** (4-6h)
   - Day 1: Integrate trigger system

### **Week 5-6: Agent Integration**
7. **Deliverable 7: Access & Advocacy** (2-3d)
   - Week 5: Create interface, wire to orchestrator

8. **Deliverable 8: Toxicity Risk** (2-3d)
   - Week 6: Create interface, wire to orchestrator

---

## ✅ SUCCESS CRITERIA

### **Week 1 Complete:**
- ✅ Can upload VCF/PDF/MAF and extract mutations
- ✅ Drug efficacy ranking works end-to-end
- ✅ Pipeline: Upload → Extract → Rank → Display

### **Week 2 Complete:**
- ✅ Nutrition recommendations integrated
- ✅ Data validation working
- ✅ All core agents can run end-to-end

### **Week 3 Complete:**
- ✅ Error handling prevents cascade failures
- ✅ State recovery works after crash
- ✅ Orchestrator infrastructure solid

### **Week 4 Complete:**
- ✅ Trigger system detects events and triggers actions

### **Week 5-6 Complete:**
- ✅ All agent integration points defined and tested
- ✅ Orchestrator ready for plumber work

---

**See Also:**
- [13_TASK_DELEGATION.md](13_TASK_DELEGATION.md) - Complete delegation breakdown
- [14_ZOS_FOCUSED_PLAN.md](14_ZOS_FOCUSED_PLAN.md) - Focused workload
- [agent-implementation-guide.mdc](agent-implementation-guide.mdc) - Implementation patterns

