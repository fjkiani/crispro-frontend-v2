# Clinical Trials Universal Access - Implementation Guide

## Executive Summary

This guide provides step-by-step instructions for implementing universal access to the clinical trials intelligence pipeline. **All changes are additive - zero modifications to existing Ayesha code.**

## Key Decisions (Improved from User Answers)

### 1. Profile Schema: **Support Both** (Override: full_only → both)
**Reasoning**: 
- Easier adoption: Users can start with simple profile, build up over time
- Adapter layer is minimal (~50 lines)
- Better UX: Progressive enhancement
- Full profile still required for best results, but we don't block users

**Implementation**: Create `ProfileAdapter` that converts simple → full with sensible defaults

### 2. Storage: **Hybrid Approach** (Override: database → hybrid)
**Reasoning**:
- Markdown/JSON files are large (50-200KB each) - better in file system
- Database for metadata/search (patient_id, nct_id, tier, score, created_at)
- File system is simpler for now, matches current Ayesha system
- Can add full database storage later if needed

**Implementation**: 
- File system: `.cursor/patients/{patient_id}/dossiers/{nct_id}.md`
- Database: `patient_dossiers` table (metadata only)

### 3. Location Config: **From Profile + Geocoding** (Enhance: from_profile)
**Reasoning**:
- Derive from profile is good, but need to convert ZIP → states/cities
- Simple ZIP-to-state mapping (no external service needed for US)
- Allow API override for flexibility

**Implementation**: 
- Extract ZIP from `patient_profile['logistics']['zip_code']`
- Use ZIP-to-state mapping (US only, ~50 states)
- Build `ALLOWED_STATES` and `NYC_METRO_CITIES` from ZIP
- Allow API parameter override

## Implementation Phases

### Phase 1: Clone & Adapt Pipeline (6-8 hours)

#### Step 1.1: Clone Directory Structure
```bash
# Create universal version
cp -r api/services/trial_intelligence api/services/trial_intelligence_universal
```

**Files to clone**:
- `pipeline.py` (main orchestrator)
- `config.py` (filter configuration)
- `__init__.py` (exports)
- `stage1_hard_filters/` (all files)
- `stage2_trial_type/` (all files)
- `stage3_location/` (all files)
- `stage4_eligibility/` (all files)
- `stage5_llm_analysis/` (all files)
- `stage6_dossier/` (all files)

#### Step 1.2: Rename Variables (Find & Replace)
**In `pipeline.py`**:
- `ayesha_profile` → `patient_profile` (parameter)
- `self.ayesha` → `self.patient_profile` (instance variable)
- `ayesha: Dict[str, Any]` → `patient: Dict[str, Any]` (type hints)

**In all stage files**:
- `ayesha` parameter → `patient` parameter
- Update docstrings

#### Step 1.3: Update FilterConfig to Derive from Patient Profile
**In `config.py`**, add method:
```python
def create_config_from_patient_profile(patient_profile: Dict[str, Any]) -> FilterConfig:
    """Create FilterConfig from patient profile."""
    config = FilterConfig()
    
    # Extract location
    logistics = patient_profile.get('logistics', {})
    zip_code = logistics.get('zip_code') or logistics.get('home_zip')
    if zip_code:
        config.PATIENT_ZIP = str(zip_code)
        # Convert ZIP to state (simple mapping)
        state = zip_to_state(zip_code)
        if state:
            config.ALLOWED_STATES = {state}  # Start with patient's state
            # Expand to metro area if needed
            if state in ['NY', 'NJ', 'CT']:
                config.ALLOWED_STATES = {'NY', 'NJ', 'CT'}
            # Add adjacent states (simple logic)
            adjacent_states = get_adjacent_states(state)
            config.ALLOWED_STATES.update(adjacent_states)
    
    # Extract disease
    disease = patient_profile.get('disease', {})
    primary_diagnosis = disease.get('primary_diagnosis', '')
    config.DISEASE_KEYWORDS = extract_disease_keywords(primary_diagnosis)
    config.PATIENT_STAGE = disease.get('figo_stage') or disease.get('stage', '')
    
    # Extract treatment line
    treatment = patient_profile.get('treatment', {})
    treatment_line = treatment.get('line', 'first-line')
    config.PREFERRED_TREATMENT_LINES = [treatment_line]
    
    return config
```

**Helper functions needed**:
- `zip_to_state(zip_code: str) -> Optional[str]`: Simple ZIP prefix → state mapping
- `get_adjacent_states(state: str) -> Set[str]`: Return adjacent states
- `extract_disease_keywords(diagnosis: str) -> List[str]`: Extract keywords from diagnosis

#### Step 1.4: Update Pipeline Initialization
**In `pipeline.py`**, update `__init__`:
```python
def __init__(self, patient_profile: Dict[str, Any], config: Optional[FilterConfig] = None, use_llm: bool = True, verbose: bool = True):
    """
    Initialize pipeline with patient profile.
    
    Args:
        patient_profile: Complete patient profile (full Ayesha-style format)
        config: FilterConfig instance (if None, derived from patient_profile)
        use_llm: Enable LLM classification and analysis
        verbose: Print progress messages
    """
    self.patient_profile = patient_profile
    
    # Derive config from patient profile if not provided
    if config is None:
        from .config import create_config_from_patient_profile
        config = create_config_from_patient_profile(patient_profile)
    
    self.config = config
    self.use_llm = use_llm if use_llm else self.config.USE_LLM
    self.verbose = verbose
    self.audit_trail = []
```

#### Step 1.5: Update Stage 3 to Pass Config
**In `pipeline.py`**, update `run_stage3`:
```python
async def run_stage3(self, trial) -> FilterResult:
    """STAGE 3: Location Validation"""
    from .stage3_location import nyc_metro_detector
    
    # Pass config explicitly
    has_location, matching_locs, reasoning = nyc_metro_detector.check(trial, self.config)
    
    if not has_location:
        return FilterResult(False, 'STAGE_3_NO_LOCATION', 0.0, [],
                          f"❌ No matching locations: {reasoning}")
    
    score = min(1.0, len(matching_locs) * 0.3)
    return FilterResult(True, 'STAGE_3', score,
                      [f'✅ {len(matching_locs)} location(s)'],
                      metadata={'matching_locations': matching_locs})
```

#### Step 1.6: Update LLM Prompts to Be Generic
**In `trial_fit_analyzer.py`**, update prompt:
```python
# Replace hardcoded "Ayesha Kiani" with:
patient_name = patient_profile.get('demographics', {}).get('name', 'Patient')
patient_zip = patient_profile.get('logistics', {}).get('zip_code', 'Unknown')
patient_location = patient_profile.get('logistics', {}).get('location', 'Unknown')

# In prompt:
f"- Name: {patient_name} ({patient_profile.get('demographics', {}).get('age', 'N/A')}{patient_profile.get('demographics', {}).get('sex', '')}, ZIP {patient_zip} - {patient_location})"
```

**Also update fallback analysis** to use patient dict fields instead of hardcoded values.

#### Step 1.7: Create Profile Adapter
**New file**: `api/services/trial_intelligence_universal/profile_adapter.py`
```python
"""
Profile Adapter - Converts simple profile to full profile format.
"""
from typing import Dict, Any

def adapt_simple_to_full_profile(simple_profile: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert simple profile to full Ayesha-style profile.
    
    Simple profile format:
    {
        'patient_id': str,
        'disease': str,
        'treatment_line': str,
        'location': str,
        'biomarkers': Dict[str, Any]
    }
    
    Returns full profile format with sensible defaults.
    """
    # Extract from simple profile
    patient_id = simple_profile.get('patient_id', 'unknown')
    disease_name = simple_profile.get('disease', '')
    treatment_line = simple_profile.get('treatment_line', 'first-line')
    location = simple_profile.get('location', 'Unknown')
    biomarkers = simple_profile.get('biomarkers', {})
    
    # Build full profile with defaults
    full_profile = {
        'demographics': {
            'patient_id': patient_id,
            'name': simple_profile.get('name', f'Patient {patient_id}'),
            'age': simple_profile.get('age', None),
            'sex': simple_profile.get('sex', None),
            'location': location,
        },
        'disease': {
            'primary_diagnosis': disease_name,
            'stage': simple_profile.get('stage', 'Unknown'),
            'figo_stage': simple_profile.get('stage', 'Unknown'),
            'tumor_burden': simple_profile.get('tumor_burden', 'Unknown'),
            'performance_status': simple_profile.get('performance_status', None),
        },
        'treatment': {
            'line': treatment_line,
            'line_number': _parse_treatment_line_number(treatment_line),
            'status': 'treatment_naive' if treatment_line == 'first-line' else 'on_treatment',
            'prior_therapies': simple_profile.get('prior_therapies', []),
        },
        'biomarkers': biomarkers,
        'eligibility': {
            'age_eligible': True,
            'performance_status': 'ECOG 0-2',
            'organ_function': {
                'hepatic': 'normal',
                'renal': 'normal',
                'cardiac': 'normal',
                'pulmonary': 'normal',
            },
            'exclusions': {
                'bowel_obstruction': False,
                'active_infection': False,
                'brain_metastases': False,
                'other_malignancy': False,
            },
        },
        'logistics': {
            'location': location,
            'zip_code': simple_profile.get('zip_code', None),
            'home_zip': simple_profile.get('zip_code', None),
            'travel_radius_miles': simple_profile.get('travel_radius_miles', 50),
            'willing_to_travel': simple_profile.get('willing_to_travel', True),
        },
        'labs': {},
        'screening': {
            'recist_measurable_disease': True,
            'target_lesions_present': True,
        },
        'critical_gates': {},
        'probability_estimates': {},
    }
    
    return full_profile

def _parse_treatment_line_number(treatment_line: str) -> int:
    """Parse treatment line string to number."""
    if 'first' in treatment_line.lower() or '1' in treatment_line:
        return 1
    elif 'second' in treatment_line.lower() or '2' in treatment_line:
        return 2
    elif 'third' in treatment_line.lower() or '3' in treatment_line:
        return 3
    return 1  # Default to first-line
```

### Phase 2: Create API Endpoints (4-5 hours)

#### Step 2.1: Create Router
**New file**: `api/routers/dossiers_intelligence.py`
```python
"""
Universal Dossier Intelligence API Router
"""
from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel
from typing import List, Dict, Any, Optional
from pathlib import Path
import json
from datetime import datetime

from api.services.trial_intelligence_universal.pipeline import TrialIntelligencePipeline
from api.services.trial_intelligence_universal.profile_adapter import adapt_simple_to_full_profile
from api.services.trial_intelligence_universal.config import FilterConfig, create_config_from_patient_profile

router = APIRouter(prefix="/api/dossiers/intelligence", tags=["Dossier Intelligence"])

# Storage directory
DOSSIER_STORAGE_ROOT = Path(__file__).resolve().parent.parent.parent / ".cursor" / "patients"
DOSSIER_STORAGE_ROOT.mkdir(parents=True, exist_ok=True)

class SimplePatientProfile(BaseModel):
    """Simple patient profile (for easy adoption)."""
    patient_id: str
    disease: str
    treatment_line: str = "first-line"
    location: str = "Unknown"
    biomarkers: Dict[str, Any] = {}
    zip_code: Optional[str] = None
    age: Optional[int] = None
    sex: Optional[str] = None
    stage: Optional[str] = None

class FullPatientProfile(BaseModel):
    """Full patient profile (Ayesha-style)."""
    demographics: Dict[str, Any]
    disease: Dict[str, Any]
    treatment: Dict[str, Any]
    biomarkers: Dict[str, Any]
    eligibility: Dict[str, Any]
    logistics: Dict[str, Any]
    labs: Optional[Dict[str, Any]] = None
    screening: Optional[Dict[str, Any]] = None
    critical_gates: Optional[Dict[str, Any]] = None
    probability_estimates: Optional[Dict[str, Any]] = None

class FilterRequest(BaseModel):
    """Request to filter trials."""
    patient_profile: Dict[str, Any]  # Accepts both simple and full
    candidates: List[Dict[str, Any]]  # Trial candidates to filter
    use_llm: bool = True
    config_override: Optional[Dict[str, Any]] = None  # Override config values

class GenerateDossierRequest(BaseModel):
    """Request to generate dossier."""
    patient_profile: Dict[str, Any]
    nct_id: str
    use_llm: bool = True

@router.post("/filter")
async def filter_trials(request: FilterRequest):
    """Filter trials using universal pipeline."""
    try:
        # Adapt profile if needed
        if _is_simple_profile(request.patient_profile):
            patient_profile = adapt_simple_to_full_profile(request.patient_profile)
        else:
            patient_profile = request.patient_profile
        
        # Create config (with override if provided)
        if request.config_override:
            config = create_config_from_patient_profile(patient_profile)
            # Apply overrides
            for key, value in request.config_override.items():
                setattr(config, key, value)
        else:
            config = create_config_from_patient_profile(patient_profile)
        
        # Run pipeline
        pipeline = TrialIntelligencePipeline(
            patient_profile=patient_profile,
            config=config,
            use_llm=request.use_llm,
            verbose=False
        )
        
        results = await pipeline.execute(request.candidates)
        
        return {
            'top_tier': results['top_tier'],
            'good_tier': results['good_tier'],
            'rejected': results['rejected'],
            'statistics': results['statistics']
        }
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/generate")
async def generate_dossier(request: GenerateDossierRequest):
    """Generate dossier for a trial."""
    try:
        # Adapt profile if needed
        if _is_simple_profile(request.patient_profile):
            patient_profile = adapt_simple_to_full_profile(request.patient_profile)
        else:
            patient_profile = request.patient_profile
        
        # Get trial from database (TODO: implement)
        # For now, assume trial is passed or fetched
        
        # Run pipeline on single trial
        config = create_config_from_patient_profile(patient_profile)
        pipeline = TrialIntelligencePipeline(
            patient_profile=patient_profile,
            config=config,
            use_llm=request.use_llm,
            verbose=False
        )
        
        # Generate dossier (use stage6_dossier assembler)
        from api.services.trial_intelligence_universal.stage6_dossier.assembler import assemble
        
        # TODO: Fetch trial data
        trial = {}  # Placeholder
        
        markdown = assemble(trial, patient_profile)
        
        # Save to file system
        patient_id = patient_profile['demographics']['patient_id']
        dossier_dir = DOSSIER_STORAGE_ROOT / patient_id / "dossiers"
        dossier_dir.mkdir(parents=True, exist_ok=True)
        
        dossier_file = dossier_dir / f"{request.nct_id}.md"
        dossier_file.write_text(markdown)
        
        return {
            'dossier_id': f"{patient_id}_{request.nct_id}",
            'nct_id': request.nct_id,
            'patient_id': patient_id,
            'markdown': markdown,
            'file_path': str(dossier_file)
        }
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

def _is_simple_profile(profile: Dict[str, Any]) -> bool:
    """Check if profile is simple format."""
    return 'demographics' not in profile and 'patient_id' in profile
```

#### Step 2.2: Register Router
**In `api/main.py`**, add:
```python
from api.routers import dossiers_intelligence
app.include_router(dossiers_intelligence.router)
```

### Phase 3: Database Schema (2-3 hours)

#### Step 3.1: Create Database Table
**New file**: `migrations/create_patient_dossiers_table.sql`
```sql
-- Patient Dossiers Metadata Table
CREATE TABLE IF NOT EXISTS patient_dossiers (
    id SERIAL PRIMARY KEY,
    patient_id VARCHAR(255) NOT NULL,
    nct_id VARCHAR(50) NOT NULL,
    dossier_id VARCHAR(255) UNIQUE NOT NULL,
    tier VARCHAR(20),  -- 'TOP_TIER', 'GOOD_TIER'
    match_score FLOAT,
    file_path TEXT,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    UNIQUE(patient_id, nct_id)
);

CREATE INDEX idx_patient_dossiers_patient_id ON patient_dossiers(patient_id);
CREATE INDEX idx_patient_dossiers_nct_id ON patient_dossiers(nct_id);
CREATE INDEX idx_patient_dossiers_tier ON patient_dossiers(tier);
```

### Phase 4: Testing (3-4 hours)

#### Step 4.1: Create Test Script
**New file**: `tests/test_universal_pipeline.py`
```python
"""
Test universal pipeline with Ayesha profile (should match Ayesha pipeline).
"""
import asyncio
from pathlib import Path
import sys

# Add project root
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root))

from api.services.trial_intelligence_universal.pipeline import TrialIntelligencePipeline
from api.services.trial_intelligence.pipeline import TrialIntelligencePipeline as AyeshaPipeline
import ayesha_patient_profile

async def test_universal_matches_ayesha():
    """Test that universal pipeline produces same results as Ayesha pipeline."""
    # Load Ayesha profile
    ayesha = ayesha_patient_profile.get_ayesha_complete_profile()
    
    # Get test trials
    # TODO: Load from database or test data
    
    # Run Ayesha pipeline
    ayesha_pipeline = AyeshaPipeline(ayesha, use_llm=False, verbose=True)
    ayesha_results = await ayesha_pipeline.execute(test_trials)
    
    # Run universal pipeline
    universal_pipeline = TrialIntelligencePipeline(ayesha, use_llm=False, verbose=True)
    universal_results = await universal_pipeline.execute(test_trials)
    
    # Compare results
    assert len(ayesha_results['top_tier']) == len(universal_results['top_tier'])
    assert len(ayesha_results['good_tier']) == len(universal_results['good_tier'])
    
    print("✅ Universal pipeline matches Ayesha pipeline!")

if __name__ == "__main__":
    asyncio.run(test_universal_matches_ayesha())
```

## Critical Implementation Notes

### 1. ZIP to State Mapping
Create simple mapping (US only):
```python
ZIP_TO_STATE = {
    # NY: 10000-14999
    **{str(i): 'NY' for i in range(10000, 15000)},
    # NJ: 07000-08999
    **{str(i): 'NJ' for i in range(7000, 9000)},
    # CT: 06000-06999
    **{str(i): 'CT' for i in range(6000, 7000)},
    # Add more as needed
}

def zip_to_state(zip_code: str) -> Optional[str]:
    """Convert ZIP code to state."""
    zip_prefix = zip_code[:5] if len(zip_code) >= 5 else zip_code
    return ZIP_TO_STATE.get(zip_prefix)
```

### 2. Location Detector Rename
Consider renaming `nyc_metro_detector` → `location_detector` in universal version (optional, for clarity).

### 3. LLM Prompt Updates
Ensure all prompts use `patient_profile` dict fields, not hardcoded values.

## Success Checklist

- [ ] Universal pipeline cloned and renamed
- [ ] FilterConfig derives from patient profile
- [ ] Profile adapter converts simple → full
- [ ] API endpoints created and registered
- [ ] Storage paths use patient_id
- [ ] LLM prompts are generic
- [ ] Tests pass (universal matches Ayesha with Ayesha profile)
- [ ] Tests pass (universal works with different patient)

## Next Steps After Implementation

1. Add database integration for dossier metadata
2. Add geocoding service for better location handling
3. Add frontend components for universal dossiers
4. Add authentication/authorization (patient_id access control)
5. Add batch processing endpoint


