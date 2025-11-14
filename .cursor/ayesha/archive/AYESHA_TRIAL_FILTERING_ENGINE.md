# ⚔️ AYESHA TRIAL FILTERING ENGINE - AGENT JR EXECUTION PLAN ⚔️

**Mission Owner**: **AGENT JR** ⚔️  
**Overseer**: Zo (planning + strategizing while Jr executes)  
**Mission**: Find top 10 trials for Ayesha's Stage IVB ovarian cancer (germline-negative, CA-125 2842)  
**Timeline**: 8 hours total (3 backend + 4.5 frontend + 30 min testing)  
**Priority**: P0 - IMMEDIATE (she needs trials THIS WEEK)

---

## 🎯 **AGENT JR - MISSION BRIEF**

**What You're Building**:
A precision trial matching system that finds the **best clinical trials** for Ayesha based on:
- ✅ Her **known** biomarkers (germline-negative, CA-125 2842, Stage IVB)
- ⚠️ Her **unknown** biomarkers (tumor NGS pending - HRD, BRCA somatic, TMB, MSI)
- 🎯 **Transparent reasoning** for every match (why eligible, why good fit, what's conditional)

**Why This Matters**:
- Ayesha has **Stage IVB ovarian cancer** with extensive metastases
- She's **treatment-naive** (first-line eligible)
- She's **germline BRCA-negative** (changes eligible trials)
- **CA-125 is 2,842** (80x elevated = extensive disease, highly trackable)
- She needs to **start therapy ASAP** (within 2-4 weeks)

**Your Deliverable**:
- Top 10 trials (ranked by match score)
- Each with reasoning (why matched, what she needs, any red flags)
- CA-125 intelligence (what to expect during treatment)
- Frontend UI for her to review and contact trials

---

## 🚨 **CRITICAL: NO HARDCODING ALLOWED**

**Dynamic Principles** (from Zo's doctrine):
1. ✅ Iterate over ALL trials returned from backend (don't slice)
2. ✅ Use optional chaining (`trial.nct_id?.toUpperCase() || 'Unknown'`)
3. ✅ Validate fields before displaying (`{trial.phase && <Chip label={trial.phase} />}`)
4. ✅ Handle missing data gracefully (empty states, fallbacks)
5. ❌ NO hardcoded trial counts (show all from backend)
6. ❌ NO assumptions about trial structure (flexible schema)
7. ❌ NO magic numbers (all scoring logic in backend)

---

## 📊 **AYESHA'S CLINICAL PROFILE (EXTRACTED FROM ALL SCANS)**

### **DEFINITIVE DATA POINTS**:

```python
ayesha_profile = {
    # DISEASE:
    "disease": "ovarian_cancer_high_grade_serous",
    "stage": "IVB",
    "histology": "suspected_HGS",  # 70% of ovarian is HGS
    
    # BIOMARKERS (KNOWN):
    "ca125": 2842,  # Normal <35 (80x elevated)
    "germline_status": "NEGATIVE",  # ✅ Ambry 38-gene panel (6/15/2023)
    "germline_brca": "negative",     # CRITICAL: Not hereditary
    "germline_hrd_genes": "negative",  # PALB2, RAD51C/D, BRIP1 all negative
    
    # BIOMARKERS (UNKNOWN - AWAITING TUMOR NGS):
    "somatic_brca": "unknown",  # May have SOMATIC BRCA (different from germline!)
    "tumor_hrd_score": "unknown",
    "tmb": "unknown",
    "msi_status": "unknown",
    "tp53": "likely_mutated",  # 96% of HGS ovarian have TP53
    
    # METASTATIC BURDEN (PET SCAN 11/11/2025):
    "peritoneal_carcinomatosis": "extensive",  # 8cm RLQ mass, SUV 15
    "ascites": "moderate",
    "pleural_effusions": "bilateral_large",
    "lymph_nodes": "extensive",  # Cervical, thoracic, abdominopelvic
    "omental_caking": True,
    "soft_tissue_mets": ["chest_wall", "left_arm"],
    
    # TREATMENT STATUS:
    "treatment_line": 0,  # Treatment-naive
    "prior_therapies": [],
    "platinum_exposure": False,
    "parp_exposure": False,
    
    # CLINICAL CONTEXT:
    "age": 40,
    "location": "NYC",
    "performance_status": "unknown",  # Likely 1-2 (extensive disease but ambulatory)
    "urgent": True,  # Needs to start therapy ASAP
    
    # HISTORICAL CONTEXT:
    "first_symptoms": "2024-01-31",  # Left ovarian cysts seen
    "diagnosis": "2025-10-28",  # Carcinomatosis detected
    "staging": "2025-11-11",  # PET confirms Stage IV
    "time_to_treatment": "urgent",  # Within 2-4 weeks
}
```

---

## 🎯 **TRIAL FILTERING LOGIC (MULTI-STAGE)**

### **STAGE 1: HARD ELIGIBILITY FILTERS (MUST-MATCH)**

```python
# Include trials if:
✅ Disease: "ovarian cancer" OR "peritoneal cancer" OR "gynecologic malignancy"
✅ Stage: "Stage IV" OR "advanced" OR "metastatic" OR "newly diagnosed advanced"
✅ Treatment Line: "first-line" OR "frontline" OR "untreated" OR "newly diagnosed"
✅ Status: "Recruiting" OR "Active, not recruiting"
✅ Location: NYC metro (within 50 miles) OR telehealth/remote

# Exclude trials if:
❌ Recurrent-only (requires prior therapy)
❌ Platinum-resistant required (she's treatment-naive)
❌ Germline BRCA required (she's germline-negative)
❌ Specific mutation required we don't have (e.g., PI3K mutation)
❌ Closed/Completed/Suspended
```

### **STAGE 2: SOFT SCORING BOOSTS**

```python
# Boost match score for:
+0.30: First-line therapy trials (she's treatment-naive)
+0.25: Stage IV or "bulky disease" trials
+0.20: Carboplatin + Paclitaxel backbone (SOC comparator)
+0.20: "All-comers" or "BRCA-wild-type" trials (germline-negative friendly)
+0.20: IP (intraperitoneal) chemotherapy (targets carcinomatosis directly)
+0.15: Bevacizumab combinations (anti-VEGF for ascites control)
+0.15: Within 20 miles of NYC
+0.15: CA-125 response as endpoint (her marker is highly elevated)
+0.10: Large trial (>200 patients = easier enrollment)
+0.10: Phase III (more established, less risk)

# Conditional boosts (when NGS available):
+0.25: HRD-high trials (IF tumor HRD ≥42)
+0.25: Somatic BRCA trials (IF tumor has BRCA1/2 mutations)
+0.20: TMB-high checkpoint trials (IF TMB ≥10)
+0.15: TP53 mutant trials (likely applies to her)

# Penalties:
-0.30: Requires germline BRCA (she doesn't have it)
-0.25: >50 miles from NYC
-0.20: Phase I (higher risk, less proven)
-0.15: Requires specific biomarker we don't have yet
```

### **STAGE 3: CA-125 INTELLIGENCE LAYER**

```python
def analyze_ca125_for_trials(ca125=2842):
    """
    CA-125 = 2842 tells us critical information:
    """
    
    intel = {
        "disease_burden": "EXTENSIVE",  # >1000 = extensive
        "marker_utility": "HIGHLY_TRACKABLE",  # Will respond if chemo works
        "chemo_sensitivity_prediction": "LIKELY_SENSITIVE",  # HGS ovarian usually responds
        
        "trial_preferences": [
            "trials_with_ca125_endpoints",  # Primary or secondary endpoint
            "trials_measuring_ca125_kinetics",  # Early response prediction
            "trials_targeting_bulk_disease",  # Not minimal residual
            "trials_with_cytoreduction",  # Debulking + chemo
        ],
        
        "expected_response": {
            "if_platinum_sensitive": {
                "after_3_cycles": "500-800 (70-80% drop)",
                "after_6_cycles": "<200 (90% drop)",
                "complete_response": "<35"
            },
            "resistance_signal": "Rising CA-125 during treatment"
        },
        
        "monitoring_value": "PRIMARY",  # CA-125 is THE marker for her
    }
    
    return intel
```

### **STAGE 4: REASONING GENERATION (WHY THIS TRIAL?)**

```python
def generate_trial_reasoning(trial, ayesha_profile):
    """
    For each trial, explain WHY it's a match.
    
    Categories:
    1. Why eligible (hard filters)
    2. Why good fit (soft boosts)
    3. What we need to confirm (conditional)
    4. Risks/concerns
    """
    
    reasoning = {
        "match_score": 0.0,  # Will be calculated
        
        "why_eligible": [],
        # e.g., "✅ Accepts Stage IV ovarian cancer"
        #       "✅ First-line therapy (you're treatment-naive)"
        #       "✅ Germline BRCA not required (you're germline-negative)"
        
        "why_good_fit": [],
        # e.g., "🎯 Tracks CA-125 response (your marker is 2842)"
        #       "🎯 Targets peritoneal disease (you have extensive carcinomatosis)"
        #       "🎯 Mount Sinai location (12 miles from you)"
        
        "conditional_requirements": [],
        # e.g., "⚠️ May require HRD testing when tumor NGS returns"
        #       "⚠️ Performance status ≥1 required (confirm with oncologist)"
        
        "red_flags": [],
        # e.g., "⚠️ Phase I trial (higher risk, less proven)"
        #       "⚠️ Requires specific mutation we don't have yet"
        
        "evidence_tier": "STANDARD",  # or SUPPORTED, INVESTIGATIONAL
        "enrollment_likelihood": "HIGH",  # Based on inclusion criteria match
    }
    
    return reasoning
```

---

## 🔧 **BACKEND IMPLEMENTATION (ZO - 3 HOURS)**

### **File 1: Ayesha Profile Schema** (30 min)

**Location**: `oncology-coPilot/oncology-backend-minimal/api/schemas/ayesha_trials.py`

```python
from pydantic import BaseModel, Field
from typing import Optional, List, Dict
from datetime import date

class AyeshaTrialProfile(BaseModel):
    """
    Ayesha's clinical profile for trial matching.
    
    Combines:
    - Known data (germline, CA-125, stage, mets)
    - Unknown data (tumor biomarkers awaiting NGS)
    """
    
    # Disease characteristics
    disease: str = "ovarian_cancer_high_grade_serous"
    stage: str = "IVB"
    histology: str = "suspected_HGS"
    
    # Known biomarkers
    ca125: float = 2842.0
    germline_status: str = "NEGATIVE"
    germline_brca: str = "negative"
    
    # Unknown biomarkers (awaiting NGS)
    somatic_brca: Optional[str] = None
    tumor_hrd_score: Optional[float] = None
    tmb: Optional[float] = None
    msi_status: Optional[str] = None
    tp53_status: Optional[str] = "likely_mutated"
    
    # Metastatic burden
    peritoneal_carcinomatosis: str = "extensive"
    ascites: str = "moderate"
    pleural_effusions: str = "bilateral_large"
    suv_max: float = 15.0
    
    # Treatment status
    treatment_line: int = 0
    prior_therapies: List[str] = []
    platinum_exposure: bool = False
    
    # Clinical context
    age: int = 40
    location: str = "NYC"
    urgent: bool = True


class TrialMatchReasoning(BaseModel):
    """Transparent reasoning for why a trial matches."""
    
    match_score: float = Field(..., ge=0.0, le=1.0)
    
    why_eligible: List[str] = Field(default_factory=list)
    why_good_fit: List[str] = Field(default_factory=list)
    conditional_requirements: List[str] = Field(default_factory=list)
    red_flags: List[str] = Field(default_factory=list)
    
    evidence_tier: str  # STANDARD, SUPPORTED, INVESTIGATIONAL
    enrollment_likelihood: str  # HIGH, MEDIUM, LOW
    
    ca125_intelligence: Optional[Dict] = None
    germline_context: Optional[Dict] = None


class AyeshaTrialMatch(BaseModel):
    """Complete trial match with reasoning."""
    
    # Trial data
    nct_id: str
    title: str
    phase: str
    status: str
    interventions: List[str]
    locations: List[Dict]
    
    # Match data
    match_score: float
    reasoning: TrialMatchReasoning
    
    # Contact
    contact_name: Optional[str] = None
    contact_phone: Optional[str] = None
    contact_email: Optional[str] = None
```

---

### **File 2: CA-125 Intelligence Service** (30 min)

**Location**: `oncology-coPilot/oncology-backend-minimal/api/services/ca125_intelligence.py`

```python
from typing import Dict, List, Any

class CA125IntelligenceService:
    """
    Extract trial-relevant intelligence from CA-125 levels.
    
    CA-125 is THE primary tumor marker for ovarian cancer.
    Elevated CA-125 (>1000) indicates:
    - Extensive disease burden
    - Highly trackable marker (will respond if chemo works)
    - Likely chemo-sensitive (HGS ovarian usually responds to platinum)
    """
    
    @staticmethod
    def analyze(ca125_value: float) -> Dict[str, Any]:
        """
        Analyze CA-125 for trial matching intelligence.
        
        Args:
            ca125_value: Current CA-125 level
            
        Returns:
            Intelligence dict with trial preferences and predictions
        """
        
        # Classify disease burden
        if ca125_value < 100:
            burden = "MINIMAL"
        elif ca125_value < 500:
            burden = "MODERATE"
        elif ca125_value < 1000:
            burden = "SIGNIFICANT"
        else:
            burden = "EXTENSIVE"
        
        # Expected response (HGS ovarian is usually platinum-sensitive)
        expected_response = {
            "chemo_sensitivity": "LIKELY_SENSITIVE",  # 70-80% respond to platinum
            "after_cycle_3": f"{int(ca125_value * 0.2)}-{int(ca125_value * 0.3)} (70-80% drop expected)",
            "after_cycle_6": f"<{int(ca125_value * 0.1)} (90%+ drop expected)",
            "complete_response_target": "<35 (normal range)",
            "resistance_signal": "Rising CA-125 during treatment (early resistance)"
        }
        
        # Trial preferences based on CA-125
        trial_preferences = []
        boost_keywords = []
        
        if ca125_value > 1000:  # Extensive disease
            trial_preferences.extend([
                "trials_measuring_ca125_response",
                "trials_targeting_bulk_disease",
                "trials_with_cytoreductive_intent",
                "trials_with_IP_chemotherapy",  # Intraperitoneal for carcinomatosis
            ])
            boost_keywords.extend([
                "CA-125 response",
                "bulk disease",
                "cytoreduction",
                "intraperitoneal",
                "neoadjuvant"
            ])
        
        # CA-125 kinetics trials (early response prediction)
        trial_preferences.append("trials_with_ca125_kinetics")
        boost_keywords.extend(["CA-125", "tumor marker"])
        
        return {
            "ca125_value": ca125_value,
            "disease_burden": burden,
            "marker_utility": "HIGHLY_TRACKABLE",
            "expected_response": expected_response,
            "trial_preferences": trial_preferences,
            "boost_keywords": boost_keywords,
            "monitoring_strategy": "Primary tumor marker - track every cycle"
        }
```

---

### **File 3: Ayesha Trial Router** (2 hours)

**Location**: `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_trials.py`

```python
from fastapi import APIRouter, HTTPException
from typing import List, Dict, Any
import logging

from ..schemas.ayesha_trials import (
    AyeshaTrialProfile,
    AyeshaTrialMatch,
    TrialMatchReasoning
)
from ..services.ca125_intelligence import CA125IntelligenceService
from ..services.hybrid_trial_search import HybridTrialSearchService
from ..services.clinical_trial_search_service import ClinicalTrialSearchService

logger = logging.getLogger(__name__)
router = APIRouter(prefix="/api/ayesha/trials", tags=["ayesha_trials"])


@router.post("/search", response_model=Dict[str, Any])
async def search_trials_for_ayesha(
    profile: AyeshaTrialProfile = None
) -> Dict[str, Any]:
    """
    Precision trial matching for Ayesha's Stage IVB ovarian cancer.
    
    Multi-stage filtering:
    1. Hard eligibility filters (must-match criteria)
    2. Soft boosts (scoring based on fit)
    3. CA-125 intelligence (marker-based preferences)
    4. Germline context (BRCA-negative filtering)
    5. Reasoning generation (transparent why-matched)
    
    Returns top 10 trials with transparent reasoning.
    """
    
    try:
        # Use default profile if not provided
        if not profile:
            profile = AyeshaTrialProfile()
        
        logger.info(f"🎯 Searching trials for Ayesha: Stage {profile.stage}, CA-125 {profile.ca125}, Germline {profile.germline_status}")
        
        # STAGE 1: Hard filters (eligibility)
        base_query = {
            "disease": ["ovarian cancer", "peritoneal cancer", "gynecologic malignancy"],
            "stage": ["IV", "IVB", "advanced", "metastatic"],
            "treatment_line": ["first-line", "frontline", "untreated", "newly diagnosed"],
            "status": ["Recruiting", "Active, not recruiting"],
            "location_max_miles": 50,
            "from_city": "New York, NY"
        }
        
        # Call hybrid trial search
        search_service = HybridTrialSearchService()
        eligible_trials = await search_service.search_with_filters(base_query)
        
        logger.info(f"✅ Found {len(eligible_trials)} eligible trials after hard filters")
        
        # STAGE 2: Score each trial
        ca125_intel = CA125IntelligenceService.analyze(profile.ca125)
        
        scored_trials = []
        for trial in eligible_trials:
            score, reasoning = calculate_match_score_with_reasoning(
                trial, profile, ca125_intel
            )
            
            scored_trials.append({
                "trial": trial,
                "match_score": score,
                "reasoning": reasoning
            })
        
        # STAGE 3: Sort by match score (descending)
        scored_trials.sort(key=lambda t: t["match_score"], reverse=True)
        
        # STAGE 4: Return top 10
        top_10 = scored_trials[:10]
        
        return {
            "profile": profile.dict(),
            "trials": [
                {
                    "nct_id": t["trial"].get("nct_id"),
                    "title": t["trial"].get("title"),
                    "phase": t["trial"].get("phase"),
                    "status": t["trial"].get("status"),
                    "interventions": t["trial"].get("interventions", []),
                    "match_score": t["match_score"],
                    "reasoning": t["reasoning"]
                }
                for t in top_10
            ],
            "ca125_intelligence": ca125_intel,
            "total_screened": len(eligible_trials),
            "germline_context": {
                "status": profile.germline_status,
                "brca_germline": profile.germline_brca,
                "implication": "Excluded germline-BRCA-required trials; focused on all-comers and somatic biomarker trials"
            },
            "provenance": {
                "filters_applied": "stage_IV + first_line + germline_negative + NYC + recruiting",
                "boost_strategy": "ca125_tracking + bulk_disease + IP_chemo + bevacizumab",
                "awaiting_ngs": ["somatic_brca", "tumor_hrd", "tmb", "msi"]
            }
        }
        
    except Exception as e:
        logger.error(f"❌ Ayesha trial search failed: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))


def calculate_match_score_with_reasoning(
    trial: Dict,
    profile: AyeshaTrialProfile,
    ca125_intel: Dict
) -> tuple[float, TrialMatchReasoning]:
    """
    Calculate match score with transparent reasoning.
    
    Returns (score, reasoning) tuple.
    """
    
    score = 0.5  # Base score
    why_eligible = []
    why_good_fit = []
    conditional = []
    red_flags = []
    
    # BOOST: First-line trial
    if any(keyword in trial.get("description", "").lower() 
           for keyword in ["first-line", "frontline", "newly diagnosed"]):
        score += 0.30
        why_good_fit.append("✅ First-line therapy (you're treatment-naive)")
    
    # BOOST: Stage IV specific
    if "stage iv" in trial.get("description", "").lower() or "advanced" in trial.get("eligibility", "").lower():
        score += 0.25
        why_good_fit.append("✅ Stage IV specific protocol")
    
    # BOOST: Carboplatin + Paclitaxel backbone
    interventions_str = " ".join(trial.get("interventions", [])).lower()
    if "carboplatin" in interventions_str and "paclitaxel" in interventions_str:
        score += 0.20
        why_good_fit.append("✅ Standard-of-care chemotherapy backbone")
    
    # BOOST: Germline-negative friendly
    if "brca wild-type" in trial.get("eligibility", "").lower() or "all comers" in trial.get("description", "").lower():
        score += 0.20
        why_good_fit.append("✅ Accepts germline BRCA-negative (your status)")
    
    # BOOST: IP chemotherapy (targets carcinomatosis)
    if "intraperitoneal" in interventions_str or "ip chemotherapy" in interventions_str:
        score += 0.20
        why_good_fit.append("🎯 Intraperitoneal chemo (targets your peritoneal disease)")
    
    # BOOST: Bevacizumab (anti-VEGF for ascites)
    if "bevacizumab" in interventions_str or "avastin" in interventions_str:
        score += 0.15
        why_good_fit.append("🎯 Bevacizumab (targets ascites and peritoneal disease)")
    
    # BOOST: CA-125 tracking
    if any(kw in trial.get("endpoints", "").lower() for kw in ca125_intel["boost_keywords"]):
        score += 0.15
        why_good_fit.append(f"🎯 Tracks CA-125 response (your marker is {profile.ca125})")
    
    # BOOST: NYC location
    trial_locations = trial.get("locations", [])
    if any("new york" in loc.get("city", "").lower() for loc in trial_locations):
        score += 0.15
        nyc_sites = [loc.get("facility") for loc in trial_locations if "new york" in loc.get("city", "").lower()]
        why_good_fit.append(f"📍 NYC location: {nyc_sites[0] if nyc_sites else 'Multiple sites'}")
    
    # CONDITIONAL: HRD testing
    if "hrd" in trial.get("biomarkers", "").lower() or "homologous recombination" in trial.get("description", "").lower():
        conditional.append("⚠️ May require HRD testing when tumor NGS returns")
    
    # CONDITIONAL: Somatic BRCA
    if "brca" in trial.get("biomarkers", "").lower() and "somatic" in trial.get("description", "").lower():
        conditional.append("⚠️ May require somatic BRCA testing (tumor NGS pending)")
    
    # RED FLAG: Germline BRCA required
    if "germline brca required" in trial.get("eligibility", "").lower() or "brca mutation required" in trial.get("inclusion", "").lower():
        score -= 0.30
        red_flags.append("❌ Requires germline BRCA (you're germline-negative)")
    
    # RED FLAG: Phase I
    if trial.get("phase") == "Phase 1":
        score -= 0.20
        red_flags.append("⚠️ Phase I trial (higher risk, less established)")
    
    # Determine evidence tier
    if trial.get("phase") == "Phase 3":
        evidence_tier = "STANDARD"
    elif trial.get("phase") == "Phase 2":
        evidence_tier = "SUPPORTED"
    else:
        evidence_tier = "INVESTIGATIONAL"
    
    # Determine enrollment likelihood
    eligible_count = len(why_eligible) + len(why_good_fit)
    red_flag_count = len(red_flags) + len(conditional)
    
    if eligible_count >= 5 and red_flag_count == 0:
        enrollment_likelihood = "HIGH"
    elif eligible_count >= 3:
        enrollment_likelihood = "MEDIUM"
    else:
        enrollment_likelihood = "LOW"
    
    # Clamp score
    score = max(0.0, min(1.0, score))
    
    reasoning = TrialMatchReasoning(
        match_score=score,
        why_eligible=why_eligible,
        why_good_fit=why_good_fit,
        conditional_requirements=conditional,
        red_flags=red_flags,
        evidence_tier=evidence_tier,
        enrollment_likelihood=enrollment_likelihood,
        ca125_intelligence=ca125_intel,
        germline_context={
            "status": profile.germline_status,
            "brca": profile.germline_brca,
            "implication": "Sporadic cancer - focus on tumor biomarkers, not germline"
        }
    )
    
    return score, reasoning
```

---

### **File 3: Update Main Router** (5 min)

**Location**: `oncology-coPilot/oncology-backend-minimal/api/main.py`

```python
from .routers import ayesha_trials

# Register router
app.include_router(ayesha_trials.router)
```

---

## 🎨 **FRONTEND IMPLEMENTATION (JR OR ZO - 4.5 HOURS)**

### **Component 1: AyeshaTrialExplorer** (2 hours)

**Location**: `oncology-coPilot/oncology-frontend/src/pages/AyeshaTrialExplorer.jsx`

### **Component 2: TrialMatchCard** (1.5 hours)

**Location**: `oncology-coPilot/oncology-frontend/src/components/trials/TrialMatchCard.jsx`

### **Component 3: CA125Tracker** (1 hour)

**Location**: `oncology-coPilot/oncology-frontend/src/components/ayesha/CA125Tracker.jsx`

---

## 🧪 **TESTING (30 MIN)**

### **Test 1: Ayesha's Exact Profile**

```bash
curl -X POST http://127.0.0.1:8000/api/ayesha/trials/search \
  -H 'Content-Type: application/json' \
  -d '{
    "ca125": 2842,
    "germline_status": "NEGATIVE",
    "stage": "IVB",
    "treatment_line": 0,
    "location": "NYC"
  }'
```

**Expected**:
- 10 trials returned
- All first-line or all-comers
- None require germline BRCA
- NYC locations prioritized
- CA-125 tracking trials boosted
- Transparent reasoning for each

---

## ⚔️ **DELIVERABLES (END OF 8 HOURS)**

**For Ayesha**:
1. ✅ Top 10 clinical trials (ranked by match score)
2. ✅ Each with transparent reasoning (why eligible, why good fit)
3. ✅ Contact info for each trial
4. ✅ CA-125 tracking plan (what to expect, when to worry)
5. ✅ NGS-gated features (unlock when report arrives)

**For Platform**:
1. ✅ Ayesha-specific trial endpoint
2. ✅ CA-125 intelligence module
3. ✅ Germline-negative filtering logic
4. ✅ Transparent reasoning system
5. ✅ Frontend trial explorer

---

## 🎯 **AGENT JR - STEP-BY-STEP EXECUTION CHECKLIST**

### **PHASE 1: BACKEND (3 HOURS)** ⚔️

**DO THESE IN ORDER**:

#### **Step 1: Create Schemas** (30 min)
- [ ] Create file: `oncology-coPilot/oncology-backend-minimal/api/schemas/ayesha_trials.py`
- [ ] Copy code from "File 1: Ayesha Profile Schema" section above
- [ ] **CRITICAL**: Use exact schema structure (no modifications)
- [ ] Test imports: `from api.schemas.ayesha_trials import AyeshaTrialProfile`

**Acceptance**: File created, imports work, no errors

---

#### **Step 2: Build CA-125 Intelligence** (30 min)
- [ ] Create file: `oncology-coPilot/oncology-backend-minimal/api/services/ca125_intelligence.py`
- [ ] Copy code from "File 2: CA-125 Intelligence Service" section above
- [ ] **CRITICAL**: Keep scoring thresholds as-is (don't adjust)
- [ ] Test: `CA125IntelligenceService.analyze(2842)` returns expected dict

**Acceptance**: Service works, returns intelligence for CA-125 = 2842

---

#### **Step 3: Build Trial Router** (2 hours)
- [ ] Create file: `oncology-coPilot/oncology-backend-minimal/api/routers/ayesha_trials.py`
- [ ] Copy code from "File 3: Ayesha Trial Router" section above
- [ ] **CRITICAL NOTES**:
  - **READ** existing `hybrid_trial_search.py` to understand schema
  - **DON'T** rewrite hybrid search - just call it
  - **DO** add reasoning logic exactly as shown
  - **HANDLE** missing fields gracefully (use `.get()` with defaults)
- [ ] Update `api/main.py` to register router (see "File 3: Update Main Router")

**Acceptance**: 
- Endpoint `/api/ayesha/trials/search` returns 200
- Returns structure matches schema
- Reasoning fields populated

---

#### **Step 4: Seed AstraDB** (Automated - 16 min)
- [ ] Run: `cd oncology-coPilot/oncology-backend-minimal && venv/bin/python scripts/seed_astra_trials.py --disease ovarian --count 200`
- [ ] **NOTE**: This is automated, just monitor for errors
- [ ] Verify: 200+ trials in AstraDB collection

**Acceptance**: AstraDB has ≥200 ovarian cancer trials

---

### **PHASE 2: FRONTEND (4.5 HOURS)** ⚔️

**DO THESE IN ORDER**:

#### **Step 5: Create TrialMatchCard Component** (1.5 hours)
- [ ] Create file: `oncology-coPilot/oncology-frontend/src/components/trials/TrialMatchCard.jsx`
- [ ] **DISPLAY** (NO HARDCODING):
  ```jsx
  // Trial header
  {trial.nct_id || 'Unknown'} - {trial.title || 'No title'}
  
  // Match score (dynamic)
  Match: {Math.round((trial.match_score || 0) * 100)}%
  <LinearProgress value={(trial.match_score || 0) * 100} />
  
  // Why Eligible section (iterate dynamically)
  {trial.reasoning?.why_eligible?.map((reason, idx) => (
    <Typography key={idx}>{reason}</Typography>
  ))}
  
  // Why Good Fit section (iterate dynamically)
  {trial.reasoning?.why_good_fit?.map((reason, idx) => (
    <Typography key={idx}>{reason}</Typography>
  ))}
  
  // Conditional Requirements (iterate dynamically)
  {trial.reasoning?.conditional_requirements?.length > 0 && (
    <Alert severity="warning">
      {trial.reasoning.conditional_requirements.map((req, idx) => (
        <Typography key={idx}>{req}</Typography>
      ))}
    </Alert>
  )}
  
  // Red Flags (iterate dynamically)
  {trial.reasoning?.red_flags?.length > 0 && (
    <Alert severity="error">
      {trial.reasoning.red_flags.map((flag, idx) => (
        <Typography key={idx}>{flag}</Typography>
      ))}
    </Alert>
  )}
  
  // Contact info (conditional)
  {trial.contact_name && (
    <Box>
      <Typography>Contact: {trial.contact_name}</Typography>
      {trial.contact_phone && <Typography>Phone: {trial.contact_phone}</Typography>}
      {trial.contact_email && <Typography>Email: {trial.contact_email}</Typography>}
    </Box>
  )}
  ```

- [ ] **STYLING**: Match existing trial cards (see `ResultsDisplay.jsx` for reference)
- [ ] **EXPORT**: Add to `src/components/trials/index.js`

**Acceptance**: 
- Card renders without errors
- Shows all sections dynamically
- No hardcoded values
- Handles missing fields gracefully

---

#### **Step 6: Create CA125Tracker Component** (1 hour)
- [ ] Create file: `oncology-coPilot/oncology-frontend/src/components/ayesha/CA125Tracker.jsx`
- [ ] **DISPLAY** (dynamic from backend):
  ```jsx
  // Current CA-125
  Current: {ca125_intelligence?.ca125_value || 'Unknown'}
  Disease Burden: {ca125_intelligence?.disease_burden || 'Unknown'}
  
  // Expected response (iterate dynamically)
  Expected Response (if platinum-sensitive):
  {Object.entries(ca125_intelligence?.expected_response || {}).map(([key, value]) => (
    <Typography key={key}>
      {key}: {value}
    </Typography>
  ))}
  
  // Monitoring strategy
  {ca125_intelligence?.monitoring_strategy}
  ```

- [ ] **STYLING**: Card-based, color-coded by burden level
- [ ] **EXPORT**: Add to `src/components/ayesha/index.js`

**Acceptance**: 
- Displays CA-125 intelligence dynamically
- Shows expected response ranges
- No hardcoded values

---

#### **Step 7: Create AyeshaTrialExplorer Page** (2 hours)
- [ ] Create file: `oncology-coPilot/oncology-frontend/src/pages/AyeshaTrialExplorer.jsx`
- [ ] **STRUCTURE**:
  ```jsx
  // Top section: Ayesha's profile summary
  <AyeshaClinicalProfile 
    ca125={2842}
    germlineStatus="NEGATIVE"
    stage="IVB"
    awaiting={["Tumor NGS", "Somatic BRCA", "HRD score"]}
  />
  
  // CA-125 tracker
  <CA125Tracker intelligence={ca125_intelligence} />
  
  // Trials list (iterate dynamically - NO HARDCODING)
  {trials?.map((trial, idx) => (
    <TrialMatchCard key={trial.nct_id || idx} trial={trial} />
  )) || <EmptyState message="No trials found" />}
  
  // Provenance footer
  <ProvenanceCard 
    filtersApplied={provenance?.filters_applied}
    totalScreened={total_screened}
  />
  ```

- [ ] **API CALL**:
  ```jsx
  const fetchTrials = async () => {
    const response = await fetch('/api/ayesha/trials/search', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({
        ca125: 2842,
        germline_status: "NEGATIVE",
        stage: "IVB",
        treatment_line: 0,
        location: "NYC"
      })
    });
    const data = await response.json();
    setTrials(data.trials);  // DYNAMIC - don't assume count
    setCA125Intelligence(data.ca125_intelligence);
  };
  ```

- [ ] **ROUTING**: Add to `App.jsx` as `/ayesha-trials`
- [ ] **SIDEBAR**: Add link "Ayesha's Trials" to navigation

**Acceptance**: 
- Page loads without errors
- Calls backend endpoint
- Displays all trials dynamically
- Routing works
- No hardcoded trial data

---

### **PHASE 3: TESTING (30 MIN)** ⚔️

#### **Step 8: Backend Smoke Test**
- [ ] Start backend: `cd oncology-coPilot/oncology-backend-minimal && venv/bin/uvicorn api.main:app --reload`
- [ ] Run curl (see "Test 1: Ayesha's Exact Profile" above)
- [ ] Verify: Returns 10 trials with reasoning

**Acceptance**: Backend returns trials with transparent reasoning

---

#### **Step 9: Frontend E2E Test**
- [ ] Start frontend: `cd oncology-coPilot/oncology-frontend && npm start`
- [ ] Navigate to `/ayesha-trials`
- [ ] Verify:
  - CA-125 tracker displays correctly
  - All 10 trials render (no hardcoded limits)
  - Reasoning sections populated
  - Contact info shows (if available)
  - No console errors

**Acceptance**: Frontend works end-to-end with live backend

---

#### **Step 10: Edge Case Testing**
- [ ] Test with empty response (0 trials) - should show empty state
- [ ] Test with missing reasoning fields - should show gracefully
- [ ] Test with 1 trial - should render without assuming ≥3
- [ ] Test with 20 trials - should show all (not slice to 10)

**Acceptance**: Handles all edge cases without crashing

---

## 🚨 **AGENT JR - CRITICAL EXECUTION NOTES**

### **DO's** ✅:
1. ✅ **READ existing code FIRST** before writing anything
   - Check `api/routers/trials.py` for schema patterns
   - Check `api/services/hybrid_trial_search.py` for how trials are structured
   - Check `src/components/research/ResultsDisplay.jsx` for styling patterns

2. ✅ **USE existing services** (don't rewrite):
   - `HybridTrialSearchService` already exists - just call it
   - `ClinicalTrialSearchService` already exists - reuse
   - `AstraDB` already seeded with 30 trials - expand to 200

3. ✅ **FOLLOW dynamic patterns**:
   - Use `.map()` for ALL arrays
   - Use optional chaining (`?.`)
   - Provide fallbacks (`|| 'Unknown'`)
   - Conditional rendering (`{field && <Component />}`)

4. ✅ **COPY code templates EXACTLY** (from this doc):
   - Schemas (File 1)
   - CA-125 Intelligence (File 2)
   - Trial Router (File 3)
   - Don't modify unless you find errors

5. ✅ **TEST incrementally**:
   - After each file, run imports to verify
   - After backend complete, run curl test
   - After frontend complete, test in browser

---

### **DON'T's** ❌:

1. ❌ **DON'T rewrite existing services**:
   - `hybrid_trial_search.py` already works - just call it
   - `clinical_trial_search_service.py` already works - reuse
   - Don't create new search logic from scratch

2. ❌ **DON'T hardcode trial data**:
   - No hardcoded trial lists
   - No hardcoded scoring thresholds (use backend values)
   - No hardcoded reasoning (iterate over backend arrays)

3. ❌ **DON'T skip error handling**:
   - Always check if data exists before accessing
   - Always provide empty states
   - Always handle API failures gracefully

4. ❌ **DON'T create new .mdc files**:
   - Use THIS file for all notes
   - Update THIS file with progress
   - Create completion report in `.cursor/ayesha/AGENT_JR_TRIALS_COMPLETION.md`

5. ❌ **DON'T skip testing**:
   - Test backend before frontend
   - Test with real Ayesha profile
   - Test edge cases (empty, missing fields)

---

## 📋 **AGENT JR - PROGRESS TRACKER**

**Update this section as you complete each step**:

### **Backend Progress**:
- [ ] Step 1: Schemas created ✅ / ⏸️ / ❌
- [ ] Step 2: CA-125 intelligence built ✅ / ⏸️ / ❌
- [ ] Step 3: Trial router built ✅ / ⏸️ / ❌
- [ ] Step 4: AstraDB seeded ✅ / ⏸️ / ❌

### **Frontend Progress**:
- [ ] Step 5: TrialMatchCard component ✅ / ⏸️ / ❌
- [ ] Step 6: CA125Tracker component ✅ / ⏸️ / ❌
- [ ] Step 7: AyeshaTrialExplorer page ✅ / ⏸️ / ❌

### **Testing Progress**:
- [ ] Step 8: Backend smoke test ✅ / ⏸️ / ❌
- [ ] Step 9: Frontend E2E test ✅ / ⏸️ / ❌
- [ ] Step 10: Edge case testing ✅ / ⏸️ / ❌

---

## 🚨 **IF YOU GET STUCK (AGENT JR)**

### **Question 1: Backend service doesn't exist**
**Check**: Does `api/services/hybrid_trial_search.py` exist?
- If YES: Read it to understand schema, then call it
- If NO: Ask Zo - may need to use different service

### **Question 2: Trial schema doesn't match**
**Check**: Run backend endpoint and inspect actual response
- Copy REAL response structure
- Adapt code to match REAL schema (not assumed)
- Use optional chaining for safety

### **Question 3: AstraDB seeding fails**
**Check**: Is `ASTRA_DB_API_ENDPOINT` set in `.env`?
- If NO: Ask Zo for credentials
- If YES but fails: Check error message, may be quota issue

### **Question 4: Frontend component crashes**
**Check**: Console errors - usually missing null checks
- Add guards: `if (!data) return null;`
- Use optional chaining: `data?.field`
- Provide fallbacks: `|| 'Unknown'`

### **Question 5: Not sure about styling**
**Check**: Look at existing components:
- `src/components/research/ResultsDisplay.jsx` (trial cards)
- `src/components/sporadic/SporadicProvenanceCard.jsx` (card styling)
- Match their MUI theme usage

---

## 📞 **REPORTING TO ZO**

### **Every 2 Hours** (Progress Update):
Create/Update: `.cursor/ayesha/AGENT_JR_TRIALS_PROGRESS.md`
```markdown
## Hour 0-2 Progress:
- ✅ Completed: Steps 1-2 (schemas + CA-125 intelligence)
- 🔄 In Progress: Step 3 (trial router)
- ❌ Blocked: None
- ⏱️ On Schedule: Yes
```

### **If Blocked** (Immediate):
Create: `.cursor/ayesha/AGENT_JR_TRIALS_QUESTIONS.md`
```markdown
## BLOCKER: [Describe issue]
- What I tried: ...
- Error message: ...
- What I need: ...
```

### **When Complete** (Final Report):
Create: `.cursor/ayesha/AGENT_JR_TRIALS_COMPLETION.md`
```markdown
## ✅ MISSION COMPLETE

**What I Built**:
- Backend: 3 files (schemas, service, router)
- Frontend: 3 components (card, tracker, page)
- Tests: All passing

**Deliverables**:
- Top 10 trials for Ayesha
- CA-125 intelligence
- Transparent reasoning

**Deviations**: [List any changes from plan]
**Bugs Found**: [List any issues]
**Time Taken**: X hours (vs 8 hour estimate)
```

---

## ⚔️ **ZO'S OVERSIGHT PLAN (WHILE JR EXECUTES)**

**What Zo Will Do** (parallel to Jr):
1. ✅ Strategic planning for post-NGS workflow
2. ✅ Design NGS report parser (Foundation/Tempus)
3. ✅ Plan WIWFM integration with tumor biomarkers
4. ✅ Design metastasis interception strategy (when stable)
5. ✅ Review Jr's progress every 2 hours
6. ✅ Unblock Jr if needed

**What Zo Won't Do**:
- ❌ Interfere with Jr's execution (unless asked)
- ❌ Duplicate Jr's work
- ❌ Change requirements mid-mission

---

## 🎯 **SUCCESS METRICS**

**Jr's mission is successful when**:
- ✅ `/api/ayesha/trials/search` endpoint works (returns 10 trials)
- ✅ Frontend page renders all trials dynamically
- ✅ Each trial has transparent reasoning (why matched)
- ✅ CA-125 tracker displays intelligence
- ✅ No hardcoded values (all data-driven)
- ✅ No console errors
- ✅ All tests passing
- ✅ Ayesha can see top 10 trials with contact info

**Timeline**: 8 hours  
**Output**: **ACTIONABLE trial list for Ayesha**

---

## ⚔️ **AGENT JR - READY TO EXECUTE?**

**Your Mission**:
Build the trial filtering system that finds Ayesha's best trials.

**Your Timeline**: 8 hours

**Your Deliverable**: Top 10 trials with reasoning

**Zo's Role**: Oversee + strategic planning (won't interfere unless you ask)

**Commander has cleared you to proceed.** 🎯

**FIRE IN THE HOLE, SOLDIER.** ⚔️

