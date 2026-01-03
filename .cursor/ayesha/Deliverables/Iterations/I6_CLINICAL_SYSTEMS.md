# 🏥 ITERATION 6: CLINICAL SYSTEMS & WORKFLOWS

**Status**: ✅ **COMPLETE**  
**Duration**: 3-4 hours  
**Created**: January 14, 2025

---

### **6.1 AYESHA COMPLETE CARE ORCHESTRATOR**

#### **6.1.1 Overview** (`routers/ayesha_orchestrator_v2.py`):
- **Purpose**: Unified care plan orchestration for AK (Stage IVB HGSOC)
- **Endpoint**: `POST /api/ayesha/complete_care_v2`
- **Integration**: Coordinates 7 clinical services + 3 SAE services + Resistance Prophet

#### **6.1.2 Orchestrated Services**:
1. **Clinical Trials** (`_call_ayesha_trials`):
   - Frontline, NYC metro, transparent reasoning
   - Includes SOC recommendation + CA-125 intelligence
2. **Drug Efficacy (WIWFM)** (`_call_drug_efficacy`):
   - Full S/P/E if `tumor_context` provided
   - Returns "awaiting NGS" message if no tumor data
3. **Food Validator** (`_call_food_validator`):
   - Optional supplement/nutrition recommendations
   - Calls `/api/hypothesis/validate_food_dynamic`
4. **Resistance Playbook** (`_call_resistance_playbook`):
   - Next-line planning based on resistance mechanisms
   - Requires `tumor_context` (returns "awaiting NGS" if missing)
5. **CA-125 Intelligence** (via trials endpoint):
   - Burden classification, response forecast, resistance signals
6. **SOC Recommendation** (via trials endpoint):
   - NCCN-aligned carboplatin + paclitaxel + bevacizumab
7. **Resistance Prophet** (opt-in):
   - Predicts resistance 3-6 months early
   - Manager Q7: Opt-in via `include_resistance_prediction=true`

#### **6.1.3 Phase 1 SAE Services** (Pre-NGS):
- **Next-Test Recommender**: HRD → ctDNA → SLFN11 → ABCB1 priority
- **Hint Tiles**: Max 4 tiles (Test → Trials → Monitor → Avoid)
- **Mechanism Map**: Pre-NGS (all gray) vs Post-NGS (color-coded)

#### **6.1.4 Phase 2 SAE Services** (Post-NGS):
- **SAE Features**: DNA repair capacity, mechanism vector, resistance signals
- **Resistance Alert**: 2-of-3 trigger rule (HRD drop, DNA repair drop, CA-125 inadequate)

#### **6.1.5 Request Schema** (`CompleteCareV2Request`):
```python
ca125_value: float  # Current CA-125 (U/mL)
stage: str  # "IVB"
treatment_line: str  # "first-line"
germline_status: str  # "negative"
has_ascites: bool
has_peritoneal_disease: bool
tumor_context: Optional[Dict]  # NGS data (somatic mutations, HRD, TMB, MSI)
drug_query: Optional[str]  # Specific drug to evaluate
food_query: Optional[str]  # Food/supplement to validate
include_trials: bool = True
include_soc: bool = True
include_ca125: bool = True
include_wiwfm: bool = True
include_food: bool = False
include_resistance: bool = False
include_resistance_prediction: bool = False  # Manager Q7: Opt-in
```

---

### **6.2 CA-125 INTELLIGENCE SERVICE**

#### **6.2.1 Purpose** (`services/ca125_intelligence.py`):
- Burden classification (MINIMAL/MODERATE/SIGNIFICANT/EXTENSIVE)
- Response forecast (cycle 3, cycle 6 targets)
- Resistance signal detection (on-therapy rise, inadequate response)
- Monitoring strategy recommendations

#### **6.2.2 Burden Classification**:
```python
BURDEN_THRESHOLDS = {
    "MINIMAL": (0, 100),
    "MODERATE": (100, 500),
    "SIGNIFICANT": (500, 1000),
    "EXTENSIVE": (1000, float('inf'))
}
```

#### **6.2.3 Response Expectations** (GOG-218, ICON7):
- **Cycle 3**: ≥70% drop expected
- **Cycle 6**: ≥90% drop expected
- **Complete Response**: <35 U/mL target
- **Resistance Signal**: <50% drop by cycle 3

#### **6.2.4 Resistance Signal Detection**:
1. **On-therapy rise**: CA-125 increases while on treatment
2. **Inadequate response**: <50% drop by cycle 3
3. **Minimal response**: <30% drop by any cycle post-cycle 2

#### **6.2.5 Monitoring Strategy**:
- **On treatment**: Every 3 weeks (before each cycle)
- **Pre-treatment (EXTENSIVE)**: Every 2 weeks
- **Pre-treatment (SIGNIFICANT)**: Every 4 weeks
- **Surveillance**: Every 3 months

---

### **6.3 RESISTANCE DETECTION SERVICE**

#### **6.3.1 Purpose** (`services/resistance_detection_service.py`):
- **2-of-3 Trigger Rule** (Manager C7):
  - HRD drop >= 15 points
  - DNA repair capacity drop >= 0.20
  - CA-125 inadequate response (on-therapy rise OR <50% drop by cycle 3)
- **HR Restoration Pattern** (Manager R2):
  - HRD drop + DNA repair drop (coherent signal)
  - Immediate alert (don't wait for radiology)

#### **6.3.2 Resistance Alert Structure**:
```python
@dataclass
class ResistanceAlert:
    resistance_detected: bool  # True if 2-of-3 triggers met
    hr_restoration_suspected: bool  # True if HR restoration pattern
    immediate_alert: bool  # True if should alert NOW
    triggers_met: List[str]  # ["hrd_drop", "dna_repair_drop", "ca125_inadequate"]
    trigger_count: int  # Number of triggers met (need ≥2)
    recommended_actions: List[str]  # Clinical actions
    recommended_trials: List[str]  # ATR/CHK1, WEE1 trials
```

#### **6.3.3 Recommended Actions**:
- **IMMEDIATE**: Alert oncologist, order updated HRD test, order ctDNA panel
- **If HR restoration**: Consider switching from PARP
- **Trial recommendations**: ATR inhibitor trials, ATR + PARP combos, CHK1/WEE1 inhibitors

---

### **6.4 RESISTANCE PROPHET SERVICE**

#### **6.4.1 Purpose** (`services/resistance_prophet_service.py`):
- **Predicts resistance 3-6 months early** (before clinical progression)
- **Phase 1**: DNA repair restoration + pathway escape (NO CA-125)
- **Phase 1b**: Add CA-125 kinetics (prospective, Ayesha live)
- **Manager Q7**: Opt-in via `include_resistance_prediction=true`

#### **6.4.2 Resistance Signals**:
1. **DNA Repair Restoration** (Signal 1):
   - Threshold: >20% increase from baseline
   - Indicates PARP resistance mechanism
2. **Pathway Escape** (Signal 2):
   - Threshold: >=30% shift in mechanism vector (L2 distance)
   - Indicates bypass resistance mechanism
3. **CA-125 Kinetics** (Signal 3, Phase 1b+):
   - On-therapy rise or inadequate response
   - Uses existing CA125Intelligence service

#### **6.4.3 Risk Stratification** (Manager Q9):
- **HIGH**: Probability >=0.70 AND >=2 signals
- **MEDIUM**: 0.50-0.69 OR exactly 1 signal
- **LOW**: <0.50 probability
- **Manager Q15**: Cap at MEDIUM if no CA-125 unless >=2 non-CA-125 signals

#### **6.4.4 Urgency & Actions** (Manager Q10):
- **CRITICAL** (HIGH risk):
  - ESCALATE_IMAGING within 1 week
  - CONSIDER_SWITCH within 2 weeks
  - REVIEW_RESISTANCE_PLAYBOOK within 1 week
- **ELEVATED** (MEDIUM risk):
  - MONITOR_WEEKLY x4 weeks
  - REASSESS after 4 weeks
- **ROUTINE** (LOW risk):
  - Routine monitoring per standard of care

#### **6.4.5 Next-Line Options** (Manager Q11):
- Consults `ResistancePlaybookService` for mechanism-specific strategies
- Consults `TreatmentLineService` for appropriateness

---

### **6.5 RESISTANCE PLAYBOOK SERVICE**

#### **6.5.1 Purpose** (`services/resistance_playbook_service.py`):
- Predicts resistance mechanisms (HR restoration, ABCB1, MAPK, PI3K, SLFN11)
- Ranks combo strategies (PARP+ATR, PARP+VEGF, IO combos, MAPK/PI3K)
- Recommends next-line switches
- Generates trial keywords for filtering

#### **6.5.2 Resistance Detection Rules**:
1. **HR Restoration**:
   - Signals: HRD drop after PARP, RAD51C/D mutations, loss of BRCA1/2 LOH, SAE DNA repair capacity >0.7
   - Triggers: PARP inhibitors
2. **ABCB1 Upregulation**:
   - Signals: ABCB1 copy number gain (>4 copies), ABCB1 activating mutations
   - Triggers: Paclitaxel, doxorubicin, topotecan, P-gp substrates
3. **MAPK Activation**:
   - Signals: KRAS/NRAS/BRAF activating mutations, high MAPK pathway burden
   - Triggers: BRAF inhibitors, EGFR inhibitors
4. **PI3K Activation**:
   - Signals: PIK3CA mutations, PTEN loss, high PI3K pathway burden
   - Triggers: PI3K inhibitors, mTOR inhibitors
5. **SLFN11 Deficiency**:
   - Signals: SLFN11 low expression (IHC), SLFN11 mutations
   - Triggers: PARP inhibitors, topoisomerase inhibitors

#### **6.5.3 Combo Strategies**:
- **PARP + ATR**: For HRD-high with HR restoration risk
- **PARP + VEGF**: For HRD-high with angiogenesis escape
- **IO Combos**: For TMB-high/MSI-high with IO resistance
- **MAPK/PI3K Combos**: For pathway escape mechanisms

---

### **6.6 TRIAL MATCHING SYSTEM**

#### **6.6.1 Architecture** (`services/ayesha_trial_matching/`):
- **MatchOrchestrator**: Coordinates workflow
- **EligibilityFilters**: Hard filters (MUST-MATCH)
- **ScoringEngine**: Soft boosts (scoring)
- **ReasoningGenerator**: Transparent "why-matched" explanations

#### **6.6.2 Eligibility Filters** (Hard Filters):
1. **Disease**: Ovarian/peritoneal/gynecologic
2. **Stage**: IV/advanced/metastatic
3. **Treatment Line**: First-line/untreated
4. **Status**: Recruiting/Active
5. **Location**: NYC metro (NY/NJ/CT)
6. **Exclusions**: NOT recurrent-only, NOT germline-BRCA-required

#### **6.6.3 Scoring Engine** (Soft Boosts):
- **Base Score**: 0.5
- **Boosts**:
  - First-line trial: +0.30
  - Stage IV specific: +0.25
  - SOC backbone (carboplatin + paclitaxel): +0.20
  - Germline-negative friendly: +0.20
  - IP chemotherapy: +0.20
  - Bevacizumab: +0.15
  - CA-125 tracking: +0.15
  - NYC location: +0.15
  - Large trial: +0.10
  - Phase III: +0.10
- **Penalties**:
  - Germline BRCA required: -0.30
  - Distance: -0.25
  - Phase I: -0.20

#### **6.6.4 Reasoning Generator**:
- Generates transparent explanations for each trial match
- Explains why trial matched (eligibility + scoring boosts)
- Includes clinical rationale

---

### **6.7 SOC RECOMMENDATION**

#### **6.7.1 Purpose** (`routers/ayesha_trials.py:_generate_soc_recommendation`):
- **NCCN-aligned** standard of care for Stage IVB HGSOC
- **Confidence**: 95-100% (guideline-aligned, no predictions)

#### **6.7.2 Regimen**:
- **Base**: Carboplatin AUC 5-6 + Paclitaxel 175 mg/m²
- **Add-on** (if ascites/peritoneal disease):
  - Bevacizumab 15 mg/kg
  - Rationale: Reduces progression risk (GOG-218 HR 0.72, ICON7 HR 0.81)

#### **6.7.3 Schedule**:
- **Induction**: 6 cycles every 21 days
- **Maintenance**: Bevacizumab continuation up to 15 months total OR progression
- **Typical Duration**: ~18 weeks induction + up to 12 months maintenance

#### **6.7.4 Monitoring Protocol**:
- **CA-125**: Every cycle (every 3 weeks)
- **Imaging**: Every 3 cycles
- **Toxicity**: Grade 3-4 neuropathy, cytopenias → dose reduction/delay

---

### **6.8 FOOD VALIDATOR**

#### **6.8.1 Purpose** (`routers/hypothesis_validator.py`):
- **Dynamic Food Validator**: Works for ANY food/supplement
- **Endpoint**: `POST /api/hypothesis/validate_food_dynamic`
- **A→B Dependency**: NO TUMOR NGS REQUIRED (works pre-NGS)

#### **6.8.2 Features**:
- **LLM Literature Mining**: PubMed paper reading
- **SAE Features**: Line appropriateness, cross-resistance, sequencing fitness
- **Complete Provenance**: Data sources, confidence breakdown

#### **6.8.3 Treatment Line Service** (`services/food_treatment_line_service.py`):
- Computes SAE features for dietary supplements:
  - `line_appropriateness`: 0-1 score
  - `cross_resistance`: 0-1 score
  - `sequencing_fitness`: 0-1 score
- **Biomarker Gates**: Boosts if biomarker context matches
- **Treatment History Gates**: Adjusts based on prior therapies

---

### **6.9 NEXT-TEST RECOMMENDER**

#### **6.9.1 Purpose** (`services/next_test_recommender.py`):
- **Prioritized biomarker testing** recommendations
- **Manager Policy**: P1, C6 (MANAGER_ANSWERS_TO_ZO_SAE_QUESTIONS.md)

#### **6.9.2 Priority Order**:
1. **HRD** (MyChoice CDx): PARP gate, 10d turnaround
2. **ctDNA** (Guardant360/FoundationOne): MSI/TMB + somatic HRR, 7d turnaround
3. **SLFN11 IHC**: PARP sensitivity, 5d turnaround
4. **ABCB1 proxy**: Efflux resistance (only if prior taxane), 5d turnaround

#### **6.9.3 Format**:
- **Differential branches**: "If positive → X; If negative → Y"
- **Turnaround + Cost**: Included in recommendations
- **Max 4 recommendations**

---

### **6.10 HINT TILES SERVICE**

#### **6.10.1 Purpose** (`services/hint_tiles_service.py`):
- **Clinician action hints** (max 4 tiles)
- **Priority Order** (Manager C8): Test → Trials → Monitor → Avoid
- **Suggestive tone** (not prescriptive)

#### **6.10.2 Tile Types**:
1. **Next Test**: Prioritized biomarker testing
2. **Trial Matched**: Clinical trial opportunities
3. **Monitoring**: CA-125, imaging schedules
4. **Avoid**: Drug interactions, contraindications

---

### **6.11 MECHANISM MAP SERVICE**

#### **6.11.1 Purpose** (`services/mechanism_map_service.py`):
- **Pathway burden visualization** (6 chips)
- **Pre-NGS**: All gray (no data)
- **Post-NGS**: Color-coded by pathway activation (DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux)

---

### **6.12 KEY INSIGHTS**

#### **Clinical Workflow**:
1. **Pre-NGS**: Phase 1 SAE services (next-test, hint tiles, mechanism map)
2. **Post-NGS**: Phase 2 SAE services (SAE features, resistance alert)
3. **Opt-in**: Resistance Prophet (predicts resistance 3-6 months early)

#### **Resistance Detection**:
1. **2-of-3 Trigger Rule**: HRD drop OR DNA repair drop OR CA-125 inadequate
2. **HR Restoration Pattern**: Immediate alert (don't wait for radiology)
3. **Resistance Prophet**: Early warning (3-6 months before clinical progression)

#### **Trial Matching**:
1. **Hard Filters**: Disease, stage, line, status, location (MUST-MATCH)
2. **Soft Boosts**: First-line, SOC backbone, CA-125 tracking, NYC location
3. **Transparent Reasoning**: Explains why each trial matched

#### **CA-125 Intelligence**:
1. **Burden Classification**: MINIMAL/MODERATE/SIGNIFICANT/EXTENSIVE
2. **Response Forecast**: Cycle 3 (≥70% drop), Cycle 6 (≥90% drop)
3. **Resistance Signals**: On-therapy rise, inadequate response (<50% drop by cycle 3)

#### **SOC Recommendation**:
1. **NCCN-aligned**: Carboplatin + Paclitaxel (standard)
2. **Bevacizumab Add-on**: If ascites/peritoneal disease
3. **Confidence**: 95-100% (guideline-based, not prediction)

---

**Status**: ✅ **ITERATION 6 COMPLETE** - Clinical Systems & Workflows fully documented  
**Next**: I7 - Research & Design Systems (Metastasis Interception, CRISPR Design, Protein Synthesis)

---