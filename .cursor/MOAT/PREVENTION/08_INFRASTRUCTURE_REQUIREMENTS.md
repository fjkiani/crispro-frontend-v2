# 🛠️ INFRASTRUCTURE REQUIREMENTS

**Date:** December 24, 2025  
**Status:** 📋 **BUILD PLAN**  
**Location:** `.cursor/MOAT/PREVENTION/08_INFRASTRUCTURE_REQUIREMENTS.md`

---

## 🎯 WHAT WE NEED TO BUILD (2025-2027)

To enable DDR_bin-guided intervention, we need 4 core infrastructure components:

| # | Component | Purpose | Priority |
|---|-----------|---------|----------|
| 1 | **Real-Time DDR_bin Dashboard** | Monitor pathway bins over time | 🔴 HIGH |
| 2 | **Clinical Decision Support** | Recommend interventions | 🔴 HIGH |
| 3 | **Longitudinal Biospecimen** | Track clonal evolution | 🟡 MEDIUM |
| 4 | **Adaptive Trial Platform** | DDR_bin-guided clinical trials | 🟡 MEDIUM |

---

## 1️⃣ REAL-TIME DDR_BIN MONITORING DASHBOARD

### **Patient View Mockup:**

```
┌─────────────────────────────────────────────────────────────────────────┐
│ 🩺 PATIENT: AK (AK_001) - Ovarian Cancer, Stage IV           │
├─────────────────────────────────────────────────────────────────────────┤
│ Current Treatment: Olaparib 300mg BID (Cycle 12)                        │
│ Days on therapy: 270                                                     │
│                                                                          │
│ DDR_bin Score: 0.82 ⚠️ (Baseline: 0.88)                                 │
│ Trend: ↓ Declining (6% drop over 9 months)                              │
│                                                                          │
│ 🚨 ALERT: HR restoration detected                                       │
│    Resistant clone estimated: 5-8%                                       │
│    Recommendation: Consider combination therapy                          │
│                                                                          │
│ Pathway Status:                                                          │
│  DDR_bin:  0.82 ⚠️ [████████████░░░░░░░░] 82%                          │
│  MAPK_bin: 0.15 ✅ [███░░░░░░░░░░░░░░░░] 15%                           │
│  PI3K_bin: 0.22 ✅ [████░░░░░░░░░░░░░░░] 22%                           │
│  IO_bin:   0.05 ✅ [█░░░░░░░░░░░░░░░░░░] 5%                            │
│  Efflux:   0.08 ✅ [█░░░░░░░░░░░░░░░░░░] 8%                            │
│                                                                          │
│ Feature Breakdown (9 DDR diamond features):                              │
│  Feature 27607: 0.87 (↓ from 0.90) - TP53 dominant                      │
│  Feature 16337: 0.78 (↓ from 0.82) - TP53 dominant                      │
│  Feature 26220: 0.84 (↓ from 0.88) - TP53 dominant                      │
│  ...                                                                     │
│                                                                          │
│ Historical DDR_bin:                                                      │
│  Month 0:  0.88 ████████████████████████████████████████ Baseline       │
│  Month 3:  0.87 ███████████████████████████████████████░ Stable         │
│  Month 6:  0.86 ██████████████████████████████████████░░ Stable         │
│  Month 9:  0.82 ████████████████████████████████░░░░░░░░ ⚠️ ALERT      │
│                                                                          │
│ Next Steps:                                                              │
│  1. Order ctDNA for RAD51C/BRCA1/PALB2 sequencing                       │
│  2. Review combination therapy options                                   │
│  3. Schedule oncology discussion within 2 weeks                          │
└─────────────────────────────────────────────────────────────────────────┘
```

### **Technical Requirements:**

| Component | Technology | Status |
|-----------|------------|--------|
| DDR_bin calculation | TRUE SAE extraction | ✅ VALIDATED |
| Historical tracking | State management | ✅ EXISTS |
| Trend analysis | Time-series analysis | ⏳ TO BUILD |
| Alert thresholds | Rule engine | ⏳ TO BUILD |
| Dashboard UI | React/Mission Control | ⏳ TO BUILD |

### **Plumber Tasks:**

1. **Add DDR_bin to patient state** - Store DDR_bin with each ctDNA sample
2. **Build trend analysis** - Calculate slope, detect decline
3. **Implement alert thresholds** - DDR_bin < 0.85 = WARNING, < 0.80 = ALERT
4. **Create dashboard component** - Real-time DDR_bin visualization

---

## 2️⃣ CLINICAL DECISION SUPPORT SYSTEM

### **Intervention Recommendation Engine:**

```python
# Pseudocode for recommendation engine

def recommend_intervention(ddr_bin, mapk_bin, pi3k_bin, current_therapy, cycle):
    """
    Generate tiered intervention recommendations based on pathway bins.
    """
    
    # Calculate resistant clone estimate
    resistant_clone = estimate_resistant_clone(ddr_bin, baseline=0.88)
    
    # Decision logic
    if ddr_bin >= 0.85:
        return {
            "recommendation": "CONTINUE_CURRENT",
            "rationale": "DDR_bin stable, no resistance signal",
            "urgency": "ROUTINE"
        }
    
    elif ddr_bin >= 0.80:
        return {
            "recommendation": "TIER_1_COMBINATION",
            "action": "ADD carboplatin AUC 2-3 every 4 weeks",
            "evidence": "PAOLA-1 trial (combination extends PFS)",
            "rationale": f"Resistant clone ~{resistant_clone}%, suppress early",
            "urgency": "HIGH",
            
            "tier_2": {
                "action": "ENROLL in RAD51 inhibitor trial",
                "evidence": "Preclinical synthetic lethality"
            },
            
            "tier_3": {
                "action": "SWITCH to ATR inhibitor (evolutionary steering)",
                "evidence": "Game theory models (Gatenby lab)"
            }
        }
    
    elif ddr_bin >= 0.70:
        return {
            "recommendation": "TIER_2_ESCALATE",
            "action": "ATR inhibitor OR full-dose platinum",
            "rationale": f"Resistant clone ~{resistant_clone}%, major escalation needed",
            "urgency": "CRITICAL"
        }
    
    else:  # ddr_bin < 0.70
        return {
            "recommendation": "TIER_3_SWITCH",
            "action": "Switch to platinum-based therapy",
            "rationale": "Resistant clone dominant, abandon PARP",
            "urgency": "CRITICAL"
        }
```

### **Output Schema:**

```json
{
  "patient_id": "AK_001",
  "ddr_bin": 0.82,
  "mapk_bin": 0.15,
  "pi3k_bin": 0.22,
  "resistant_clone_estimate": "5-8%",
  "time_to_progression_estimate": "3-6 months",
  
  "recommendations": [
    {
      "tier": 1,
      "action": "ADD carboplatin AUC 2-3 q4weeks",
      "evidence_level": "Level A (RCT)",
      "source": "PAOLA-1 trial",
      "expected_benefit": "+6-12 months PFS",
      "toxicity_risk": "Moderate (Grade 2-3 neutropenia)"
    },
    {
      "tier": 2,
      "action": "ENROLL in CYT-0851 + Olaparib trial",
      "evidence_level": "Level C (Preclinical)",
      "source": "Cyteir Therapeutics Phase 1",
      "expected_benefit": "Potentially indefinite",
      "toxicity_risk": "Unknown"
    },
    {
      "tier": 3,
      "action": "SWITCH to ATR inhibitor",
      "evidence_level": "Level B (Game theory)",
      "source": "Gatenby lab prostate data",
      "expected_benefit": "+12-18 months if cycling",
      "toxicity_risk": "Moderate"
    }
  ],
  
  "not_recommended": [
    {
      "action": "Continue Olaparib monotherapy",
      "reason": "Resistance will progress, DDR_bin declining"
    },
    {
      "action": "Switch to platinum now",
      "reason": "Too early, resistant clone <10%"
    }
  ],
  
  "next_tests": [
    {
      "test": "ctDNA for RAD51C/BRCA1/PALB2",
      "rationale": "Confirm HR restoration mechanism",
      "cost": "$500",
      "turnaround": "3-5 days"
    }
  ],
  
  "provenance": {
    "model_version": "prevention_v1.0",
    "ddr_bin_source": "TRUE_SAE_29_features",
    "recommendation_engine": "rule_based_v1"
  }
}
```

### **Plumber Tasks:**

1. **Create recommendation engine** - `api/services/prevention_recommendation_service.py`
2. **Integrate with complete_care_v2** - Add to response when `include_prevention=true`
3. **Add evidence levels** - Link recommendations to clinical trials
4. **Create UI component** - Display tiered recommendations in Mission Control

---

## 3️⃣ LONGITUDINAL BIOSPECIMEN REPOSITORY

### **Sample Collection Protocol:**

```
MONTH 0 (Baseline):
  - Tumor tissue (FFPE) → Whole exome sequencing
  - Blood (ctDNA) → Baseline mutations
  - Germline DNA → MBD4, BRCA1/2, RAD51C germline status
  → Compute baseline DDR_bin = 0.88

MONTH 3, 6, 9, 12, 15, 18:
  - Blood (ctDNA) → Serial liquid biopsies
  - Compute DDR_bin at each timepoint
  - Track clonal evolution over time

PROGRESSION TIMEPOINT:
  - Tumor biopsy (if accessible) → Compare to baseline
  - Blood (ctDNA) → Resistance mutation confirmation
  - Archive samples for research
```

### **Data Schema:**

```json
{
  "patient_id": "AK_001",
  "biospecimens": [
    {
      "sample_id": "AK_001_T0",
      "timepoint": "baseline",
      "sample_type": "tumor_tissue",
      "collection_date": "2025-01-15",
      "ddr_bin": 0.88,
      "mutations": ["MBD4 c.1239delA", "TP53 p.R175H"]
    },
    {
      "sample_id": "AK_001_B3",
      "timepoint": "month_3",
      "sample_type": "ctdna",
      "collection_date": "2025-04-15",
      "ddr_bin": 0.87,
      "mutations": ["MBD4 c.1239delA (VAF 12%)", "TP53 p.R175H (VAF 8%)"]
    },
    {
      "sample_id": "AK_001_B9",
      "timepoint": "month_9",
      "sample_type": "ctdna",
      "collection_date": "2025-10-15",
      "ddr_bin": 0.82,
      "mutations": ["RAD51C c.905A>G (VAF 2.8%) NEW 🚨"]
    }
  ]
}
```

### **Plumber Tasks:**

1. **Extend patient state** - Add biospecimen tracking
2. **Add ctDNA import** - Parse liquid biopsy results
3. **Track VAF over time** - Detect rising clones
4. **Alert on new mutations** - Flag new resistance mutations

---

## 4️⃣ ADAPTIVE CLINICAL TRIAL PLATFORM

### **DDR_bin-Guided Umbrella Trial Design:**

```
ENROLLMENT: 500 patients starting PARP inhibitor therapy
BASELINE: All patients get DDR_bin computed at Month 0

MONITORING PHASE (Month 0-12):
  Every 3 months: Compute DDR_bin
  
  IF DDR_bin stable (>0.80):
    → Continue current therapy (no intervention)
  
  IF DDR_bin drops to 0.75-0.82: 🚨
    → RANDOMIZE to intervention arm:
       ARM A: Adaptive therapy (dose reduction)
       ARM B: Combination therapy (add carboplatin)
       ARM C: RAD51 inhibitor (add CYT-0851)
       ARM D: Evolutionary steering (switch to ATR inhibitor)
       ARM E: Standard care (continue until clinical progression)
  
  TRACK OUTCOMES:
    Primary: Time to progression (imaging + CA-125)
    Secondary: DDR_bin trajectory, resistant clone burden
    Exploratory: Overall survival, quality of life

EXPECTED RESULTS (2027-2029):
  ARM E (standard): Median PFS 15 months (control)
  ARM A (adaptive): Median PFS 22 months (p<0.01)
  ARM B (combo): Median PFS 26 months (p<0.001)
  ARM C (RAD51i): Median PFS 30 months (p<0.001)
  ARM D (steering): Median PFS 32 months (p<0.001)
```

### **FDA Pathway:**

1. **Companion Diagnostic** - DDR_bin as biomarker for intervention
2. **Breakthrough Device** - Early detection of resistance
3. **De Novo 510(k)** - Novel resistance monitoring device

---

## 📊 IMPLEMENTATION TIMELINE

| Phase | Timeline | Deliverables |
|-------|----------|--------------|
| **Phase 1** | Q1 2025 | DDR_bin dashboard MVP, recommendation engine v1 |
| **Phase 2** | Q2 2025 | Longitudinal tracking, ctDNA import, Mission Control integration |
| **Phase 3** | Q3-Q4 2025 | Clinical trial design, FDA pre-submission |
| **Phase 4** | 2026-2027 | Clinical trial launch, validation |
| **Phase 5** | 2028-2029 | FDA clearance, clinical deployment |

---

## 🛠️ PLUMBER IMMEDIATE TASKS

### **Week 1-2:**
1. Add DDR_bin to patient state schema
2. Create `prevention_recommendation_service.py`
3. Add DDR_bin trend analysis

### **Week 3-4:**
1. Build dashboard component (DDR_bin over time)
2. Implement alert thresholds
3. Integrate with Mission Control

### **Month 2:**
1. Create ctDNA import pipeline
2. Add biospecimen tracking
3. Implement recommendation tiers

---

*Document Owner: Zo*  
*Last Updated: December 24, 2025*

