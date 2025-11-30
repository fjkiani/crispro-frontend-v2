# Resistance Prophet: Service Architecture

## 📋 **PHASE 1: RESISTANCE PROPHET PREDICTION ENGINE (WEEK 1-2)**

### **Objective:** Build CA-125 kinetics-based resistance prediction (MVP)

### **Why CA-125 First (Manager-Approved):**
- ✅ **Service exists:** `ca125_intelligence.py` already operational
- ✅ **Data available:** CA-125 tracking already implemented
- ✅ **Proven signal:** GOG-218 shows CA-125 kinetics predict progression 3-6 weeks early
- ✅ **Fast to ship:** 2 weeks vs 4+ weeks for full multi-signal

---

## **Architecture: Resistance Prophet Service**

**File:** `oncology-backend-minimal/api/services/resistance_prophet_service.py`

### **Core Capabilities:**
1. CA-125 kinetics analysis (trend detection)
2. DNA repair capacity monitoring (SAE mechanism vector)
3. Pathway activation shift detection (escape mechanisms)
4. Multi-signal integration (2-of-3 confidence rule)

---

## **Service Dependencies:**

```python
"""
Resistance Prophet Service - Predicts treatment resistance 3-6 months early.

Methodology:
- CA-125 kinetics (rising trend detection)
- SAE mechanism vector (pathway activation shifts)
- Treatment line context (cross-resistance patterns)
- 2-of-3 signal detection (high confidence threshold)

Integration:
- Uses existing CA125Intelligence service
- Uses existing SAE feature extraction
- Uses existing TreatmentLineIntelligence framework
"""

from typing import Dict, List, Optional
from datetime import datetime, timedelta
import numpy as np

class ResistanceProphetService:
    """
    Predicts treatment resistance before clinical progression.
    """
    
    def __init__(self, ca125_service, sae_service, treatment_line_service):
        self.ca125_service = ca125_service
        self.sae_service = sae_service
        self.treatment_line_service = treatment_line_service
```

---

## **Main API Contract:**

### **Input:**
```python
{
    "patient_id": "PAT12345",
    "current_therapy": "carboplatin_paclitaxel_bevacizumab",
    "treatment_start_date": "2024-07-01",
    "ca125_history": [
        {"date": "2024-07-01", "value": 2842},
        {"date": "2024-07-22", "value": 1200},
        {"date": "2024-08-12", "value": 450},
        {"date": "2024-09-02", "value": 600}  # Rising (resistance signal)
    ],
    "current_mutations": [
        {"gene": "BRCA1", "hgvs_p": "C64R", "chrom": "17", "pos": 43044295}
    ],
    "baseline_sae_features": {
        "dna_repair_capacity": 0.30,
        "mechanism_vector": [0.70, 0.30, 0.20, 0.15, 0.10, 0.25, 0.15]
    }
}
```

### **Output:**
```python
{
    "resistance_probability": 0.75,  # 0-1 scale
    "confidence": "HIGH",  # HIGH/MEDIUM/LOW
    "predicted_progression_date": "2025-04-15",
    "signals_detected": [
        {
            "type": "CA125_RISE",
            "confidence": 0.85,
            "mechanism": "On-therapy CA-125 rise detected",
            "evidence": "Cycle 4: 600 U/mL (↑33% from Cycle 3: 450 U/mL)"
        },
        {
            "type": "DNA_REPAIR_RESTORATION",
            "confidence": 0.70,
            "mechanism": "DNA repair capacity increased",
            "evidence": "Baseline: 0.30 → Current: 0.52 (↑73%)"
        }
    ],
    "recommendations": [
        {
            "action": "ESCALATE_IMAGING",
            "urgency": "IMMEDIATE",
            "rationale": "CA-125 rise suggests progression 3-6 weeks before imaging"
        },
        {
            "action": "CONSIDER_SWITCH",
            "urgency": "HIGH",
            "rationale": "2-of-3 resistance signals detected (confidence >70%)"
        },
        {
            "action": "REVIEW_RESISTANCE_PLAYBOOK",
            "urgency": "HIGH",
            "rationale": "Prepare next-line therapy (Niraparib + Bevacizumab 96.6%)"
        }
    ],
    "next_line_options": [
        {
            "therapy": "Niraparib + Bevacizumab",
            "rank": 0.966,
            "rationale": "Dual mechanism (PARP + anti-angio) overcomes HR restoration"
        }
    ]
}
```

---

## **Signal Detection Framework:**

### **Signal 1: CA-125 Kinetics**
**Confidence Levels:**
- On-therapy rise (ANY cycle) → 85% confidence
- <50% drop by Cycle 3 → 65% confidence
- <70% drop by Cycle 6 → 70% confidence

**Evidence Base:**
- GOG-218
- ICON7
- SOLO-1 post-hoc analyses

---

### **Signal 2: DNA Repair Restoration**
**Confidence Levels:**
- DNA repair capacity increase >0.20 → 75% confidence
- BRCA1/2 reversion mutations → 80% confidence
- RAD51C/D restoration → 70% confidence

**Evidence Base:**
- PMID:28588062
- PMID:27775701

**Relevance:**
- Only for PARP/platinum therapies
- Indicates homologous recombination repair (HRR) restoration

---

### **Signal 3: Pathway Activation Shift**
**Confidence Levels:**
- MAPK activation (>0.60) → 65% confidence
- PI3K/AKT activation (>0.60) → 60% confidence
- Efflux pathway (>0.50) → 70% confidence

**Evidence Base:**
- PMID:29945933
- PMID:30952680

**Mechanism:**
- Indicates escape via alternative pathways
- Common in acquired resistance to targeted therapies

---

## **Probability Calculation:**

### **Method:**
Weighted average of signal confidences

### **Rules:**
- 0 signals: 0.30 (low baseline risk)
- 1 signal: Average confidence (MEDIUM risk)
- 2+ signals: Average confidence × 1.1 (HIGH risk if >0.70)

### **Confidence Levels:**
- **HIGH**: probability ≥0.70 AND 2+ signals
- **MEDIUM**: probability 0.50-0.70 OR 1 signal
- **LOW**: probability <0.50

---

## **Progression Timeline Prediction:**

### **Method:**
Based on detected signals

### **Timeline Estimates:**
- CA-125 rise: 3-6 weeks (use midpoint: 4.5 weeks)
- Inadequate response: 8-12 weeks (use midpoint: 10 weeks)
- Pathway-based only: 12-16 weeks (use midpoint: 14 weeks)

---

## **Recommendations Framework:**

### **Priority Levels:**

1. **IMMEDIATE:**
   - Trigger: CA-125 rise detected
   - Action: Escalate imaging (order CT/PET scan within 1 week)

2. **HIGH:**
   - Trigger: 2+ signals OR probability ≥0.70
   - Actions:
     - Consider therapy switch (discuss within 2 weeks)
     - Review resistance playbook (prepare next-line options within 1 week)

3. **MEDIUM:**
   - Trigger: 1 signal
   - Action: Monitor more frequently (weekly CA-125 for next 4 weeks)

---

## **Integration Points:**

### **Existing Services Used:**
1. `CA125Intelligence` - Burden classification, forecast
2. `SAEFeatureService` - Mechanism vector, DNA repair capacity
3. `TreatmentLineIntelligence` - Next-line options, cross-resistance

### **Data Requirements:**
1. CA-125 history (time series)
2. Mutation data (baseline + current)
3. Treatment context (therapy, start date, line)
4. Baseline SAE features

### **Output Integration:**
- Enhanced `/api/ayesha/complete_care_v2` endpoint
- Mission Control dashboard (Critical Alerts, Resistance Playbook, Next Best Actions)

