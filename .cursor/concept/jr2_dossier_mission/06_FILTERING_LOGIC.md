# ⚔️ FILTERING LOGIC - REPLICATE ZO'S "1 IN 700" ⚔️

**Purpose**: Filter 50 candidates to find top 5-10 trials (like Zo found 1 in 700)

**Reference**: `oncology-backend-minimal/api/routers/ayesha_trials.py` lines 180-250

---

## 📋 **FILTERING FUNCTION**

```python
def filter_50_candidates(trials: List[Dict], patient: Dict) -> Dict:
    """
    Filter 50 candidates to find top 5-10 (like Zo found 1 in 700).
    
    Returns:
        {
            'top_tier': [...],      # 5-10 trials (pass ALL filters)
            'good_tier': [...],     # 10-15 trials (pass MOST filters)
            'ok_tier': [...],       # 15-20 trials (interesting but not immediate)
            'rejected': [...]       # 10-20 trials (not eligible)
        }
    """
    top_tier = []
    good_tier = []
    ok_tier = []
    rejected = []
    
    for trial in trials:
        # FILTER 1: Stage IV allowed
        stage_match = check_stage_match(trial, patient)  # True if Stage IV allowed
        
        # FILTER 2: Treatment line (first-line OR maintenance)
        line_match = check_treatment_line(trial, patient)  # True if first-line or maintenance
        
        # FILTER 3: Recruiting status
        recruiting = trial['status'] in ['RECRUITING', 'NOT_YET_RECRUITING']
        
        # FILTER 4: Geography (USA)
        usa_location = check_usa_location(trial)  # True if any USA site
        
        # FILTER 5: Biomarker gates (critical gates)
        biomarker_gates = check_biomarker_requirements(trial, patient)
        # Returns: {'her2': 'PENDING', 'brca': 'PASS', 'hrd': 'PENDING'}
        
        # Scoring
        if stage_match and line_match and recruiting and usa_location:
            # Check biomarker gates
            pending_gates = [k for k, v in biomarker_gates.items() if v == 'PENDING']
            failed_gates = [k for k, v in biomarker_gates.items() if v == 'FAIL']
            
            if failed_gates:
                rejected.append({**trial, 'reason': f"Biomarker gate failed: {failed_gates}"})
            elif len(pending_gates) == 0:
                top_tier.append(trial)  # Perfect match!
            elif len(pending_gates) <= 2:
                good_tier.append(trial)  # Good match, pending tests
            else:
                ok_tier.append(trial)  # Too many unknowns
        elif stage_match and usa_location:
            # Partial match (wrong line or not recruiting)
            ok_tier.append(trial)
        else:
            rejected.append({**trial, 'reason': 'Stage or geography mismatch'})
    
    return {
        'top_tier': sorted(top_tier, key=lambda t: t.get('match_score', 0), reverse=True),
        'good_tier': sorted(good_tier, key=lambda t: t.get('match_score', 0), reverse=True),
        'ok_tier': sorted(ok_tier, key=lambda t: t.get('match_score', 0), reverse=True),
        'rejected': rejected
    }
```

---

## 🎯 **FILTER CRITERIA**

### **TOP-TIER (Priority 1)** - Must meet ALL:
- ✅ Stage IV allowed
- ✅ First-line OR Maintenance (post-frontline)
- ✅ Recruiting OR Not Yet Recruiting (start date < 6 months)
- ✅ USA location
- ✅ 0-2 pending biomarker gates

### **GOOD-TIER (Priority 2)** - Must meet MOST:
- ✅ Stage III/IV allowed
- ⚠️ Maintenance OR Recurrent platinum-sensitive
- ✅ Recruiting (or recently completed with extension)
- ✅ USA or nearby (Canada/Mexico)
- ⚠️ 2-3 pending biomarker gates

### **OK-TIER (Priority 3)** - Interesting but not immediate:
- ⚠️ Stage III only
- ⚠️ Recurrent platinum-resistant
- ⚠️ Not yet recruiting (start date > 6 months)
- ⚠️ International only
- ⚠️ 3+ pending biomarker gates

---

## 🔍 **HELPER FUNCTIONS**

```python
def check_stage_match(trial: Dict, patient: Dict) -> bool:
    """Check if trial allows patient's stage."""
    eligibility = trial.get('eligibility_text', '').lower()
    patient_stage = patient.get('stage', '').upper()
    
    # Check for Stage IV
    if 'stage iv' in eligibility or 'stage 4' in eligibility:
        return True
    if 'stage iii or iv' in eligibility or 'stage 3 or 4' in eligibility:
        return True
    if patient_stage.startswith('IV'):
        return 'stage iv' in eligibility or 'stage 4' in eligibility
    
    return False

def check_treatment_line(trial: Dict, patient: Dict) -> bool:
    """Check if trial matches patient's treatment line."""
    eligibility = trial.get('eligibility_text', '').lower()
    description = trial.get('description_text', '').lower()
    combined = eligibility + ' ' + description
    
    patient_line = patient.get('treatment_line', '').lower()
    
    # First-line keywords
    if patient_line == 'first-line':
        return any(kw in combined for kw in ['first-line', 'frontline', 'newly diagnosed', 'treatment-naive'])
    
    # Maintenance keywords
    if patient_line == 'maintenance':
        return 'maintenance' in combined
    
    return False

def check_usa_location(trial: Dict) -> bool:
    """Check if trial has any USA location."""
    locations = trial.get('locations_data', [])
    return any(loc.get('country', '').lower() == 'united states' for loc in locations)

def check_biomarker_requirements(trial: Dict, patient: Dict) -> Dict:
    """Check biomarker gates (HER2, BRCA, HRD, etc.)."""
    eligibility = trial.get('eligibility_text', '').lower()
    gates = {}
    
    # HER2 gate
    if 'her2' in eligibility:
        patient_her2 = patient.get('her2_status', 'UNKNOWN')
        if patient_her2 == 'UNKNOWN':
            gates['her2'] = 'PENDING'
        elif 'her2-expressing' in eligibility or 'her2 ihc 3+/2+/1+' in eligibility:
            if patient_her2 in ['IHC_3+', 'IHC_2+', 'IHC_1+']:
                gates['her2'] = 'PASS'
            else:
                gates['her2'] = 'FAIL'
    
    # BRCA gate
    if 'brca' in eligibility:
        patient_brca = patient.get('brca_status', 'UNKNOWN')
        if patient_brca == 'UNKNOWN':
            gates['brca'] = 'PENDING'
        elif 'brca wildtype' in eligibility or 'without brca' in eligibility:
            if patient_brca == 'WILDTYPE':
                gates['brca'] = 'PASS'
            else:
                gates['brca'] = 'FAIL'
    
    # HRD gate
    if 'hrd' in eligibility or 'homologous recombination' in eligibility:
        patient_hrd = patient.get('hrd_score', None)
        if patient_hrd is None:
            gates['hrd'] = 'PENDING'
        elif 'hrd-positive' in eligibility or 'hrd-high' in eligibility:
            if patient_hrd >= 42:
                gates['hrd'] = 'PASS'
            else:
                gates['hrd'] = 'FAIL'
    
    return gates
```

---

**Your Goal**: Find 5-10 "top_tier" trials from Zo's 50 candidates (like Zo found 1 in 700)

