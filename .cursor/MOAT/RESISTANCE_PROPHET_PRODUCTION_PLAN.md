# 🛠️ RESISTANCE PROPHET PRODUCTION PLAN

**Version:** 1.0  
**Date:** January 28, 2025  
**Status:** EXECUTION READY  
**Source of Truth:** `RESISTANCE_PROPHET_PRODUCTION_AUDIT.md`

---

## 🎯 PRODUCTION DEFINITION (No CA-125 Required)

### What We Ship Now (Baseline-Only Mode)

| **Capability** | **Evidence** | **Status** |
|---------------|--------------|------------|
| NF1 mutation → platinum resistance (OV) | RR=2.10, TCGA-OV validated | ✅ READY |
| DIS3 mutation → mortality risk (MM) | RR=2.08, p=0.0145, MMRF validated | ✅ READY |
| Treatment alternatives + monitoring | Rule-based playbook | ✅ READY |
| Input completeness caps (L0/L1/L2) | Confidence limiting | ✅ READY |
| Evidence tier tracking | Provenance in response | ✅ READY |

### What We DON'T Ship (Until CA-125 Kinetics)

| **Claim** | **Issue** | **Status** |
|-----------|-----------|------------|
| "Predicts resistance 3-6 months early" | No longitudinal validation | ❌ REMOVE |
| "AUROC ≥ 0.70" | Actual: 0.464 | ❌ REMOVE |
| "MAPK RR = 1.97" | Hard-coded, not computed | ❌ REVALIDATE |
| Resistance probability score | Zero sensitivity | ❌ RUO ONLY |

---

## 📋 PHASE 1: REMOVE HARD-CODED CLAIMS (1 day)

### Task 1.1: Downgrade MAPK/KRAS Evidence Level

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/resistance_playbook_service.py`

**Location:** Line ~476

**Current:**
```python
"KRAS": {
    "relative_risk": 1.97,  # Part of MAPK pathway
    "p_value": 0.05,
    "evidence_level": EvidenceLevel.VALIDATED,
    ...
}
```

**Change To:**
```python
"KRAS": {
    "relative_risk": None,  # PENDING REVALIDATION - was hard-coded
    "p_value": None,
    "evidence_level": EvidenceLevel.LITERATURE_BASED,  # Downgrade until validated
    "validation_source": "PENDING_REVALIDATION",
    ...
}
```

**Acceptance Criteria:**
- [ ] `relative_risk` set to `None` for KRAS, NRAS, BRAF (MAPK genes)
- [ ] `evidence_level` changed to `LITERATURE_BASED`
- [ ] `validation_source` changed to `"PENDING_REVALIDATION"`
- [ ] Unit test confirms no numeric RR returned for MAPK genes

### Task 1.2: Add Baseline-Only Mode Flag

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/resistance_prophet_service.py`

**Add to response:**
```python
# In predict_resistance() response
{
    "mode": "baseline_only" if not has_ca125_kinetics else "kinetics_enabled",
    "mode_disclaimer": "Baseline risk stratification only. Early detection requires CA-125 kinetics.",
    "claims_disabled": ["early_detection", "lead_time_analysis"] if mode == "baseline_only" else [],
    ...
}
```

**Acceptance Criteria:**
- [ ] Response includes `mode` field
- [ ] Response includes `mode_disclaimer` when baseline-only
- [ ] Response includes `claims_disabled` array

### Task 1.3: Add RUO Disclaimer to All Resistance Outputs

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/resistance_prophet_service.py`

**Add:**
```python
RUO_DISCLAIMER = (
    "Research Use Only (RUO). Resistance risk stratification based on baseline genetics. "
    "Not validated for clinical decision-making without CA-125 kinetics. "
    "Validated markers: NF1 (OV, RR=2.10), DIS3 (MM, RR=2.08, p=0.0145)."
)

# Add to every response
response["ruo_disclaimer"] = RUO_DISCLAIMER
```

**Acceptance Criteria:**
- [ ] Every `/api/resistance/predict` response includes `ruo_disclaimer`
- [ ] Every `/api/care/resistance_playbook_v2` response includes `ruo_disclaimer`

---

## 📋 PHASE 2: RECEIPT-FIRST REVALIDATION (2 days)

### Task 2.1: Create MAPK Validation Script

**File:** `oncology-coPilot/oncology-backend-minimal/scripts/validation/validate_mapk_ov_platinum.py`

**Content:**
```python
#!/usr/bin/env python3
"""
MAPK Pathway Validation for OV Platinum Resistance
Computes actual RR from TCGA-OV cohort (not hard-coded)
"""

import json
import numpy as np
from scipy import stats
from pathlib import Path

# MAPK gene set (from audit)
MAPK_GENES = {"KRAS", "NRAS", "BRAF", "NF1", "MAP2K1", "MAP2K2"}

def load_tcga_ov_cohort(cohort_path: str) -> list:
    """Load TCGA-OV cohort with platinum response labels"""
    with open(cohort_path) as f:
        return json.load(f)

def has_mapk_mutation(patient: dict) -> bool:
    """Check if patient has any MAPK pathway mutation"""
    mutations = patient.get("mutations", [])
    patient_genes = {m.get("gene", "").upper() for m in mutations}
    return bool(patient_genes & MAPK_GENES)

def compute_contingency_table(cohort: list) -> dict:
    """Compute 2x2 contingency table"""
    mapk_resistant = 0
    mapk_sensitive = 0
    wt_resistant = 0
    wt_sensitive = 0
    
    for patient in cohort:
        is_resistant = patient.get("platinum_response") in ["resistant", "refractory"]
        is_mapk = has_mapk_mutation(patient)
        
        if is_mapk and is_resistant:
            mapk_resistant += 1
        elif is_mapk and not is_resistant:
            mapk_sensitive += 1
        elif not is_mapk and is_resistant:
            wt_resistant += 1
        else:
            wt_sensitive += 1
    
    return {
        "mapk_resistant": mapk_resistant,
        "mapk_sensitive": mapk_sensitive,
        "wt_resistant": wt_resistant,
        "wt_sensitive": wt_sensitive
    }

def compute_relative_risk(table: dict) -> dict:
    """Compute RR with 95% CI"""
    a = table["mapk_resistant"]
    b = table["mapk_sensitive"]
    c = table["wt_resistant"]
    d = table["wt_sensitive"]
    
    n_mapk = a + b
    n_wt = c + d
    
    if n_mapk == 0 or n_wt == 0:
        return {"error": "Insufficient data"}
    
    risk_mapk = a / n_mapk
    risk_wt = c / n_wt
    
    if risk_wt == 0:
        return {"error": "Zero risk in wildtype"}
    
    rr = risk_mapk / risk_wt
    
    # 95% CI using log method
    log_rr = np.log(rr)
    se_log_rr = np.sqrt((1/a - 1/n_mapk) + (1/c - 1/n_wt)) if a > 0 and c > 0 else np.inf
    ci_lower = np.exp(log_rr - 1.96 * se_log_rr)
    ci_upper = np.exp(log_rr + 1.96 * se_log_rr)
    
    # Chi-square test
    contingency = [[a, b], [c, d]]
    chi2, p_value, dof, expected = stats.chi2_contingency(contingency)
    
    return {
        "relative_risk": round(rr, 3),
        "ci_lower": round(ci_lower, 3),
        "ci_upper": round(ci_upper, 3),
        "p_value": round(p_value, 4),
        "chi2": round(chi2, 3),
        "n_mapk": n_mapk,
        "n_wt": n_wt,
        "risk_mapk": round(risk_mapk, 4),
        "risk_wt": round(risk_wt, 4),
        "significant": p_value < 0.05
    }

def main():
    cohort_path = "data/validation/tcga_ov_platinum_response_with_genomics.json"
    output_dir = Path("scripts/validation/out/mapk_ov_platinum")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    cohort = load_tcga_ov_cohort(cohort_path)
    table = compute_contingency_table(cohort)
    result = compute_relative_risk(table)
    
    report = {
        "endpoint_contract": "ov_platinum_response_tcga_style",
        "gene_set": list(MAPK_GENES),
        "contingency_table": table,
        "result": result,
        "acceptance": {
            "target_rr": 1.5,
            "target_p": 0.05,
            "passed": result.get("significant", False) and result.get("relative_risk", 0) >= 1.5
        },
        "action": "COHORT_VALIDATED" if result.get("significant") else "LITERATURE_BASED"
    }
    
    with open(output_dir / "report.json", "w") as f:
        json.dump(report, f, indent=2)
    
    print(f"MAPK Validation Report:")
    print(f"  RR: {result.get('relative_risk')} [{result.get('ci_lower')}, {result.get('ci_upper')}]")
    print(f"  p-value: {result.get('p_value')}")
    print(f"  Significant: {result.get('significant')}")
    print(f"  Action: {report['action']}")

if __name__ == "__main__":
    main()
```

**Acceptance Criteria:**
- [ ] Script runs without errors
- [ ] Produces `report.json` with actual RR computation
- [ ] RR is computed from data, not hard-coded
- [ ] If RR < 1.5 or p ≥ 0.05, action = "LITERATURE_BASED"

### Task 2.2: Create PI3K Validation Script

**File:** `oncology-coPilot/oncology-backend-minimal/scripts/validation/validate_pi3k_ov_platinum.py`

**Same pattern as MAPK, with:**
```python
PI3K_GENES = {"PIK3CA", "PIK3CB", "PIK3R1", "AKT1", "AKT2", "PTEN"}
```

**Acceptance Criteria:**
- [ ] Script runs without errors
- [ ] Produces `pi3k_ov_report.json`
- [ ] If significant, add PI3K to playbook code
- [ ] If not significant, remove PI3K claims from docs

### Task 2.3: Confirm NF1 Validation (Re-run Existing)

**File:** `oncology-coPilot/oncology-backend-minimal/scripts/validation/validate_ov_nf1_playbook.py`

**Command:**
```bash
python3 scripts/validation/validate_ov_nf1_playbook.py
```

**Acceptance Criteria:**
- [ ] Report confirms RR ≥ 2.0
- [ ] Report confirms p < 0.05
- [ ] Output matches `resistance_playbook_service.py:427`

### Task 2.4: Update Playbook With Validated RRs

**After Phase 2.1-2.3 complete:**

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/resistance_playbook_service.py`

**Actions:**
1. If MAPK validated: Update RR from report.json
2. If MAPK not validated: Keep `relative_risk: None, evidence_level: LITERATURE_BASED`
3. If PI3K validated: Add to OV_RESISTANCE_MARKERS
4. If PI3K not validated: Remove PI3K claims from docs

---

## 📋 PHASE 3: PRODUCTION-READY UX (1 day)

### Task 3.1: Update ResistancePanel.jsx Copy

**File:** `oncology-coPilot/oncology-frontend/src/components/resistance/ResistancePanel.jsx`

**Changes:**

1. **Remove early detection claim:**
```jsx
// OLD
<h2>Resistance Prediction (3-6 months early)</h2>

// NEW
<h2>Baseline Resistance Risk Stratification</h2>
<span className="ruo-badge">Research Use Only</span>
```

2. **Add evidence tier badges:**
```jsx
const EvidenceBadge = ({ level }) => {
  const badges = {
    COHORT_VALIDATED: { color: "green", label: "Validated" },
    COHORT_TREND: { color: "yellow", label: "Trend" },
    LITERATURE_BASED: { color: "gray", label: "Literature" },
    PENDING_REVALIDATION: { color: "red", label: "Pending" }
  };
  const badge = badges[level] || badges.LITERATURE_BASED;
  return <span className={`evidence-badge ${badge.color}`}>{badge.label}</span>;
};
```

3. **Add missing kinetics banner:**
```jsx
{mode === "baseline_only" && (
  <div className="kinetics-banner warning">
    <AlertTriangle size={16} />
    <span>CA-125 kinetics not provided. Early detection layer disabled.</span>
  </div>
)}
```

**Acceptance Criteria:**
- [ ] No "3-6 months early" language anywhere
- [ ] RUO badge displayed
- [ ] Evidence tier badges per signal
- [ ] Missing kinetics banner when baseline-only

### Task 3.2: Update API Response Contract

**File:** `oncology-coPilot/oncology-backend-minimal/api/contracts/resistance_contract.py`

**Add fields:**
```python
class ResistanceContract(BaseModel):
    # Existing fields...
    
    # New baseline-only fields
    mode: Literal["baseline_only", "kinetics_enabled"] = "baseline_only"
    mode_disclaimer: str = Field(
        default="Baseline risk stratification only. Early detection requires CA-125 kinetics."
    )
    claims_disabled: List[str] = Field(default_factory=list)
    ruo_disclaimer: str = Field(
        default="Research Use Only. Not validated for clinical decision-making without CA-125 kinetics."
    )
    
    # Evidence summary
    evidence_summary: Dict[str, int] = Field(
        default_factory=lambda: {"COHORT_VALIDATED": 0, "COHORT_TREND": 0, "LITERATURE_BASED": 0}
    )
```

**Acceptance Criteria:**
- [ ] Contract includes all new fields
- [ ] Default mode is `baseline_only`
- [ ] Evidence summary counts signals by tier

---

## 📋 PHASE 4: OPERATIONAL HARDENING (1 day)

### Task 4.1: Add Ring-1 Validation Gate

**File:** `oncology-coPilot/oncology-backend-minimal/scripts/validation/run_resistance_validation_suite.py`

**Add checks:**
```python
def check_no_hard_coded_rr():
    """Ensure no MAPK RR is surfaced without validation receipt"""
    playbook_path = "api/services/resistance_playbook_service.py"
    
    # Parse playbook and check KRAS/NRAS/BRAF entries
    # If relative_risk is numeric and evidence_level != COHORT_VALIDATED
    # AND no report.json exists → FAIL
    
    mapk_genes = ["KRAS", "NRAS", "BRAF"]
    for gene in mapk_genes:
        report_path = f"scripts/validation/out/{gene.lower()}_ov_platinum/report.json"
        if not Path(report_path).exists():
            # Check if RR is numeric in playbook
            # If yes → FAIL
            pass
    
    return True  # or False

def check_baseline_mode_enforced():
    """Ensure baseline-only mode never claims early detection"""
    # Make test request with no CA-125
    # Assert mode == "baseline_only"
    # Assert "early_detection" in claims_disabled
    pass
```

**Acceptance Criteria:**
- [ ] Ring-1 fails if hard-coded RR without receipt
- [ ] Ring-1 fails if baseline-only mode claims early detection

### Task 4.2: Add Deterministic Fixtures

**File:** `oncology-coPilot/oncology-backend-minimal/scripts/validation/fixtures/resistance_fixtures.json`

**Content:**
```json
{
  "fixtures": [
    {
      "id": "nf1_positive_ov",
      "disease": "ovarian",
      "mutations": [{"gene": "NF1", "hgvs_p": "p.R1276*"}],
      "expected_risk_level": "HIGH",
      "expected_evidence": "COHORT_VALIDATED"
    },
    {
      "id": "dis3_positive_mm",
      "disease": "myeloma",
      "mutations": [{"gene": "DIS3", "hgvs_p": "p.C562Y"}],
      "expected_risk_level": "HIGH",
      "expected_evidence": "COHORT_VALIDATED"
    },
    {
      "id": "kras_pending_ov",
      "disease": "ovarian",
      "mutations": [{"gene": "KRAS", "hgvs_p": "p.G12D"}],
      "expected_risk_level": "MEDIUM",
      "expected_evidence": "LITERATURE_BASED"
    },
    {
      "id": "no_signal_baseline",
      "disease": "ovarian",
      "mutations": [],
      "expected_risk_level": "LOW",
      "expected_mode": "baseline_only"
    }
  ]
}
```

**Acceptance Criteria:**
- [ ] Fixtures exist for validated markers (NF1, DIS3)
- [ ] Fixtures exist for pending markers (KRAS)
- [ ] Fixtures exist for no-signal baseline

### Task 4.3: Add Observability Logging

**File:** `oncology-coPilot/oncology-backend-minimal/api/services/resistance_prophet_service.py`

**Add:**
```python
import structlog
logger = structlog.get_logger()

# In predict_resistance()
logger.info(
    "resistance_prediction",
    mode=response["mode"],
    input_completeness=input_completeness_level,
    evidence_tiers=response["evidence_summary"],
    signals_triggered=[s["type"] for s in response.get("signals_detected", [])],
    patient_id=request.patient_id
)
```

**Acceptance Criteria:**
- [ ] Every prediction logs mode + evidence tiers
- [ ] Logs are structured (JSON-parseable)

---

## 📋 PHASE 5: DOCUMENTATION CLEANUP (1 day)

### Task 5.1: Update 02_RESISTANCE_PREDICTION.md

**File:** `.cursor/MOAT/ADVANCED_CARE_PLAN/02_RESISTANCE_PREDICTION.md`

**Changes:**
1. Remove MAPK row from "WHAT WE ACTUALLY VALIDATED" table (or mark as PENDING)
2. Add note: "MAPK RR was identified as hard-coded; revalidation in progress"
3. Update status table to match `RESISTANCE_PROPHET_PRODUCTION_AUDIT.md`

### Task 5.2: Update AYESHA Resistance Prediction Doc

**File:** `.cursor/ayesha/ADVANCED_CARE_PLAN_RESISTANCE_PREDICTION.md`

**Changes:**
1. Remove all "3-6 months early" claims
2. Remove MAPK RR claims until validated
3. Add RUO disclaimer section
4. Update "What We Can Say" to match audit

### Task 5.3: Update RESISTANCE_VALIDATION_PLAN.md

**File:** `.cursor/MOAT/RESISTANCE_VALIDATION_PLAN.md`

**Add section:**
```markdown
## Baseline-Only Production Spec

### What's Shipped (Mode: baseline_only)
- NF1 OV (COHORT_VALIDATED)
- DIS3 MM (COHORT_VALIDATED)
- Playbook alternatives (rule-based)
- Input completeness caps
- Evidence tier tracking

### What's NOT Shipped (Until kinetics_enabled)
- Resistance probability score
- Early detection claims
- MAPK RR (pending revalidation)
- Lead-time analysis
```

### Task 5.4: Create Single Source of Truth Link

**All docs must reference:**
```markdown
**Source of Truth:** `.cursor/MOAT/RESISTANCE_PROPHET_PRODUCTION_AUDIT.md`
```

---

## 📋 PHASE 6: FINAL CHECKLIST

### Pre-Production Gate

| **Check** | **Verification** | **Owner** |
|-----------|------------------|-----------|
| MAPK RR downgraded | `grep "KRAS.*LITERATURE" resistance_playbook_service.py` | Plumber |
| Baseline mode added | Response includes `mode: "baseline_only"` | Plumber |
| RUO disclaimer added | Response includes `ruo_disclaimer` | Plumber |
| MAPK validator exists | `ls scripts/validation/validate_mapk_ov_platinum.py` | Plumber |
| Ring-1 passes | `python3 run_resistance_validation_suite.py --ring1` | CI |
| Frontend updated | No "3-6 months" text | Plumber |
| Docs aligned | All docs reference audit | Manager |

### Post-Validation Gate (After MAPK/PI3K revalidation)

| **Check** | **Verification** | **Owner** |
|-----------|------------------|-----------|
| MAPK report exists | `cat scripts/validation/out/mapk_ov_platinum/report.json` | Plumber |
| PI3K report exists | `cat scripts/validation/out/pi3k_ov_platinum/report.json` | Plumber |
| Playbook updated | RRs match reports | Plumber |
| Evidence levels correct | All markers have correct tier | Plumber |

---

## 🎯 DELIVERABLE SUMMARY

| **Phase** | **Duration** | **Key Deliverables** |
|-----------|--------------|----------------------|
| Phase 1 | 1 day | MAPK downgraded, baseline mode added, RUO disclaimer |
| Phase 2 | 2 days | MAPK/PI3K validators, NF1 confirmed, playbook updated |
| Phase 3 | 1 day | Frontend copy updated, evidence badges, kinetics banner |
| Phase 4 | 1 day | Ring-1 gates, fixtures, observability |
| Phase 5 | 1 day | All docs aligned to audit |

**Total: 6 days to production-ready (without CA-125)**

---

## 🔥 HONEST CLAIMS (What We Can Ship)

> "Our Resistance Playbook provides **evidence-tiered baseline resistance risk stratification** for Ovarian Cancer (NF1 RR=2.10, COHORT_VALIDATED) and Multiple Myeloma (DIS3 RR=2.08, p=0.0145, COHORT_VALIDATED), with **treatment alternatives**, **monitoring recommendations**, and **confidence caps** based on input completeness. Research Use Only until CA-125 kinetics validation complete."

---

**⚔️ PLAN COMPLETE. EXECUTE PHASES 1-5 IN ORDER. ⚔️**

