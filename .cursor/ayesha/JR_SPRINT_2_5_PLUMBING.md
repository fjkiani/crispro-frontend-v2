# 🔧 JR AGENT PLUMBING TASKS - Sprint 2.5

**Assigned to:** JR Agent  
**Tracked in:** `.cursor/ayesha/JR_AGENT_TASKS.md` (Task 7)  
**Parent Sprint:** `.cursor/MOAT/SPRINT_2_5_SAE_RESISTANCE_VALIDATION.md`  
**Timeline:** Days 1-3  
**Status:** 🚀 **READY TO START**

---

## 🎯 Goal

Create fixture files, validator scaffolds, and Ring-1 integration for Sprint 2.5 SAE resistance validation.

**Zo implements validation logic. JR does plumbing.**

---

## Task 7.1: Create Fixture Files (Day 1)

**Deliverable:** 4 JSON fixture files in `scripts/validation/fixtures/sae_resistance/`

### Commands:
```bash
# Create directory
mkdir -p oncology-coPilot/oncology-backend-minimal/scripts/validation/fixtures/sae_resistance

# Create fixture files
cd oncology-coPilot/oncology-backend-minimal/scripts/validation/fixtures/sae_resistance
touch hr_restoration_fixtures.json
touch dna_repair_formula_fixtures.json
touch ca125_inadequate_fixtures.json
touch 2_of_3_triggers_fixtures.json
```

### Fixture 1: `hr_restoration_fixtures.json`

```json
[
  {
    "name": "BRCA1_reversion_hr_restoration",
    "patient_id": "FIXTURE_HR_RESTORATION",
    "baseline": {
      "hrd_score": 58,
      "dna_repair_capacity": 0.75,
      "pathway_burden_ddr": 0.70
    },
    "current": {
      "hrd_score": 43,
      "dna_repair_capacity": 0.50,
      "pathway_burden_ddr": 0.45
    },
    "expected": {
      "hr_restoration_detected": true,
      "repair_change": -0.25,
      "risk_level": "HIGH"
    }
  },
  {
    "name": "stable_disease_no_restoration",
    "patient_id": "FIXTURE_STABLE",
    "baseline": {
      "hrd_score": 52,
      "dna_repair_capacity": 0.70,
      "pathway_burden_ddr": 0.65
    },
    "current": {
      "hrd_score": 50,
      "dna_repair_capacity": 0.68,
      "pathway_burden_ddr": 0.63
    },
    "expected": {
      "hr_restoration_detected": false,
      "repair_change": -0.02,
      "risk_level": "LOW"
    }
  },
  {
    "name": "treatment_naive_no_baseline",
    "patient_id": "FIXTURE_TREATMENT_NAIVE",
    "baseline": null,
    "current": {
      "dna_repair_capacity": 0.60
    },
    "expected": {
      "hr_restoration_detected": false,
      "reason": "no_baseline",
      "false_positive": false
    }
  }
]
```

### Fixture 2: `dna_repair_formula_fixtures.json`

```json
[
  {
    "name": "BRCA1_biallelic",
    "patient_id": "FIXTURE_BRCA1_BIALLELIC",
    "inputs": {
      "pathway_burden_ddr": 0.70,
      "essentiality_hrr": 0.85,
      "exon_disruption": 0.90
    },
    "expected": {
      "dna_repair_capacity": 0.82,
      "formula": "0.6×0.70 + 0.2×0.85 + 0.2×0.90 = 0.82"
    }
  },
  {
    "name": "MBD4_frameshift",
    "patient_id": "FIXTURE_MBD4_FRAMESHIFT",
    "inputs": {
      "pathway_burden_ddr": 1.0,
      "essentiality_hrr": 0.0,
      "exon_disruption": 0.0
    },
    "expected": {
      "dna_repair_capacity": 0.60,
      "formula": "0.6×1.0 + 0.2×0.0 + 0.2×0.0 = 0.60"
    }
  },
  {
    "name": "wild_type",
    "patient_id": "FIXTURE_WILD_TYPE",
    "inputs": {
      "pathway_burden_ddr": 0.0,
      "essentiality_hrr": 0.0,
      "exon_disruption": 0.0
    },
    "expected": {
      "dna_repair_capacity": 0.0,
      "formula": "0.6×0.0 + 0.2×0.0 + 0.2×0.0 = 0.0"
    }
  }
]
```

### Fixture 3: `ca125_inadequate_fixtures.json`

```json
[
  {
    "name": "inadequate_response_47_percent",
    "patient_id": "FIXTURE_CA125_INADEQUATE",
    "ca125_history": [
      {"cycle": 0, "value": 2842, "date": "2025-11-17"},
      {"cycle": 3, "value": 1500, "date": "2026-02-17"}
    ],
    "expected": {
      "drop_percent": 47.2,
      "inadequate_response": true,
      "trigger": true
    }
  },
  {
    "name": "adequate_response_then_rise",
    "patient_id": "FIXTURE_CA125_RISE",
    "ca125_history": [
      {"cycle": 0, "value": 2842, "date": "2025-11-17"},
      {"cycle": 3, "value": 1000, "date": "2026-02-17"},
      {"cycle": 6, "value": 1200, "date": "2026-05-17"}
    ],
    "expected": {
      "drop_percent": 64.8,
      "inadequate_response": false,
      "on_therapy_rise": true,
      "trigger": true
    }
  },
  {
    "name": "adequate_response_68_percent",
    "patient_id": "FIXTURE_CA125_ADEQUATE",
    "ca125_history": [
      {"cycle": 0, "value": 2842, "date": "2025-11-17"},
      {"cycle": 3, "value": 900, "date": "2026-02-17"}
    ],
    "expected": {
      "drop_percent": 68.3,
      "inadequate_response": false,
      "trigger": false
    }
  }
]
```

### Fixture 4: `2_of_3_triggers_fixtures.json`

```json
[
  {
    "name": "high_risk_2_of_3",
    "patient_id": "FIXTURE_HIGH_RISK",
    "signals": {
      "dna_repair_restoration": true,
      "pathway_escape": false,
      "ca125_inadequate": true
    },
    "probability": 0.75,
    "expected": {
      "risk_level": "HIGH",
      "signals_detected": 2,
      "trigger": true
    }
  },
  {
    "name": "medium_risk_1_of_3",
    "patient_id": "FIXTURE_MEDIUM_RISK",
    "signals": {
      "dna_repair_restoration": true,
      "pathway_escape": false,
      "ca125_inadequate": false
    },
    "probability": 0.60,
    "expected": {
      "risk_level": "MEDIUM",
      "signals_detected": 1,
      "trigger": false
    }
  },
  {
    "name": "low_risk_no_signals",
    "patient_id": "FIXTURE_LOW_RISK",
    "signals": {
      "dna_repair_restoration": false,
      "pathway_escape": false,
      "ca125_inadequate": false
    },
    "probability": 0.30,
    "expected": {
      "risk_level": "LOW",
      "signals_detected": 0,
      "trigger": false
    }
  }
]
```

**Acceptance:**
- [ ] 4 fixture files created
- [ ] Each file has 3 test cases
- [ ] JSON is valid (can parse)

---

## Task 7.2: Create Validator Scaffolds (Day 2)

**Deliverable:** 4 Python validator files with boilerplate

### Commands:
```bash
cd oncology-coPilot/oncology-backend-minimal/scripts/validation
touch validate_hr_restoration_pattern.py
touch validate_dna_repair_formula.py
touch validate_ca125_inadequate_response.py
touch validate_2_of_3_triggers.py
chmod +x validate_*.py
```

### Boilerplate Template (All 4 Files):

```python
#!/usr/bin/env python3
"""
Validator: [NAME]
Purpose: [PURPOSE]
Fixtures: scripts/validation/fixtures/sae_resistance/[NAME]_fixtures.json
Output: scripts/validation/out/[NAME]/report.json
Owner: Zo (implements validate_fixture logic)
Plumber: JR (creates scaffold)
"""
import json
from pathlib import Path
from datetime import datetime

def load_fixtures():
    """Load test fixtures"""
    fixture_path = Path(__file__).parent / "fixtures/sae_resistance/[NAME]_fixtures.json"
    with open(fixture_path) as f:
        return json.load(f)

def validate_fixture(fixture):
    """
    Validate single fixture
    
    ZO IMPLEMENTS THIS - Validation logic goes here
    
    Args:
        fixture: Test case from fixtures file
    
    Returns:
        dict: {"passed": bool, "name": str, "expected": ..., "actual": ..., "error": str or None}
    """
    # TODO: Zo implements validation logic
    return {
        "passed": False,
        "name": fixture["name"],
        "expected": fixture.get("expected"),
        "actual": None,
        "error": "NOT_IMPLEMENTED"
    }

def main():
    """Run all fixtures and generate report"""
    fixtures = load_fixtures()
    results = {
        "validator": "[NAME]",
        "status": "PASS",
        "fixtures_tested": len(fixtures),
        "fixtures_passed": 0,
        "fixtures_failed": 0,
        "fixtures": [],
        "metrics": {
            "sensitivity": 0.0,
            "specificity": 0.0,
            "false_positive_rate": 0.0
        },
        "provenance": {
            "code_version": "sprint_2.5",
            "date": datetime.now().isoformat()
        }
    }
    
    for fixture in fixtures:
        result = validate_fixture(fixture)
        results["fixtures"].append(result)
        if result["passed"]:
            results["fixtures_passed"] += 1
        else:
            results["fixtures_failed"] += 1
    
    # Update status
    if results["fixtures_failed"] > 0:
        results["status"] = "FAIL"
    
    # Write report.json
    output_dir = Path(__file__).parent / "out/[NAME]"
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "report.json", "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"{'✅' if results['status'] == 'PASS' else '❌'} {results['fixtures_passed']}/{len(fixtures)} fixtures passed")
    return results["status"] == "PASS"

if __name__ == "__main__":
    exit(0 if main() else 1)
```

**Replace `[NAME]` with:**
- `hr_restoration_pattern`
- `dna_repair_formula`
- `ca125_inadequate_response`
- `2_of_3_triggers`

**Acceptance:**
- [ ] 4 validator files created
- [ ] All have boilerplate structure
- [ ] Can run (will fail with NOT_IMPLEMENTED until Zo adds logic)

---

## Task 7.3: Update Ring-1 Suite (Day 3)

**File:** `oncology-coPilot/oncology-backend-minimal/scripts/validation/run_resistance_validation_suite.py`

**Find the RING1_VALIDATORS list and add:**

```python
RING1_VALIDATORS = [
    # Existing validators
    "validate_ov_nf1_playbook.py",
    "validate_synthetic_lethality_pilot_benchmark.py",
    
    # NEW - Sprint 2.5 (SAE Resistance Detection)
    "validate_hr_restoration_pattern.py",
    "validate_dna_repair_formula.py",
    "validate_ca125_inadequate_response.py",
    "validate_2_of_3_triggers.py",
]
```

**Test command:**
```bash
cd oncology-coPilot/oncology-backend-minimal
python3 scripts/validation/run_resistance_validation_suite.py --ring1
```

**Expected:** All validators run (will fail with NOT_IMPLEMENTED until Zo adds logic)

**Acceptance:**
- [ ] Ring-1 suite updated
- [ ] Can run suite without errors
- [ ] 4 new validators appear in output

---

## ✅ JR Completion Checklist

**Task 7 is COMPLETE when:**
- [ ] 4 fixture files exist with correct JSON data
- [ ] 4 validator files exist with boilerplate
- [ ] Ring-1 suite includes 4 new validators
- [ ] Can run `python3 scripts/validation/run_resistance_validation_suite.py --ring1`
- [ ] All files committed to git

**Then:** Hand off to Zo to implement validation logic in `validate_fixture()` functions.

---

## 🔗 Links

- **Parent Sprint:** `.cursor/MOAT/SPRINT_2_5_SAE_RESISTANCE_VALIDATION.md`
- **JR Task Tracker:** `.cursor/ayesha/JR_AGENT_TASKS.md`
- **Validation Plan:** `.cursor/MOAT/RESISTANCE_VALIDATION_PLAN.md`

---

**Status:** 🚀 **READY FOR JR** - All fixtures defined, boilerplate specified, clear acceptance criteria.
