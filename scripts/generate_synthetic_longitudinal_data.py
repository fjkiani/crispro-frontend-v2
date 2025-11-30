"""
Generate SYNTHETIC longitudinal data for Resistance Prophet FULL capability demonstration.

This creates realistic time-series data showing HOW Resistance Prophet would work
in a real prospective monitoring scenario (like Ayesha's case).

Synthetic data models:
1. Sensitive patients: Stable DNA repair, declining CA-125, no bypass activation
2. Resistant patients: DNA repair restoration, rising CA-125, bypass pathway activation

Author: Zo
Date: January 14, 2025
For: Demonstrating FULL Resistance Prophet capability
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime, timedelta
import random
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Resistance mechanisms
DDR_GENES = ["BRCA1", "BRCA2", "RAD51", "RAD51C", "PALB2", "ATM"]
BYPASS_GENES = ["KRAS", "NRAS", "PIK3CA", "AKT1", "MTOR"]


def generate_sensitive_patient(patient_id: str, baseline_ca125: float) -> dict:
    # Baseline mutations (pre-treatment)
    baseline_mutations = {
        "BRCA1": "pathogenic",
        "TP53": "pathogenic",
    }

    # CA-125 time-series (declining)
    ca125_timepoints = []
    current = baseline_ca125
    for month in range(0, 13, 2):
        if month == 0:
            value = current
        else:
            value = current * (0.70 ** (month / 2)) * random.uniform(0.9, 1.1)
        ca125_timepoints.append({
            "month": month,
            "value": max(value, 5.0),
            "date": (datetime.now() - timedelta(days=365 - month * 30)).isoformat(),
        })

    progression_mutations = baseline_mutations.copy()

    return {
        "patient_id": patient_id,
        "baseline_mutations": list(baseline_mutations.keys()),
        "progression_mutations": list(progression_mutations.keys()),
        "new_mutations": [],
        "ca125_timeseries": ca125_timepoints,
        "platinum_response": "sensitive",
        "resistant_binary": 0,
        "ground_truth_resistance_month": None,
    }


def generate_resistant_patient(patient_id: str, baseline_ca125: float, resistance_month: int = 6) -> dict:
    # Baseline mutations (pre-treatment)
    baseline_mutations = {
        "BRCA2": "pathogenic",
        "TP53": "pathogenic",
    }

    # CA-125 time-series: decline then rise starting ~2 months before resistance
    ca125_timepoints = []
    current = baseline_ca125
    for month in range(0, 13, 2):
        if month == 0:
            value = current
        elif month < resistance_month - 2:
            value = current * (0.70 ** (month / 2)) * random.uniform(0.9, 1.1)
        else:
            value = current * (1.40 ** ((month - (resistance_month - 2)) / 2)) * random.uniform(0.9, 1.1)
        ca125_timepoints.append({
            "month": month,
            "value": max(value, 5.0),
            "date": (datetime.now() - timedelta(days=365 - month * 30)).isoformat(),
        })

    # Resistance mechanisms appear by resistance month
    progression_mutations = baseline_mutations.copy()
    progression_mutations["RAD51"] = "gain_of_function"  # HR restoration
    progression_mutations[random.choice(BYPASS_GENES)] = "activating"  # bypass

    new_mutations = [m for m in progression_mutations if m not in baseline_mutations]

    return {
        "patient_id": patient_id,
        "baseline_mutations": list(baseline_mutations.keys()),
        "progression_mutations": list(progression_mutations.keys()),
        "new_mutations": new_mutations,
        "ca125_timeseries": ca125_timepoints,
        "platinum_response": "resistant",
        "resistant_binary": 1,
        "ground_truth_resistance_month": resistance_month,
    }


def generate_synthetic_cohort(n_sensitive: int = 80, n_resistant: int = 40) -> list:
    patients = []
    for i in range(n_sensitive):
        pid = f"SYN_SENS_{i+1:03d}"
        baseline = random.uniform(500, 3000)
        patients.append(generate_sensitive_patient(pid, baseline))
    for i in range(n_resistant):
        pid = f"SYN_RESIST_{i+1:03d}"
        baseline = random.uniform(500, 3000)
        resistance_month = random.choice([4, 6, 8, 10])
        patients.append(generate_resistant_patient(pid, baseline, resistance_month))
    return patients


def main():
    output_dir = Path("data/validation")
    output_dir.mkdir(parents=True, exist_ok=True)

    synthetic_patients = generate_synthetic_cohort(80, 40)

    output_file = output_dir / "synthetic_longitudinal_resistance_dataset.json"
    payload = {
        "metadata": {
            "source": "Synthetic longitudinal data for Resistance Prophet demonstration",
            "generation_date": datetime.utcnow().isoformat(),
            "total_patients": len(synthetic_patients),
            "sensitive": sum(1 for p in synthetic_patients if p["platinum_response"] == "sensitive"),
            "resistant": sum(1 for p in synthetic_patients if p["platinum_response"] == "resistant"),
            "has_ca125_timeseries": len(synthetic_patients),
        },
        "patients": synthetic_patients,
    }

    with open(output_file, "w") as f:
        json.dump(payload, f, indent=2)

    logger.info(f"✅ SYNTHETIC DATASET SAVED: {output_file}")
    logger.info(
        f"✅ {payload['metadata']['sensitive']} sensitive / {payload['metadata']['resistant']} resistant with full longitudinal data"
    )

    return output_file


if __name__ == "__main__":
    path = main()
    print("\n🔥 SYNTHETIC DATA GENERATED - DEMONSTRATES FULL CAPABILITY")
    print("🎯 NEXT STEP: Run validation with longitudinal data")
    print("python scripts/validate_resistance_prophet_longitudinal.py")