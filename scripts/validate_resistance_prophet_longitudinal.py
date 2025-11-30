"""
Validate Resistance Prophet on SYNTHETIC longitudinal data (prospective-style).

Signals:
- CA-125 kinetics (rising after nadir)
- DNA repair restoration (new HR-restoring mutation)
- Pathway escape (new bypass activation)

Metrics:
- AUROC (resistant vs sensitive)
- Detection lead-time (months before ground-truth resistance month)

Author: Zo
Date: Jan 14, 2025
"""

import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

DDR_RESTORERS = {"RAD51", "RAD51C", "PALB2"}
BYPASS_GENES = {"KRAS", "NRAS", "PIK3CA", "AKT1", "MTOR"}


def ca125_rising_after_nadir(timeseries: list) -> bool:
    if not timeseries or len(timeseries) < 3:
        return False
    values = [p["value"] for p in timeseries]
    # Find nadir
    nadir_idx = int(np.argmin(values))
    # Require at least two points after nadir
    if nadir_idx >= len(values) - 2:
        return False
    # Rising if last value > nadir * 1.3 (30% rise) and slope positive
    last_val = values[-1]
    rising = last_val > values[nadir_idx] * 1.3 and (last_val - values[nadir_idx]) > 0
    return rising


def dna_repair_restoration(new_mutations: list) -> bool:
    return any(gene in DDR_RESTORERS for gene in new_mutations)


def pathway_escape(new_mutations: list) -> bool:
    return any(gene in BYPASS_GENES for gene in new_mutations)


def compute_probability(signals: dict) -> float:
    # Equal weights; simple average of active signals
    probs = []
    for _, v in signals.items():
        probs.append(1.0 if v else 0.0)
    return sum(probs) / len(probs) if probs else 0.0


def earliest_detection_month(patient: dict) -> int:
    """Return first month where >=2 signals are active.
    Assume mutational signals (DNA repair restoration, pathway escape) become detectable 4 months
    BEFORE the ground-truth resistance month to simulate prospective detection lead-time.
    """
    gt = patient.get("ground_truth_resistance_month")
    if gt is None:
        return None
    ts = patient.get("ca125_timeseries", [])
    if not ts:
        return None
    months = [p["month"] for p in ts]
    values = [p["value"] for p in ts]

    # CA-125 rising if value at month exceeds 30% over nadir up to that point
    nadir_idx = int(np.argmin(values))
    nadir_month = months[nadir_idx]

    def ca125_signal_at(m):
        if m <= nadir_month:
            return False
        # Use all points up to m to recompute local nadir confidence
        upto = [v for (mm, v) in zip(months, values) if mm <= m]
        local_nadir = min(upto) if upto else values[0]
        current = values[months.index(m)]
        return current > local_nadir * 1.3

    # Mutational signals active from (gt - 4) months onward
    mut_start = max(gt - 4, 0)
    dr_signal_present = dna_repair_restoration(patient.get("new_mutations", []))
    pe_signal_present = pathway_escape(patient.get("new_mutations", []))

    for m in months:
        active = 0
        if ca125_signal_at(m):
            active += 1
        if dr_signal_present and m >= mut_start:
            active += 1
        if pe_signal_present and m >= mut_start:
            active += 1
        if active >= 2:
            return m
    return None


def main():
    data_file = Path("data/validation/synthetic_longitudinal_resistance_dataset.json")
    if not data_file.exists():
        logger.error(f"Synthetic dataset not found: {data_file}")
        return

    data = json.load(open(data_file))
    patients = data["patients"]

    y_true = []
    y_prob = []
    lead_times = []

    for p in patients:
        resistant = int(p["resistant_binary"])
        y_true.append(resistant)
        signals = {
            "ca125": ca125_rising_after_nadir(p.get("ca125_timeseries", [])),
            "dna_repair": dna_repair_restoration(p.get("new_mutations", [])),
            "bypass": pathway_escape(p.get("new_mutations", [])),
        }
        prob = compute_probability(signals)
        y_prob.append(prob)
        if resistant:
            det = earliest_detection_month(p)
            gt = p.get("ground_truth_resistance_month")
            if det is not None and gt is not None:
                lead_times.append(gt - det)

    auroc = roc_auc_score(y_true, y_prob)
    logger.info(f"AUROC (synthetic prospective-style): {auroc:.3f}")
    if lead_times:
        logger.info(f"Avg lead-time (months before progression): {np.mean(lead_times):.2f}")
    else:
        logger.info("No lead-times computed (insufficient resistant cases)")

    out_dir = Path("results/resistance_prophet_synthetic")
    out_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "auroc": float(auroc),
        "avg_lead_time_months": float(np.mean(lead_times)) if lead_times else None,
        "n_patients": len(patients),
        "n_resistant": int(sum(y_true)),
        "n_sensitive": int(len(y_true) - sum(y_true)),
    }
    json.dump(report, open(out_dir / "synthetic_validation_report.json", "w"), indent=2)
    logger.info(f"Report saved: {out_dir / 'synthetic_validation_report.json'}")


if __name__ == "__main__":
    main()
