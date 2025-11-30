#!/usr/bin/env python3
"""
Merge HRD Scores - Day 2 Task 2.3
Agent: Jr | Date: January 28, 2025
"""

import json
import sys
from pathlib import Path

DATA_DIR = Path("data/validation")
VALIDATED_FILE = DATA_DIR / "tcga_ov_469_validated.json"
HRD_FILE = DATA_DIR / "tcga_ov_hrd_scores_raw.json"
OUTPUT_FILE = DATA_DIR / "tcga_ov_469_with_hrd.json"

def main():
    print("🧬 Merge HRD Scores - Day 2")
    
    with open(VALIDATED_FILE, 'r') as f:
        patients = json.load(f)
    print(f"✅ Loaded {len(patients)} patients")
    
    with open(HRD_FILE, 'r') as f:
        hrd_data = json.load(f)
    hrd_scores = hrd_data.get("patients", {})
    print(f"✅ Loaded HRD for {len(hrd_scores)} patients")
    
    merged = 0
    for p in patients:
        pid = p.get("patient_id", "")
        hrd = hrd_scores.get(pid, {})
        if hrd:
            if "biomarkers" not in p:
                p["biomarkers"] = {}
            p["biomarkers"]["hrd_score"] = hrd.get("hrd_sum")
            p["biomarkers"]["hrd_status"] = hrd.get("hrd_status")
            p["biomarkers"]["loh_score"] = hrd.get("loh_score")
            p["biomarkers"]["tai_score"] = hrd.get("tai_score")
            p["biomarkers"]["lst_score"] = hrd.get("lst_score")
            merged += 1
    
    print(f"✅ Merged: {merged}/{len(patients)}")
    
    hrd_plus = sum(1 for p in patients if p.get("biomarkers", {}).get("hrd_status") == "HRD+")
    hrd_minus = sum(1 for p in patients if p.get("biomarkers", {}).get("hrd_status") == "HRD-")
    print(f"📊 HRD+: {hrd_plus} | HRD-: {hrd_minus}")
    
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(patients, f, indent=2)
    
    summary = {"total": len(patients), "with_hrd": merged, "hrd_plus": hrd_plus, "hrd_minus": hrd_minus}
    with open(DATA_DIR / "tcga_ov_hrd_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"✅ Saved to {OUTPUT_FILE}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
