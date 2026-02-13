#!/usr/bin/env python3
"""
TIMING ENGINE - RIGOROUS AUDIT PACKAGE GENERATOR
Collects Git state, Environment info, Data Provenance, and Validation Logs.
"""

import subprocess
import json
import hashlib
from pathlib import Path
from datetime import datetime

# CONFIG
REPO_ROOT = Path.cwd()
DATA_FILE = REPO_ROOT / "oncology-coPilot/oncology-backend-minimal/data/TCGA-OV/ds_cci.17.00096-1.xlsx"
VALIDATION_SCRIPT = "scripts/cohorts/validate_timing_engine_REAL.py"
OUTPUT_DIR = REPO_ROOT / "artifacts/audit_pkg"

def run_cmd(cmd, cwd=None):
    try:
        if cwd:
            cwd = str(cwd)
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd)
        return {
            "cmd": cmd,
            "returncode": result.returncode,
            "stdout": result.stdout.strip(),
            "stderr": result.stderr.strip()
        }
    except Exception as e:
        return {"cmd": cmd, "error": str(e)}

def sha256_file(filepath):
    try:
        h = hashlib.sha256()
        with open(filepath, "rb") as f:
            for b in iter(lambda: f.read(4096), b""):
                h.update(b)
        return h.hexdigest()
    except FileNotFoundError:
        return "FILE_NOT_FOUND"

def generate_audit():
    timestamp = datetime.now().isoformat()
    audit_data = {
        "timestamp": timestamp,
        "repo_state": {},
        "environment": {},
        "data_provenance": {},
        "validation_run": {}
    }
    
    print(f"Starting Audit Generation... {timestamp}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # 1. Repo State
    print("Capturing Git State...")
    audit_data["repo_state"]["rev_parse_head"] = run_cmd("git rev-parse HEAD")
    audit_data["repo_state"]["status"] = run_cmd("git status --porcelain")
    audit_data["repo_state"]["diff_stat"] = run_cmd("git diff --stat")
    
    # Save full diff
    diff_res = run_cmd("git diff")
    with open(OUTPUT_DIR / "current.diff", "w") as f:
        f.write(diff_res["stdout"])
        
    # 2. Environment
    print("Capturing Environment...")
    audit_data["environment"]["python_version"] = run_cmd("python3 --version")
    
    # Save pip freeze
    pip_res = run_cmd("pip freeze")
    with open(OUTPUT_DIR / "pip_freeze.txt", "w") as f:
        f.write(pip_res["stdout"])
        
    # 3. Data Provenance
    print("Hashing Data Files...")
    audit_data["data_provenance"]["input_data"] = {
        str(DATA_FILE): sha256_file(DATA_FILE)
    }
    
    # 4. Validation Run
    print("Executing Validation Script...")
    val_res = run_cmd(f"python3 {VALIDATION_SCRIPT}")
    audit_data["validation_run"]["exit_code"] = val_res["returncode"]
    
    # Save full logs
    with open(OUTPUT_DIR / "validation_stdout.txt", "w") as f:
        f.write(val_res["stdout"])
    with open(OUTPUT_DIR / "validation_stderr.txt", "w") as f:
        f.write(val_res["stderr"])
        
    # Check for Circularity (Grep check)
    print("Checking for Circular Imports/Data Leakage...")
    engine_path = "oncology-coPilot/oncology-backend-minimal/api/services/resistance/biomarkers/therapeutic/timing_chemo_features.py"
    leak_check = run_cmd(f"grep 'ds_cci' {engine_path}") # Ensure engine doesn't read the dataset directly
    audit_data["data_provenance"]["leakage_check_cmd"] = leak_check["cmd"]
    audit_data["data_provenance"]["leakage_check_result"] = "PASS" if not leak_check["stdout"] else "FAIL (Filename found in engine)"

    # Write Manifest
    manifest_path = OUTPUT_DIR / "audit_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(audit_data, f, indent=2)
        
    print(f"Audit Complete. Artifacts saved to {OUTPUT_DIR}")
    print(f"Manifest: {manifest_path}")

if __name__ == "__main__":
    generate_audit()
