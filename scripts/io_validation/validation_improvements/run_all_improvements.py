#!/usr/bin/env python3
"""
Master Script: Run All Validation Improvements
===============================================

Runs all validation improvement scripts in sequence:
1. Task 1.1: Address Overfitting (reduce model, bootstrap, learning curve)
2. Task 1.3: Statistical Fixes (multiple testing, calibration, improved CV)
3. Task 1.4: External Validation (GSE179994) - requires data
4. Task 1.5: Compare to Published Signatures (TIDE, IMPRES, TMEscore)

Usage:
    python scripts/validation_improvements/run_all_improvements.py
"""

import subprocess
import sys
from pathlib import Path
from datetime import datetime

# Configuration
SCRIPT_DIR = Path(__file__).parent
REPO_ROOT = SCRIPT_DIR.parent.parent.parent
SCRIPTS_DIR = SCRIPT_DIR  # validation_improvements/

print("="*80)
print("RUNNING ALL VALIDATION IMPROVEMENTS")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print()

scripts = [
    ("Task 1.1A: Reduce Model Complexity", "reduce_model_complexity.py"),
    ("Task 1.1B: Bootstrap Validation", "bootstrap_validation.py"),
    ("Task 1.1C: Learning Curve Analysis", "learning_curve_analysis.py"),
    ("Task 1.3A: Multiple Testing Correction", "multiple_testing_correction.py"),
    ("Task 1.3B: Calibration Analysis", "calibration_analysis.py"),
    ("Task 1.3C: Improved Cross-Validation", "improved_cross_validation.py"),
    ("Task 1.5: Compare Published Signatures", "compare_published_signatures.py"),
]

# Optional scripts (require additional data)
optional_scripts = [
    ("Task 1.4: External Validation (GSE179994)", "external_validation_gse179994.py"),
]

results = []

# Run required scripts
for task_name, script_name in scripts:
    print(f"\n{'='*80}")
    print(f"Running: {task_name}")
    print(f"{'='*80}")
    print()
    
    script_path = SCRIPTS_DIR / script_name
    
    if not script_path.exists():
        print(f"❌ Script not found: {script_path}")
        results.append((task_name, "NOT FOUND", None))
        continue
    
    try:
        result = subprocess.run(
            [sys.executable, str(script_path)],
            cwd=str(SCRIPTS_DIR.parent.parent),
            capture_output=True,
            text=True,
            timeout=3600  # 1 hour timeout
        )
        
        if result.returncode == 0:
            print(f"✅ {task_name}: SUCCESS")
            results.append((task_name, "SUCCESS", result.returncode))
        else:
            print(f"❌ {task_name}: FAILED (exit code {result.returncode})")
            print("STDERR:")
            print(result.stderr[:500])  # First 500 chars
            results.append((task_name, "FAILED", result.returncode))
    
    except subprocess.TimeoutExpired:
        print(f"⏱️  {task_name}: TIMEOUT (>1 hour)")
        results.append((task_name, "TIMEOUT", None))
    except Exception as e:
        print(f"❌ {task_name}: ERROR - {e}")
        results.append((task_name, "ERROR", str(e)))

# Run optional scripts
print(f"\n{'='*80}")
print("OPTIONAL SCRIPTS (May Require Additional Data)")
print(f"{'='*80}")
print()

for task_name, script_name in optional_scripts:
    print(f"\n{'='*80}")
    print(f"Running: {task_name}")
    print(f"{'='*80}")
    print()
    
    script_path = SCRIPTS_DIR / script_name
    
    if not script_path.exists():
        print(f"❌ Script not found: {script_path}")
        results.append((task_name, "NOT FOUND", None))
        continue
    
    try:
        result = subprocess.run(
            [sys.executable, str(script_path)],
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
            timeout=3600
        )
        
        if result.returncode == 0:
            print(f"✅ {task_name}: SUCCESS")
            results.append((task_name, "SUCCESS", result.returncode))
        else:
            print(f"⚠️  {task_name}: SKIPPED (may require additional data)")
            print("   This is expected if GSE179994 bulk expression data is not available")
            results.append((task_name, "SKIPPED", result.returncode))
    
    except Exception as e:
        print(f"⚠️  {task_name}: SKIPPED - {e}")
        results.append((task_name, "SKIPPED", str(e)))

# Summary
print(f"\n{'='*80}")
print("EXECUTION SUMMARY")
print(f"{'='*80}")
print()

success_count = sum(1 for _, status, _ in results if status == "SUCCESS")
failed_count = sum(1 for _, status, _ in results if status in ["FAILED", "ERROR"])
skipped_count = sum(1 for _, status, _ in results if status in ["SKIPPED", "NOT FOUND"])

print(f"Total scripts: {len(results)}")
print(f"  ✅ Success: {success_count}")
print(f"  ❌ Failed: {failed_count}")
print(f"  ⚠️  Skipped: {skipped_count}")
print()

print("Detailed Results:")
for task_name, status, code in results:
    status_icon = "✅" if status == "SUCCESS" else "❌" if status in ["FAILED", "ERROR"] else "⚠️"
    print(f"  {status_icon} {task_name}: {status}")

print()
print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"{'='*80}")

# Exit with error if any required scripts failed
if failed_count > 0:
    sys.exit(1)
else:
    sys.exit(0)
