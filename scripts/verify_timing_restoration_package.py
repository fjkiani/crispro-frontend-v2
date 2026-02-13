#!/usr/bin/env python3
"""
TIMING ENGINE RESTORATION - FINAL VERIFICATION PACKAGE
Generates the evidence package required for "MISSION COMPLETE" sign-off.

Outputs:
1. File Manifest (Paths + SHA256 Hashes)
2. API Contract Check (Import + Schema Validation)
3. Reproducible Validation Run (Executes validate_timing_engine_REAL.py)
"""

import sys
import hashlib
import json
import subprocess
from pathlib import Path
from datetime import datetime

# Define base path
REPO_ROOT = Path("oncology-coPilot/oncology-backend-minimal").absolute()

def calculate_sha256(filepath):
    """Calculate SHA256 hash of a file."""
    try:
        sha256_hash = hashlib.sha256()
        with open(filepath, "rb") as f:
            for byte_block in iter(lambda: f.read(4096), b""):
                sha256_hash.update(byte_block)
        return sha256_hash.hexdigest()
    except FileNotFoundError:
        return "FILE_NOT_FOUND"

def verify_files():
    """Verify existence and integrity of restored files."""
    files_to_check = [
        "api/services/resistance/biomarkers/therapeutic/timing_chemo_features.py",
        "api/services/resistance/biomarkers/therapeutic/ca125_kelim_ovarian.py",
        "api/routers/resistance.py"
    ]
    
    manifest = {}
    print("="*80)
    print("1. FILE MANIFEST")
    print("="*80)
    
    for rel_path in files_to_check:
        full_path = REPO_ROOT / rel_path
        status = "✅ FOUND" if full_path.exists() else "❌ MISSING"
        file_hash = calculate_sha256(full_path)
        manifest[rel_path] = file_hash
        print(f"{status} {rel_path}")
        print(f"   SHA256: {file_hash}")
        
    return manifest

def verify_api_contract():
    """Verify API router import and schema definitions."""
    print("\n" + "="*80)
    print("2. API CONTRACT CHECK")
    print("="*80)
    
    try:
        sys.path.insert(0, str(REPO_ROOT))
        from api.routers.resistance import router, TimingChemoFeaturesRequest, TimingChemoFeaturesResponse
        
        print("✅ Router Import Successful")
        print(f"✅ Request Model: {TimingChemoFeaturesRequest.__name__}")
        print(f"✅ Response Model: {TimingChemoFeaturesResponse.__name__}")
        
        # Check specific fields
        req_fields = TimingChemoFeaturesRequest.model_fields.keys()
        res_fields = TimingChemoFeaturesResponse.model_fields.keys()
        
        print(f"   Request Fields: {', '.join(req_fields)}")
        print(f"   Response Fields: {', '.join(res_fields)}")
        
        if 'regimen_table' in req_fields and 'timing_features_table' in res_fields:
            print("✅ Schema Compliant with Frontend Contract")
        else:
            print("❌ Schema Mismatch Check Failed")
            
    except Exception as e:
        print(f"❌ API Contract Check Failed: {str(e)}")

def run_validation_script():
    """Execute the cohort validation script and capture output."""
    print("\n" + "="*80)
    print("3. REPRODUCIBLE VALIDATION RUN")
    print("="*80)
    
    script_path = "scripts/cohorts/validate_timing_engine_REAL.py"
    cmd = ["python3", script_path]
    
    print(f"Executing: {' '.join(cmd)}")
    print("-" * 40)
    
    try:
        result = subprocess.run(
            cmd,
            cwd=str(Path.cwd()),
            capture_output=True,
            text=True
        )
        print(result.stdout)
        if result.stderr:
            print("STDERR:")
            print(result.stderr)
            
        if result.returncode == 0:
            print(f"\n✅ Validation Script Exit Code: 0")
        else:
            print(f"\n❌ Validation Script Failed (Exit Code: {result.returncode})")
            
    except Exception as e:
        print(f"❌ Execution Failed: {str(e)}")

if __name__ == "__main__":
    print(f"TIMING ENGINE RESTORATION AUDIT - {datetime.now().isoformat()}")
    print("-" * 80)
    
    verify_files()
    verify_api_contract()
    run_validation_script()
