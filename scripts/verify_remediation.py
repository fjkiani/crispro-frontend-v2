
import sys
import os

# Add the project root to sys.path to allow imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# Try to import the modified constants
try:
    from oncology_coPilot.oncology_backend_minimal.api.services.synthetic_lethality.constants import SYNTHETIC_LETHALITY_MAP
    
    print("=== SYNTHETIC_LETHALITY_MAP VERIFICATION ===")
    
    # 1. Check for BER key
    if "BER" in SYNTHETIC_LETHALITY_MAP:
        print("FAIL: 'BER' key still present in map!")
        print(f"Content: {SYNTHETIC_LETHALITY_MAP['BER']}")
        sys.exit(1)
    else:
        print("PASS: 'BER' key successfully removed.")
        
    # 2. Check for PARP key (should remain for HR)
    if "HR" in SYNTHETIC_LETHALITY_MAP:
        hr_deps = SYNTHETIC_LETHALITY_MAP["HR"]
        has_parp = any(d["pathway_id"] == "PARP" for d in hr_deps)
        if has_parp:
            print("PASS: 'HR' -> 'PARP' link preserved (Correct).")
        else:
            print("FAIL: 'HR' -> 'PARP' link missing!")
            sys.exit(1)
            
    print("\nRemediation Verified Successful.")
    
except ImportError as e:
    # Handle the complex import path issue locally
    # Just read the file content directly as fallback
    print(f"Import failed ({e}), checking file content directly...")
    
    target_file = "oncology-coPilot/oncology-backend-minimal/api/services/synthetic_lethality/constants.py"
    if not os.path.exists(target_file):
        print(f"FAIL: File not found at {target_file}")
        sys.exit(1)
        
    with open(target_file, "r") as f:
        content = f.read()
        
    if "'BER': [" in content:
        print("FAIL: 'BER': [ found in file content!")
        sys.exit(1)
    else:
        print("PASS: 'BER': [ not found in file content.")

if __name__ == "__main__":
    pass
