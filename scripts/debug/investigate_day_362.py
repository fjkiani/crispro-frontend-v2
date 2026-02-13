
import sys
import os
import pandas as pd
from datetime import datetime

# Ensure backend modules are importable
sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from api.services.ca125_intelligence_service import CA125IntelligenceService

def investigate_day_362():
    print("=== Investigation: OV-045 Day 362 (Lead Time Opportunity) ===")
    
    # 1. Setup Data
    history = [
        {"date": "Day 159", "value": 11.0},
        {"date": "Day 173", "value": 11.0},
        {"date": "Day 265", "value": 9.0},  # Nadir
        {"date": "Day 362", "value": 14.0}, # Target
    ]
    
    # 2. Run Current Logic (Noise Floor 20)
    service = CA125IntelligenceService()
    # Mocking sort order manually
    result_std = None
    
    # We need to run it async? No, logic is pure data? 
    # analyze_kinetics is async.
    import asyncio
    
    current_val = 14.0
    nadir = 9.0
    ratio = current_val / nadir
    
    print(f"Data: Current={current_val}, Nadir={nadir}, Ratio={ratio:.2f}x")
    
    # Check 1: Standard
    print("\n--- Standard Logic (Noise Floor 20.0) ---")
    if current_val > 20.0:
        print("Noise Floor Passed (>20). Checking Kinetics...")
    else:
        print("❌ Blocked by Noise Floor (14.0 < 20.0). Signal Suppressed.")
        
    # Check 2: Hypothetical Lower Floor
    print("\n--- Optimized Logic (Noise Floor 10.0) ---")
    if current_val > 10.0:
        print("✅ Noise Floor Passed (>10). Checking Kinetics...")
        if ratio >= 1.5:
            print(f"✅ RISING DETECTED (Ratio {ratio:.2f} >= 1.5). Lead Time Gained.")
        else:
            print(f"❌ Ratio insufficient ({ratio:.2f} < 1.5).")
            
    # Check 3: ctDNA
    print("\n--- ctDNA Signal (Day 362) ---")
    tf = 0.01344 # From CSV
    print(f"Tumor Fraction: {tf:.5f} (1.3%)")
    limit = 0.005 # 0.5%
    if tf > limit:
        print(f"✅ MOLECULAR RECURRENCE (TF > {limit}). Lead Time Gained.")
    else:
        print("❌ TF below limit.")
        
    print("\n=== Conclusion ===")
    print("Current Phase 1c logic confirms 'Low Risk' due to Noise Floor.")
    print("Phase 2 Logic (Floor=10 or ctDNA) would detect at Day 362 (109 Days Lead Time).")

if __name__ == "__main__":
    investigate_day_362()
