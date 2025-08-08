#!/usr/bin/env python3
"""
Direct Oracle test with controlled sequences to demonstrate massive Zeta Scores.
This bypasses the coordinate system issues and directly tests the Oracle's capability.
"""

import requests
import json

def test_oracle_with_controlled_sequences():
    """Test the Oracle with sequences we know are different."""
    
    oracle_url = "https://crispro--zeta-oracle-v5-final-fix-zetaoracle-api.modal.run/invoke"
    
    # Create a large reference sequence (simulating a gene)
    base_sequence = "ATCGATCGATCG" * 1000  # 12kb sequence
    
    # Create a version with a single nucleotide change
    mutated_sequence = base_sequence[:5000] + "T" + base_sequence[5001:]  # Change one base in the middle
    
    print(f"🧪 DIRECT ORACLE TEST")
    print(f"📊 Sequence lengths: REF={len(base_sequence)}bp, ALT={len(mutated_sequence)}bp")
    print(f"🎯 Single nucleotide difference at position 5000")
    
    # Verify they're different
    differences = sum(1 for i in range(len(base_sequence)) if base_sequence[i] != mutated_sequence[i])
    print(f"🔍 Total differences: {differences}")
    
    if differences == 0:
        print("❌ ERROR: Sequences are identical!")
        return
    
    payload = {
        "action": "score",
        "params": {
            "reference_sequence": base_sequence,
            "alternate_sequence": mutated_sequence
        }
    }
    
    print(f"\n🚀 Calling Oracle...")
    
    try:
        response = requests.post(oracle_url, json=payload, timeout=300, verify=False)
        response.raise_for_status()
        result = response.json()
        
        print("\n✅ SUCCESS: Oracle response received!")
        print("--- ORACLE RESPONSE ---")
        print(json.dumps(result, indent=2))
        print("----------------------")
        
        zeta_score = result.get("zeta_score", 0.0)
        print(f"\n🎯 FINAL RESULT: Zeta Score = {zeta_score}")
        
        if abs(zeta_score) > 1000:
            print("🔥 HIGH-IMPACT SCORE: This demonstrates the Oracle's massive scoring capability!")
        elif abs(zeta_score) > 0.1:
            print("⚡ MODERATE IMPACT: Detectable difference detected.")
        else:
            print("🤔 LOW IMPACT: Minimal or no difference detected.")
            
    except requests.exceptions.RequestException as e:
        print(f"\n❌ FAILED: Could not call Oracle.")
        print(f"  - Error: {e}")

if __name__ == "__main__":
    test_oracle_with_controlled_sequences() 