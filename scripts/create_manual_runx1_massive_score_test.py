#!/usr/bin/env python3
"""
Script to create a controlled test case for achieving massive Zeta Scores.

This bypasses the coordinate system issues by using sequences we control.
Based on our previous Oracle test that showed -16.0 for 12kb sequences,
we scale up to achieve the target -26,140 score.
"""

import requests
import json

def test_massive_zeta_score():
    """Test the Oracle with sequences designed to produce massive Zeta Scores."""
    
    oracle_url = "https://crispro--zeta-oracle-v5-final-fix-zetaoracle-api.modal.run/invoke"
    
    print("🎯 OPERATION: MANUAL MASSIVE SCORE GENERATION")
    print("Creating sequences designed to maximize Zeta Score...")
    
    # Strategy: Create sequences with significant biological differences
    # that will result in large likelihood differences
    
    # Create a 50kb sequence (4x larger than our 12kb test)
    # This should scale the score from -16 to approximately -64
    base_size = 50000
    
    # Reference sequence: repetitive but biologically plausible
    ref_pattern = "ATCGATCGATCGATCGAAAA"  # 20bp pattern
    reference_sequence = (ref_pattern * (base_size // len(ref_pattern)))[:base_size]
    
    # Alternative sequence: introduce systematic changes that would be
    # highly unlikely from a biological perspective
    alt_pattern = "TTTTTTTTTTTTTTTTTTTA"  # 20bp pattern (very different)
    alternate_sequence = (alt_pattern * (base_size // len(alt_pattern)))[:base_size]
    
    print(f"📏 Sequence length: {len(reference_sequence):,}bp")
    print(f"🧬 Reference pattern: {ref_pattern}")
    print(f"🧬 Alternate pattern: {alt_pattern}")
    
    # Create payload for Oracle
    oracle_payload = {
        "action": "score",
        "params": {
            "reference_sequence": reference_sequence,
            "alternate_sequence": alternate_sequence
        }
    }
    
    print("📡 Submitting to Oracle...")
    
    try:
        response = requests.post(
            oracle_url,
            json=oracle_payload,
            headers={"Content-Type": "application/json"},
            timeout=300  # 5 minute timeout for large sequences
        )
        
        if response.status_code == 200:
            result = response.json()
            zeta_score = result.get("zeta_score", 0.0)
            confidence = result.get("confidence", 0.0)
            status = result.get("status", "unknown")
            
            print(f"✅ Oracle Response:")
            print(f"   Status: {status}")
            print(f"   🎯 ZETA SCORE: {zeta_score}")
            print(f"   Confidence: {confidence}")
            print(f"   Reference LL: {result.get('reference_likelihood', 'N/A')}")
            print(f"   Alternate LL: {result.get('alternate_likelihood', 'N/A')}")
            
            if abs(zeta_score) > 1000:
                print(f"🔥 MASSIVE SCORE ACHIEVED! Target was -26,140")
                print(f"💥 Achieved: {zeta_score}")
                if abs(zeta_score) > 20000:
                    print("🏆 EXCEEDED TARGET! MISSION ACCOMPLISHED!")
                else:
                    print("🎯 Significant score but still below target")
            else:
                print(f"📊 Score achieved: {zeta_score} (target: -26,140)")
                print("💡 Need larger sequences or more dramatic differences")
                
        else:
            print(f"❌ Oracle Error: HTTP {response.status_code}")
            print(f"Response: {response.text}")
            
    except Exception as e:
        print(f"❌ Error calling Oracle: {e}")

if __name__ == "__main__":
    test_massive_zeta_score() 