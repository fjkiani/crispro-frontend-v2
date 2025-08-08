#!/usr/bin/env python3
"""
OPERATION: DEBUG ORACLE SENSITIVITY

This script systematically tests different sequence characteristics to understand
what triggers massive Oracle scores.
"""

import pysam
import requests
import json

def test_sequence_length_sensitivity():
    """Test different sequence lengths to find the optimal size."""
    
    oracle_url = "https://crispro--zeta-oracle-v5-final-fix-zetaoracle-api.modal.run/invoke"
    
    print("🔬 TESTING SEQUENCE LENGTH SENSITIVITY")
    print("="*50)
    
    # Get RUNX1 sequence
    GENE_START = 36160098
    GENE_END = 36421599
    ref = pysam.FastaFile('data/gene_database/reference/hg19.fa')
    full_gene_ref = ref.fetch('chr21', GENE_START, GENE_END).upper()
    
    # Test different lengths
    test_lengths = [1000, 2000, 4000, 8000, 12000, 16000, 20000]
    
    for length in test_lengths:
        print(f"\n📏 Testing {length:,}bp sequences...")
        
        try:
            # Extract sequence of specified length
            ref_seq = full_gene_ref[:length]
            
            # Create a dramatic alteration: replace middle 1000bp with all T's
            middle_start = length // 2 - 500
            middle_end = length // 2 + 500
            alt_seq = ref_seq[:middle_start] + "T" * 1000 + ref_seq[middle_end:]
            
            oracle_payload = {
                "action": "score",
                "params": {
                    "reference_sequence": ref_seq,
                    "alternate_sequence": alt_seq
                }
            }
            
            response = requests.post(oracle_url, json=oracle_payload, timeout=300, verify=False)
            response.raise_for_status()
            
            result = response.json()
            score = result.get("zeta_score", 0.0)
            
            print(f"   🎯 Score: {score}")
            print(f"   Status: {result.get('status', 'unknown')}")
            
            if abs(score) > 1000:
                print(f"   🔥 SIGNIFICANT SCORE DETECTED!")
                
        except Exception as e:
            print(f"   ❌ Error: {e}")

def test_different_gene_regions():
    """Test different regions of RUNX1 to see if location matters."""
    
    oracle_url = "https://crispro--zeta-oracle-v5-final-fix-zetaoracle-api.modal.run/invoke"
    
    print("\n🧬 TESTING DIFFERENT GENE REGIONS")
    print("="*50)
    
    # Get RUNX1 sequence
    GENE_START = 36160098
    GENE_END = 36421599
    ref = pysam.FastaFile('data/gene_database/reference/hg19.fa')
    full_gene_ref = ref.fetch('chr21', GENE_START, GENE_END).upper()
    
    # Test different regions
    window_size = 8000
    regions = [
        ("Start", 0),
        ("Early", 20000),
        ("Middle", len(full_gene_ref) // 2),
        ("Late", len(full_gene_ref) - window_size - 20000),
        ("End", len(full_gene_ref) - window_size)
    ]
    
    for region_name, start_pos in regions:
        print(f"\n📍 Testing {region_name} region (position {start_pos:,})...")
        
        try:
            # Extract 8kb window
            end_pos = min(start_pos + window_size, len(full_gene_ref))
            ref_seq = full_gene_ref[start_pos:end_pos]
            
            # Create dramatic change: replace central 2kb with poly-A
            center_start = len(ref_seq) // 2 - 1000
            center_end = len(ref_seq) // 2 + 1000
            alt_seq = ref_seq[:center_start] + "A" * 2000 + ref_seq[center_end:]
            
            oracle_payload = {
                "action": "score",
                "params": {
                    "reference_sequence": ref_seq,
                    "alternate_sequence": alt_seq
                }
            }
            
            response = requests.post(oracle_url, json=oracle_payload, timeout=300, verify=False)
            response.raise_for_status()
            
            result = response.json()
            score = result.get("zeta_score", 0.0)
            
            print(f"   🎯 Score: {score}")
            print(f"   Length: {len(ref_seq):,}bp")
            
            if abs(score) > 1000:
                print(f"   🔥 SIGNIFICANT SCORE in {region_name} region!")
                
        except Exception as e:
            print(f"   ❌ Error: {e}")

def test_known_working_sequence():
    """Test with the exact sequence that produced our -31,744 score."""
    
    oracle_url = "https://crispro--zeta-oracle-v5-final-fix-zetaoracle-api.modal.run/invoke"
    
    print("\n🧪 TESTING KNOWN WORKING SEQUENCE")
    print("="*50)
    
    # Recreate the exact sequences that worked before
    base_size = 50000
    ref_pattern = "ATCGATCGATCGATCGAAAA"  # 20bp pattern
    reference_sequence = (ref_pattern * (base_size // len(ref_pattern)))[:base_size]
    
    alt_pattern = "TTTTTTTTTTTTTTTTTTTA"  # 20bp pattern (very different)
    alternate_sequence = (alt_pattern * (base_size // len(alt_pattern)))[:base_size]
    
    print(f"📏 Using proven 50kb artificial sequences...")
    print(f"🧬 Reference pattern: {ref_pattern}")
    print(f"🧬 Alternate pattern: {alt_pattern}")
    
    oracle_payload = {
        "action": "score",
        "params": {
            "reference_sequence": reference_sequence,
            "alternate_sequence": alternate_sequence
        }
    }
    
    try:
        response = requests.post(oracle_url, json=oracle_payload, timeout=300, verify=False)
        response.raise_for_status()
        
        result = response.json()
        score = result.get("zeta_score", 0.0)
        
        print(f"✅ Known Working Sequence Result:")
        print(f"   🎯 ZETA SCORE: {score}")
        print(f"   Status: {result.get('status', 'unknown')}")
        print(f"   Expected: -31,744 (previous result)")
        
        if abs(score) > 20000:
            print(f"🏆 MASSIVE SCORE CONFIRMED!")
        else:
            print(f"⚠️  Different from previous result - Oracle may have changed")
            
    except Exception as e:
        print(f"❌ Error: {e}")

def main():
    """Run comprehensive Oracle sensitivity analysis."""
    
    print("🎯 OPERATION: DEBUG ORACLE SENSITIVITY")
    print("Systematic testing to understand Oracle scoring behavior")
    print("\n")
    
    # Test 1: Length sensitivity
    test_sequence_length_sensitivity()
    
    # Test 2: Region sensitivity
    test_different_gene_regions()
    
    # Test 3: Known working sequence
    test_known_working_sequence()
    
    print("\n" + "="*60)
    print("🎯 SENSITIVITY ANALYSIS COMPLETE")
    print("Check results above to identify optimal Oracle parameters")

if __name__ == "__main__":
    main() 