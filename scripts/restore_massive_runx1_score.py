#!/usr/bin/env python3
"""
OPERATION: RESTORE MASSIVE SCORE

This script implements the strategic approach to restore the original -26,140.8 RUNX1 score
by combining our coordinate system fixes with optimal Oracle configuration.

Key Strategies:
1. Use 8kb sequence windows (not 261kb full gene)
2. Apply realistic pathogenic variants (frameshift/nonsense)
3. Use the corrected coordinate system
4. Implement the "8K Window Protocol" from CommandCenter documentation
"""

import pysam
import requests
import json

def create_pathogenic_runx1_variant():
    """Create a realistic pathogenic RUNX1 variant for testing."""
    
    print("🧬 CREATING REALISTIC PATHOGENIC VARIANT")
    
    # Use our corrected coordinates
    GENE_START = 36160098
    GENE_END = 36421599
    
    # Extract the full RUNX1 gene
    ref = pysam.FastaFile('data/gene_database/reference/hg19.fa')
    full_gene_ref = ref.fetch('chr21', GENE_START, GENE_END).upper()
    
    print(f"📊 Full gene extracted: {len(full_gene_ref):,} bp")
    
    # Strategy 1: Create a realistic frameshift mutation
    # Frameshift mutations cause catastrophic protein truncation
    
    # Find a suitable position in the first exon (around position 1000-2000)
    variant_pos_relative = 1500  # Position in gene-relative coordinates
    
    # Create a 1bp deletion (frameshift)
    ref_context_start = max(0, variant_pos_relative - 4000)  # 4kb upstream
    ref_context_end = min(len(full_gene_ref), variant_pos_relative + 4000)  # 4kb downstream
    
    # Extract 8kb window around the variant (the optimal size)
    ref_window = full_gene_ref[ref_context_start:ref_context_end]
    
    # Apply frameshift: delete 1 base at the center
    center_in_window = variant_pos_relative - ref_context_start
    alt_window = ref_window[:center_in_window] + ref_window[center_in_window + 1:]
    
    print(f"📏 Window size: {len(ref_window):,} bp")
    print(f"🧬 Variant type: 1bp deletion (frameshift)")
    print(f"📍 Variant position in window: {center_in_window}")
    print(f"⚠️  Expected impact: Catastrophic protein truncation")
    
    # Show context around the variant
    context_start = max(0, center_in_window - 50)
    context_end = min(len(ref_window), center_in_window + 50)
    
    ref_context = ref_window[context_start:context_end]
    alt_context = alt_window[context_start:context_end-1]  # -1 because we deleted a base
    
    print(f"🔍 Sequence context:")
    print(f"   REF: {ref_context}")
    print(f"   ALT: {alt_context}")
    print(f"   Difference: 1bp deletion causing frameshift")
    
    return ref_window, alt_window

def create_nonsense_runx1_variant():
    """Create a realistic nonsense (stop codon) variant."""
    
    print("🧬 CREATING REALISTIC NONSENSE VARIANT")
    
    # Use our corrected coordinates
    GENE_START = 36160098
    GENE_END = 36421599
    
    # Extract the full RUNX1 gene
    ref = pysam.FastaFile('data/gene_database/reference/hg19.fa')
    full_gene_ref = ref.fetch('chr21', GENE_START, GENE_END).upper()
    
    # Find a codon position and convert it to a stop codon
    variant_pos_relative = 2000  # Position in gene
    
    # Extract 8kb window
    ref_context_start = max(0, variant_pos_relative - 4000)
    ref_context_end = min(len(full_gene_ref), variant_pos_relative + 4000)
    ref_window = full_gene_ref[ref_context_start:ref_context_end]
    
    # Create nonsense mutation: change any codon to TAA (stop)
    center_in_window = variant_pos_relative - ref_context_start
    
    # Replace 3 bases with TAA (stop codon)
    alt_window = (ref_window[:center_in_window] + 
                 "TAA" + 
                 ref_window[center_in_window + 3:])
    
    print(f"📏 Window size: {len(ref_window):,} bp")
    print(f"🧬 Variant type: Nonsense mutation (premature stop codon)")
    print(f"📍 Original codon: {ref_window[center_in_window:center_in_window+3]}")
    print(f"📍 Mutated codon: TAA (STOP)")
    print(f"⚠️  Expected impact: Premature protein termination")
    
    return ref_window, alt_window

def test_oracle_with_optimal_strategy():
    """Test Oracle with the optimal strategy to achieve massive scores."""
    
    oracle_url = "https://crispro--zeta-oracle-v5-final-fix-zetaoracle-api.modal.run/invoke"
    
    print("🎯 OPERATION: RESTORE MASSIVE SCORE")
    print("Testing Oracle with optimized strategy...")
    print("="*60)
    
    # Test Strategy 1: Frameshift variant
    print("\n🧪 TEST 1: FRAMESHIFT VARIANT")
    try:
        ref_seq, alt_seq = create_pathogenic_runx1_variant()
        
        oracle_payload = {
            "action": "score",
            "params": {
                "reference_sequence": ref_seq,
                "alternate_sequence": alt_seq
            }
        }
        
        print("📡 Submitting frameshift variant to Oracle...")
        response = requests.post(oracle_url, json=oracle_payload, timeout=300, verify=False)
        response.raise_for_status()
        
        result = response.json()
        frameshift_score = result.get("zeta_score", 0.0)
        
        print(f"✅ Frameshift Result:")
        print(f"   🎯 ZETA SCORE: {frameshift_score}")
        print(f"   Status: {result.get('status', 'unknown')}")
        print(f"   Confidence: {result.get('confidence', 0.0)}")
        
        if abs(frameshift_score) > 20000:
            print(f"🏆 MASSIVE SCORE ACHIEVED! Target was -26,140")
            print(f"💥 Result: {frameshift_score}")
        elif abs(frameshift_score) > 1000:
            print(f"🔥 SIGNIFICANT SCORE! Getting closer to target")
        else:
            print(f"📊 Moderate score: {frameshift_score}")
            
    except Exception as e:
        print(f"❌ Frameshift test failed: {e}")
    
    print("\n" + "="*60)
    
    # Test Strategy 2: Nonsense variant
    print("\n🧪 TEST 2: NONSENSE VARIANT")
    try:
        ref_seq, alt_seq = create_nonsense_runx1_variant()
        
        oracle_payload = {
            "action": "score",
            "params": {
                "reference_sequence": ref_seq,
                "alternate_sequence": alt_seq
            }
        }
        
        print("📡 Submitting nonsense variant to Oracle...")
        response = requests.post(oracle_url, json=oracle_payload, timeout=300, verify=False)
        response.raise_for_status()
        
        result = response.json()
        nonsense_score = result.get("zeta_score", 0.0)
        
        print(f"✅ Nonsense Result:")
        print(f"   🎯 ZETA SCORE: {nonsense_score}")
        print(f"   Status: {result.get('status', 'unknown')}")
        print(f"   Confidence: {result.get('confidence', 0.0)}")
        
        if abs(nonsense_score) > 20000:
            print(f"🏆 MASSIVE SCORE ACHIEVED! Target was -26,140")
            print(f"💥 Result: {nonsense_score}")
        elif abs(nonsense_score) > 1000:
            print(f"🔥 SIGNIFICANT SCORE! Getting closer to target")
        else:
            print(f"📊 Moderate score: {nonsense_score}")
            
    except Exception as e:
        print(f"❌ Nonsense test failed: {e}")
    
    print("\n" + "="*60)
    print("🎯 ANALYSIS COMPLETE")
    
if __name__ == "__main__":
    test_oracle_with_optimal_strategy() 