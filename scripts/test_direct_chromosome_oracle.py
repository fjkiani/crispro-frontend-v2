#!/usr/bin/env python3
"""
Direct chromosome-based Oracle test for RUNX1 variant.
This bypasses the reference file issues by working directly with the chromosome.
"""

import pysam
import requests
import json

def test_runx1_oracle_direct():
    """Test Oracle with direct chromosome extraction."""
    
    print("🎯 OPERATION: DIRECT CHROMOSOME ORACLE TEST")
    print("Working directly with chromosome coordinates...")
    
    # NCBI verified coordinates
    GENE_START = 36160098
    GENE_END = 36421599
    VCF_POS = 36250941
    
    # Direct chromosome access
    ref = pysam.FastaFile('data/gene_database/reference/hg19.fa')
    
    # Extract the full RUNX1 gene from chromosome
    print(f"📊 Extracting RUNX1 from chr21:{GENE_START}-{GENE_END}")
    full_gene_ref = ref.fetch('chr21', GENE_START, GENE_END).upper()
    
    # Calculate relative position for the variant
    relative_pos = VCF_POS - GENE_START - 1  # Convert to 0-based indexing
    
    print(f"🧬 VARIANT DETAILS:")
    print(f"   VCF position: {VCF_POS}")
    print(f"   Relative position: {relative_pos}")
    
    # Get the actual base at this position
    actual_ref_base = full_gene_ref[relative_pos]
    print(f"   Actual ref base: '{actual_ref_base}'")
    
    # Create a T->C mutation (actual reference is T)
    if actual_ref_base == 'T':
        alt_allele = 'C'
        print(f"   Creating mutation: T->C")
    elif actual_ref_base == 'C':
        alt_allele = 'T'
        print(f"   Creating mutation: C->T")
    else:
        # Default to creating an A->G mutation
        alt_allele = 'G' if actual_ref_base == 'A' else 'A'
        print(f"   Creating mutation: {actual_ref_base}->{alt_allele}")
    
    # Create mutated sequence
    full_gene_alt = (
        full_gene_ref[:relative_pos] + 
        alt_allele + 
        full_gene_ref[relative_pos + 1:]
    )
    
    print(f"📏 SEQUENCES:")
    print(f"   Reference length: {len(full_gene_ref):,} bp")
    print(f"   Mutated length: {len(full_gene_alt):,} bp")
    print(f"   Differences: {sum(1 for i in range(len(full_gene_ref)) if full_gene_ref[i] != full_gene_alt[i])}")
    
    # Show context
    context_start = max(0, relative_pos - 50)
    context_end = min(len(full_gene_ref), relative_pos + 50)
    ref_context = full_gene_ref[context_start:context_end]
    alt_context = full_gene_alt[context_start:context_end]
    
    print(f"   Context (+/- 50bp):")
    print(f"   REF: {ref_context}")
    print(f"   ALT: {alt_context}")
    
    # Test with Oracle
    oracle_url = "https://crispro--zeta-oracle-v5-final-fix-zetaoracle-api.modal.run/invoke"
    
    oracle_payload = {
        "action": "score",
        "params": {
            "reference_sequence": full_gene_ref,
            "alternate_sequence": full_gene_alt
        }
    }
    
    print(f"🚀 Submitting to Oracle...")
    try:
        response = requests.post(oracle_url, json=oracle_payload, timeout=1200, verify=False)
        response.raise_for_status()
        
        result = response.json()
        print(f"✅ Oracle Response:")
        print(json.dumps(result, indent=2))
        
        if 'zeta_score' in result:
            score = result['zeta_score']
            print(f"\n🎯 FINAL SCORE: {score}")
            if abs(score) > 1000:
                print(f"🏆 MASSIVE SCORE ACHIEVED!")
            else:
                print(f"📊 Moderate score - may need larger changes")
        
    except Exception as e:
        print(f"❌ Oracle request failed: {e}")
        return False
    
    return True

if __name__ == "__main__":
    test_runx1_oracle_direct() 