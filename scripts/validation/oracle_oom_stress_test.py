import os
import sys
import httpx
import pysam
from loguru import logger

# --- DOCTRINAL SETUP ---
# Ensure we can import from the project root
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

def setup_config():
    """Setup configuration for the stress test"""
    print("--- [STRESS TEST CONFIG] Initializing CUDA OOM Stress Test ---\n")
    
    config = {
        "REFERENCE_GENOME_PATH": "data/gene_database/reference/hg19.fa",
        "ORACLE_URL": "https://crispro--zeta-oracle-v3-log-fix-zetaoracle-api.modal.run/invoke",
        "RUNX1_CHROM": "chr21",
        "RUNX1_GENE_REGION_START": 36160098,
        "RUNX1_GENE_REGION_END": 36421599,
    }
    
    print("✅ Configuration loaded for OOM stress test.")
    print(f"   Oracle Target: {config['ORACLE_URL']}")
    print(f"   Gene Region: {config['RUNX1_CHROM']}:{config['RUNX1_GENE_REGION_START']}-{config['RUNX1_GENE_REGION_END']}")
    print()
    
    return config

def acquire_full_sequence(config):
    """Acquire the complete RUNX1 gene sequence"""
    print("--- [SEQUENCE ACQUISITION] Loading full RUNX1 gene sequence ---\n")
    
    try:
        fasta_handle = pysam.FastaFile(config["REFERENCE_GENOME_PATH"])
        sequence = fasta_handle.fetch(
            config["RUNX1_CHROM"], 
            config["RUNX1_GENE_REGION_START"], 
            config["RUNX1_GENE_REGION_END"]
        )
        fasta_handle.close()
        
        print(f"✅ Successfully loaded RUNX1 gene sequence")
        print(f"   Length: {len(sequence):,} base pairs")
        print(f"   Size: ~{len(sequence) * 2 / 1024:.1f} KB (as UTF-8 strings)")
        print(f"   Expected GPU Memory: ~{len(sequence) * 32 / 1024 / 1024:.1f} MB (float32 embeddings)")
        print()
        
        return sequence.upper()
    except Exception as e:
        print(f"💥 ERROR: Failed to fetch sequence from FASTA file")
        print(f"   Path: {config['REFERENCE_GENOME_PATH']}")
        print(f"   Details: {e}")
        return None

def forge_mutation(wild_type_dna):
    """Create a single-base mutation in the middle of the sequence"""
    print("--- [MUTATION FORGE] Creating single-base mutation ---\n")
    
    # Place mutation at the center of the sequence
    mutation_pos = len(wild_type_dna) // 2
    original_allele = wild_type_dna[mutation_pos]
    alt_allele = 'A' if original_allele != 'A' else 'G'
    
    # Create mutated sequence
    mutated_dna_list = list(wild_type_dna)
    mutated_dna_list[mutation_pos] = alt_allele
    mutated_dna = "".join(mutated_dna_list)
    
    print(f"✅ Mutation forged at position {mutation_pos:,}")
    print(f"   Change: {original_allele} -> {alt_allele}")
    print(f"   Mutated sequence length: {len(mutated_dna):,}bp")
    print()
    
    return mutated_dna, mutation_pos, original_allele, alt_allele

def invoke_oracle_stress_test(config, ref_sequence, alt_sequence, mutation_info):
    """Send the full sequences to the Oracle for stress testing"""
    print("--- [ORACLE STRESS TEST] Launching full-sequence payload ---\n")
    
    mutation_pos, original_allele, alt_allele = mutation_info
    
    print(f"🚀 Preparing to fire massive payload at Oracle:")
    print(f"   Reference sequence: {len(ref_sequence):,}bp")
    print(f"   Alternate sequence: {len(alt_sequence):,}bp") 
    print(f"   Total payload size: ~{(len(ref_sequence) + len(alt_sequence)) * 2 / 1024:.1f} KB")
    print(f"   Mutation: pos {mutation_pos:,} ({original_allele} -> {alt_allele})")
    print()
    
    payload = {
        "action": "score",
        "params": {
            "reference_sequence": ref_sequence,
            "alternate_sequence": alt_sequence
        }
    }
    
    print("💥 FIRING MASSIVE PAYLOAD AT ORACLE...")
    print("   (This will likely trigger CUDA out-of-memory errors)")
    print()
    
    try:
        with httpx.Client(timeout=600.0, follow_redirects=True) as client:  # Extended timeout for large payload
            response = client.post(config["ORACLE_URL"], json=payload)
            response.raise_for_status()
            result = response.json()
            delta_score = result.get("zeta_score", 0.0)
            
            print("🎉 IMPOSSIBLE SUCCESS!")
            print(f"   Oracle processed the full sequence without OOM!")
            print(f"   Delta Score: {delta_score:.6f}")
            print(f"   This indicates the Oracle has been successfully optimized.")
            return delta_score
            
    except httpx.HTTPStatusError as e:
        print(f"💥 HTTP ERROR: Status {e.response.status_code}")
        response_text = e.response.text
        print(f"   Response: {response_text}")
        
        if "out of memory" in response_text.lower() or "oom" in response_text.lower():
            print("✅ EXPECTED FAILURE: CUDA out-of-memory error detected.")
            print("   This confirms the Oracle requires memory optimization for large sequences.")
        elif "tuple indices must be integers" in response_text:
            print("✅ TUPLE ERROR: The original bug is still present.")
            print("   This indicates the Oracle deployment is not using the patched code.")
        else:
            print("❓ UNEXPECTED ERROR: Different failure mode detected.")
            
        return None
        
    except Exception as e:
        print(f"💥 CLIENT ERROR: {e}")
        return None

def main():
    """Execute the full CUDA OOM stress test"""
    print("=" * 70)
    print("🧬 ZETA ORACLE CUDA OOM STRESS TEST")
    print("=" * 70)
    print("This test sends the complete RUNX1 gene sequence (~261kb) to the Oracle")
    print("to trigger CUDA out-of-memory conditions and test memory optimization.")
    print("=" * 70)
    print()
    
    # Phase 1: Setup
    config = setup_config()
    
    # Phase 2: Sequence acquisition
    wild_type_dna = acquire_full_sequence(config)
    if not wild_type_dna:
        print("💥 STRESS TEST ABORTED: Failed to acquire sequence.")
        return
    
    # Phase 3: Mutation forging
    mutated_dna, mutation_pos, original_allele, alt_allele = forge_mutation(wild_type_dna)
    mutation_info = (mutation_pos, original_allele, alt_allele)
    
    # Phase 4: Oracle stress test
    result = invoke_oracle_stress_test(config, wild_type_dna, mutated_dna, mutation_info)
    
    # Phase 5: Results
    print("=" * 70)
    print("🧬 STRESS TEST COMPLETE")
    print("=" * 70)
    if result is not None:
        print(f"Final Result: SUCCESS (Delta Score: {result:.6f})")
        print("The Oracle successfully processed the full sequence without memory errors.")
    else:
        print("Final Result: EXPECTED FAILURE")
        print("The Oracle failed as expected due to memory constraints or other issues.")
        print("This data can be used to optimize the Oracle's memory handling.")
    print("=" * 70)

if __name__ == "__main__":
    main() 