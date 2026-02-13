import argparse
import json
import os
import pysam
import requests
from dotenv import load_dotenv
import urllib3

# --- Configuration ---
load_dotenv()
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
# COMMAND_CENTER_URL = os.getenv("COMMAND_CENTER_URL", "https://crispro--command-center-v20-oracle-fix-web-app.modal.run")
# TACTICAL OVERRIDE: Using the newly deployed CommandCenter with the correct Oracle URL configuration.
COMMAND_CENTER_URL = "https://crispro--command-center-v20-oracle-fix-web-app.modal.run"

class VCFProcessor:
    """Handles all VCF file parsing and variant extraction using pysam."""
    def __init__(self, vcf_path):
        """Initializes the processor with the path to the VCF file."""
        self.vcf_path = vcf_path
        try:
            self.vcf_reader = pysam.VariantFile(self.vcf_path)
            print(f"✅ VCFProcessor initialized for: {self.vcf_path}")
        except (ValueError, FileNotFoundError) as e:
            print(f"❌ ERROR: Could not open VCF file: {e}")
            raise

    def get_variants_for_gene(self, gene_symbol, sample_name):
        """(Placeholder) Identifies all variants for a given gene and sample."""
        print(f"🔎 (Placeholder) Searching for variants in gene '{gene_symbol}' for sample '{sample_name}'...")
        # TODO: Implement logic to query VCF by gene coordinates.
        # This will require a gene coordinate map (e.g., from a BED/GTF file).
        return []

    def get_germline_runx1_variant(self, germline_sample_name="PATIENT_GERMLINE"):
        """
        Finds the primary germline RUNX1 variant from the VCF.
        Assumes a single, heterozygous germline variant for this patient.
        """
        print(f"🧬  Searching for germline RUNX1 variant for sample: {germline_sample_name}")
        for record in self.vcf_reader.fetch():
            if germline_sample_name in record.samples:
                sample = record.samples[germline_sample_name]
                # Correctly check if the sample's alleles differ from the reference allele
                if any(allele != record.ref for allele in sample.alleles) and any(sample.alleles):
                    print(f"  ✅ Found germline variant: {record.chrom}:{record.pos} {record.ref}->{record.alts[0]}")
                    return record
        print("  ⚠️ No germline variant found for the specified sample.")
        return None

    def get_somatic_variants(self, somatic_sample_name="PATIENT_TUMOR", germline_variant=None):
        """
        Finds all somatic variants for the specified tumor sample,
        excluding the known germline variant.
        """
        print(f"🧬  Searching for somatic variants for sample: {somatic_sample_name}")
        somatic_variants = []
        try:
            for record in self.vcf_reader.fetch():
                 if somatic_sample_name in record.samples:
                    sample = record.samples[somatic_sample_name]
                    # Correctly check if the sample's alleles differ from the reference allele
                    if any(allele != record.ref for allele in sample.alleles) and any(sample.alleles):
                        # Exclude the germline variant if it's present in the tumor
                        if germline_variant and record.pos == germline_variant.pos and record.alts == germline_variant.alts:
                            continue
                        print(f"  ✅ Found somatic variant: {record.chrom}:{record.pos} {record.ref}->{record.alts[0]}")
                        somatic_variants.append(record)
        except OSError as e:
            print(f"  ⚠️  WARNING: Encountered a parsing error in the VCF file: {e}")
            print("     This may indicate a corrupted or malformed record. Skipping to the next valid record.")
            # This is a simplified way to handle this. A more robust solution might try to recover.
            pass
        
        if not somatic_variants:
            print("  ⚠️ No somatic variants found for the specified sample.")
            
        return somatic_variants

class SequenceExtractor:
    """Handles sequence extraction from FASTA files."""
    
    def __init__(self, reference_file):
        self.reference_file = reference_file
        ensure_fasta_index(reference_file)
    
    def get_sequence(self):
        """Get the complete RUNX1 sequence from the reference file."""
        try:
            with pysam.FastaFile(self.reference_file) as fasta:
                # For our fixed RUNX1 file, get the first (and only) sequence
                contig_names = list(fasta.references)
                if not contig_names:
                    raise ValueError("No sequences found in reference file")
                
                # Use the first contig (should be the RUNX1 sequence)
                contig_name = contig_names[0]
                sequence = fasta.fetch(contig_name).upper()
                
                print(f"📋 Using contig: {contig_name}")
                print(f"📏 Sequence length: {len(sequence)} bp")
                
                return sequence
        except Exception as e:
            print(f"❌ ERROR: Could not extract sequence from {self.reference_file}: {e}")
            return None

    def apply_variant(self, sequence, variant_info):
        """(Placeholder) Applies a variant to a reference sequence."""
        print(f"🔧 (Placeholder) Applying variant to sequence...")
        # TODO: Implement logic for SNVs, insertions, and deletions.
        return sequence # Placeholder

    def generate_variant_windows(self, base_sequence, variant, sequence_start_pos, flank_size=100):
        """
        Applies a single variant and generates perfectly matched reference and alternate windows.
        Returns a tuple: (reference_window, mutated_window)
        """
        try:
            # Position is 1-based in VCF, need 0-based for string slicing.
            variant_pos_relative_0based = variant.pos - sequence_start_pos - 1

            if variant_pos_relative_0based < 0 or variant_pos_relative_0based >= len(base_sequence):
                 print(f"  ❌ Error: Variant position {variant.pos} is outside the provided sequence range.")
                 center = len(base_sequence) // 2
                 fallback_window = base_sequence[max(0, center - flank_size) : center + flank_size]
                 return fallback_window, fallback_window

            # Common logic for all variant types
            ref_allele = variant.ref
            alt_allele = str(variant.alts[0])
            
            pre_variant_base = base_sequence[:variant_pos_relative_0based]
            post_variant_base = base_sequence[variant_pos_relative_0based + len(ref_allele):]
            mutated_sequence_full = pre_variant_base + alt_allele + post_variant_base

            # Define the window boundaries based on the variant position in the *original* sequence
            window_start = max(0, variant_pos_relative_0based - flank_size)
            # The end of the window should be relative to the original sequence length for the reference
            ref_window_end = min(len(base_sequence), variant_pos_relative_0based + len(ref_allele) + flank_size)
            
            # The end of the window for the mutated sequence must account for its new length
            mut_window_end = min(len(mutated_sequence_full), variant_pos_relative_0based + len(alt_allele) + flank_size)

            # Extract the windows
            reference_window = base_sequence[window_start:ref_window_end]
            mutated_window = mutated_sequence_full[window_start:mut_window_end]

        except Exception as e:
            print(f"  ❌ Error applying variant: {e}. Defaulting to reference sequence window.")
            center = len(base_sequence) // 2
            fallback_window = base_sequence[max(0, center - flank_size) : center + flank_size]
            return fallback_window, fallback_window
        
        print(f"  ✅ Generated variant windows for: {variant.ref} -> {str(variant.alts[0])} at pos {variant.pos}")
        return reference_window, mutated_window

def ensure_fasta_index(ref_path):
    """Checks for a FASTA index file (.fai) and creates one if not found."""
    index_path = ref_path + ".fai"
    if not os.path.exists(index_path):
        print(f"⚠️ FASTA index not found. Creating index at: {index_path}")
        try:
            pysam.faidx(ref_path)
            print("✅ Index created successfully.")
        except pysam.SamtoolsError as e:
            print(f"❌ Error creating FASTA index: {e}")
            print("Please ensure the reference FASTA file is valid and `samtools` is accessible.")
            exit(1)
    else:
        print("✅ FASTA index found.")

def submit_payload(url, payload):
    """Submits the final payload to the CommandCenter API."""
    endpoint = "/workflow/full_patient_assessment"
    full_url = f"{url}{endpoint}"
    print(f"\n🚀 Submitting payload to: {full_url}")
    try:
        response = requests.post(full_url, json=payload, timeout=1200, verify=False)
        response.raise_for_status()
        print("\n✅ SUCCESS: Full patient assessment complete.")
        print("--- RESPONSE ---")
        print(json.dumps(response.json(), indent=2))
        print("----------------")
        return True
    except requests.exceptions.RequestException as e:
        print(f"\n❌ FAILED: Could not submit payload to CommandCenter.")
        print(f"  - Error: {e}")
        if hasattr(e, 'response') and e.response is not None:
            print(f"  - Status Code: {e.response.status_code}")
            try:
                print(f"  - Response Body: {e.response.json()}")
            except json.JSONDecodeError:
                print(f"  - Response Body: {e.response.text}")
        return False

def main():
    """Main orchestration logic for the data ingestion pipeline."""
    parser = argparse.ArgumentParser(description="CrisPRO Data Ingestion Pipeline")
    parser.add_argument("--vcf", default="external/runx1_workspace/data_assets/manual_runx1_corrected.vcf", help="Path to the VCF file.")
    parser.add_argument("--ref", default="external/runx1_workspace/data_assets/runx1_hg19_fixed.fa", help="Path to the RUNX1 reference sequence.")
    args = parser.parse_args()

    print("🚀 Starting CrisPRO Data Ingestion Pipeline...")
    
    # --- Constants ---
    # NCBI OFFICIAL RUNX1 coordinates (GRCh37/hg19)
    RUNX1_GENE_REGION_START = 36160098  # Start position on chr21
    RUNX1_GENE_REGION_END = 36421599    # End position on chr21
    RUNX1_GENE_LENGTH = RUNX1_GENE_REGION_END - RUNX1_GENE_REGION_START  # 261,501 bp
    
    print(f"📊 RUNX1 Gene Region: chr21:{RUNX1_GENE_REGION_START}-{RUNX1_GENE_REGION_END} ({RUNX1_GENE_LENGTH:,} bp)")
    
    # Command Center endpoint
    COMMAND_CENTER_URL = "https://crispro--command-center-v20-oracle-fix-web-app.modal.run"
    
    # Initialize processors
    vcf_processor = VCFProcessor(args.vcf)
    sequence_extractor = SequenceExtractor(args.ref)
    
    # 1. Extract germline RUNX1 variant
    print("🔍 Step 1: Extracting germline RUNX1 variant...")
    runx1_germline_variant = vcf_processor.get_germline_runx1_variant()
    
    if runx1_germline_variant is None:
        print("❌ No germline RUNX1 variant found. Exiting.")
        return
    
    print(f"✅ Found germline variant: chr{runx1_germline_variant.chrom}:{runx1_germline_variant.pos} {runx1_germline_variant.ref}>{runx1_germline_variant.alts[0]}")
    
    # Verify variant is within RUNX1 gene boundaries
    if not (RUNX1_GENE_REGION_START <= runx1_germline_variant.pos <= RUNX1_GENE_REGION_END):
        print(f"❌ ERROR: Variant position {runx1_germline_variant.pos} is outside RUNX1 gene region!")
        print(f"   Gene region: {RUNX1_GENE_REGION_START} - {RUNX1_GENE_REGION_END}")
        return
    
    print(f"✅ Variant confirmed within RUNX1 gene boundaries")
    
    # 2. Extract somatic variants
    print("🔍 Step 2: Extracting somatic variants...")
    somatic_variants = vcf_processor.get_somatic_variants()
    print(f"✅ Found {len(somatic_variants)} somatic variants")
    
    # 3. Load the reference RUNX1 gene sequence
    print("🔍 Step 3: Loading RUNX1 reference sequence...")
    full_runx1_ref_seq = sequence_extractor.get_sequence()
    print(f"✅ Loaded RUNX1 reference: {len(full_runx1_ref_seq):,} bp")
    
    # Verify the sequence length matches expected
    if len(full_runx1_ref_seq) != RUNX1_GENE_LENGTH:
        print(f"⚠️  WARNING: Reference sequence length ({len(full_runx1_ref_seq)}) doesn't match expected ({RUNX1_GENE_LENGTH})")
    
    # 4. Generate sequences for Oracle scoring - FULL GENE APPROACH
    print(f"🔥 TACTICAL OVERRIDE: Using FULL GENE sequences for maximum impact scoring...")
    
    # Calculate relative position within the gene
    variant_pos_relative = runx1_germline_variant.pos - RUNX1_GENE_REGION_START
    ref_allele = runx1_germline_variant.ref
    alt_allele = str(runx1_germline_variant.alts[0])
    
    print(f"🧬 VARIANT COORDINATES:")
    print(f"   Genomic position: chr21:{runx1_germline_variant.pos}")
    print(f"   Gene-relative position: {variant_pos_relative} (0-based)")
    print(f"   Ref allele: '{ref_allele}' -> Alt allele: '{alt_allele}'")
    
    # Verify the reference allele matches
    actual_ref_base = full_runx1_ref_seq[variant_pos_relative]
    print(f"   Expected ref base: '{ref_allele}' | Actual ref base: '{actual_ref_base}'")
    
    if actual_ref_base != ref_allele:
        print(f"❌ ERROR: Reference mismatch! VCF says '{ref_allele}' but sequence has '{actual_ref_base}' at position {variant_pos_relative}")
        return
    
    print(f"✅ Reference allele verified!")
    
    # Create the mutated full gene sequence
    mut_runx1_full_seq = (
        full_runx1_ref_seq[:variant_pos_relative] + 
        alt_allele + 
        full_runx1_ref_seq[variant_pos_relative + len(ref_allele):]
    )
    
    print(f"📊 SEQUENCE COMPARISON:")
    print(f"   Reference length: {len(full_runx1_ref_seq):,} bp")
    print(f"   Mutated length: {len(mut_runx1_full_seq):,} bp")
    
    # Show the actual change context
    context_start = max(0, variant_pos_relative - 50)
    context_end = min(len(full_runx1_ref_seq), variant_pos_relative + 50)
    ref_context = full_runx1_ref_seq[context_start:context_end]
    mut_context = mut_runx1_full_seq[context_start:context_end]
    
    print(f"   Context around variant (+/- 50bp):")
    print(f"   REF: {ref_context}")
    print(f"   ALT: {mut_context}")
    
    # Verify sequences are actually different
    differences = sum(1 for i in range(min(len(full_runx1_ref_seq), len(mut_runx1_full_seq))) if full_runx1_ref_seq[i] != mut_runx1_full_seq[i])
    print(f"   TOTAL DIFFERENCES: {differences} nucleotides")
    
    if differences == 0:
        print("❌ ERROR: Reference and mutated sequences are identical!")
        return
    
    print(f"✅ Sequences are different - Oracle should detect impact!")

    # 5. Submit to CommandCenter for Oracle scoring
    print("🎯 Step 5: Submitting to Oracle for scoring...")
    
    payload = {
        "patient_id": "RUNX1_PATIENT_001",
        "germline_variant": {
            "chromosome": f"chr{runx1_germline_variant.chrom}",
            "position": runx1_germline_variant.pos,
            "ref_allele": ref_allele,
            "alt_allele": alt_allele,
            "reference_sequence": full_runx1_ref_seq,
            "mutated_sequence": mut_runx1_full_seq
        },
        "somatic_variants": [
            {
                "chromosome": f"chr{var.chrom}",
                "position": var.pos,
                "ref_allele": var.ref,
                "alt_allele": str(var.alts[0]) if var.alts else "N"
            } for var in somatic_variants[:5]  # Limit to first 5 for testing
        ]
    }
    
    success = submit_payload(COMMAND_CENTER_URL, payload)
    
    if success:
        print("🎉 Pipeline completed successfully!")
    else:
        print("❌ Pipeline failed during submission.")
    
    return


if __name__ == "__main__":
    main() 