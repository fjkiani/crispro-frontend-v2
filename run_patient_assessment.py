import argparse
import json
import os
import vcf # Use PyVCF3, which installs as 'vcf'
import pysam
import requests
from dotenv import load_dotenv
import urllib3

# Placeholder for PyVCF and Pysam. We will add these once we install them.
# import vcf
# import pysam

# --- Configuration ---
load_dotenv()
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
COMMAND_CENTER_URL = os.getenv("COMMAND_CENTER_URL", "https://crispro--command-center-commandcenter-api.modal.run")

class VCFProcessor:
    """Handles all VCF file parsing and variant extraction."""
    def __init__(self, vcf_path):
        """Initializes the processor with the path to the VCF file."""
        self.vcf_path = vcf_path
        self.vcf_reader = vcf.Reader(filename=self.vcf_path)
        print(f"✅ VCFProcessor initialized for: {self.vcf_path}")

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
        for record in self.vcf_reader:
            # We are assuming the first variant for the germline sample is the one we want.
            # In a real scenario, we'd have more sophisticated filtering (e.g., by impact).
            sample = record.genotype(germline_sample_name)
            if sample and sample.is_variant:
                print(f"  ✅ Found germline variant: {record.CHROM}:{record.POS} {record.REF}->{record.ALT[0]}")
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
        for record in self.vcf_reader:
            sample = record.genotype(somatic_sample_name)
            if sample and sample.is_variant:
                # Exclude the germline variant if it's present in the tumor
                if germline_variant and record.POS == germline_variant.POS and record.ALT == germline_variant.ALT:
                    continue
                print(f"  ✅ Found somatic variant: {record.CHROM}:{record.POS} {record.REF}->{record.ALT[0]}")
                somatic_variants.append(record)
        
        if not somatic_variants:
            print("  ⚠️ No somatic variants found for the specified sample.")
            
        return somatic_variants

class SequenceExtractor:
    """Fetches reference DNA sequences and constructs mutated sequences."""
    def __init__(self, ref_genome_path):
        """Initializes the extractor with the path to the reference FASTA file."""
        self.ref_genome = pysam.FastaFile(ref_genome_path)
        print(f"✅ SequenceExtractor initialized with reference: {ref_genome_path}")

    def get_reference_sequence(self, chrom, start, end):
        """Fetches a DNA sequence from the reference genome for a given region."""
        print(f"🧬  Fetching reference sequence for {chrom}:{start}-{end}")
        try:
            return self.ref_genome.fetch(chrom, start, end).upper()
        except KeyError:
            print(f"  ❌ Error: Chromosome '{chrom}' not found in reference FASTA.")
            return None
        except Exception as e:
            print(f"  ❌ Error fetching reference sequence: {e}")
            return None

    def apply_variant(self, sequence, variant_info):
        """(Placeholder) Applies a variant to a reference sequence."""
        print(f"🔧 (Placeholder) Applying variant to sequence...")
        # TODO: Implement logic for SNVs, insertions, and deletions.
        return sequence # Placeholder

    def apply_variant_to_sequence(self, base_sequence, variant, sequence_start_pos, flank_size=100):
        """
        Applies a single variant (SNP or INDEL) to a base sequence to create the mutated version.
        Includes a flanking region around the mutation.
        
        `sequence_start_pos` is the 0-based genomic coordinate of the start of `base_sequence`.
        
        NOTE: This implementation is simplified for demonstration.
        """
        try:
            # Position is 1-based in VCF, need 0-based for string slicing.
            # Calculate the variant's position relative to the start of the provided sequence.
            variant_pos_relative_0based = variant.POS - sequence_start_pos - 1

            if variant_pos_relative_0based < 0 or variant_pos_relative_0based >= len(base_sequence):
                 print(f"  ❌ Error: Variant position {variant.POS} is outside the provided sequence range.")
                 # Fallback: return the central part of the original sequence
                 center = len(base_sequence) // 2
                 return base_sequence[max(0, center - flank_size) : center + flank_size]

            if variant.is_snp:
                # The end position needs to account for the length of the reference allele
                ref_allele = variant.REF
                alt_allele = str(variant.ALT[0])

                # Construct the mutated sequence segment
                pre_variant = base_sequence[:variant_pos_relative_0based]
                post_variant = base_sequence[variant_pos_relative_0based + len(ref_allele):]
                mutated_sequence = pre_variant + alt_allele + post_variant
                
                # Extract the final segment with flanking regions centered on the mutation
                final_start = max(0, variant_pos_relative_0based - flank_size)
                final_end = min(len(mutated_sequence), variant_pos_relative_0based + len(alt_allele) + flank_size)
                mutated_segment = mutated_sequence[final_start:final_end]

            elif variant.is_indel:
                ref_allele = variant.REF
                alt_allele = str(variant.ALT[0])
                
                # Construct the mutated sequence with the indel
                pre_variant = base_sequence[:variant_pos_relative_0based]
                post_variant = base_sequence[variant_pos_relative_0based + len(ref_allele):]
                mutated_sequence = pre_variant + alt_allele + post_variant

                # Extract the final segment with flanking regions
                final_start = max(0, variant_pos_relative_0based - flank_size)
                final_end = min(len(mutated_sequence), variant_pos_relative_0based + len(alt_allele) + flank_size)
                mutated_segment = mutated_sequence[final_start:final_end]

            else:
                print(f"  ⚠️ Skipping non-SNP/non-INDEL variant at {variant.POS}. Type: {variant.var_type}")
                center = len(base_sequence) // 2
                return base_sequence[max(0, center - flank_size) : center + flank_size]

        except Exception as e:
            print(f"  ❌ Error applying variant: {e}. Defaulting to reference sequence.")
            # Fallback for any other errors
            center = len(base_sequence) // 2
            return base_sequence[max(0, center - flank_size) : center + flank_size]
        
        print(f"  ✅ Applied variant: {variant.REF} -> {str(variant.ALT[0])} at position {variant.POS}")
        return mutated_segment

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

def submit_payload(payload):
    """Submits the final payload to the CommandCenter API."""
    endpoint = "/workflow/full_patient_assessment"
    url = f"{COMMAND_CENTER_URL}{endpoint}"
    print(f"\n🚀 Submitting payload to: {url}")
    try:
        response = requests.post(url, json=payload, timeout=1200, verify=False)
        response.raise_for_status()
        print("\n✅ SUCCESS: Full patient assessment complete.")
        print("--- RESPONSE ---")
        print(json.dumps(response.json(), indent=2))
        print("----------------")
    except requests.exceptions.RequestException as e:
        print(f"\n❌ FAILED: Could not submit payload to CommandCenter.")
        print(f"  - Error: {e}")
        if e.response is not None:
            print(f"  - Status Code: {e.response.status_code}")
            try:
                print(f"  - Response Body: {e.response.json()}")
            except json.JSONDecodeError:
                print(f"  - Response Body: {e.response.text}")

def main():
    """Main orchestration logic for the data ingestion pipeline."""
    parser = argparse.ArgumentParser(description="CrisPRO Data Ingestion Pipeline")
    parser.add_argument("--vcf", default="data/vcf/manual_runx1.vcf", help="Path to the VCF file.")
    parser.add_argument("--ref", default="data/reference/hg19.fa", help="Path to the reference genome FASTA.")
    args = parser.parse_args()

    print("🚀 Starting CrisPRO Data Ingestion Pipeline...")
    
    # --- Constants ---
    RUNX1_GENE_REGION_START = 36160098
    RUNX1_GENE_REGION_END = 36421599
    SOMATIC_FLANK_SIZE = 200 # Provides a 401bp window for gRNA design

    # Ensure the reference genome is indexed.
    ensure_fasta_index(args.ref)

    # 1. Initialize processors
    vcf_processor = VCFProcessor(args.vcf)
    sequence_extractor = SequenceExtractor(args.ref)

    # 2. Extract variants from VCF
    runx1_germline_variant = vcf_processor.get_germline_runx1_variant()
    if not runx1_germline_variant:
        print("❌ Critical error: No germline RUNX1 variant found. Aborting.")
        return

    somatic_variants = vcf_processor.get_somatic_variants(germline_variant=runx1_germline_variant)

    # 3. Get the full reference sequence for the RUNX1 gene
    full_runx1_ref_seq = sequence_extractor.get_reference_sequence(
        runx1_germline_variant.CHROM,
        RUNX1_GENE_REGION_START,
        RUNX1_GENE_REGION_END
    )
    if not full_runx1_ref_seq:
        print("❌ CRITICAL: Could not fetch RUNX1 reference sequence. Aborting.")
        exit(1)

    # 4. Generate the mutated sequence for the germline variant
    mutated_runx1_seq = sequence_extractor.apply_variant_to_sequence(
        full_runx1_ref_seq,
        runx1_germline_variant,
        sequence_start_pos=RUNX1_GENE_REGION_START,
        flank_size=200 # Use a large flank for the main gene correction
    )
    if not mutated_runx1_seq:
        print("❌ CRITICAL: Could not generate mutated RUNX1 sequence. Aborting.")
        exit(1)

    # 5. Process somatic variants
    somatic_targets = []
    for record in somatic_variants:
        # For each somatic variant, get a large flanking region for context
        somatic_region_start = record.POS - SOMATIC_FLANK_SIZE - 1
        somatic_region_end = record.POS + SOMATIC_FLANK_SIZE
        somatic_ref_seq = sequence_extractor.get_reference_sequence(
            record.CHROM,
            somatic_region_start,
            somatic_region_end
        )
        if not somatic_ref_seq:
            print(f"  ⚠️ Could not fetch reference for somatic variant at {record.POS}. Skipping.")
            continue

        mutated_somatic_seq = sequence_extractor.apply_variant_to_sequence(
            somatic_ref_seq,
            record,
            sequence_start_pos=somatic_region_start,
            flank_size=100 # This flank is for the final payload segment
        )
        
        if mutated_somatic_seq:
            somatic_targets.append({
                "id": f"{record.CHROM}:{record.POS}",
                "ref_seq": somatic_ref_seq,
                "mut_seq": mutated_somatic_seq
            })

    # 6. Assemble final API payload
    api_payload = {
        "patient_id": "FPD_001",
        "genes_of_interest": ["RUNX1", "GATA2", "FLT3"], # Example
        "germline_target": {
            "gene": "RUNX1",
            "variant": f"{runx1_germline_variant.CHROM}:{runx1_germline_variant.POS}",
            "reference_sequence": full_runx1_ref_seq,
            "mutated_sequence": mutated_runx1_seq,
        },
        "somatic_targets": somatic_targets
    }

    print("\n\n✅ Pipeline Complete. Final API Payload:")
    print("="*40)
    print(json.dumps(api_payload, indent=2))
    print("="*40)

    submit_payload(api_payload)


if __name__ == "__main__":
    main() 