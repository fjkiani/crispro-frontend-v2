import argparse
import json
import os
import vcf # Use PyVCF3, which installs as 'vcf'
import pysam

# Placeholder for PyVCF and Pysam. We will add these once we install them.
# import vcf
# import pysam

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
            print(f"  ❌ An unexpected error occurred: {e}")
            return None

    def apply_variant(self, sequence, variant_info):
        """(Placeholder) Applies a variant to a reference sequence."""
        print(f"🔧 (Placeholder) Applying variant to sequence...")
        # TODO: Implement logic for SNVs, insertions, and deletions.
        return sequence # Placeholder

    def apply_variant_to_sequence(self, base_sequence, variant, flank_size=100):
        """
        Applies a single variant (SNP) to a base sequence to create the mutated version.
        Includes a flanking region around the mutation.
        
        NOTE: This initial implementation only handles SNPs.
        """
        if not variant.is_snp:
            print(f"  ⚠️ Skipping non-SNP variant at {variant.POS}. Implementation needed.")
            return None

        # Position is 1-based in VCF, need 0-based for string slicing
        variant_pos_0based = variant.POS - 1
        
        # Ensure we don't go out of bounds with the flanking regions
        start = max(0, variant_pos_0based - flank_size)
        end = min(len(base_sequence), variant_pos_0based + flank_size + 1)

        # Reconstruct the sequence with the variant
        ref_allele = variant.REF
        alt_allele = str(variant.ALT[0])

        # Basic check
        if base_sequence[variant_pos_0based:variant_pos_0based+len(ref_allele)] != ref_allele:
            print(f"  ❌ Mismatch! Reference allele in FASTA does not match VCF at {variant.POS}.")
            print(f"    FASTA: {base_sequence[variant_pos_0based:variant_pos_0based+len(ref_allele)]}")
            print(f"    VCF:   {ref_allele}")
            return None

        # Construct the mutated sequence segment
        mutated_segment = base_sequence[start:variant_pos_0based] + alt_allele + base_sequence[variant_pos_0based + len(ref_allele):end]
        
        print(f"  ✅ Applied SNP: {ref_allele} -> {alt_allele} at position {variant.POS}")
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

def main():
    """Main orchestration logic for the data ingestion pipeline."""
    parser = argparse.ArgumentParser(description="CrisPRO Data Ingestion Pipeline")
    parser.add_argument("--vcf", default="data/vcf/manual_runx1.vcf", help="Path to the VCF file.")
    parser.add_argument("--ref", default="data/reference/hg19.fa", help="Path to the reference genome FASTA.")
    args = parser.parse_args()

    print("🚀 Starting CrisPRO Data Ingestion Pipeline...")
    
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
    # Using coordinates from our previous search
    runx1_ref_seq = sequence_extractor.get_reference_sequence("chr21", 36160098, 36421599)

    # 4. Generate the mutated sequence for the primary germline variant
    # We pass the full gene sequence here
    runx1_mut_seq_germline = sequence_extractor.apply_variant_to_sequence(runx1_ref_seq, runx1_germline_variant)
    
    # 5. Generate mutated sequences for all somatic variants
    somatic_mut_sequences = []
    for svar in somatic_variants:
        # For somatic variants, we get a smaller window around them
        flank_size = 150
        somatic_ref_seq = sequence_extractor.get_reference_sequence(
            svar.CHROM, 
            svar.POS - flank_size - 1, # -1 for 0-based
            svar.POS + flank_size
        )
        mut_seq = sequence_extractor.apply_variant_to_sequence(somatic_ref_seq, svar, flank_size=flank_size)
        if mut_seq:
            somatic_mut_sequences.append({
                "id": f"{svar.CHROM}:{svar.POS}",
                "ref_seq": somatic_ref_seq,
                "mut_seq": mut_seq
            })


    # 6. Assemble the final payload for CommandCenter
    api_payload = {
        "patient_id": "FPD_001",
        "genes_of_interest": ["RUNX1", "GATA2", "FLT3"], # Example
        "germline_target": {
            "gene": "RUNX1",
            "variant": f"{runx1_germline_variant.CHROM}:{runx1_germline_variant.POS}",
            "reference_sequence": runx1_ref_seq,
            "mutated_sequence": runx1_mut_seq_germline,
        },
        "somatic_targets": somatic_mut_sequences
    }

    print("\n\n✅ Pipeline Complete. Final API Payload:")
    print("="*40)
    print(json.dumps(api_payload, indent=2))
    print("="*40)


if __name__ == "__main__":
    main() 