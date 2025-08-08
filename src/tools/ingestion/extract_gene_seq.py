import pysam
import argparse

def extract_sequence(ref_path, region, output_path):
    """
    Extracts a specific genomic region from a FASTA file.
    """
    try:
        with pysam.FastaFile(ref_path) as fasta:
            seq = fasta.fetch(region=region)
            with open(output_path, "w") as out_file:
                out_file.write(f">{region}\\n")
                out_file.write(seq + "\\n")
        print(f"✅ Successfully extracted region '{region}' to '{output_path}'")
    except Exception as e:
        print(f"❌ Error extracting sequence: {e}")
        exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract a genomic sequence from a reference FASTA.")
    parser.add_argument("--ref", required=True, help="Path to the reference FASTA file.")
    parser.add_argument("--region", required=True, help="Genomic region to extract (e.g., 'chr21:1000-2000').")
    parser.add_argument("--out", required=True, help="Path to the output FASTA file.")
    args = parser.parse_args()

    extract_sequence(args.ref, args.region, args.out) 