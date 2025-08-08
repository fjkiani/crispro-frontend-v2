import pandas as pd
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_vcf(vcf_path):
    """
    Parses a VCF file to extract variant information.

    This function reads a VCF file, skipping the header lines, and returns a 
    structured list of variants. It currently extracts basic information 
    (chromosome, position, ID, reference and alternate alleles).

    Args:
        vcf_path (str): The file path to the VCF file.

    Returns:
        list: A list of dictionaries, where each dictionary represents a variant.
              Returns an empty list if the file cannot be parsed.
    """
    variants = []
    try:
        with open(vcf_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\\t')
                if len(fields) >= 5:
                    variant_info = {
                        'chr': fields[0],
                        'pos': int(fields[1]),
                        'id': fields[2],
                        'ref': fields[3],
                        'alt': fields[4]
                        # Additional fields like INFO can be parsed here as needed
                    }
                    variants.append(variant_info)
        logging.info(f"Successfully parsed {len(variants)} variants from {vcf_path}")
        return variants
    except FileNotFoundError:
        logging.error(f"Error: The file was not found at {vcf_path}")
        return []
    except Exception as e:
        logging.error(f"An error occurred while parsing the VCF file: {e}")
        return []

# Example Usage
if __name__ == '__main__':
    # To run this example, create a dummy VCF file named 'sample.vcf'
    # with content like:
    # ##fileformat=VCFv4.2
    # #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    # chr1	12345	rs123	A	G	.	.	.
    # chr2	54321	.	C	T	.	.	.
    
    # Create a dummy VCF for testing purposes
    try:
        with open("sample.vcf", "w") as f:
            f.write("##fileformat=VCFv4.2\\n")
            f.write("#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n")
            f.write("chr1\\t1000\\trs123\\tA\\tG\\t.\\t.\\t.\\n")
            f.write("chrX\\t2000\\trs456\\tC\\tT,A\\t.\\t.\\t.\\n")

        parsed_variants = parse_vcf('sample.vcf')
        if parsed_variants:
            print("\\nParsed Variants:")
            for variant in parsed_variants:
                print(variant)
    except Exception as e:
        logging.error(f"Failed to create or parse dummy VCF file: {e}") 