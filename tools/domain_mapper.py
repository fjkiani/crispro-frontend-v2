import sys
import os

# Ensure the parent directory is in the Python path to allow for relative imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tools.coordinate_handler import GeneCoordinateConverter

def map_protein_domain_to_genomic_coords(gene_symbol, transcript_id, domain_name, domain_aa_start, domain_aa_end, converter):
    """
    Maps a protein domain's amino acid coordinates to genomic coordinates.

    Args:
        gene_symbol (str): The symbol of the gene (e.g., "ASXL1").
        transcript_id (str): The RefSeq ID of the specific transcript to use.
        domain_name (str): The name of the protein domain.
        domain_aa_start (int): The start amino acid position of the domain.
        domain_aa_end (int): The end amino acid position of the domain.
        converter (GeneCoordinateConverter): An instance of the coordinate converter.

    Returns:
        A dictionary with the genomic coordinates or an error message.
    """
    exons = converter.get_exons(gene_symbol, transcript_id=transcript_id)
    if not exons:
        return {"error": f"Could not retrieve exons for {gene_symbol} transcript {transcript_id}."}

    # The canonical transcript for ASXL1 (NP_056153.2) is on the reverse strand.
    # We must sort the exons by their genomic start position.
    # For a gene on the reverse strand, transcription proceeds from higher to lower coordinates.
    # The CDS start is at the end of the last exon, and the CDS end is at the start of the first exon.
    strand = exons[0].strand
    if strand == '-':
        exons.sort(key=lambda e: e.start, reverse=True)
    else:
        exons.sort(key=lambda e: e.start)

    cds_len_bp = 0
    genomic_pos_map = {}  # Map from cumulative CDS position (in bp) to genomic position

    for exon in exons:
        for i in range(exon.length):
            cds_len_bp += 1
            if strand == '-':
                # On the reverse strand, we count down from the exon end
                genomic_pos = exon.end - i
            else:
                # On the forward strand, we count up from the exon start
                genomic_pos = exon.start + i
            genomic_pos_map[cds_len_bp] = genomic_pos

    # Convert amino acid positions to 0-indexed base positions in the CDS
    start_bp_in_cds = (domain_aa_start - 1) * 3 + 1
    end_bp_in_cds = domain_aa_end * 3

    # Find the corresponding genomic positions
    start_genomic_pos = genomic_pos_map.get(start_bp_in_cds)
    end_genomic_pos = genomic_pos_map.get(end_bp_in_cds)

    if start_genomic_pos is None or end_genomic_pos is None:
        return {"error": f"Could not map domain {domain_name} to genomic coordinates. Check AA coordinates."}

    # The "start" of the domain on the genome might be a larger number if on the reverse strand
    final_start = min(start_genomic_pos, end_genomic_pos)
    final_end = max(start_genomic_pos, end_genomic_pos)

    return {
        "domain": domain_name,
        "chromosome": exons[0].chromosome,
        "start": final_start,
        "end": final_end,
        "strand": strand,
    }

if __name__ == "__main__":
    # --- Configuration ---
    REFGENE_PATH = "data/refGene.txt"
    GENOME_ASSEMBLY = "hg19"
    GENE_SYMBOL = "ASXL1"
    # Canonical transcript for ASXL1 corresponding to NP_056153.2 is NM_015338
    TRANSCRIPT_ID = "NM_015338" 

    # ASXL1 functional domains (from UniProt/NCBI for NP_056153.2)
    DOMAINS = [
        {"name": "HARE-HTH", "start_aa": 11, "end_aa": 83},
        {"name": "ASXH", "start_aa": 236, "end_aa": 359},
        {"name": "PHD_3", "start_aa": 1506, "end_aa": 1539},
    ]

    # --- Execution ---
    print(f"Loading gene annotations from: {REFGENE_PATH}")
    converter = GeneCoordinateConverter(refgene_file=REFGENE_PATH, assembly=GENOME_ASSEMBLY)

    gene_region = converter.gene_to_region(GENE_SYMBOL, transcript_id=TRANSCRIPT_ID)
    if gene_region:
        print(f"\\n--- Gene Region for {GENE_SYMBOL} (Transcript: {TRANSCRIPT_ID}) ---")
        print(f"Chromosome: {gene_region.chromosome}, Start: {gene_region.start}, End: {gene_region.end}, Strand: {gene_region.strand}")
    else:
        print(f"Could not find region for gene {GENE_SYMBOL}")

    print(f"\\n--- Domain Mapping for {GENE_SYMBOL} (Transcript: {TRANSCRIPT_ID}) ---")
    for domain in DOMAINS:
        result = map_protein_domain_to_genomic_coords(
            GENE_SYMBOL,
            TRANSCRIPT_ID,
            domain["name"],
            domain["start_aa"],
            domain["end_aa"],
            converter
        )
        print(result) 