from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable

def apply_protein_change(dna_sequence: str, protein_change: str) -> str | None:
    """
    Applies a simple protein change (e.g., "R175H") to a DNA sequence.

    NOTE: This is a simplified implementation for demonstration. It finds the first
    matching codon and replaces it. It does not handle all edge cases.
    """
    if not protein_change or len(protein_change) < 3:
        return None

    ref_aa = protein_change[0]
    position = int(protein_change[1:-1])
    alt_aa = protein_change[-1]

    dna_seq_obj = Seq(dna_sequence)
    protein_seq = dna_seq_obj.translate(to_stop=True)

    if position > len(protein_seq):
        print(f"Error: Position {position} is out of bounds for protein of length {len(protein_seq)}")
        return None

    if protein_seq[position - 1] != ref_aa:
        print(f"Error: Reference AA mismatch at position {position}. Expected {ref_aa}, found {protein_seq[position-1]}")
        return None

    # Find the codon for the target amino acid
    codon_start = (position - 1) * 3
    codon_end = codon_start + 3
    original_codon = dna_sequence[codon_start:codon_end]

    # Find a new codon for the alternate amino acid
    standard_table = CodonTable.unambiguous_dna_by_id[1]
    alt_codons = [k for k, v in standard_table.forward_table.items() if v == alt_aa]
    
    if not alt_codons:
        print(f"Error: No codon found for alternate AA: {alt_aa}")
        return None

    # To be simple, we just pick the first possible codon.
    new_codon = alt_codons[0]

    # Create the new DNA sequence
    mutated_dna = dna_sequence[:codon_start] + new_codon + dna_sequence[codon_end:]
    
    return mutated_dna 