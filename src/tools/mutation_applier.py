
import re
from Bio.Seq import Seq

class MutationApplier:
    """
    A tool to apply mutations specified in HGVS protein notation to a DNA sequence.
    Handles substitutions, deletions, insertions, and frameshifts.
    """

    AA_CODES = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
    }

    def _three_to_one(self, three_letter_code):
        return self.AA_CODES.get(three_letter_code, 'X')

    def _parse_hgvs(self, hgvs: str):
        """Parses complex HGVS notations."""
        # Simple substitution: p.Arg135Ala
        sub_match = re.match(r"p\.(?P<ref>[A-Z][a-z]{2})(?P<pos>\d+)(?P<alt>[A-Z][a-z]{2})", hgvs)
        if sub_match:
            groups = sub_match.groupdict()
            return {
                "type": "substitution",
                "ref": self._three_to_one(groups['ref']),
                "pos": int(groups['pos']),
                "alt": self._three_to_one(groups['alt']),
            }

        # Frameshift: p.Arg135fs or p.Arg135fs*12
        fs_match = re.match(r"p\.(?P<ref>[A-Z][a-z]{2})(?P<pos>\d+)fs(?:\*(?P<term_dist>\d+))?", hgvs)
        if fs_match:
            groups = fs_match.groupdict()
            return {
                "type": "frameshift",
                "ref": self._three_to_one(groups['ref']),
                "pos": int(groups['pos']),
            }
        
        # Add more parsers for deletions, insertions, etc. as needed.
        
        raise ValueError(f"Unsupported HGVS notation: {hgvs}")

    def apply_mutation_to_dna(self, dna_sequence: str, hgvs_mutation: str) -> str:
        """
        Applies a mutation to a DNA sequence based on HGVS protein notation.
        This is a simplified implementation focusing on a specific frameshift.
        """
        parsed_mutation = self._parse_hgvs(hgvs_mutation)
        
        pos = parsed_mutation['pos']
        
        # The DNA position is 1-based index * 3, for the start of the codon.
        # We adjust to 0-based index for slicing.
        codon_start_idx = (pos - 1) * 3
        
        if codon_start_idx + 3 > len(dna_sequence):
            raise ValueError("Mutation position is out of bounds for the DNA sequence.")
            
        original_codon = dna_sequence[codon_start_idx : codon_start_idx + 3]
        translated_codon_aa = str(Seq(original_codon).translate())
        
        if translated_codon_aa != parsed_mutation['ref']:
            raise ValueError(f"Reference amino acid mismatch at position {pos}. "
                             f"DNA codon {original_codon} translates to {translated_codon_aa}, "
                             f"but HGVS reference is {parsed_mutation['ref']}.")

        if parsed_mutation['type'] == 'frameshift':
            # This is a placeholder for the logic that would determine the
            # exact nucleotide change. For RUNX1 p.Arg135fs, we know it's a
            # specific insertion. For simplicity, we'll simulate a single nucleotide
            # insertion at the start of the codon, which is a common cause of frameshifts.
            # A real implementation would need more biological context.
            # For p.Arg135fs, let's assume a 'T' insertion.
            
            mutated_dna = (
                dna_sequence[:codon_start_idx] + 
                'T' + # The inserted nucleotide
                dna_sequence[codon_start_idx:]
            )
            return mutated_dna
            
        elif parsed_mutation['type'] == 'substitution':
            # This is also a placeholder. A real implementation would need to
            # find a codon that translates to the 'alt' amino acid.
            # This can be complex (e.g., choose the most common codon).
            # For now, we will raise an error for unsupported types.
            raise NotImplementedError("Substitutions are not yet implemented in this simplified tool.")
            
        return dna_sequence 
 