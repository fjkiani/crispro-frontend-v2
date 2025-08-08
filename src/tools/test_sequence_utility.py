import unittest
import os
from tools.sequence_utility import apply_protein_change

class TestSequenceUtility(unittest.TestCase):

    def setUp(self):
        """Set up the test case by reading the official FASTA file."""
        fasta_path = "data/reference/tp53_sequence.fasta"
        self.assertTrue(os.path.exists(fasta_path), f"FASTA file not found at {fasta_path}")
        
        with open(fasta_path, "r") as f:
            lines = f.readlines()
            # Skip the header line (e.g., '>tp53_cds') and join the rest
            self.tp53_cds = "".join(line.strip() for line in lines[1:])

        self.assertGreater(len(self.tp53_cds), 0, "Failed to read a non-empty sequence from FASTA file.")

    def test_r175h_mutation_creates_different_sequence(self):
        """
        Tests if applying the R175H mutation results in a different DNA sequence.
        This is the most critical test to validate the core logic.
        """
        protein_change = "R175H"
        mutated_sequence = apply_protein_change(self.tp53_cds, protein_change)

        self.assertIsNotNone(mutated_sequence, "Mutated sequence should not be None.")
        self.assertNotEqual(self.tp53_cds, mutated_sequence, 
                            f"Applying mutation {protein_change} did not change the original DNA sequence. This is the root cause of the failure.")
        
        print(f"\\n✅ Test passed: Applying {protein_change} successfully altered the DNA sequence.")

    def test_benign_mutation_also_changes_sequence(self):
        """
        Tests a different mutation to ensure the logic isn't hardcoded.
        """
        protein_change = "V272G"
        mutated_sequence = apply_protein_change(self.tp53_cds, protein_change)

        self.assertIsNotNone(mutated_sequence)
        self.assertNotEqual(self.tp53_cds, mutated_sequence)
        print(f"✅ Test passed: Applying {protein_change} successfully altered the DNA sequence.")

    def test_invalid_position_returns_none(self):
        """
        Tests that a position outside the bounds of the protein returns None.
        """
        protein_change = "R999H" # Position 999 is out of bounds
        mutated_sequence = apply_protein_change(self.tp53_cds, protein_change)
        self.assertIsNone(mutated_sequence, "Should return None for an out-of-bounds mutation.")
        print("✅ Test passed: Out-of-bounds mutation correctly handled.")

if __name__ == '__main__':
    unittest.main() 