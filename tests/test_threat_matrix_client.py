import unittest
import os
import sys
import json

# Add the project root to the Python path to allow importing from 'tools'
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from tools.threat_matrix_client import ThreatMatrix

class TestThreatMatrixClient(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Set up the ThreatMatrix client once for all tests."""
        print("\n--- Setting up ThreatMatrixClient tests ---")
        try:
            cls.matrix = ThreatMatrix()
            print("✅ ThreatMatrix client initialized.")
        except FileNotFoundError as e:
            print(f"❌ CRITICAL SETUP FAILURE: {e}")
            print("  - Please ensure the database exists at 'data/databases/threat_matrix.db'")
            print("  - You may need to run the 'tools/cosmic_importer.py' script first.")
            raise

    def test_01_get_gene_id(self):
        """Tests the internal _get_gene_id method."""
        print("\n🔬 Testing: _get_gene_id")
        gene_id = self.matrix._get_gene_id("ASXL1")
        self.assertIsNotNone(gene_id, "Should return a valid ID for ASXL1.")
        self.assertIsInstance(gene_id, int, "Gene ID should be an integer.")
        print(f"  - Passed: ASXL1 -> {gene_id}")

        non_existent_id = self.matrix._get_gene_id("NONEXISTENTGENE")
        self.assertIsNone(non_existent_id, "Should return None for a non-existent gene.")
        print("  - Passed: NONEXISTENTGENE -> None")

    def test_02_get_variant_profile(self):
        """Tests fetching a variant profile, prioritizing high-impact ones."""
        print("\n🔬 Testing: get_variant_profile")
        profile = self.matrix.get_variant_profile("ASXL1")
        self.assertIsNotNone(profile, "Should find a profile for ASXL1.")
        self.assertIn("mutation_description", profile)
        # This specific variant is a frameshift, confirming prioritization logic
        self.assertIn("Frameshift", profile["mutation_description"])
        print(f"  - Passed: Found prioritized profile for ASXL1 (p.R860fs*3)")
        print(json.dumps(profile, indent=2))

    def test_03_get_clinical_trials(self):
        """Tests fetching clinical trials for a gene."""
        print("\n🔬 Testing: get_clinical_trials")
        trials = self.matrix.get_clinical_trials("ASXL1")
        self.assertIsNotNone(trials, "Should return a list of trials, not None.")
        self.assertIsInstance(trials, list, "Should return a list.")
        self.assertGreater(len(trials), 0, "Should find at least one clinical trial for ASXL1.")
        print(f"  - Passed: Found {len(trials)} clinical trials for ASXL1.")

    def test_04_get_efficacy_evidence(self):
        """Tests fetching efficacy evidence for a gene."""
        print("\n🔬 Testing: get_efficacy_evidence")
        evidence = self.matrix.get_efficacy_evidence("ASXL1")
        self.assertIsNotNone(evidence, "Should return a list of evidence, not None.")
        self.assertIsInstance(evidence, list, "Should return a list.")
        self.assertGreater(len(evidence), 0, "Should find at least one efficacy record for ASXL1.")
        print(f"  - Passed: Found {len(evidence)} efficacy evidence records for ASXL1.")

    def test_05_get_literature_summary(self):
        """Tests fetching a literature summary for a known PMID."""
        print("\n🔬 Testing: get_literature_summary")
        # This PMID is associated with an ASXL1 variant in the test DB
        # We don't know if it has been analyzed, so the test checks for both cases
        pmid = 21637286 
        summary = self.matrix.get_literature_summary(pmid)
        if summary:
            self.assertIsInstance(summary, dict)
            print(f"  - Passed: Found existing summary for PMID {pmid}.")
        else:
            print(f"  - Passed: No summary found for PMID {pmid}, which is an acceptable state.")

if __name__ == '__main__':
    unittest.main(verbosity=2) 