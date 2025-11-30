"""
Test Suite for Hotspot Mutation Detector (P0 Fix #3)
=====================================================
Tests COSMIC hotspot detection for KRAS/BRAF/NRAS mutations.

Manager Policy: C2 (MANAGER_ANSWERS_TO_ZO_SAE_QUESTIONS.md)
Owner: Zo
Date: January 13, 2025
"""

import sys
import os
import pytest

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../oncology-coPilot/oncology-backend-minimal'))

from api.services.hotspot_detector import HotspotDetector, detect_hotspot_mutation


class TestHotspotDetector:
    """Test COSMIC hotspot detection service."""
    
    def test_kras_g12d_hotspot(self):
        """Test KRAS G12D (highly recurrent hotspot)"""
        detector = HotspotDetector()
        result = detector.detect_hotspot("KRAS", "p.G12D")
        
        assert result.is_hotspot == True, "KRAS G12D should be a hotspot"
        assert result.gene == "KRAS"
        assert result.mutation == "G12D"
        assert result.pathway == "MAPK"
        assert result.evidence in ["highly_recurrent", "recurrent"]
    
    def test_kras_g12c_hotspot(self):
        """Test KRAS G12C (highly recurrent hotspot)"""
        detector = HotspotDetector()
        result = detector.detect_hotspot("KRAS", "p.G12C")
        
        assert result.is_hotspot == True, "KRAS G12C should be a hotspot"
        assert result.mutation == "G12C"
        assert result.pathway == "MAPK"
    
    def test_kras_g12v_hotspot(self):
        """Test KRAS G12V (highly recurrent hotspot)"""
        detector = HotspotDetector()
        result = detector.detect_hotspot("KRAS", "G12V")  # Without "p." prefix
        
        assert result.is_hotspot == True, "KRAS G12V should be a hotspot"
        assert result.mutation == "G12V"
    
    def test_braf_v600e_hotspot(self):
        """Test BRAF V600E (highly recurrent hotspot)"""
        detector = HotspotDetector()
        result = detector.detect_hotspot("BRAF", "p.V600E")
        
        assert result.is_hotspot == True, "BRAF V600E should be a hotspot"
        assert result.gene == "BRAF"
        assert result.mutation == "V600E"
        assert result.pathway == "MAPK"
        assert result.frequency > 0.5, "BRAF V600E is very common in melanoma"
    
    def test_nras_q61k_hotspot(self):
        """Test NRAS Q61K (highly recurrent hotspot)"""
        detector = HotspotDetector()
        result = detector.detect_hotspot("NRAS", "p.Q61K")
        
        assert result.is_hotspot == True, "NRAS Q61K should be a hotspot"
        assert result.gene == "NRAS"
        assert result.mutation == "Q61K"
        assert result.pathway == "MAPK"
    
    def test_nras_q61r_hotspot(self):
        """Test NRAS Q61R (highly recurrent hotspot)"""
        detector = HotspotDetector()
        result = detector.detect_hotspot("NRAS", "Q61R")  # Without "p." prefix
        
        assert result.is_hotspot == True, "NRAS Q61R should be a hotspot"
        assert result.mutation == "Q61R"
    
    def test_kras_g12a_non_highly_recurrent(self):
        """Test KRAS G12A (recurrent but not highly recurrent)"""
        detector = HotspotDetector()
        result = detector.detect_hotspot("KRAS", "p.G12A")
        
        assert result.is_hotspot == True, "KRAS G12A should be a hotspot (recurrent)"
        assert result.evidence == "recurrent"  # Not "highly_recurrent"
    
    def test_kras_non_hotspot(self):
        """Test KRAS non-hotspot mutation"""
        detector = HotspotDetector()
        result = detector.detect_hotspot("KRAS", "p.A146T")  # Not in COSMIC hotspots
        
        assert result.is_hotspot == False, "KRAS A146T should not be a hotspot"
        assert result.gene is None
        assert result.mutation is None
    
    def test_brca1_not_in_hotspot_database(self):
        """Test gene not in COSMIC hotspot database"""
        detector = HotspotDetector()
        result = detector.detect_hotspot("BRCA1", "p.C61G")
        
        assert result.is_hotspot == False, "BRCA1 not in COSMIC hotspot database"
    
    def test_invalid_hgvs(self):
        """Test invalid HGVS format"""
        detector = HotspotDetector()
        result = detector.detect_hotspot("KRAS", None)
        
        assert result.is_hotspot == False, "None HGVS should return False"
    
    def test_case_insensitive_gene(self):
        """Test case-insensitive gene matching"""
        detector = HotspotDetector()
        result1 = detector.detect_hotspot("kras", "p.G12D")  # lowercase
        result2 = detector.detect_hotspot("KRAS", "p.G12D")  # uppercase
        
        assert result1.is_hotspot == result2.is_hotspot == True
        assert result1.mutation == result2.mutation == "G12D"
    
    def test_batch_detection(self):
        """Test batch hotspot detection"""
        detector = HotspotDetector()
        
        mutations = [
            {"gene": "KRAS", "hgvs_p": "p.G12D"},
            {"gene": "BRAF", "hgvs_p": "p.V600E"},
            {"gene": "NRAS", "hgvs_p": "p.Q61K"},
            {"gene": "BRCA1", "hgvs_p": "p.C61G"}  # Not a hotspot
        ]
        
        results = detector.detect_batch(mutations)
        
        assert len(results) == 4
        assert results["KRAS:p.G12D"].is_hotspot == True
        assert results["BRAF:p.V600E"].is_hotspot == True
        assert results["NRAS:p.Q61K"].is_hotspot == True
        assert results["BRCA1:p.C61G"].is_hotspot == False
    
    def test_convenience_function(self):
        """Test module-level convenience function"""
        result = detect_hotspot_mutation("KRAS", "p.G12D")
        
        assert isinstance(result, dict)
        assert result["is_hotspot"] == True
        assert result["gene"] == "KRAS"
        assert result["mutation"] == "G12D"
        assert result["pathway"] == "MAPK"
        assert "cosmic_id" in result
        assert "evidence" in result
    
    def test_cosmic_provenance(self):
        """Test COSMIC provenance tracking"""
        detector = HotspotDetector()
        result = detector.detect_hotspot("KRAS", "p.G12D")
        
        assert result.cosmic_id is not None, "Should have COSMIC ID"
        assert result.source == "COSMIC v98"
        assert result.frequency is not None
        assert result.cancers is not None
        assert isinstance(result.cancers, list)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])







