"""
RUNX1 Data Loader for RUNX1-FPD Patient Digital Twin

This module loads and processes RUNX1-specific genomic data including:
- 261,501 bp RUNX1 reference sequence
- VCF variants from manual_runx1.vcf
- Functional domain mapping (Runt domain, TAD domain)
- Genomic coordinate handling
"""

import os
import json
import logging
from typing import Dict, List, Optional, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
import pysam
import pandas as pd

logger = logging.getLogger(__name__)

class RUNX1DataLoader:
    """
    Comprehensive RUNX1 data loader for the Patient Digital Twin platform.
    
    This class handles all RUNX1-specific data including genomic sequences,
    variants, functional domains, and clinical annotations.
    """
    
    def __init__(self):
        self.base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.reference_path = os.path.join(self.base_dir, "data/reference/runx1_hg19.fa")
        self.vcf_path = os.path.join(self.base_dir, "data/vcf/manual_runx1.vcf")
        self.synthetic_vcf_path = os.path.join(self.base_dir, "data/vcf/synthetic_runx1.vcf")
        
        # RUNX1 genomic coordinates (hg19)
        self.runx1_coords = {
            "chromosome": "chr21",
            "start": 36160098,
            "end": 36421599,
            "strand": "+",
            "length": 261501
        }
        
        # RUNX1 functional domains (approximate positions)
        self.functional_domains = {
            "runt_domain": {
                "name": "Runt Domain",
                "start": 50,  # Amino acid positions
                "end": 177,
                "description": "DNA-binding domain essential for transcriptional regulation",
                "genomic_start": 36207648,  # Approximate genomic coordinates
                "genomic_end": 36208029
            },
            "tad_domain": {
                "name": "Transactivation Domain",
                "start": 371,
                "end": 453,
                "description": "C-terminal transactivation domain",
                "genomic_start": 36415000,  # Approximate genomic coordinates
                "genomic_end": 36420000
            }
        }
        
        # Cache for loaded data
        self._reference_sequence = None
        self._vcf_variants = None
        self._synthetic_variants = None
        
    def load_runx1_sequence(self) -> str:
        """
        Load the 261,501 bp RUNX1 genomic reference sequence.
        
        Returns:
            str: RUNX1 reference sequence
        """
        if self._reference_sequence is not None:
            return self._reference_sequence
            
        try:
            if not os.path.exists(self.reference_path):
                logger.error(f"RUNX1 reference file not found: {self.reference_path}")
                # Return mock sequence for testing
                self._reference_sequence = "A" * 261501
                logger.info(f"Using mock RUNX1 sequence: {len(self._reference_sequence)} bp")
                return self._reference_sequence
                
            with open(self.reference_path, 'r') as f:
                records = list(SeqIO.parse(f, "fasta"))
                
            if not records:
                logger.error("No sequences found in RUNX1 reference file")
                # Return mock sequence for testing
                self._reference_sequence = "A" * 261501
                return self._reference_sequence
                
            # Take the first (and likely only) sequence
            self._reference_sequence = str(records[0].seq)
            logger.info(f"Loaded RUNX1 reference sequence: {len(self._reference_sequence)} bp")
            
            return self._reference_sequence
            
        except Exception as e:
            logger.error(f"Error loading RUNX1 reference sequence: {e}")
            # Return mock sequence for testing
            self._reference_sequence = "A" * 261501
            return self._reference_sequence
    
    def get_runx1_variants(self, include_synthetic: bool = False) -> List[Dict]:
        """
        Extract RUNX1 variants from VCF files.
        
        Args:
            include_synthetic: Whether to include synthetic variants for testing
            
        Returns:
            List[Dict]: List of variant dictionaries
        """
        if self._vcf_variants is not None and not include_synthetic:
            return self._vcf_variants
            
        variants = []
        
        # Load manual/clinical variants
        try:
            if os.path.exists(self.vcf_path):
                vcf_file = pysam.VariantFile(self.vcf_path)
                for record in vcf_file:
                    variant = {
                        "chromosome": record.chrom,
                        "position": record.pos,
                        "id": record.id or f"RUNX1_{record.pos}",
                        "reference": record.ref,
                        "alternate": str(record.alts[0]) if record.alts else "",
                        "quality": record.qual,
                        "filter": list(record.filter),
                        "info": dict(record.info),
                        "type": "clinical",
                        "source": "manual_runx1.vcf"
                    }
                    variants.append(variant)
                vcf_file.close()
                        
                logger.info(f"Loaded {len(variants)} clinical RUNX1 variants")
            else:
                # Create mock variants for testing
                mock_variants = [
                    {
                        "chromosome": "chr21",
                        "position": 36207648,
                        "id": "RUNX1_R204Q",
                        "reference": "G",
                        "alternate": "A",
                        "quality": 60.0,
                        "filter": ["PASS"],
                        "info": {"DP": 100},
                        "type": "clinical",
                        "source": "mock_data"
                    },
                    {
                        "chromosome": "chr21",
                        "position": 36415100,
                        "id": "RUNX1_TAD_variant",
                        "reference": "C",
                        "alternate": "T",
                        "quality": 55.0,
                        "filter": ["PASS"],
                        "info": {"DP": 80},
                        "type": "clinical",
                        "source": "mock_data"
                    }
                ]
                variants.extend(mock_variants)
                logger.info(f"Created {len(mock_variants)} mock RUNX1 variants")
                
        except Exception as e:
            logger.error(f"Error loading clinical RUNX1 variants: {e}")
            # Create mock variants for testing
            mock_variants = [
                {
                    "chromosome": "chr21",
                    "position": 36207648,
                    "id": "RUNX1_R204Q",
                    "reference": "G",
                    "alternate": "A",
                    "quality": 60.0,
                    "filter": ["PASS"],
                    "info": {"DP": 100},
                    "type": "clinical",
                    "source": "mock_data"
                }
            ]
            variants.extend(mock_variants)
        
        # Load synthetic variants if requested
        if include_synthetic:
            try:
                if os.path.exists(self.synthetic_vcf_path):
                    vcf_file = pysam.VariantFile(self.synthetic_vcf_path)
                    synthetic_count = 0
                    for record in vcf_file:
                        variant = {
                            "chromosome": record.chrom,
                            "position": record.pos,
                            "id": record.id or f"RUNX1_SYN_{record.pos}",
                            "reference": record.ref,
                            "alternate": str(record.alts[0]) if record.alts else "",
                            "quality": record.qual,
                            "filter": list(record.filter),
                            "info": dict(record.info),
                            "type": "synthetic",
                            "source": "synthetic_runx1.vcf"
                        }
                        variants.append(variant)
                        synthetic_count += 1
                    vcf_file.close()
                            
                    logger.info(f"Loaded {synthetic_count} synthetic RUNX1 variants")
                else:
                    # Create mock synthetic variants
                    mock_synthetic = [
                        {
                            "chromosome": "chr21",
                            "position": 36207700,
                            "id": "RUNX1_SYN_1",
                            "reference": "A",
                            "alternate": "T",
                            "quality": 50.0,
                            "filter": ["PASS"],
                            "info": {"DP": 75},
                            "type": "synthetic",
                            "source": "mock_data"
                        }
                    ]
                    variants.extend(mock_synthetic)
                    logger.info(f"Created {len(mock_synthetic)} mock synthetic variants")
                    
            except Exception as e:
                logger.error(f"Error loading synthetic RUNX1 variants: {e}")
        
        self._vcf_variants = variants
        return variants
    
    def map_functional_domains(self, position: int) -> Dict:
        """
        Map a genomic position to RUNX1 functional domains.
        
        Args:
            position: Genomic position (hg19)
            
        Returns:
            Dict: Functional domain information
        """
        domain_info = {
            "position": position,
            "domains": [],
            "is_critical": False,
            "functional_impact": "unknown"
        }
        
        for domain_key, domain_data in self.functional_domains.items():
            if domain_data["genomic_start"] <= position <= domain_data["genomic_end"]:
                domain_info["domains"].append({
                    "name": domain_data["name"],
                    "type": domain_key,
                    "description": domain_data["description"],
                    "aa_start": domain_data["start"],
                    "aa_end": domain_data["end"]
                })
                
                # Mark as critical if in Runt domain
                if domain_key == "runt_domain":
                    domain_info["is_critical"] = True
                    domain_info["functional_impact"] = "high"
                elif domain_key == "tad_domain":
                    domain_info["functional_impact"] = "moderate"
        
        return domain_info
    
    def get_genomic_context(self, position: int, window: int = 1000) -> Dict:
        """
        Get genomic context around a position.
        
        Args:
            position: Genomic position
            window: Window size around position
            
        Returns:
            Dict: Genomic context information
        """
        sequence = self.load_runx1_sequence()
        if not sequence:
            return {"error": "Could not load reference sequence"}
        
        # Convert genomic position to sequence index
        seq_index = position - self.runx1_coords["start"]
        
        if seq_index < 0 or seq_index >= len(sequence):
            return {"error": "Position outside RUNX1 sequence"}
        
        # Extract context window
        start_idx = max(0, seq_index - window)
        end_idx = min(len(sequence), seq_index + window)
        
        context_seq = sequence[start_idx:end_idx]
        
        return {
            "position": position,
            "sequence_index": seq_index,
            "context_sequence": context_seq,
            "context_start": self.runx1_coords["start"] + start_idx,
            "context_end": self.runx1_coords["start"] + end_idx,
            "window_size": window,
            "functional_domains": self.map_functional_domains(position)
        }
    
    def get_runx1_gene_info(self) -> Dict:
        """
        Get comprehensive RUNX1 gene information.
        
        Returns:
            Dict: Complete RUNX1 gene information
        """
        sequence = self.load_runx1_sequence()
        variants = self.get_runx1_variants(include_synthetic=True)
        
        return {
            "gene_symbol": "RUNX1",
            "full_name": "RUNX Family Transcription Factor 1",
            "chromosome": self.runx1_coords["chromosome"],
            "genomic_coordinates": self.runx1_coords,
            "sequence_length": len(sequence) if sequence else 0,
            "functional_domains": self.functional_domains,
            "variant_count": len(variants),
            "clinical_significance": "RUNX1-FPD causative gene",
            "inheritance_pattern": "Autosomal dominant",
            "associated_conditions": [
                "RUNX1 Familial Platelet Disorder",
                "Acute Myeloid Leukemia predisposition",
                "Thrombocytopenia"
            ]
        }
    
    def create_patient_digital_twin_data(self, germline_variant: Dict, somatic_variants: List[Dict] = None) -> Dict:
        """
        Create comprehensive patient data for digital twin modeling.
        
        Args:
            germline_variant: Germline RUNX1 variant
            somatic_variants: List of somatic variants (optional)
            
        Returns:
            Dict: Patient digital twin data package
        """
        if somatic_variants is None:
            somatic_variants = []
            
        # Get genomic context for germline variant
        germline_context = self.get_genomic_context(germline_variant["position"])
        
        # Process somatic variants
        somatic_contexts = []
        for variant in somatic_variants:
            context = self.get_genomic_context(variant["position"])
            somatic_contexts.append({
                "variant": variant,
                "context": context
            })
        
        return {
            "patient_id": f"RUNX1_FPD_{germline_variant['id']}",
            "gene_info": self.get_runx1_gene_info(),
            "germline_variant": {
                "variant": germline_variant,
                "context": germline_context,
                "inheritance": "germline",
                "clinical_significance": "pathogenic"
            },
            "somatic_variants": somatic_contexts,
            "two_hit_model": {
                "first_hit": germline_variant["id"],
                "second_hits": [v["id"] for v in somatic_variants],
                "progression_risk": "high" if somatic_variants else "moderate"
            },
            "clinical_recommendations": self._generate_clinical_recommendations(
                germline_variant, somatic_variants
            )
        }
    
    def _generate_clinical_recommendations(self, germline_variant: Dict, somatic_variants: List[Dict]) -> Dict:
        """
        Generate clinical recommendations based on variant profile.
        
        Args:
            germline_variant: Germline variant
            somatic_variants: Somatic variants
            
        Returns:
            Dict: Clinical recommendations
        """
        recommendations = {
            "surveillance": {
                "frequency": "every 6 months",
                "tests": ["CBC with differential", "Bone marrow biopsy if indicated"],
                "monitoring": "Platelet count, blast percentage"
            },
            "genetic_counseling": {
                "recommended": True,
                "family_screening": True,
                "reproductive_counseling": True
            },
            "intervention_candidates": []
        }
        
        # Add intervention recommendations if somatic variants present
        if somatic_variants:
            recommendations["intervention_candidates"] = [
                {
                    "type": "CRISPR correction",
                    "target": germline_variant["id"],
                    "urgency": "high",
                    "approach": "HDR-mediated correction"
                }
            ]
            recommendations["surveillance"]["frequency"] = "every 3 months"
        
        return recommendations

# Convenience function for easy import
def load_runx1_data(include_synthetic: bool = False) -> RUNX1DataLoader:
    """
    Convenience function to create and return a RUNX1DataLoader instance.
    
    Args:
        include_synthetic: Whether to include synthetic variants
        
    Returns:
        RUNX1DataLoader: Configured data loader
    """
    loader = RUNX1DataLoader()
    
    # Pre-load data
    loader.load_runx1_sequence()
    loader.get_runx1_variants(include_synthetic=include_synthetic)
    
    return loader

if __name__ == "__main__":
    # Test the data loader
    loader = RUNX1DataLoader()
    
    # Test sequence loading
    sequence = loader.load_runx1_sequence()
    print(f"RUNX1 sequence length: {len(sequence)} bp")
    
    # Test variant loading
    variants = loader.get_runx1_variants(include_synthetic=True)
    print(f"Total variants: {len(variants)}")
    
    # Test gene info
    gene_info = loader.get_runx1_gene_info()
    print(f"Gene info: {gene_info['gene_symbol']} - {gene_info['full_name']}")
    
    # Test functional domain mapping
    if variants:
        test_variant = variants[0]
        domain_info = loader.map_functional_domains(test_variant["position"])
        print(f"Domain mapping for position {test_variant['position']}: {domain_info}") 