"""
ChopChop Integration for RUNX1-FPD Platform

This module integrates the professional ChopChop guide RNA design and scoring
suite with the RUNX1 Patient Digital Twin workflow, providing comprehensive
guide validation and safety analysis.
"""

import os
import json
import subprocess
import tempfile
import logging
import pandas as pd
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import shutil

logger = logging.getLogger(__name__)

class ChopChopIntegration:
    """
    Professional ChopChop integration for RUNX1-FPD guide design and validation.
    
    This class provides comprehensive guide RNA scoring, off-target analysis,
    and safety validation using the ChopChop suite.
    """
    
    def __init__(self):
        self.base_dir = Path(__file__).parent.parent
        self.chopchop_dir = self.base_dir / "tools" / "chopchop"
        self.chopchop_script = self.chopchop_dir / "chopchop.py"
        self.temp_dir = None
        
        # Verify ChopChop installation
        if not self.chopchop_script.exists():
            raise FileNotFoundError(f"ChopChop script not found at {self.chopchop_script}")
        
        # Default ChopChop parameters for RUNX1
        self.default_params = {
            "genome": "hg19",
            "pam": "NGG",
            "guide_size": 20,
            "max_mismatches": 3,
            "max_offtargets": 300,
            "scoring_method": "DOENCH_2016",
            "filter_gc_min": 40,
            "filter_gc_max": 70,
            "target_type": "CODING"
        }
        
    def setup_temp_directory(self) -> str:
        """Create temporary directory for ChopChop operations."""
        if self.temp_dir is None:
            self.temp_dir = tempfile.mkdtemp(prefix="chopchop_runx1_")
        return self.temp_dir
    
    def cleanup_temp_directory(self):
        """Clean up temporary directory."""
        if self.temp_dir and os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
            self.temp_dir = None
    
    def score_guides_for_runx1(self, guides: List[Dict], target_sequence: str = None) -> pd.DataFrame:
        """
        Score guide RNAs specifically for RUNX1 targeting.
        
        Args:
            guides: List of guide dictionaries with 'sequence' field
            target_sequence: Optional target sequence context
            
        Returns:
            pd.DataFrame: Scored guides with ChopChop metrics
        """
        if not guides:
            return pd.DataFrame()
        
        temp_dir = self.setup_temp_directory()
        
        try:
            # Create input file for ChopChop
            input_file = os.path.join(temp_dir, "runx1_target.txt")
            
            if target_sequence:
                with open(input_file, 'w') as f:
                    f.write(f">RUNX1_target\n{target_sequence}\n")
            else:
                # Use RUNX1 gene name for ChopChop lookup
                input_file = "RUNX1"
            
            # Run ChopChop for each guide
            scored_guides = []
            
            for i, guide in enumerate(guides):
                guide_seq = guide.get('sequence', guide.get('seq', ''))
                if not guide_seq:
                    continue
                
                try:
                    # Run ChopChop scoring
                    chopchop_result = self._run_chopchop_analysis(
                        input_file, guide_seq, temp_dir, guide_id=f"guide_{i}"
                    )
                    
                    # Merge guide data with ChopChop results
                    scored_guide = {**guide, **chopchop_result}
                    scored_guides.append(scored_guide)
                    
                except Exception as e:
                    logger.error(f"Error scoring guide {i}: {e}")
                    # Add guide with error status
                    scored_guide = {
                        **guide,
                        "chopchop_error": str(e),
                        "efficiency_score": 0.0,
                        "off_target_count": 999,
                        "safety_score": 0.0
                    }
                    scored_guides.append(scored_guide)
            
            return pd.DataFrame(scored_guides)
            
        finally:
            self.cleanup_temp_directory()
    
    def _run_chopchop_analysis(self, target_input: str, guide_sequence: str, 
                              temp_dir: str, guide_id: str) -> Dict:
        """
        Run ChopChop analysis for a single guide.
        
        Args:
            target_input: Target sequence file or gene name
            guide_sequence: Guide RNA sequence
            temp_dir: Temporary directory
            guide_id: Unique guide identifier
            
        Returns:
            Dict: ChopChop analysis results
        """
        output_dir = os.path.join(temp_dir, f"output_{guide_id}")
        os.makedirs(output_dir, exist_ok=True)
        
        # Build ChopChop command
        cmd = [
            "python2.7",  # ChopChop requires Python 2.7
            str(self.chopchop_script),
            "-G", self.default_params["genome"],
            "-o", output_dir,
            "-M", self.default_params["pam"],
            "--maxMismatches", str(self.default_params["max_mismatches"]),
            "-g", str(self.default_params["guide_size"]),
            "--scoringMethod", self.default_params["scoring_method"],
            "--filterGCmin", str(self.default_params["filter_gc_min"]),
            "--filterGCmax", str(self.default_params["filter_gc_max"]),
            "-t", self.default_params["target_type"],
            target_input
        ]
        
        try:
            # Run ChopChop
            result = subprocess.run(
                cmd, 
                cwd=self.chopchop_dir,
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )
            
            if result.returncode != 0:
                logger.error(f"ChopChop failed: {result.stderr}")
                return {"chopchop_error": result.stderr}
            
            # Parse ChopChop output
            return self._parse_chopchop_output(output_dir, guide_sequence)
            
        except subprocess.TimeoutExpired:
            logger.error(f"ChopChop timeout for guide {guide_id}")
            return {"chopchop_error": "Timeout"}
        except Exception as e:
            logger.error(f"ChopChop execution error: {e}")
            return {"chopchop_error": str(e)}
    
    def _parse_chopchop_output(self, output_dir: str, guide_sequence: str) -> Dict:
        """
        Parse ChopChop output files.
        
        Args:
            output_dir: ChopChop output directory
            guide_sequence: Guide sequence to match
            
        Returns:
            Dict: Parsed ChopChop results
        """
        results = {
            "efficiency_score": 0.0,
            "off_target_count": 0,
            "safety_score": 0.0,
            "gc_content": 0.0,
            "self_complementarity": 0.0,
            "chopchop_rank": 999
        }
        
        try:
            # Look for ChopChop results file
            results_file = None
            for file in os.listdir(output_dir):
                if file.endswith('.txt') and 'results' in file.lower():
                    results_file = os.path.join(output_dir, file)
                    break
            
            if not results_file:
                logger.warning(f"No ChopChop results file found in {output_dir}")
                return results
            
            # Parse results file
            with open(results_file, 'r') as f:
                lines = f.readlines()
            
            # Find matching guide in results
            for line in lines:
                if guide_sequence in line:
                    # Parse tab-separated values
                    fields = line.strip().split('\t')
                    if len(fields) >= 10:
                        try:
                            results["efficiency_score"] = float(fields[5]) if fields[5] != 'N/A' else 0.0
                            results["off_target_count"] = int(fields[6]) if fields[6] != 'N/A' else 0
                            results["gc_content"] = float(fields[7]) if fields[7] != 'N/A' else 0.0
                            results["self_complementarity"] = float(fields[8]) if fields[8] != 'N/A' else 0.0
                            results["chopchop_rank"] = int(fields[1]) if fields[1] != 'N/A' else 999
                            
                            # Calculate safety score (inverse of off-targets, scaled by efficiency)
                            if results["off_target_count"] == 0:
                                results["safety_score"] = results["efficiency_score"]
                            else:
                                results["safety_score"] = results["efficiency_score"] / (1 + results["off_target_count"])
                            
                            break
                        except (ValueError, IndexError) as e:
                            logger.error(f"Error parsing ChopChop results: {e}")
            
            return results
            
        except Exception as e:
            logger.error(f"Error parsing ChopChop output: {e}")
            return results
    
    def validate_guide_safety(self, guide_sequence: str, target_gene: str = "RUNX1") -> Dict:
        """
        Comprehensive safety analysis for a single guide.
        
        Args:
            guide_sequence: Guide RNA sequence
            target_gene: Target gene name
            
        Returns:
            Dict: Comprehensive safety assessment
        """
        temp_dir = self.setup_temp_directory()
        
        try:
            # Run comprehensive ChopChop analysis
            safety_results = self._run_chopchop_analysis(
                target_gene, guide_sequence, temp_dir, "safety_check"
            )
            
            # Calculate comprehensive safety score
            safety_assessment = self._calculate_safety_score(safety_results)
            
            return {
                "guide_sequence": guide_sequence,
                "target_gene": target_gene,
                "safety_assessment": safety_assessment,
                "chopchop_results": safety_results,
                "recommendations": self._generate_safety_recommendations(safety_assessment)
            }
            
        finally:
            self.cleanup_temp_directory()
    
    def _calculate_safety_score(self, chopchop_results: Dict) -> Dict:
        """
        Calculate comprehensive safety score from ChopChop results.
        
        Args:
            chopchop_results: ChopChop analysis results
            
        Returns:
            Dict: Safety assessment
        """
        efficiency = chopchop_results.get("efficiency_score", 0.0)
        off_targets = chopchop_results.get("off_target_count", 999)
        gc_content = chopchop_results.get("gc_content", 0.0)
        self_comp = chopchop_results.get("self_complementarity", 0.0)
        
        # Calculate individual safety metrics
        efficiency_score = min(efficiency / 0.8, 1.0)  # Normalize to 0.8 as excellent
        specificity_score = 1.0 / (1.0 + off_targets)  # Inverse of off-targets
        gc_score = 1.0 - abs(gc_content - 50.0) / 50.0  # Optimal GC ~50%
        structure_score = max(0.0, 1.0 - self_comp)  # Lower self-complementarity is better
        
        # Calculate composite safety score
        composite_score = (
            efficiency_score * 0.4 +
            specificity_score * 0.3 +
            gc_score * 0.2 +
            structure_score * 0.1
        )
        
        # Determine safety category
        if composite_score >= 0.8:
            safety_category = "excellent"
        elif composite_score >= 0.6:
            safety_category = "good"
        elif composite_score >= 0.4:
            safety_category = "acceptable"
        else:
            safety_category = "poor"
        
        return {
            "composite_score": composite_score,
            "safety_category": safety_category,
            "efficiency_score": efficiency_score,
            "specificity_score": specificity_score,
            "gc_score": gc_score,
            "structure_score": structure_score,
            "individual_metrics": {
                "efficiency": efficiency,
                "off_targets": off_targets,
                "gc_content": gc_content,
                "self_complementarity": self_comp
            }
        }
    
    def _generate_safety_recommendations(self, safety_assessment: Dict) -> List[str]:
        """
        Generate safety recommendations based on assessment.
        
        Args:
            safety_assessment: Safety assessment results
            
        Returns:
            List[str]: Safety recommendations
        """
        recommendations = []
        
        category = safety_assessment["safety_category"]
        metrics = safety_assessment["individual_metrics"]
        
        if category == "excellent":
            recommendations.append("Guide shows excellent safety profile - recommended for use")
        elif category == "good":
            recommendations.append("Guide shows good safety profile - suitable for most applications")
        elif category == "acceptable":
            recommendations.append("Guide shows acceptable safety - consider additional validation")
        else:
            recommendations.append("Guide shows poor safety profile - not recommended")
        
        # Specific recommendations based on metrics
        if metrics["efficiency"] < 0.3:
            recommendations.append("Low efficiency score - consider alternative guides")
        
        if metrics["off_targets"] > 10:
            recommendations.append("High off-target count - perform additional off-target validation")
        
        if metrics["gc_content"] < 30 or metrics["gc_content"] > 70:
            recommendations.append("Suboptimal GC content - may affect guide stability")
        
        if metrics["self_complementarity"] > 0.5:
            recommendations.append("High self-complementarity - may form secondary structures")
        
        return recommendations
    
    def batch_guide_analysis(self, guides: List[Dict], target_sequence: str = None) -> Dict:
        """
        Comprehensive batch analysis of multiple guides.
        
        Args:
            guides: List of guide dictionaries
            target_sequence: Optional target sequence
            
        Returns:
            Dict: Comprehensive batch analysis results
        """
        # Score all guides
        scored_guides_df = self.score_guides_for_runx1(guides, target_sequence)
        
        if scored_guides_df.empty:
            return {"error": "No guides could be analyzed"}
        
        # Rank guides by safety score
        scored_guides_df = scored_guides_df.sort_values('safety_score', ascending=False)
        
        # Generate summary statistics
        summary_stats = {
            "total_guides": len(scored_guides_df),
            "excellent_guides": len(scored_guides_df[scored_guides_df['safety_score'] >= 0.8]),
            "good_guides": len(scored_guides_df[scored_guides_df['safety_score'] >= 0.6]),
            "acceptable_guides": len(scored_guides_df[scored_guides_df['safety_score'] >= 0.4]),
            "poor_guides": len(scored_guides_df[scored_guides_df['safety_score'] < 0.4]),
            "mean_efficiency": scored_guides_df['efficiency_score'].mean(),
            "mean_off_targets": scored_guides_df['off_target_count'].mean(),
            "top_guide_score": scored_guides_df['safety_score'].max()
        }
        
        # Get top recommendations
        top_guides = scored_guides_df.head(5).to_dict('records')
        
        return {
            "summary_statistics": summary_stats,
            "top_guides": top_guides,
            "all_guides": scored_guides_df.to_dict('records'),
            "recommendations": self._generate_batch_recommendations(summary_stats, top_guides)
        }
    
    def _generate_batch_recommendations(self, summary_stats: Dict, top_guides: List[Dict]) -> List[str]:
        """
        Generate recommendations for batch analysis.
        
        Args:
            summary_stats: Summary statistics
            top_guides: Top-ranked guides
            
        Returns:
            List[str]: Batch recommendations
        """
        recommendations = []
        
        if summary_stats["excellent_guides"] > 0:
            recommendations.append(f"Found {summary_stats['excellent_guides']} excellent guides - proceed with top candidates")
        elif summary_stats["good_guides"] > 0:
            recommendations.append(f"Found {summary_stats['good_guides']} good guides - suitable for most applications")
        elif summary_stats["acceptable_guides"] > 0:
            recommendations.append(f"Found {summary_stats['acceptable_guides']} acceptable guides - consider additional validation")
        else:
            recommendations.append("No high-quality guides found - consider redesigning or expanding search")
        
        if summary_stats["mean_off_targets"] > 5:
            recommendations.append("High average off-target count - recommend additional specificity screening")
        
        if summary_stats["mean_efficiency"] < 0.4:
            recommendations.append("Low average efficiency - consider alternative target sites")
        
        if top_guides:
            best_guide = top_guides[0]
            recommendations.append(f"Top recommendation: Guide with {best_guide['safety_score']:.3f} safety score")
        
        return recommendations

# Convenience functions for easy integration
def score_runx1_guides(guides: List[Dict], target_sequence: str = None) -> pd.DataFrame:
    """
    Convenience function to score RUNX1 guides.
    
    Args:
        guides: List of guide dictionaries
        target_sequence: Optional target sequence
        
    Returns:
        pd.DataFrame: Scored guides
    """
    integrator = ChopChopIntegration()
    return integrator.score_guides_for_runx1(guides, target_sequence)

def validate_runx1_guide_safety(guide_sequence: str) -> Dict:
    """
    Convenience function to validate guide safety.
    
    Args:
        guide_sequence: Guide RNA sequence
        
    Returns:
        Dict: Safety assessment
    """
    integrator = ChopChopIntegration()
    return integrator.validate_guide_safety(guide_sequence, "RUNX1")

if __name__ == "__main__":
    # Test the ChopChop integration
    integrator = ChopChopIntegration()
    
    # Test with sample guides
    test_guides = [
        {"sequence": "GCTGAAGTTCTTATCGTCGT", "id": "test_guide_1"},
        {"sequence": "AGTCGATCGATCGATCGATC", "id": "test_guide_2"}
    ]
    
    print("Testing ChopChop integration...")
    
    try:
        # Test batch analysis
        results = integrator.batch_guide_analysis(test_guides)
        print(f"Batch analysis completed: {results['summary_statistics']['total_guides']} guides analyzed")
        
        # Test individual safety validation
        if test_guides:
            safety_results = integrator.validate_guide_safety(test_guides[0]["sequence"])
            print(f"Safety validation: {safety_results['safety_assessment']['safety_category']}")
            
    except Exception as e:
        print(f"Test failed: {e}")
        logger.error(f"ChopChop integration test failed: {e}") 