"""
RUNX1 Integration Plan - Complete Implementation Guide

This module provides the exact integration plan for connecting all RUNX1 components
to the Patient Digital Twin, with specific code examples and implementation steps.
"""

import streamlit as st
from typing import Dict, List, Optional
import logging
import pandas as pd # Added missing import for pandas

# Import all the RUNX1 components we've built - now local imports
from runx1_data_loader import RUNX1DataLoader
from runx1_progression_modeler import RUNX1ProgressionModeler
from chopchop_integration import ChopChopIntegration
# Import from tools directory for components that remain there
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'tools'))
from intelligent_guide_finder import find_intelligent_guides
from experiment_advisor import generate_protocol, estimate_experiment_timeline, recommend_reagents
from command_center_client import formulate_battle_plan, get_all_patients, initiate_germline_correction
from ui_enhancements import AppleStyleResultsDisplay, LEAPGrantAligner

logger = logging.getLogger(__name__)

class RUNX1DigitalTwinIntegrator:
    """
    Complete integration of RUNX1 tools into the Patient Digital Twin workflow.
    
    This class orchestrates all RUNX1-specific components and provides a unified
    interface for the Patient Digital Twin page.
    """
    
    def __init__(self):
        self.data_loader = RUNX1DataLoader()
        self.progression_modeler = RUNX1ProgressionModeler()
        self.chopchop_integration = ChopChopIntegration()
        
        # Pre-load RUNX1 data
        self.runx1_sequence = self.data_loader.load_runx1_sequence()
        self.runx1_variants = self.data_loader.get_runx1_variants(include_synthetic=True)
        self.runx1_gene_info = self.data_loader.get_runx1_gene_info()
        
    def enhanced_two_hit_analysis(self, germline_variant: Dict, somatic_variant: Dict = None, 
                                patient_age: float = 35.0) -> Dict:
        """
        Enhanced two-hit analysis with progression modeling and intervention planning.
        
        This replaces the basic two-hit analysis in the Patient Digital Twin with
        sophisticated progression modeling and clinical decision support.
        
        Args:
            germline_variant: Germline RUNX1 variant data
            somatic_variant: Somatic variant data (optional)
            patient_age: Patient age
            
        Returns:
            Dict: Complete two-hit analysis with progression modeling
        """
        logger.info("Running enhanced RUNX1 two-hit analysis")
        
        # Step 1: Load genomic context for variants
        germline_context = self.data_loader.get_genomic_context(
            germline_variant.get("position", 36207648)  # Default to Runt domain
        )
        
        somatic_context = None
        if somatic_variant:
            somatic_context = self.data_loader.get_genomic_context(
                somatic_variant.get("position", 31022441)  # Default ASXL1 position
            )
        
        # Step 2: Run progression modeling
        somatic_variants = [somatic_variant] if somatic_variant else []
        progression_model = self.progression_modeler.model_clonal_evolution(
            germline_variant, somatic_variants, patient_age
        )
        
        # Step 3: Generate clinical recommendations
        clinical_recommendations = self._generate_clinical_recommendations(
            progression_model, germline_context, somatic_context
        )
        
        # Step 4: Identify intervention opportunities
        intervention_opportunities = self.progression_modeler.identify_intervention_opportunities(
            progression_model
        )
        
        return {
            "germline_analysis": {
                "variant": germline_variant,
                "genomic_context": germline_context,
                "functional_impact": self._assess_functional_impact(germline_variant, germline_context)
            },
            "somatic_analysis": {
                "variant": somatic_variant,
                "genomic_context": somatic_context,
                "functional_impact": self._assess_functional_impact(somatic_variant, somatic_context) if somatic_variant else None
            },
            "progression_model": progression_model,
            "clinical_recommendations": clinical_recommendations,
            "intervention_opportunities": intervention_opportunities,
            "risk_stratification": self._stratify_risk(progression_model),
            "next_steps": self._generate_next_steps(intervention_opportunities)
        }
    
    def design_precision_intervention(self, target_variant: Dict, intervention_type: str = "crispr_correction") -> Dict:
        """
        Design precision intervention for identified target variant.
        
        This integrates guide design, safety analysis, and experimental protocols
        into a complete therapeutic blueprint.
        
        Args:
            target_variant: Target variant for intervention
            intervention_type: Type of intervention (crispr_correction, base_editing, etc.)
            
        Returns:
            Dict: Complete intervention design
        """
        logger.info(f"Designing precision intervention for {target_variant.get('id', 'unknown')}")
        
        # Step 1: Extract genomic locus for guide design
        genomic_position = target_variant.get("position", 36207648)
        locus = f"chr21:{genomic_position}"
        
        # Step 2: Design intelligent guides
        guides = find_intelligent_guides(locus, "hg19", "NGG")
        
        if not guides:
            return {"error": "No suitable guides found for target locus"}
        
        # Step 3: Score guides with ChopChop
        try:
            scored_guides_df = self.chopchop_integration.score_guides_for_runx1(guides)
        except Exception as e:
            logger.error(f"ChopChop scoring failed: {e}")
            scored_guides_df = pd.DataFrame()
        
        # Step 4: Select top guide
        if not scored_guides_df.empty:
            top_guide = scored_guides_df.iloc[0].to_dict()
        else:
            top_guide = guides[0]  # Fallback to first guide
        
        # Step 5: Generate experimental protocol
        protocol = generate_protocol(
            experiment_type="knock-in" if intervention_type == "crispr_correction" else "knockout",
            target_gene="RUNX1",
            guide_info=top_guide,
            cell_type="CD34+ HSCs"
        )
        
        # Step 6: Safety assessment
        try:
            safety_assessment = self.chopchop_integration.validate_guide_safety(
                top_guide.get("sequence", top_guide.get("seq", ""))
            )
        except Exception as e:
            logger.error(f"Safety assessment failed: {e}")
            safety_assessment = {"error": str(e)}
        
        # Step 7: Try to use Command Center for germline correction design
        try:
            if intervention_type == "crispr_correction":
                # Get genomic context for the target
                context = self.data_loader.get_genomic_context(genomic_position, window=2000)
                corrected_sequence = context.get("context_sequence", "")
                target_site_context = f"chr21:{genomic_position-2000}-{genomic_position+2000}"
                
                # Call Command Center for germline correction
                correction_blueprint = initiate_germline_correction(
                    corrected_sequence, target_site_context, arm_length=500
                )
                
                if correction_blueprint and len(correction_blueprint) == 2:
                    guides_df, blueprint_data = correction_blueprint
                    if not guides_df.empty:
                        top_guide = guides_df.iloc[0].to_dict()
        except Exception as e:
            logger.error(f"Command Center integration failed: {e}")
        
        return {
            "target_variant": target_variant,
            "intervention_type": intervention_type,
            "guide_design": {
                "top_guide": top_guide,
                "all_guides": guides,
                "chopchop_scores": scored_guides_df.to_dict('records') if not scored_guides_df.empty else []
            },
            "safety_assessment": safety_assessment,
            "experimental_protocol": protocol,
            "success_probability": self._calculate_success_probability(top_guide, safety_assessment),
            "timeline_estimate": estimate_experiment_timeline("knock-in" if intervention_type == "crispr_correction" else "knockout"),
            "resource_requirements": recommend_reagents("knock-in" if intervention_type == "crispr_correction" else "knockout", "CD34+ HSCs")
        }
    
    def create_genomic_browser_data(self, variant_position: int, window: int = 5000) -> Dict:
        """
        Create data for interactive genomic browser visualization.
        
        Args:
            variant_position: Genomic position to center on
            window: Window size around position
            
        Returns:
            Dict: Genomic browser data
        """
        # Get genomic context
        context = self.data_loader.get_genomic_context(variant_position, window)
        
        # Map functional domains
        domain_mapping = self.data_loader.map_functional_domains(variant_position)
        
        # Get nearby variants
        nearby_variants = []
        for variant in self.runx1_variants:
            if abs(variant["position"] - variant_position) <= window:
                nearby_variants.append(variant)
        
        return {
            "center_position": variant_position,
            "window_size": window,
            "sequence_context": context,
            "functional_domains": domain_mapping,
            "nearby_variants": nearby_variants,
            "gene_structure": self._get_gene_structure_data(),
            "conservation_scores": self._get_conservation_scores(variant_position, window)
        }
    
    def generate_clinical_report(self, analysis_results: Dict, patient_info: Dict = None) -> Dict:
        """
        Generate comprehensive clinical report for physicians.
        
        Args:
            analysis_results: Complete analysis results
            patient_info: Patient demographic information
            
        Returns:
            Dict: Clinical report
        """
        progression_model = analysis_results.get("progression_model", {})
        risk_stratification = analysis_results.get("risk_stratification", {})
        
        return {
            "patient_summary": {
                "patient_id": patient_info.get("id", "unknown") if patient_info else "unknown",
                "age": patient_info.get("age", 35) if patient_info else 35,
                "diagnosis": "RUNX1 Familial Platelet Disorder",
                "analysis_date": st.session_state.get("analysis_date", "today")
            },
            "genetic_profile": {
                "germline_mutation": analysis_results["germline_analysis"]["variant"],
                "somatic_mutations": [analysis_results["somatic_analysis"]["variant"]] if analysis_results["somatic_analysis"]["variant"] else [],
                "functional_impact": analysis_results["germline_analysis"]["functional_impact"]
            },
            "risk_assessment": {
                "transformation_probability": progression_model.get("transformation_probability", {}),
                "risk_tier": risk_stratification.get("risk_tier", "unknown"),
                "surveillance_frequency": risk_stratification.get("surveillance_frequency", "every 6 months")
            },
            "clinical_recommendations": analysis_results.get("clinical_recommendations", {}),
            "intervention_options": analysis_results.get("intervention_opportunities", []),
            "surveillance_plan": self._generate_surveillance_plan(risk_stratification),
            "family_counseling": self._generate_family_counseling_recommendations(analysis_results)
        }
    
    def _assess_functional_impact(self, variant: Dict, genomic_context: Dict) -> Dict:
        """Assess functional impact of variant."""
        if not variant:
            return {"impact": "none", "description": "No variant provided"}
        
        # Check if variant is in functional domain
        functional_domains = genomic_context.get("functional_domains", {})
        
        if functional_domains.get("is_critical", False):
            impact_level = "high"
            description = f"Variant in critical {functional_domains['domains'][0]['name']}"
        elif functional_domains.get("domains", []):
            impact_level = "moderate"
            description = f"Variant in {functional_domains['domains'][0]['name']}"
        else:
            impact_level = "low"
            description = "Variant outside known functional domains"
        
        return {
            "impact": impact_level,
            "description": description,
            "functional_domains": functional_domains.get("domains", []),
            "molecular_consequence": variant.get("consequence", "unknown")
        }
    
    def _generate_clinical_recommendations(self, progression_model: Dict, 
                                        germline_context: Dict, somatic_context: Dict) -> Dict:
        """Generate clinical recommendations based on progression model."""
        transformation_prob = progression_model.get("transformation_probability", {})
        current_prob = transformation_prob.get("current_probability", 0.0)
        
        recommendations = {
            "immediate_actions": [],
            "surveillance_plan": {},
            "intervention_considerations": [],
            "family_screening": {}
        }
        
        # Immediate actions based on risk
        if current_prob > 0.5:
            recommendations["immediate_actions"] = [
                "Urgent hematology consultation",
                "Bone marrow biopsy within 2 weeks",
                "Consider clinical trial enrollment",
                "Multidisciplinary team consultation"
            ]
        elif current_prob > 0.3:
            recommendations["immediate_actions"] = [
                "Hematology consultation within 1 month",
                "Enhanced surveillance protocol",
                "Consider intervention planning"
            ]
        else:
            recommendations["immediate_actions"] = [
                "Routine hematology follow-up",
                "Standard surveillance protocol",
                "Genetic counseling"
            ]
        
        # Surveillance plan
        if current_prob > 0.3:
            surveillance_frequency = "every 3 months"
        elif current_prob > 0.1:
            surveillance_frequency = "every 6 months"
        else:
            surveillance_frequency = "annually"
        
        recommendations["surveillance_plan"] = {
            "frequency": surveillance_frequency,
            "tests": [
                "Complete blood count with differential",
                "Peripheral blood smear",
                "Flow cytometry if indicated",
                "Bone marrow biopsy if blasts >5%"
            ],
            "monitoring_parameters": [
                "Platelet count",
                "Blast percentage",
                "Cytogenetic analysis"
            ]
        }
        
        # Intervention considerations
        if current_prob > 0.1:
            recommendations["intervention_considerations"] = [
                "CRISPR-based germline correction",
                "Targeted therapy for somatic mutations",
                "Clinical trial participation"
            ]
        
        # Family screening
        recommendations["family_screening"] = {
            "recommended": True,
            "priority": "high" if current_prob > 0.3 else "moderate",
            "tests": [
                "RUNX1 sequencing",
                "Platelet count assessment",
                "Family history documentation"
            ]
        }
        
        return recommendations
    
    def _stratify_risk(self, progression_model: Dict) -> Dict:
        """Stratify patient risk based on progression model."""
        transformation_prob = progression_model.get("transformation_probability", {})
        current_prob = transformation_prob.get("current_probability", 0.0)
        
        if current_prob > 0.5:
            return {
                "risk_tier": "very_high",
                "surveillance_frequency": "monthly",
                "intervention_urgency": "immediate",
                "description": "Very high risk of malignant transformation"
            }
        elif current_prob > 0.3:
            return {
                "risk_tier": "high",
                "surveillance_frequency": "quarterly",
                "intervention_urgency": "urgent",
                "description": "High risk of malignant transformation"
            }
        elif current_prob > 0.1:
            return {
                "risk_tier": "moderate",
                "surveillance_frequency": "every 6 months",
                "intervention_urgency": "high",
                "description": "Moderate risk of malignant transformation"
            }
        else:
            return {
                "risk_tier": "low",
                "surveillance_frequency": "annually",
                "intervention_urgency": "moderate",
                "description": "Low risk of malignant transformation"
            }
    
    def _generate_next_steps(self, intervention_opportunities: List[Dict]) -> List[Dict]:
        """Generate prioritized next steps."""
        next_steps = []
        
        # Sort opportunities by priority
        sorted_opportunities = sorted(
            intervention_opportunities, 
            key=lambda x: x.get("priority_rank", 999)
        )
        
        for i, opportunity in enumerate(sorted_opportunities[:3]):  # Top 3
            next_steps.append({
                "step": i + 1,
                "action": f"Plan {opportunity['intervention_type']} intervention",
                "target": opportunity["target_mutation"],
                "timeline": f"Within {opportunity['time_window'][1] - opportunity['time_window'][0]:.1f} years",
                "urgency": opportunity["urgency"],
                "description": opportunity["description"]
            })
        
        return next_steps
    
    def _calculate_success_probability(self, guide: Dict, safety_assessment: Dict) -> float:
        """Calculate overall success probability for intervention."""
        # Base success probability for CRISPR correction
        base_success = 0.6
        
        # Modify based on guide quality
        guide_quality = guide.get("safety_score", 0.5)
        safety_score = safety_assessment.get("safety_assessment", {}).get("composite_score", 0.5)
        
        # Combined success probability
        success_prob = base_success * guide_quality * safety_score
        
        return min(0.95, max(0.1, success_prob))
    
    def _estimate_intervention_timeline(self, protocol: Dict) -> Dict:
        """Estimate timeline for intervention."""
        timeline = protocol.get("timeline", {})
        
        total_weeks = 0
        for stage in timeline.get("stages", []):
            total_weeks += stage.get("duration_weeks", 1)
        
        return {
            "total_duration_weeks": total_weeks,
            "total_duration_months": total_weeks / 4.0,
            "key_milestones": timeline.get("stages", []),
            "critical_path": "Guide design → Cell preparation → Editing → Validation"
        }
    
    def _calculate_resource_requirements(self, protocol: Dict) -> Dict:
        """Calculate resource requirements."""
        return {
            "personnel": [
                "Molecular biologist",
                "Cell culture specialist",
                "Bioinformatician",
                "Clinical hematologist"
            ],
            "equipment": [
                "Cell culture facilities",
                "Electroporation equipment",
                "Flow cytometer",
                "NGS sequencer"
            ],
            "reagents": protocol.get("reagents", []),
            "estimated_cost": "$50,000 - $100,000",
            "timeline": "3-6 months"
        }
    
    def _get_gene_structure_data(self) -> Dict:
        """Get RUNX1 gene structure data for visualization."""
        return {
            "exons": [
                {"start": 36207648, "end": 36207848, "number": 1},
                {"start": 36208000, "end": 36208200, "number": 2},
                {"start": 36415000, "end": 36415200, "number": 8}
            ],
            "introns": [
                {"start": 36207848, "end": 36208000},
                {"start": 36208200, "end": 36415000}
            ],
            "functional_domains": self.data_loader.functional_domains
        }
    
    def _get_conservation_scores(self, position: int, window: int) -> List[Dict]:
        """Get conservation scores for genomic region."""
        # Mock conservation scores - in real implementation, 
        # this would query phyloP/phastCons databases
        scores = []
        for i in range(position - window, position + window, 100):
            scores.append({
                "position": i,
                "phylop_score": 0.5 + (i % 100) / 200.0,  # Mock score
                "phastcons_score": 0.3 + (i % 50) / 100.0  # Mock score
            })
        return scores
    
    def _generate_surveillance_plan(self, risk_stratification: Dict) -> Dict:
        """Generate detailed surveillance plan."""
        frequency = risk_stratification.get("surveillance_frequency", "every 6 months")
        
        return {
            "frequency": frequency,
            "tests": [
                "Complete blood count with differential",
                "Peripheral blood smear",
                "Flow cytometry for blast assessment",
                "Bone marrow biopsy if indicated"
            ],
            "monitoring_parameters": [
                "Platelet count (<50,000 concerning)",
                "Blast percentage (>5% concerning)",
                "Cytogenetic abnormalities",
                "Clonal evolution markers"
            ],
            "escalation_criteria": [
                "Platelet count <50,000",
                "Blast percentage >5%",
                "New cytogenetic abnormalities",
                "Rapid clonal expansion"
            ]
        }
    
    def _generate_family_counseling_recommendations(self, analysis_results: Dict) -> Dict:
        """Generate family counseling recommendations."""
        return {
            "genetic_counseling": {
                "recommended": True,
                "urgency": "high",
                "topics": [
                    "RUNX1-FPD inheritance pattern",
                    "Family screening recommendations",
                    "Reproductive counseling",
                    "Surveillance protocols"
                ]
            },
            "family_screening": {
                "recommended_relatives": [
                    "First-degree relatives",
                    "Children (if planning pregnancy)"
                ],
                "screening_tests": [
                    "RUNX1 sequencing",
                    "Platelet count",
                    "Bleeding history assessment"
                ]
            },
            "reproductive_counseling": {
                "preconception_counseling": True,
                "prenatal_testing_options": [
                    "Chorionic villus sampling",
                    "Amniocentesis",
                    "Preimplantation genetic testing"
                ]
            }
        }

# Integration functions for easy use in Patient Digital Twin

def integrate_runx1_analysis(germline_variant: Dict, somatic_variant: Dict = None, 
                           patient_age: float = 35.0) -> Dict:
    """
    Main integration function for RUNX1 analysis in Patient Digital Twin.
    
    Args:
        germline_variant: Germline RUNX1 variant
        somatic_variant: Somatic variant (optional)
        patient_age: Patient age
        
    Returns:
        Dict: Complete RUNX1 analysis results
    """
    integrator = RUNX1DigitalTwinIntegrator()
    return integrator.enhanced_two_hit_analysis(germline_variant, somatic_variant, patient_age)

def design_runx1_intervention(target_variant: Dict, intervention_type: str = "crispr_correction") -> Dict:
    """
    Design precision intervention for RUNX1 variant.
    
    Args:
        target_variant: Target variant
        intervention_type: Intervention type
        
    Returns:
        Dict: Intervention design
    """
    integrator = RUNX1DigitalTwinIntegrator()
    return integrator.design_precision_intervention(target_variant, intervention_type)

def create_runx1_genomic_browser(variant_position: int, window: int = 5000) -> Dict:
    """
    Create genomic browser data for RUNX1 visualization.
    
    Args:
        variant_position: Genomic position
        window: Window size
        
    Returns:
        Dict: Genomic browser data
    """
    integrator = RUNX1DigitalTwinIntegrator()
    return integrator.create_genomic_browser_data(variant_position, window)

def generate_runx1_clinical_report(analysis_results: Dict, patient_info: Dict = None) -> Dict:
    """
    Generate clinical report for RUNX1 analysis.
    
    Args:
        analysis_results: Analysis results
        patient_info: Patient information
        
    Returns:
        Dict: Clinical report
    """
    integrator = RUNX1DigitalTwinIntegrator()
    return integrator.generate_clinical_report(analysis_results, patient_info)

# Example usage and testing
if __name__ == "__main__":
    # Test the integration
    integrator = RUNX1DigitalTwinIntegrator()
    
    # Test data
    test_germline = {
        "id": "RUNX1_R204Q",
        "gene": "RUNX1",
        "position": 36207648,
        "consequence": "missense_variant"
    }
    
    test_somatic = {
        "id": "ASXL1_G646fs",
        "gene": "ASXL1",
        "position": 31022441,
        "consequence": "frameshift_variant"
    }
    
    print("Testing RUNX1 integration...")
    
    try:
        # Test enhanced analysis
        results = integrator.enhanced_two_hit_analysis(test_germline, test_somatic, 45.0)
        print(f"Analysis completed: {results['risk_stratification']['risk_tier']} risk")
        
        # Test intervention design
        intervention = integrator.design_precision_intervention(test_germline)
        print(f"Intervention designed: {intervention['success_probability']:.3f} success probability")
        
        # Test clinical report
        report = integrator.generate_clinical_report(results)
        print(f"Clinical report generated: {report['risk_assessment']['risk_tier']} risk tier")
        
    except Exception as e:
        print(f"Integration test failed: {e}")
        logger.error(f"RUNX1 integration test failed: {e}") 