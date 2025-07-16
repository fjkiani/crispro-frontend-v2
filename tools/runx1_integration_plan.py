"""
RUNX1 Integration Plan - Complete Implementation Guide

This module provides the exact integration plan for connecting all RUNX1 components
to the Patient Digital Twin, with specific code examples and implementation steps.
"""

import streamlit as st
from typing import Dict, List, Optional
import logging
import pandas as pd # Added missing import for pandas
import time
import requests
import gosling as gos

# Import all the RUNX1 components we've built
from tools.runx1_data_loader import RUNX1DataLoader
from tools.runx1_progression_modeler import RUNX1ProgressionModeler
# from tools.chopchop_integration import ChopChopIntegration # ARCHIVED: V3 workflow uses Evo service
# from tools.intelligent_guide_finder import find_intelligent_guides # ARCHIVED
from tools.experiment_advisor import generate_protocol, estimate_experiment_timeline, recommend_reagents
from tools.ui_enhancements import AppleStyleResultsDisplay, LEAPGrantAligner

logger = logging.getLogger(__name__)

# This now points to our live V3 deployment, corrected from deployment logs.
COMMAND_CENTER_URL = "https://crispro--crispr-assistant-command-center-v3-commandcenter-api.modal.run"


def run_live_intervention_design(target_variant_data):
    """
    Calls the full V3 Command Center workflow and polls for results.
    This is a blocking function to work with Streamlit's execution model.
    """
    try:
        # Step 1: Formulate Battle Plan
        st.info("🛰️ Contacting Command Center: Formulating Battle Plan...")
        formulate_payload = {
            "patient_identifier": target_variant_data.get("id", "runx1_digital_twin_v3"),
            "gene": target_variant_data.get("gene", "RUNX1"),
            "mutation_hgvs_p": target_variant_data.get("protein_change", "p.Arg135fs"),
            # NOTE: Sequence is currently a placeholder, a future enhancement will be to fetch this dynamically.
            "sequence_for_perplexity": "CCTGCCTGGGCTCCCCAGCCCCAGCAGCCCTGCAGCCCAGCTCTGGCCTGCCGGAGGCCACGCGC"
        }
        formulate_url = f"{COMMAND_CENTER_URL}/v3/workflow/formulate_battle_plan"
        response = requests.post(formulate_url, json=formulate_payload, timeout=60)
        response.raise_for_status()
        plan_id = response.json()["plan_id"]
        st.info(f"✅ Battle Plan {plan_id} created.")

        # Step 2: Execute Guide Design
        st.info(f"🚀 Executing Guide Design Campaign for Plan {plan_id}...")
        execute_url = f"{COMMAND_CENTER_URL}/v3/workflow/execute_guide_design_campaign/{plan_id}"
        response = requests.post(execute_url, timeout=60)
        response.raise_for_status()
        st.info("💥 Campaign initiated. Polling for results... (this may take several minutes)")

        # Step 3: Poll for Results
        status_url = f"{COMMAND_CENTER_URL}/v3/battle_plan/{plan_id}"
        for i in range(20):  # Poll for ~10 minutes
            time.sleep(30)
            st.info(f"  ...checking status (attempt {i+1}/20)")
            response = requests.get(status_url, timeout=60)
            if response.status_code == 200:
                data = response.json()
                if data.get("status") == "interventions_designed":
                    st.success("🎉 Intervention design complete!")
                    return {"guides": data.get("results", {}).get("guides", [])}
                elif data.get("status") == "design_failed":
                    error_msg = data.get('results', {}).get('error', 'Unknown error')
                    st.error(f"Intervention design failed: {error_msg}")
                    return {"error": "Design process failed."}
        
        st.warning("Polling timed out. The process is taking longer than expected.")
        return {"error": "Polling timed out."}

    except requests.exceptions.RequestException as e:
        st.error(f"API Error: Failed to connect to CommandCenter. Details: {e}")
        return {"error": str(e)}


class RUNX1DigitalTwinIntegrator:
    """
    Complete integration of RUNX1 tools into the Patient Digital Twin workflow.
    
    This class orchestrates all RUNX1-specific components and provides a unified
    interface for the Patient Digital Twin page.
    """
    
    def __init__(self):
        self.data_loader = RUNX1DataLoader()
        self.progression_modeler = RUNX1ProgressionModeler()
        # self.chopchop_integration = ChopChopIntegration() # ARCHIVED
        
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
    
    def _start_guide_design_campaign(self, battle_plan_id: int) -> bool:
        """Calls the command center to start the async guide design job."""
        url = f"{COMMAND_CENTER_URL}/v3/workflow/execute_guide_design_campaign/{battle_plan_id}"
        try:
            response = requests.post(url, timeout=30)
            response.raise_for_status()
            logger.info(f"Successfully started guide design campaign for battle plan {battle_plan_id}")
            return True
        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to start guide design campaign for BP {battle_plan_id}: {e}")
            return False

    def _poll_for_results(self, battle_plan_id: int, timeout: int = 1800) -> Optional[Dict]:
        """Polls the command center for the results of a guide design campaign."""
        url = f"{COMMAND_CENTER_URL}/v3/battle_plan/{battle_plan_id}"
        start_time = time.time()
        while time.time() - start_time < timeout:
            try:
                response = requests.get(url, timeout=30)
                response.raise_for_status()
                data = response.json()
                
                status = data.get("status")
                if status == "interventions_designed":
                    logger.info(f"Guide design complete for battle plan {battle_plan_id}")
                    return data.get("results")
                elif status == "design_failed":
                    logger.error(f"Guide design failed for battle plan {battle_plan_id}")
                    return {"error": "The guide design process failed."}
                else: # Still in progress
                    logger.info(f"Guide design for BP {battle_plan_id} is '{status}'. Waiting...")
                    time.sleep(10) # Wait for 10 seconds before polling again
            except requests.exceptions.RequestException as e:
                logger.error(f"Polling failed for BP {battle_plan_id}: {e}")
                time.sleep(15) # Wait longer on connection error
        
        logger.error(f"Polling timed out for battle plan {battle_plan_id}")
        return {"error": "Polling for results timed out."}


    def design_precision_intervention(self, target_variant: Dict, patient_id: str, intervention_type: str = "crispr_correction") -> Dict:
        """
        Designs a precision intervention by orchestrating the AI General Command Center.
        This is the new, refactored V3 workflow.
        """
        logger.info(f"V3 - Designing precision intervention for {target_variant.get('id', 'unknown')}")
        
        # Step 1: Get sequence context for the variant
        variant_pos = target_variant.get("position", 36207648)
        # Fetch a 200bp window around the variant for guide design
        context_data = self.data_loader.get_genomic_context(variant_pos, window=100)
        target_sequence = context_data.get("context_sequence")

        if not target_sequence:
            return {"error": "Could not retrieve target sequence for guide design."}

        # Step 2: Formulate a Battle Plan via the Command Center
        # The new endpoint is /v3/workflow/formulate_battle_plan
        # Construct the data packet required by the endpoint
        patient_data = {
            "patient_identifier": patient_id,
            "gene": target_variant.get("gene", "RUNX1"),
            "mutation_hgvs_c": target_variant.get("hgvs_c", "unknown"),
            "mutation_hgvs_p": target_variant.get("protein_change", "unknown"),
            "transcript_id": self.runx1_gene_info.get("transcript_id", "ENST00000546xxx"),
            "bam_file_path": "/path/to/patient/bam.bam", # This is a placeholder
            "protein_sequence": "PROTEINSEQ", # Placeholder
            "sequence_for_perplexity": target_sequence,
        }
        
        st.info("🛰️ Contacting AI General Command Center to formulate battle plan...")
        # This function needs to be updated or created to call the v3 endpoint
        plan_response = self._formulate_battle_plan_v3(patient_data)
        if not plan_response or "plan_id" not in plan_response:
            st.error("❌ Failed to formulate battle plan with Command Center.")
            return {"error": "Failed to formulate battle plan."}
        
        battle_plan_id = plan_response["plan_id"]
        st.success(f"✅ Battle Plan #{battle_plan_id} formulated. Now designing guides...")

        # Step 3: Start the asynchronous guide design campaign
        campaign_started = self._start_guide_design_campaign(battle_plan_id)
        if not campaign_started:
            st.error("❌ Failed to start the guide design campaign.")
            return {"error": "Failed to start guide design campaign."}

        # Step 4: Poll for the results
        with st.spinner(f"🔬 AI is designing and validating potential guides for Battle Plan #{battle_plan_id}... This may take several minutes."):
            results = self._poll_for_results(battle_plan_id)

        if not results or not results.get("guides"):
            st.error("❌ Guide design process did not return valid results.")
            return {"error": results.get("error", "Unknown error during guide design.")}

        # Step 5: Format results to match the UI's expectations
        st.success("🎉 Intervention design complete!")
        guides = results.get("guides", [])
        top_guide = guides[0] if guides else {}

        return {
            "target_variant": target_variant,
            "intervention_type": intervention_type,
            "guide_design": {
                "top_guide": top_guide,
                "all_guides": guides
            },
            "safety_assessment": {"summary": f"{len(guides)} guides validated. Top guide off-target score: {top_guide.get('off_target_score', 'N/A')}"},
            "experimental_protocol": generate_protocol(
                experiment_type="knock-in",
                target_gene=target_variant.get("gene", "RUNX1"),
                guide_info=top_guide,
                cell_type="CD34+ HSCs"
            ),
        }

    def _formulate_battle_plan_v3(self, patient_data: dict) -> Optional[dict]:
        """Calls the new V3 endpoint to formulate a battle plan."""
        url = f"{COMMAND_CENTER_URL}/v3/workflow/formulate_battle_plan"
        try:
            response = requests.post(url, json=patient_data, headers={"Content-Type": "application/json"}, timeout=30)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to formulate V3 battle plan: {e}")
            return None

    def create_genomic_browser_data(self, variant_position: int, window: int = 5000, guides: list = None) -> Dict:
        """
        Create data for interactive genomic browser visualization using gosling.
        """
        logger.info(f"Creating gosling genomic browser for position {variant_position}")
        
        # 1. Get sequence context and gene structure
        context = self.data_loader.get_genomic_context(variant_position, window)
        gene_structure = self.data_loader.get_runx1_gene_info()
        
        # 2. Prepare tracks for Gosling
        
        # Track 1: Gene Annotations (Exons)
        exons = [{
            "start": start, "end": end, "name": f"Exon {i+1}", "strand": "+"
        } for i, (start, end) in enumerate(gene_structure.get("exons", []))]
        
        gene_track = gos.Track(
            data=gos.Data(
                values=exons,
                type='json',
                genomicFields=["start", "end"]
            )
        ).mark_rect().encode(
            x=gos.X("start:G", axis="bottom"),
            xe="end:G",
            color=gos.Color("name:N", legend=True),
            tooltip=[
                gos.Tooltip(field="start", type="genomic", alt="Start Position"),
                gos.Tooltip(field="end", type="genomic", alt="End Position"),
                gos.Tooltip(field="name", type="nominal", alt="Annotation"),
            ]
        ).properties(title="RUNX1 Gene Structure")

        # Track 2: Variant and Guide Positions
        points_data = [{
            "position": variant_position, "label": "Patient Mutation", "value": "Mutation", "size": 20
        }]
        if guides:
            for guide in guides:
                # Placeholder for guide position calculation
                guide_seq = guide.get("guide_sequence", "")
                if guide_seq in context.get("context_sequence", ""):
                    guide_start = variant_position - (window // 2) + context.get("context_sequence", "").find(guide_seq)
                    points_data.append({
                        "position": guide_start,
                        "label": f"Guide ({guide.get('assassin_score', 0):.2f})",
                        "value": "Guide RNA",
                        "size": 10
                    })

        points_track = gos.Track(
            data=gos.Data(
                values=points_data,
                type='json',
                genomicFields=["position"]
            )
        ).mark_point().encode(
            x="position:G",
            y=gos.Y("value:N"),
            color=gos.Color("value:N", domain=["Mutation", "Guide RNA"], range=["red", "purple"]),
            size="size:Q",
            tooltip=[
                gos.Tooltip(field="position", type="genomic", alt="Position"),
                gos.Tooltip(field="label", type="nominal", alt="Label"),
            ]
        ).properties(title="Key Positions")

        # Combine tracks into a final visualization spec
        combined_view = gos.stack(gene_track, points_track).properties(
            title="RUNX1 Genomic Browser",
            centerRadius=0.8,
            layout="linear",
            spacing=0,
            xDomain=gos.Domain(chromosome="chr21", interval=[variant_position - window, variant_position + window])
        )

        return combined_view.spec()

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

# This is the new top-level function that should be called by the UI
def design_runx1_intervention_v3(target_variant: Dict, patient_id: str) -> Dict:
    """Entry point for the V3 RUNX1 intervention design workflow."""
    integrator = RUNX1DigitalTwinIntegrator()
    return integrator.design_precision_intervention(target_variant, patient_id)

def create_runx1_genomic_browser(variant_position: int, window: int = 5000, guides: list = None) -> Dict:
    """Helper function to call the integrator's browser data method."""
    integrator = RUNX1DigitalTwinIntegrator()
    return integrator.create_genomic_browser_data(variant_position, window, guides)

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
        # This part needs a patient_id for the new design_precision_intervention
        # For now, we'll pass a placeholder or skip if not available
        # Assuming a placeholder patient_id for testing
        test_patient_id = "TEST_PT_001" 
        intervention = design_runx1_intervention_v3(test_germline, test_patient_id)
        print(f"Intervention designed: {intervention['success_probability']:.3f} success probability")
        
        # Test clinical report
        report = integrator.generate_clinical_report(results)
        print(f"Clinical report generated: {report['risk_assessment']['risk_tier']} risk tier")
        
    except Exception as e:
        print(f"Integration test failed: {e}")
        logger.error(f"RUNX1 integration test failed: {e}") 