"""
RUNX1 Progression Modeler for RUNX1-FPD Patient Digital Twin

This module implements sophisticated two-hit progression modeling, clonal evolution
analysis, and intervention opportunity identification for RUNX1-FPD patients.
"""

import logging
import json
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from datetime import datetime, timedelta
import math

# Remove problematic import
# from tools.command_center_client import CommandCenterClient
from runx1_data_loader import RUNX1DataLoader

logger = logging.getLogger(__name__)

@dataclass
class ProgressionEvent:
    """Represents a single event in clonal evolution progression."""
    event_type: str
    time_point: float  # Years from germline mutation
    mutation_id: str
    clone_fraction: float
    fitness_advantage: float
    intervention_window: bool
    description: str

@dataclass
class InterventionOpportunity:
    """Represents an intervention opportunity in progression."""
    time_window: Tuple[float, float]  # Start and end time in years
    intervention_type: str
    target_mutation: str
    success_probability: float
    clinical_feasibility: str
    description: str
    urgency: str

class RUNX1ProgressionModeler:
    """
    Sophisticated progression modeling for RUNX1-FPD patients.
    
    This class implements the two-hit hypothesis with clonal evolution dynamics,
    modeling the progression from germline predisposition to malignant transformation.
    """
    
    def __init__(self):
        # Remove problematic initialization
        # self.command_center = CommandCenterClient()
        self.data_loader = RUNX1DataLoader()
        
        # Progression model parameters
        self.model_params = {
            "baseline_mutation_rate": 1e-9,  # Mutations per base per year
            "runx1_mutation_advantage": 1.5,  # Fitness advantage from RUNX1 loss
            "clonal_expansion_rate": 0.1,  # Per year
            "transformation_threshold": 0.1,  # Clone fraction for transformation
            "age_acceleration_factor": 1.02,  # Exponential age effect
            "microenvironment_factor": 1.0,  # Bone marrow microenvironment effect
        }
        
        # Known high-risk second-hit genes for RUNX1-FPD
        self.high_risk_second_hits = {
            "ASXL1": {
                "fitness_advantage": 2.0,
                "transformation_risk": 0.8,
                "typical_latency": 5.0,  # Years
                "mechanism": "Epigenetic dysregulation"
            },
            "TET2": {
                "fitness_advantage": 1.8,
                "transformation_risk": 0.7,
                "typical_latency": 7.0,
                "mechanism": "DNA demethylation disruption"
            },
            "CBL": {
                "fitness_advantage": 1.6,
                "transformation_risk": 0.6,
                "typical_latency": 6.0,
                "mechanism": "Growth factor signaling"
            },
            "NRAS": {
                "fitness_advantage": 2.2,
                "transformation_risk": 0.9,
                "typical_latency": 3.0,
                "mechanism": "RAS pathway activation"
            },
            "FLT3": {
                "fitness_advantage": 1.9,
                "transformation_risk": 0.8,
                "typical_latency": 4.0,
                "mechanism": "Receptor tyrosine kinase activation"
            }
        }
    
    def model_clonal_evolution(self, germline_variant: Dict, 
                              somatic_variants: List[Dict] = None,
                              patient_age: float = 35.0) -> Dict:
        """
        Model complete clonal evolution from germline to malignant transformation.
        
        Args:
            germline_variant: Germline RUNX1 variant
            somatic_variants: List of somatic variants (if known)
            patient_age: Current patient age
            
        Returns:
            Dict: Complete clonal evolution model
        """
        logger.info(f"Modeling clonal evolution for RUNX1 variant {germline_variant.get('id', 'unknown')}")
        
        # Initialize progression timeline
        progression_events = []
        
        # Event 1: Germline mutation (time = 0)
        germline_event = ProgressionEvent(
            event_type="germline_mutation",
            time_point=0.0,
            mutation_id=germline_variant.get("id", "RUNX1_germline"),
            clone_fraction=1.0,  # All cells carry germline mutation
            fitness_advantage=self.model_params["runx1_mutation_advantage"],
            intervention_window=True,
            description="Germline RUNX1 mutation establishes predisposition"
        )
        progression_events.append(germline_event)
        
        # Model somatic mutation acquisition
        if somatic_variants:
            for variant in somatic_variants:
                somatic_event = self._model_somatic_acquisition(
                    variant, patient_age, germline_variant
                )
                progression_events.append(somatic_event)
        else:
            # Predict likely somatic mutations
            predicted_events = self._predict_somatic_mutations(
                germline_variant, patient_age
            )
            progression_events.extend(predicted_events)
        
        # Calculate transformation probability
        transformation_prob = self._calculate_transformation_probability(
            progression_events, patient_age
        )
        
        # Identify intervention opportunities
        intervention_opportunities = self._identify_intervention_opportunities(
            progression_events, patient_age
        )
        
        return {
            "patient_age": patient_age,
            "germline_variant": germline_variant,
            "somatic_variants": somatic_variants or [],
            "progression_events": [self._event_to_dict(e) for e in progression_events],
            "transformation_probability": transformation_prob,
            "intervention_opportunities": [self._opportunity_to_dict(o) for o in intervention_opportunities],
            "model_predictions": self._generate_model_predictions(progression_events, patient_age),
            "clinical_timeline": self._generate_clinical_timeline(progression_events, patient_age)
        }
    
    def _model_somatic_acquisition(self, variant: Dict, patient_age: float, 
                                  germline_variant: Dict) -> ProgressionEvent:
        """
        Model the acquisition of a somatic mutation.
        
        Args:
            variant: Somatic variant
            patient_age: Patient age
            germline_variant: Germline variant context
            
        Returns:
            ProgressionEvent: Somatic acquisition event
        """
        # Estimate acquisition time based on variant type
        gene_symbol = variant.get("gene", "unknown")
        
        if gene_symbol in self.high_risk_second_hits:
            gene_info = self.high_risk_second_hits[gene_symbol]
            typical_latency = gene_info["typical_latency"]
            fitness_advantage = gene_info["fitness_advantage"]
        else:
            typical_latency = 10.0  # Default latency
            fitness_advantage = 1.3  # Default advantage
        
        # Add age-related acceleration
        age_factor = (patient_age / 50.0) ** 2  # Quadratic age effect
        acquisition_time = typical_latency / (1 + age_factor)
        
        # Estimate clone fraction at detection
        years_since_acquisition = max(0.1, patient_age - acquisition_time)
        clone_fraction = min(0.5, 0.01 * math.exp(
            self.model_params["clonal_expansion_rate"] * years_since_acquisition
        ))
        
        return ProgressionEvent(
            event_type="somatic_mutation",
            time_point=acquisition_time,
            mutation_id=variant.get("id", f"{gene_symbol}_somatic"),
            clone_fraction=clone_fraction,
            fitness_advantage=fitness_advantage,
            intervention_window=clone_fraction < self.model_params["transformation_threshold"],
            description=f"Somatic {gene_symbol} mutation acquired"
        )
    
    def _predict_somatic_mutations(self, germline_variant: Dict, 
                                  patient_age: float) -> List[ProgressionEvent]:
        """
        Predict likely somatic mutations based on germline context.
        
        Args:
            germline_variant: Germline variant
            patient_age: Patient age
            
        Returns:
            List[ProgressionEvent]: Predicted somatic events
        """
        predicted_events = []
        
        # Calculate risk for each potential second hit
        for gene, gene_info in self.high_risk_second_hits.items():
            # Base probability modified by age and germline context
            base_prob = 0.1  # 10% base probability
            age_modifier = (patient_age / 50.0) ** 1.5
            
            # Higher risk for certain germline mutations
            germline_modifier = 1.0
            if "frameshift" in germline_variant.get("consequence", "").lower():
                germline_modifier = 1.5
            elif "nonsense" in germline_variant.get("consequence", "").lower():
                germline_modifier = 1.4
            
            mutation_prob = base_prob * age_modifier * germline_modifier
            
            if mutation_prob > 0.05:  # 5% threshold for inclusion
                # Predict acquisition time
                acquisition_time = gene_info["typical_latency"] * (1 + np.random.normal(0, 0.3))
                acquisition_time = max(1.0, acquisition_time)  # Minimum 1 year
                
                # Only include if within reasonable timeframe
                if acquisition_time < patient_age + 20:
                    predicted_event = ProgressionEvent(
                        event_type="predicted_somatic",
                        time_point=acquisition_time,
                        mutation_id=f"{gene}_predicted",
                        clone_fraction=0.0,  # Not yet acquired
                        fitness_advantage=gene_info["fitness_advantage"],
                        intervention_window=True,
                        description=f"Predicted {gene} mutation (risk: {mutation_prob:.1%})"
                    )
                    predicted_events.append(predicted_event)
        
        return predicted_events
    
    def _calculate_transformation_probability(self, progression_events: List[ProgressionEvent], 
                                           patient_age: float) -> Dict:
        """
        Calculate malignant transformation probability.
        
        Args:
            progression_events: List of progression events
            patient_age: Patient age
            
        Returns:
            Dict: Transformation probability analysis
        """
        # Find somatic mutations
        somatic_events = [e for e in progression_events if e.event_type == "somatic_mutation"]
        
        if not somatic_events:
            # No somatic mutations detected
            base_risk = 0.1 * (patient_age / 50.0) ** 2  # Age-dependent baseline
            return {
                "current_probability": base_risk,
                "5_year_probability": base_risk * 2,
                "10_year_probability": base_risk * 4,
                "risk_category": "moderate" if base_risk > 0.1 else "low",
                "primary_risk_factors": ["germline RUNX1 mutation", "age"]
            }
        
        # Calculate combined risk from somatic mutations
        combined_fitness = 1.0
        transformation_drivers = []
        
        for event in somatic_events:
            combined_fitness *= event.fitness_advantage
            if event.clone_fraction > 0.05:  # 5% clone fraction threshold
                transformation_drivers.append(event.mutation_id)
        
        # Calculate transformation probability
        current_prob = min(0.9, combined_fitness * 0.1)
        
        # Project future risk
        year_5_prob = min(0.95, current_prob * 1.5)
        year_10_prob = min(0.98, current_prob * 2.0)
        
        # Determine risk category
        if current_prob > 0.5:
            risk_category = "very_high"
        elif current_prob > 0.3:
            risk_category = "high"
        elif current_prob > 0.1:
            risk_category = "moderate"
        else:
            risk_category = "low"
        
        return {
            "current_probability": current_prob,
            "5_year_probability": year_5_prob,
            "10_year_probability": year_10_prob,
            "risk_category": risk_category,
            "transformation_drivers": transformation_drivers,
            "combined_fitness_advantage": combined_fitness,
            "primary_risk_factors": ["germline RUNX1 mutation"] + transformation_drivers
        }
    
    def _identify_intervention_opportunities(self, progression_events: List[ProgressionEvent], 
                                          patient_age: float) -> List[InterventionOpportunity]:
        """
        Identify intervention opportunities in the progression timeline.
        
        Args:
            progression_events: Progression events
            patient_age: Patient age
            
        Returns:
            List[InterventionOpportunity]: Intervention opportunities
        """
        opportunities = []
        
        # Germline correction opportunity (always available)
        germline_opportunity = InterventionOpportunity(
            time_window=(0.0, patient_age + 50),  # Lifelong window
            intervention_type="germline_correction",
            target_mutation="RUNX1_germline",
            success_probability=0.7,  # Based on current CRISPR capabilities
            clinical_feasibility="experimental",
            description="CRISPR-mediated correction of germline RUNX1 mutation in HSCs",
            urgency="high" if patient_age > 40 else "moderate"
        )
        opportunities.append(germline_opportunity)
        
        # Somatic mutation targeting opportunities
        somatic_events = [e for e in progression_events if e.event_type == "somatic_mutation"]
        
        for event in somatic_events:
            if event.intervention_window:
                # Window for intervention before transformation
                window_start = event.time_point
                window_end = event.time_point + 5.0  # 5-year window
                
                # Success probability based on clone fraction
                success_prob = max(0.3, 0.9 - event.clone_fraction)
                
                # Determine urgency
                if event.clone_fraction > 0.05:
                    urgency = "urgent"
                elif event.clone_fraction > 0.01:
                    urgency = "high"
                else:
                    urgency = "moderate"
                
                opportunity = InterventionOpportunity(
                    time_window=(window_start, window_end),
                    intervention_type="somatic_targeting",
                    target_mutation=event.mutation_id,
                    success_probability=success_prob,
                    clinical_feasibility="investigational",
                    description=f"Target {event.mutation_id} before clonal expansion",
                    urgency=urgency
                )
                opportunities.append(opportunity)
        
        # Predicted mutation prevention opportunities
        predicted_events = [e for e in progression_events if e.event_type == "predicted_somatic"]
        
        for event in predicted_events:
            # Prevention window before acquisition
            prevention_opportunity = InterventionOpportunity(
                time_window=(patient_age, event.time_point),
                intervention_type="mutation_prevention",
                target_mutation=event.mutation_id,
                success_probability=0.5,  # Speculative
                clinical_feasibility="research",
                description=f"Prevent acquisition of {event.mutation_id}",
                urgency="low"
            )
            opportunities.append(prevention_opportunity)
        
        return opportunities
    
    def _generate_model_predictions(self, progression_events: List[ProgressionEvent], 
                                  patient_age: float) -> Dict:
        """
        Generate model predictions and insights.
        
        Args:
            progression_events: Progression events
            patient_age: Patient age
            
        Returns:
            Dict: Model predictions
        """
        somatic_events = [e for e in progression_events if e.event_type == "somatic_mutation"]
        predicted_events = [e for e in progression_events if e.event_type == "predicted_somatic"]
        
        return {
            "next_likely_mutation": predicted_events[0].mutation_id if predicted_events else None,
            "time_to_next_mutation": predicted_events[0].time_point - patient_age if predicted_events else None,
            "current_clonal_burden": sum(e.clone_fraction for e in somatic_events),
            "progression_velocity": len(somatic_events) / max(1, patient_age - 20),  # Mutations per year since age 20
            "intervention_urgency": "high" if any(e.clone_fraction > 0.05 for e in somatic_events) else "moderate",
            "optimal_intervention_target": somatic_events[0].mutation_id if somatic_events else "RUNX1_germline",
            "model_confidence": 0.8 if somatic_events else 0.6  # Higher confidence with actual data
        }
    
    def _generate_clinical_timeline(self, progression_events: List[ProgressionEvent], 
                                  patient_age: float) -> List[Dict]:
        """
        Generate clinical timeline with key milestones.
        
        Args:
            progression_events: Progression events
            patient_age: Patient age
            
        Returns:
            List[Dict]: Clinical timeline
        """
        timeline = []
        
        # Current status
        timeline.append({
            "timepoint": patient_age,
            "event_type": "current_status",
            "description": "Current patient assessment",
            "clinical_action": "Comprehensive genomic analysis and risk stratification",
            "monitoring_frequency": "every 6 months"
        })
        
        # Future monitoring points
        for i in range(1, 6):  # Next 5 years
            future_age = patient_age + i
            timeline.append({
                "timepoint": future_age,
                "event_type": "surveillance",
                "description": f"Year {i} surveillance",
                "clinical_action": "CBC, bone marrow biopsy if indicated",
                "monitoring_frequency": "annual"
            })
        
        # Intervention milestones
        somatic_events = [e for e in progression_events if e.event_type == "somatic_mutation"]
        if somatic_events:
            next_intervention = patient_age + 0.5  # 6 months
            timeline.append({
                "timepoint": next_intervention,
                "event_type": "intervention_planning",
                "description": "Plan therapeutic intervention",
                "clinical_action": "Multidisciplinary team consultation",
                "monitoring_frequency": "as needed"
            })
        
        return sorted(timeline, key=lambda x: x["timepoint"])
    
    def predict_progression_timeline(self, patient_data: Dict) -> Dict:
        """
        Predict likely progression timeline for a patient.
        
        Args:
            patient_data: Patient data including age, variants, family history
            
        Returns:
            Dict: Predicted progression timeline
        """
        germline_variant = patient_data.get("germline_variant", {})
        somatic_variants = patient_data.get("somatic_variants", [])
        patient_age = patient_data.get("age", 35.0)
        family_history = patient_data.get("family_history", {})
        
        # Run clonal evolution model
        evolution_model = self.model_clonal_evolution(
            germline_variant, somatic_variants, patient_age
        )
        
        # Add family history modifiers
        if family_history.get("leukemia_cases", 0) > 1:
            # Multiply risks by 1.5 for strong family history
            evolution_model["transformation_probability"]["current_probability"] *= 1.5
            evolution_model["transformation_probability"]["5_year_probability"] *= 1.5
            evolution_model["transformation_probability"]["10_year_probability"] *= 1.5
        
        return {
            "patient_id": patient_data.get("patient_id", "unknown"),
            "progression_model": evolution_model,
            "risk_stratification": self._stratify_patient_risk(evolution_model),
            "surveillance_recommendations": self._generate_surveillance_plan(evolution_model),
            "intervention_priorities": self._prioritize_interventions(evolution_model)
        }
    
    def identify_intervention_opportunities(self, progression_model: Dict) -> List[Dict]:
        """
        Identify and prioritize intervention opportunities.
        
        Args:
            progression_model: Progression model results
            
        Returns:
            List[Dict]: Prioritized intervention opportunities
        """
        opportunities = progression_model.get("intervention_opportunities", [])
        
        # Sort by urgency and success probability
        def priority_score(opp):
            urgency_scores = {"urgent": 4, "high": 3, "moderate": 2, "low": 1}
            urgency = urgency_scores.get(opp.get("urgency", "low"), 1)
            success_prob = opp.get("success_probability", 0.5)
            return urgency * success_prob
        
        sorted_opportunities = sorted(opportunities, key=priority_score, reverse=True)
        
        # Add priority rankings
        for i, opp in enumerate(sorted_opportunities):
            opp["priority_rank"] = i + 1
            opp["priority_score"] = priority_score(opp)
        
        return sorted_opportunities
    
    def _stratify_patient_risk(self, evolution_model: Dict) -> Dict:
        """Stratify patient risk based on evolution model."""
        transformation_prob = evolution_model["transformation_probability"]
        current_prob = transformation_prob["current_probability"]
        
        if current_prob > 0.5:
            risk_tier = "very_high"
            surveillance_frequency = "monthly"
            intervention_urgency = "immediate"
        elif current_prob > 0.3:
            risk_tier = "high"
            surveillance_frequency = "quarterly"
            intervention_urgency = "urgent"
        elif current_prob > 0.1:
            risk_tier = "moderate"
            surveillance_frequency = "every 6 months"
            intervention_urgency = "high"
        else:
            risk_tier = "low"
            surveillance_frequency = "annually"
            intervention_urgency = "moderate"
        
        return {
            "risk_tier": risk_tier,
            "surveillance_frequency": surveillance_frequency,
            "intervention_urgency": intervention_urgency,
            "transformation_probability": current_prob,
            "risk_factors": transformation_prob.get("primary_risk_factors", [])
        }
    
    def _generate_surveillance_plan(self, evolution_model: Dict) -> Dict:
        """Generate surveillance plan based on evolution model."""
        risk_stratification = self._stratify_patient_risk(evolution_model)
        
        return {
            "frequency": risk_stratification["surveillance_frequency"],
            "tests": [
                "Complete blood count with differential",
                "Peripheral blood smear",
                "Flow cytometry for blast assessment",
                "Bone marrow biopsy if blasts >5%"
            ],
            "monitoring_parameters": [
                "Platelet count",
                "Blast percentage",
                "Clonal evolution markers",
                "Cytogenetic abnormalities"
            ],
            "escalation_criteria": [
                "Platelet count <50,000",
                "Blast percentage >5%",
                "New cytogenetic abnormalities",
                "Rapid clonal expansion"
            ]
        }
    
    def _prioritize_interventions(self, evolution_model: Dict) -> List[Dict]:
        """Prioritize interventions based on evolution model."""
        opportunities = evolution_model.get("intervention_opportunities", [])
        
        prioritized = []
        for opp in opportunities:
            priority = {
                "intervention": opp["intervention_type"],
                "target": opp["target_mutation"],
                "urgency": opp["urgency"],
                "success_probability": opp["success_probability"],
                "clinical_feasibility": opp["clinical_feasibility"],
                "time_window": opp["time_window"],
                "description": opp["description"]
            }
            prioritized.append(priority)
        
        return prioritized
    
    def _event_to_dict(self, event: ProgressionEvent) -> Dict:
        """Convert ProgressionEvent to dictionary."""
        return {
            "event_type": event.event_type,
            "time_point": event.time_point,
            "mutation_id": event.mutation_id,
            "clone_fraction": event.clone_fraction,
            "fitness_advantage": event.fitness_advantage,
            "intervention_window": event.intervention_window,
            "description": event.description
        }
    
    def _opportunity_to_dict(self, opportunity: InterventionOpportunity) -> Dict:
        """Convert InterventionOpportunity to dictionary."""
        return {
            "time_window": opportunity.time_window,
            "intervention_type": opportunity.intervention_type,
            "target_mutation": opportunity.target_mutation,
            "success_probability": opportunity.success_probability,
            "clinical_feasibility": opportunity.clinical_feasibility,
            "description": opportunity.description,
            "urgency": opportunity.urgency
        }

# Convenience functions for easy integration
def model_runx1_progression(germline_variant: Dict, somatic_variants: List[Dict] = None, 
                          patient_age: float = 35.0) -> Dict:
    """
    Convenience function to model RUNX1 progression.
    
    Args:
        germline_variant: Germline RUNX1 variant
        somatic_variants: Somatic variants (optional)
        patient_age: Patient age
        
    Returns:
        Dict: Progression model results
    """
    modeler = RUNX1ProgressionModeler()
    return modeler.model_clonal_evolution(germline_variant, somatic_variants, patient_age)

def predict_runx1_timeline(patient_data: Dict) -> Dict:
    """
    Convenience function to predict RUNX1 progression timeline.
    
    Args:
        patient_data: Complete patient data
        
    Returns:
        Dict: Predicted timeline
    """
    modeler = RUNX1ProgressionModeler()
    return modeler.predict_progression_timeline(patient_data)

if __name__ == "__main__":
    # Test the progression modeler
    modeler = RUNX1ProgressionModeler()
    
    # Test with sample data
    test_germline = {
        "id": "RUNX1_R204Q",
        "gene": "RUNX1",
        "consequence": "missense_variant",
        "position": 36207648
    }
    
    test_somatic = [
        {
            "id": "ASXL1_G646fs",
            "gene": "ASXL1",
            "consequence": "frameshift_variant",
            "position": 31022441
        }
    ]
    
    print("Testing RUNX1 progression modeler...")
    
    try:
        # Test clonal evolution modeling
        evolution_model = modeler.model_clonal_evolution(
            test_germline, test_somatic, patient_age=45.0
        )
        
        print(f"Transformation probability: {evolution_model['transformation_probability']['current_probability']:.3f}")
        print(f"Intervention opportunities: {len(evolution_model['intervention_opportunities'])}")
        
        # Test timeline prediction
        patient_data = {
            "patient_id": "test_patient",
            "germline_variant": test_germline,
            "somatic_variants": test_somatic,
            "age": 45.0,
            "family_history": {"leukemia_cases": 2}
        }
        
        timeline = modeler.predict_progression_timeline(patient_data)
        print(f"Risk tier: {timeline['risk_stratification']['risk_tier']}")
        print(f"Surveillance frequency: {timeline['risk_stratification']['surveillance_frequency']}")
        
    except Exception as e:
        print(f"Test failed: {e}")
        logger.error(f"RUNX1 progression modeler test failed: {e}") 