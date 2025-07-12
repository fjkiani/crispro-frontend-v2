#!/usr/bin/env python3
"""
RUNX1-FPD Demo Scenarios
Compelling patient scenarios for LEAP Grant demonstration
"""

import streamlit as st
import pandas as pd
import json
from datetime import datetime, timedelta
from typing import Dict, List, Optional
import random
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class RUNX1DemoScenarios:
    """
    Comprehensive demo scenarios for RUNX1-FPD precision medicine platform.
    """
    
    def __init__(self):
        """Initialize demo scenarios."""
        self.scenarios = self._create_demo_scenarios()
    
    def _create_demo_scenarios(self) -> List[Dict]:
        """Create compelling demo scenarios."""
        
        scenarios = [
            {
                "id": "scenario_1_high_risk",
                "title": "🚨 High-Risk Transformation Case",
                "subtitle": "Sarah Chen - 28-year-old with rapid progression",
                "patient_info": {
                    "name": "Sarah Chen",
                    "age": 28,
                    "gender": "Female",
                    "ethnicity": "Asian",
                    "family_history": "Mother diagnosed with AML at 45",
                    "presentation": "Thrombocytopenia, easy bruising, fatigue"
                },
                "clinical_timeline": [
                    {
                        "date": "2023-01-15",
                        "event": "Initial presentation",
                        "details": "Platelet count 45,000/μL, easy bruising for 6 months"
                    },
                    {
                        "date": "2023-02-10",
                        "event": "Genetic testing",
                        "details": "RUNX1 germline mutation identified: p.Arg135fs"
                    },
                    {
                        "date": "2023-08-20",
                        "event": "Progression detected",
                        "details": "Blast count 8%, new cytogenetic abnormalities"
                    },
                    {
                        "date": "2023-09-05",
                        "event": "Somatic mutation found",
                        "details": "ASXL1 p.Gly646fs detected in bone marrow"
                    },
                    {
                        "date": "2023-09-15",
                        "event": "AI analysis completed",
                        "details": "High-risk profile confirmed, intervention recommended"
                    }
                ],
                "germline_variant": {
                    "gene": "RUNX1",
                    "protein_change": "p.Arg135fs",
                    "variant_type": "germline",
                    "id": "RUNX1_p.Arg135fs",
                    "position": 36207748,
                    "consequence": "frameshift_variant",
                    "pathogenicity": "pathogenic",
                    "functional_impact": "high"
                },
                "somatic_variant": {
                    "gene": "ASXL1",
                    "protein_change": "p.Gly646fs",
                    "variant_type": "somatic",
                    "id": "ASXL1_p.Gly646fs",
                    "position": 31022441,
                    "consequence": "frameshift_variant",
                    "pathogenicity": "pathogenic",
                    "functional_impact": "high"
                },
                "lab_values": {
                    "platelet_count": 45000,
                    "hemoglobin": 9.8,
                    "white_blood_cells": 12500,
                    "blast_percentage": 8,
                    "bone_marrow_cellularity": 85
                },
                "family_pedigree": {
                    "mother": {"status": "AML at 45", "runx1_status": "positive"},
                    "father": {"status": "healthy", "runx1_status": "negative"},
                    "sister": {"status": "healthy", "runx1_status": "untested"},
                    "maternal_grandmother": {"status": "bleeding disorder", "runx1_status": "unknown"}
                },
                "expected_analysis": {
                    "risk_tier": "high",
                    "transformation_probability": 0.75,
                    "surveillance_frequency": "monthly",
                    "intervention_urgency": "immediate"
                },
                "demo_narrative": {
                    "act_1": "Meet Sarah - a young professional with concerning symptoms",
                    "act_2": "AI discovers high-risk two-hit mutation profile",
                    "act_3": "Precision intervention designed to prevent transformation",
                    "act_4": "LEAP Grant impact - from reactive to proactive medicine"
                }
            },
            
            {
                "id": "scenario_2_moderate_risk",
                "title": "⚠️ Moderate-Risk Surveillance Case",
                "subtitle": "Michael Rodriguez - 45-year-old with stable disease",
                "patient_info": {
                    "name": "Michael Rodriguez",
                    "age": 45,
                    "gender": "Male",
                    "ethnicity": "Hispanic",
                    "family_history": "Father had thrombocytopenia",
                    "presentation": "Mild thrombocytopenia, routine screening"
                },
                "clinical_timeline": [
                    {
                        "date": "2022-06-10",
                        "event": "Family screening",
                        "details": "Tested due to daughter's RUNX1-FPD diagnosis"
                    },
                    {
                        "date": "2022-07-15",
                        "event": "RUNX1 mutation confirmed",
                        "details": "Germline RUNX1 p.Arg174Gln identified"
                    },
                    {
                        "date": "2023-01-20",
                        "event": "Surveillance bone marrow",
                        "details": "Normal cellularity, no excess blasts"
                    },
                    {
                        "date": "2023-09-10",
                        "event": "Clonal hematopoiesis detected",
                        "details": "Small TET2 clone identified (VAF 3%)"
                    },
                    {
                        "date": "2023-09-20",
                        "event": "AI risk assessment",
                        "details": "Moderate risk, enhanced surveillance recommended"
                    }
                ],
                "germline_variant": {
                    "gene": "RUNX1",
                    "protein_change": "p.Arg174Gln",
                    "variant_type": "germline",
                    "id": "RUNX1_p.Arg174Gln",
                    "position": 36207865,
                    "consequence": "missense_variant",
                    "pathogenicity": "pathogenic",
                    "functional_impact": "moderate"
                },
                "somatic_variant": {
                    "gene": "TET2",
                    "protein_change": "p.Gln1299*",
                    "variant_type": "somatic",
                    "id": "TET2_p.Gln1299*",
                    "position": 106155185,
                    "consequence": "nonsense_variant",
                    "pathogenicity": "likely_pathogenic",
                    "functional_impact": "moderate"
                },
                "lab_values": {
                    "platelet_count": 85000,
                    "hemoglobin": 12.1,
                    "white_blood_cells": 6800,
                    "blast_percentage": 2,
                    "bone_marrow_cellularity": 45
                },
                "family_pedigree": {
                    "daughter": {"status": "RUNX1-FPD", "runx1_status": "positive"},
                    "wife": {"status": "healthy", "runx1_status": "negative"},
                    "son": {"status": "healthy", "runx1_status": "negative"},
                    "father": {"status": "thrombocytopenia", "runx1_status": "unknown"}
                },
                "expected_analysis": {
                    "risk_tier": "moderate",
                    "transformation_probability": 0.25,
                    "surveillance_frequency": "every 6 months",
                    "intervention_urgency": "high"
                },
                "demo_narrative": {
                    "act_1": "Family screening reveals hidden risk",
                    "act_2": "AI detects early clonal evolution",
                    "act_3": "Precision surveillance prevents progression",
                    "act_4": "Family-based precision medicine approach"
                }
            },
            
            {
                "id": "scenario_3_pediatric",
                "title": "👶 Pediatric Early Detection Case",
                "subtitle": "Emma Johnson - 12-year-old with family history",
                "patient_info": {
                    "name": "Emma Johnson",
                    "age": 12,
                    "gender": "Female",
                    "ethnicity": "Caucasian",
                    "family_history": "Father with RUNX1-FPD, grandfather with leukemia",
                    "presentation": "Asymptomatic, family screening"
                },
                "clinical_timeline": [
                    {
                        "date": "2023-03-15",
                        "event": "Family screening",
                        "details": "Tested due to father's RUNX1-FPD diagnosis"
                    },
                    {
                        "date": "2023-04-10",
                        "event": "RUNX1 mutation confirmed",
                        "details": "Inherited father's RUNX1 p.Arg201Cys mutation"
                    },
                    {
                        "date": "2023-04-20",
                        "event": "Baseline assessment",
                        "details": "Normal blood counts, no symptoms"
                    },
                    {
                        "date": "2023-09-15",
                        "event": "AI risk modeling",
                        "details": "Long-term risk assessment for intervention planning"
                    }
                ],
                "germline_variant": {
                    "gene": "RUNX1",
                    "protein_change": "p.Arg201Cys",
                    "variant_type": "germline",
                    "id": "RUNX1_p.Arg201Cys",
                    "position": 36207946,
                    "consequence": "missense_variant",
                    "pathogenicity": "pathogenic",
                    "functional_impact": "high"
                },
                "somatic_variant": None,  # No somatic mutations yet
                "lab_values": {
                    "platelet_count": 180000,
                    "hemoglobin": 12.8,
                    "white_blood_cells": 7200,
                    "blast_percentage": 0,
                    "bone_marrow_cellularity": "not_performed"
                },
                "family_pedigree": {
                    "father": {"status": "RUNX1-FPD", "runx1_status": "positive"},
                    "mother": {"status": "healthy", "runx1_status": "negative"},
                    "brother": {"status": "healthy", "runx1_status": "negative"},
                    "paternal_grandfather": {"status": "leukemia at 62", "runx1_status": "presumed_positive"}
                },
                "expected_analysis": {
                    "risk_tier": "low_current_high_lifetime",
                    "transformation_probability": 0.05,
                    "surveillance_frequency": "annually",
                    "intervention_urgency": "preventive"
                },
                "demo_narrative": {
                    "act_1": "Early detection through family screening",
                    "act_2": "AI models lifetime risk trajectory",
                    "act_3": "Preventive intervention planning",
                    "act_4": "Next-generation precision medicine"
                }
            }
        ]
        
        return scenarios
    
    def get_scenario(self, scenario_id: str) -> Optional[Dict]:
        """Get specific scenario by ID."""
        for scenario in self.scenarios:
            if scenario["id"] == scenario_id:
                return scenario
        return None
    
    def get_all_scenarios(self) -> List[Dict]:
        """Get all available scenarios."""
        return self.scenarios
    
    def create_demo_flow(self, scenario_id: str) -> Dict:
        """Create complete demo flow for a scenario."""
        scenario = self.get_scenario(scenario_id)
        if not scenario:
            return {}
        
        return {
            "scenario": scenario,
            "demo_script": self._create_demo_script(scenario),
            "expected_results": self._create_expected_results(scenario),
            "talking_points": self._create_talking_points(scenario),
            "technical_details": self._create_technical_details(scenario)
        }
    
    def _create_demo_script(self, scenario: Dict) -> Dict:
        """Create demo script for the scenario."""
        return {
            "introduction": f"""
            Welcome to our RUNX1-FPD Precision Medicine Platform demonstration.
            
            Today we'll follow {scenario['patient_info']['name']}, a {scenario['patient_info']['age']}-year-old 
            {scenario['patient_info']['gender'].lower()} with {scenario['patient_info']['presentation']}.
            
            This case demonstrates our AI-powered approach to precision medicine for RUNX1 
            Familial Platelet Disorder - a rare genetic condition that predisposes patients to leukemia.
            """,
            
            "act_1_patient_presentation": f"""
            🎬 **Act 1: The Patient** (5 minutes)
            
            **Patient Background:**
            - Name: {scenario['patient_info']['name']}
            - Age: {scenario['patient_info']['age']}
            - Presentation: {scenario['patient_info']['presentation']}
            - Family History: {scenario['patient_info']['family_history']}
            
            **Clinical Challenge:**
            RUNX1-FPD patients appear healthy but carry a "ticking time bomb" - a germline mutation 
            that predisposes them to leukemia. The key question is: WHEN will this happen, and 
            CAN we prevent it?
            
            **Demo Action:** Show patient timeline and family pedigree
            """,
            
            "act_2_digital_twin": f"""
            🎬 **Act 2: The Digital Twin** (10 minutes)
            
            **AI-Powered Analysis:**
            Our platform creates a computational "Digital Twin" of the patient, modeling:
            1. Germline mutation impact: {scenario['germline_variant']['protein_change']}
            2. Somatic evolution: {scenario['somatic_variant']['protein_change'] if scenario['somatic_variant'] else 'None detected'}
            3. Risk stratification: {scenario['expected_analysis']['risk_tier']}
            4. Intervention opportunities
            
            **Demo Action:** 
            - Run live AI analysis
            - Show two-hit hypothesis modeling
            - Display risk assessment results
            
            **Expected Results:**
            - Risk Tier: {scenario['expected_analysis']['risk_tier']}
            - Transformation Probability: {scenario['expected_analysis']['transformation_probability']:.1%}
            - Surveillance: {scenario['expected_analysis']['surveillance_frequency']}
            """,
            
            "act_3_intervention": f"""
            🎬 **Act 3: The Intervention** (10 minutes)
            
            **Precision Therapeutic Design:**
            Based on the AI analysis, our platform designs personalized interventions:
            
            1. **Target Identification:** {scenario['somatic_variant']['gene'] if scenario['somatic_variant'] else 'Preventive approach'}
            2. **Guide RNA Design:** AI-optimized CRISPR guides
            3. **Safety Assessment:** Comprehensive off-target analysis
            4. **Protocol Generation:** Step-by-step experimental procedures
            
            **Demo Action:**
            - Design precision CRISPR intervention
            - Show guide RNA optimization
            - Display safety assessment
            - Generate experimental timeline
            
            **Clinical Impact:**
            This represents a paradigm shift from reactive cancer treatment to 
            proactive cancer interception.
            """,
            
            "act_4_leap_impact": f"""
            🎬 **Act 4: The Impact** (5 minutes)
            
            **LEAP Grant Alignment:**
            This demonstration showcases breakthrough capabilities in:
            
            **Focus Area 1 - Mechanistic Understanding:**
            ✅ Cell-autonomous mechanisms quantified
            ✅ Temporal progression modeled
            ✅ Intervention opportunities identified
            
            **Focus Area 2 - Therapeutic Development:**
            ✅ Precision targets identified
            ✅ AI-powered drug design
            ✅ Clinical trial ready protocols
            
            **Commercial Impact:**
            - Market: 2,000+ RUNX1-FPD families worldwide
            - Revenue: $50K-200K research licenses, $100K-500K clinical packages
            - Scalability: Framework applies to all germline cancer syndromes
            
            **Next Steps:**
            1. Clinical validation studies
            2. Regulatory pathway planning
            3. Commercial partnerships
            4. Platform expansion to other syndromes
            """
        }
    
    def _create_expected_results(self, scenario: Dict) -> Dict:
        """Create expected results for the scenario."""
        return {
            "germline_analysis": {
                "functional_impact": scenario['germline_variant']['functional_impact'],
                "pathogenicity": scenario['germline_variant']['pathogenicity'],
                "domain_affected": "Runt Domain" if "Arg" in scenario['germline_variant']['protein_change'] else "TAD Domain"
            },
            "somatic_analysis": {
                "clonal_evolution": "detected" if scenario['somatic_variant'] else "not_detected",
                "transformation_risk": scenario['expected_analysis']['transformation_probability'],
                "intervention_window": "immediate" if scenario['expected_analysis']['transformation_probability'] > 0.5 else "planned"
            },
            "intervention_design": {
                "guides_found": random.randint(3, 8),
                "success_probability": random.uniform(0.1, 0.4),
                "timeline_weeks": random.uniform(4, 8),
                "safety_score": random.uniform(0.7, 0.9)
            },
            "clinical_recommendations": {
                "surveillance_frequency": scenario['expected_analysis']['surveillance_frequency'],
                "intervention_urgency": scenario['expected_analysis']['intervention_urgency'],
                "family_screening": "recommended" if scenario['patient_info']['family_history'] != "none" else "optional"
            }
        }
    
    def _create_talking_points(self, scenario: Dict) -> Dict:
        """Create talking points for the demo."""
        return {
            "key_messages": [
                "AI-powered precision medicine for rare genetic diseases",
                "Proactive cancer interception vs reactive treatment",
                "Multi-AI integration for comprehensive analysis",
                "Clinical-grade results with physician-ready reports",
                "Scalable platform for all germline cancer syndromes"
            ],
            "technical_highlights": [
                "Real-time AI analysis with multiple models",
                "Comprehensive two-hit hypothesis modeling",
                "Professional-grade CRISPR guide design",
                "Evidence-based risk stratification",
                "Automated experimental protocol generation"
            ],
            "business_value": [
                "First-to-market in RUNX1-FPD precision medicine",
                "Clear revenue model with multiple customer segments",
                "Regulatory pathway well-defined",
                "Strong IP position with novel AI integration",
                "Scalable to $2-5B precision medicine market"
            ],
            "leap_grant_alignment": [
                "Directly addresses Focus Area 1 mechanistic understanding",
                "Provides foundation for Focus Area 2 therapeutic development",
                "Demonstrates AI/ML innovation in cancer research",
                "Shows clear path to clinical translation",
                "Establishes platform for future cancer interception research"
            ]
        }
    
    def _create_technical_details(self, scenario: Dict) -> Dict:
        """Create technical details for the demo."""
        return {
            "ai_models_used": [
                "ZetaOracle: Variant pathogenicity prediction",
                "Evo2: Sequence generation and optimization",
                "ChopChop: Guide RNA scoring and safety analysis",
                "Custom: Two-hit progression modeling",
                "Custom: Risk stratification algorithms"
            ],
            "data_sources": [
                "ClinVar: Variant pathogenicity database",
                "gnomAD: Population frequency data",
                "COSMIC: Somatic mutation database",
                "RUNX1-FPD: Patient registry data",
                "Literature: Curated research findings"
            ],
            "computational_pipeline": [
                "1. Variant annotation and functional impact assessment",
                "2. Two-hit hypothesis modeling with clonal evolution",
                "3. Risk stratification using machine learning",
                "4. Intervention design with AI-optimized guides",
                "5. Clinical report generation with recommendations"
            ],
            "validation_methods": [
                "Cross-validation with known pathogenic variants",
                "Comparison with clinical outcomes data",
                "Expert review by RUNX1-FPD specialists",
                "Benchmarking against existing tools",
                "Prospective validation in clinical cohorts"
            ]
        }

# Streamlit integration functions
def display_demo_scenario_selector():
    """Display demo scenario selector in Streamlit."""
    st.subheader("🎬 Demo Scenario Selection")
    
    demo_scenarios = RUNX1DemoScenarios()
    scenarios = demo_scenarios.get_all_scenarios()
    
    # Create scenario selection
    scenario_options = {
        scenario['title']: scenario['id'] for scenario in scenarios
    }
    
    selected_title = st.selectbox(
        "Choose Demo Scenario",
        options=list(scenario_options.keys()),
        help="Select a patient scenario for demonstration"
    )
    
    selected_id = scenario_options[selected_title]
    selected_scenario = demo_scenarios.get_scenario(selected_id)
    
    if selected_scenario:
        # Display scenario overview
        st.markdown(f"### {selected_scenario['title']}")
        st.markdown(f"**{selected_scenario['subtitle']}**")
        
        # Patient information
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("#### 👤 Patient Information")
            patient_info = selected_scenario['patient_info']
            st.markdown(f"**Name:** {patient_info['name']}")
            st.markdown(f"**Age:** {patient_info['age']}")
            st.markdown(f"**Gender:** {patient_info['gender']}")
            st.markdown(f"**Ethnicity:** {patient_info['ethnicity']}")
            st.markdown(f"**Presentation:** {patient_info['presentation']}")
            st.markdown(f"**Family History:** {patient_info['family_history']}")
        
        with col2:
            st.markdown("#### 🧬 Genetic Profile")
            germline = selected_scenario['germline_variant']
            st.markdown(f"**Germline Mutation:** {germline['gene']} {germline['protein_change']}")
            st.markdown(f"**Impact:** {germline['functional_impact']}")
            st.markdown(f"**Pathogenicity:** {germline['pathogenicity']}")
            
            if selected_scenario['somatic_variant']:
                somatic = selected_scenario['somatic_variant']
                st.markdown(f"**Somatic Mutation:** {somatic['gene']} {somatic['protein_change']}")
                st.markdown(f"**Impact:** {somatic['functional_impact']}")
        
        # Clinical timeline
        st.markdown("#### 📅 Clinical Timeline")
        timeline_df = pd.DataFrame(selected_scenario['clinical_timeline'])
        st.dataframe(timeline_df, use_container_width=True)
        
        # Expected analysis results
        st.markdown("#### 📊 Expected Analysis Results")
        expected = selected_scenario['expected_analysis']
        
        col3, col4, col5 = st.columns(3)
        
        with col3:
            st.metric("Risk Tier", expected['risk_tier'].title())
        
        with col4:
            st.metric("Transformation Risk", f"{expected['transformation_probability']:.1%}")
        
        with col5:
            st.metric("Surveillance", expected['surveillance_frequency'])
        
        # Demo flow
        if st.button("📋 Generate Demo Flow", type="primary"):
            demo_flow = demo_scenarios.create_demo_flow(selected_id)
            st.session_state.demo_flow = demo_flow
            st.session_state.selected_scenario = selected_scenario
            st.success("✅ Demo flow generated! Ready for presentation.")
        
        # Store selected scenario
        st.session_state.current_scenario = selected_scenario
        
        return selected_scenario
    
    return None

def display_demo_script(demo_flow: Dict):
    """Display the demo script."""
    if not demo_flow:
        return
    
    st.subheader("🎬 Demo Script")
    
    script = demo_flow['demo_script']
    
    # Introduction
    st.markdown("### 🎯 Introduction")
    st.markdown(script['introduction'])
    
    # Act 1
    st.markdown("### 🎬 Act 1: The Patient")
    st.markdown(script['act_1_patient_presentation'])
    
    # Act 2
    st.markdown("### 🎬 Act 2: The Digital Twin")
    st.markdown(script['act_2_digital_twin'])
    
    # Act 3
    st.markdown("### 🎬 Act 3: The Intervention")
    st.markdown(script['act_3_intervention'])
    
    # Act 4
    st.markdown("### 🎬 Act 4: The Impact")
    st.markdown(script['act_4_leap_impact'])

def display_talking_points(demo_flow: Dict):
    """Display talking points for the demo."""
    if not demo_flow:
        return
    
    st.subheader("💬 Talking Points")
    
    talking_points = demo_flow['talking_points']
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("#### 🎯 Key Messages")
        for point in talking_points['key_messages']:
            st.markdown(f"- {point}")
        
        st.markdown("#### 🔬 Technical Highlights")
        for point in talking_points['technical_highlights']:
            st.markdown(f"- {point}")
    
    with col2:
        st.markdown("#### 💰 Business Value")
        for point in talking_points['business_value']:
            st.markdown(f"- {point}")
        
        st.markdown("#### 🏆 LEAP Grant Alignment")
        for point in talking_points['leap_grant_alignment']:
            st.markdown(f"- {point}")

if __name__ == "__main__":
    # Test the demo scenarios
    st.title("🎬 RUNX1-FPD Demo Scenarios")
    
    scenario = display_demo_scenario_selector()
    
    if 'demo_flow' in st.session_state:
        demo_flow = st.session_state.demo_flow
        
        tab1, tab2, tab3 = st.tabs(["Demo Script", "Talking Points", "Technical Details"])
        
        with tab1:
            display_demo_script(demo_flow)
        
        with tab2:
            display_talking_points(demo_flow)
        
        with tab3:
            st.markdown("### 🔧 Technical Details")
            st.json(demo_flow['technical_details']) 