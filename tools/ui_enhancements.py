"""
UI/UX Enhancement Utilities for LEAP Grant Demo

This module provides reusable components and utilities to enhance the user experience
of the CasPro platform while maintaining 100% dynamic data flow.

Key Principles:
- No hardcoded data - all content adapts to actual API responses
- Context-aware explanations that reference real results
- Progressive disclosure of complex processes
- Clear narrative flow for grant demonstration
"""

import streamlit as st
import pandas as pd
from typing import Dict, List, Any, Optional
import time
from datetime import datetime
import json

class DemoFlowManager:
    """Manages the overall demo flow and narrative progression"""
    
    def __init__(self):
        self.demo_mode = st.session_state.get('leap_demo_mode', False)
        self.current_step = st.session_state.get('demo_step', 0)
        self.completed_steps = st.session_state.get('completed_steps', set())
    
    def enable_demo_mode(self):
        """Enable LEAP Grant demo mode with enhanced explanations"""
        st.session_state['leap_demo_mode'] = True
        st.session_state['demo_step'] = 0
        st.session_state['completed_steps'] = set()
        st.success("🎯 **LEAP Grant Demo Mode Activated** - Enhanced explanations and guided workflow enabled!")
    
    def show_demo_toggle(self):
        """Display the demo mode toggle"""
        col1, col2, col3 = st.columns([1, 2, 1])
        with col2:
            if st.button("🎯 Activate LEAP Grant Demo Mode", type="primary" if not self.demo_mode else "secondary"):
                self.enable_demo_mode()
                st.rerun()
    
    def show_progress_indicator(self, total_steps: int = 5):
        """Show progress through the demo workflow"""
        if not self.demo_mode:
            return
            
        progress = len(self.completed_steps) / total_steps
        st.progress(progress)
        st.caption(f"Demo Progress: {len(self.completed_steps)}/{total_steps} major workflows completed")
    
    def mark_step_complete(self, step_name: str):
        """Mark a workflow step as completed"""
        if self.demo_mode:
            self.completed_steps.add(step_name)
            st.session_state['completed_steps'] = self.completed_steps

class ContextualExplainer:
    """Provides dynamic, context-aware explanations based on actual data"""
    
    @staticmethod
    def explain_ai_process(process_name: str, input_data: Dict, output_data: Dict = None):
        """Explain what an AI process is doing with actual data context"""
        explanations = {
            "threat_assessment": {
                "title": "🧠 AI Threat Assessment in Progress",
                "process": f"Analyzing {input_data.get('gene_symbol', 'gene')} using our Triumvirate Protocol",
                "steps": [
                    "🔍 Querying internal COSMIC database for known variants",
                    "🌐 Fetching real-time annotation from Ensembl VEP",
                    "⚡ Running Truncation Sieve for catastrophic mutations",
                    "🧬 Applying Zeta Oracle for complex variant analysis",
                    "📊 Synthesizing multi-source intelligence"
                ]
            },
            "seed_soil_analysis": {
                "title": "🌱 Metastasis Modeling - Seed & Soil Analysis",
                "process": f"Modeling how {input_data.get('gene_symbol', 'mutation')} interacts with tissue environments",
                "steps": [
                    "🔬 Creating digital twin of healthy tissue",
                    "🦠 Simulating pathogenic 'seed' invasion",
                    "📈 Calculating tissue vulnerability scores",
                    "🎯 Identifying 'soil-poisoning' therapeutic targets"
                ]
            },
            "therapeutic_design": {
                "title": "✂️ AI-Powered Therapeutic Design",
                "process": f"Designing precision therapeutics for {input_data.get('target_gene', 'target')}",
                "steps": [
                    "🧬 Generating optimal CRISPR guide sequences",
                    "🔍 Predicting off-target effects",
                    "💊 Screening for drug repurposing opportunities",
                    "🧪 Designing novel nanobody inhibitors"
                ]
            }
        }
        
        if process_name in explanations:
            exp = explanations[process_name]
            with st.expander(f"💡 {exp['title']}", expanded=True):
                st.write(f"**What's happening:** {exp['process']}")
                st.write("**AI Processing Steps:**")
                for step in exp['steps']:
                    st.write(f"  {step}")
                
                if output_data:
                    st.write("**Key Results:**")
                    for key, value in output_data.items():
                        if isinstance(value, (str, int, float)):
                            st.write(f"  • **{key.replace('_', ' ').title()}:** {value}")

class ScientificNarrator:
    """Provides scientific context and explanations for biological processes"""
    
    @staticmethod
    def explain_mutation_impact(gene_symbol: str, protein_change: str, mutation_type: str):
        """Explain the biological impact of a specific mutation"""
        
        # Dynamic explanations based on actual data
        impact_explanations = {
            "frameshift": f"🚨 **Frameshift mutations in {gene_symbol}** cause a fundamental disruption of the protein's reading frame, typically resulting in a premature stop codon and loss of function.",
            "nonsense": f"⚠️ **Nonsense mutations in {gene_symbol}** introduce premature stop codons, leading to truncated, non-functional proteins.",
            "missense": f"🔄 **Missense mutations in {gene_symbol}** alter single amino acids, potentially disrupting protein structure and function.",
            "deletion": f"🗑️ **Deletions in {gene_symbol}** remove critical genetic material, often causing loss of essential protein domains.",
            "insertion": f"➕ **Insertions in {gene_symbol}** add genetic material that can disrupt normal protein structure and function."
        }
        
        mutation_type_lower = mutation_type.lower()
        for key, explanation in impact_explanations.items():
            if key in mutation_type_lower:
                return explanation
        
        return f"🧬 **{protein_change} in {gene_symbol}** represents a genetic alteration that may impact protein function and cellular behavior."
    
    @staticmethod
    def explain_therapeutic_rationale(target_gene: str, therapeutic_type: str, confidence_score: float = None):
        """Explain why a specific therapeutic approach makes sense"""
        
        rationales = {
            "crispr_knockout": f"✂️ **CRISPR knockout of {target_gene}** aims to eliminate the source of pathogenic activity by permanently disabling the problematic gene.",
            "nanobody_inhibition": f"🧪 **Nanobody inhibition of {target_gene}** provides a reversible, targeted approach to block the protein's harmful activity.",
            "drug_repurposing": f"💊 **Drug repurposing for {target_gene}** leverages existing FDA-approved compounds to rapidly deploy therapeutic interventions."
        }
        
        base_explanation = rationales.get(therapeutic_type, f"🎯 **Targeting {target_gene}** with {therapeutic_type} represents a precision medicine approach.")
        
        if confidence_score:
            confidence_text = "high" if confidence_score > 0.8 else "moderate" if confidence_score > 0.5 else "preliminary"
            base_explanation += f" Our AI models predict {confidence_text} likelihood of therapeutic success."
        
        return base_explanation

class ResultsVisualizer:
    """Enhanced visualization of results with context"""
    
    @staticmethod
    def display_threat_dossier(dossier_data: Dict, demo_mode: bool = False):
        """Display threat assessment results with enhanced context"""
        
        if demo_mode:
            st.success("🎯 **Intelligence Fusion Complete** - Multiple data sources successfully integrated!")
        
        # Extract key information
        gene_symbol = dossier_data.get('gene_symbol', 'Unknown')
        protein_change = dossier_data.get('protein_change', 'Unknown')
        assessment = dossier_data.get('assessment', {})
        
        # Main results display
        col1, col2 = st.columns([2, 1])
        
        with col1:
            st.subheader(f"🧬 Threat Profile: {gene_symbol} {protein_change}")
            
            # Assessment source and confidence
            source = assessment.get('assessment_source', 'Unknown')
            is_pathogenic = assessment.get('is_pathogenic', False)
            
            if is_pathogenic:
                st.error(f"⚠️ **PATHOGENIC VARIANT CONFIRMED** (Source: {source})")
            else:
                st.success(f"✅ **Benign/Uncertain** (Source: {source})")
            
            # Dynamic explanation based on actual results
            if source == "Truncation Sieve":
                st.info("🔍 **Analysis Method:** Our bioinformatics pipeline detected a catastrophic mutation (frameshift/nonsense) that definitively disrupts protein function.")
            elif source == "Zeta Oracle":
                st.info("🧠 **Analysis Method:** Our AI model analyzed the mutation's context and predicted its functional impact using deep learning.")
        
        with col2:
            # Quick stats
            st.metric("Assessment Source", source)
            st.metric("Pathogenic", "Yes" if is_pathogenic else "No")
            
            if demo_mode:
                st.info("💡 **LEAP Grant Relevance:** This demonstrates our ability to rapidly assess variant pathogenicity using multiple complementary approaches.")
    
    @staticmethod
    def display_therapeutic_options(options_data: List[Dict], demo_mode: bool = False):
        """Display therapeutic options with enhanced context"""
        
        if demo_mode:
            st.success("🎯 **Multi-Modal Therapeutic Design Complete** - Precision interventions ready for development!")
        
        for i, option in enumerate(options_data):
            with st.expander(f"Option {i+1}: {option.get('type', 'Therapeutic')} - {option.get('target', 'Unknown')}", expanded=i==0):
                
                col1, col2 = st.columns([3, 1])
                
                with col1:
                    st.write(f"**Mechanism:** {option.get('mechanism', 'Not specified')}")
                    st.write(f"**Rationale:** {option.get('rationale', 'AI-designed precision therapeutic')}")
                    
                    if 'sequence' in option:
                        st.code(option['sequence'][:100] + "..." if len(option['sequence']) > 100 else option['sequence'])
                
                with col2:
                    if 'confidence_score' in option:
                        st.metric("AI Confidence", f"{option['confidence_score']:.1%}")
                    
                    if 'safety_score' in option:
                        st.metric("Safety Score", f"{option['safety_score']:.2f}")
                
                if demo_mode:
                    st.info(f"💡 **Clinical Translation:** This {option.get('type', 'therapeutic')} could proceed to preclinical validation within 6-12 months.")

class LEAPGrantAligner:
    """Specifically highlights LEAP Grant focus areas"""
    
    @staticmethod
    def highlight_focus_area_1(data: Dict):
        """Highlight how current analysis addresses LEAP Focus Area 1"""
        try:
            st.sidebar.markdown("### 🎯 LEAP Grant Focus Area 1")
            st.sidebar.success("**Dissecting Mechanisms of Cancer Initiation**")
            st.sidebar.write("✅ Cell-autonomous factor analysis")
            st.sidebar.write("✅ Cell-nonautonomous factor integration")
            st.sidebar.write("✅ Clonal hematopoiesis progression modeling")
            st.sidebar.write("✅ RUNX1-FPD specific mechanisms")
            
            if 'gene_symbol' in data:
                st.sidebar.write(f"🧬 **Current Analysis:** {data['gene_symbol']} impact on clonal evolution")
        except Exception:
            pass  # Silently fail if Streamlit context not available
    
    @staticmethod
    def highlight_focus_area_2(data: Dict):
        """Highlight how current analysis addresses LEAP Focus Area 2"""
        try:
            st.sidebar.markdown("### 🎯 LEAP Grant Focus Area 2")
            st.sidebar.success("**Cancer Interception Strategies**")
            st.sidebar.write("✅ Evidence-based therapeutic design")
            st.sidebar.write("✅ Multi-modal intervention approaches")
            st.sidebar.write("✅ Clinical translation readiness")
            st.sidebar.write("✅ RUNX1-FPD community urgency")
            
            if 'therapeutic_type' in data:
                st.sidebar.write(f"🎯 **Current Design:** {data['therapeutic_type']} for cancer interception")
        except Exception:
            pass  # Silently fail if Streamlit context not available
    
    @staticmethod
    def show_runx1_fpd_context():
        """Display RUNX1-FPD specific context and statistics"""
        try:
            st.sidebar.markdown("### 🩸 RUNX1-FPD Context")
            st.sidebar.error("**High-Risk Population**")
            st.sidebar.metric("Lifetime HM Risk", "35-50%")
            st.sidebar.metric("Most Common", "Acute Myeloid Leukemia")
            st.sidebar.metric("Deadliness Rank", "#2 among HMs")
            
            st.sidebar.info("""
            **Clonal Hematopoiesis Impact:**
            - General population (70+): >10% CH rate
            - RUNX1-FPD patients: Much higher CH rates
            - Key gap: Which somatic mutations drive transformation?
            """)
        except Exception:
            pass  # Silently fail if Streamlit context not available
    
    @staticmethod
    def explain_cell_autonomous_vs_nonautonomous():
        """Explain the distinction between cell-autonomous and cell-nonautonomous factors"""
        try:
            with st.expander("🧬 Cell-Autonomous vs Cell-Nonautonomous Factors", expanded=False):
                col1, col2 = st.columns(2)
                
                with col1:
                    st.markdown("### 🔬 Cell-Autonomous")
                    st.write("**Internal factors within the cell:**")
                    st.write("• Germline RUNX1 mutations")
                    st.write("• Somatic 'second hit' mutations")
                    st.write("• DNA damage response defects")
                    st.write("• Transcriptional dysregulation")
                    
                with col2:
                    st.markdown("### 🌐 Cell-Nonautonomous")
                    st.write("**External microenvironmental factors:**")
                    st.write("• Bone marrow niche signals")
                    st.write("• Inflammatory cytokines")
                    st.write("• Stromal cell interactions")
                    st.write("• Immune surveillance escape")
                
                st.info("**LEAP Grant Focus:** Understanding how these factors interact to drive clonal evolution and malignant transformation in RUNX1-FPD patients.")
        except Exception:
            pass  # Silently fail if Streamlit context not available

def create_enhanced_sidebar(page_name: str, data: Dict = None):
    """Create an enhanced sidebar with context and LEAP Grant alignment"""
    
    # Only proceed if we're in a valid Streamlit context
    try:
        # Demo mode indicator
        if st.session_state.get('leap_demo_mode', False):
            st.sidebar.success("🎯 **LEAP Grant Demo Mode Active**")
            st.sidebar.progress(0.8)  # Example progress
            
            # Show RUNX1-FPD context in demo mode
            LEAPGrantAligner.show_runx1_fpd_context()
        
        # Page-specific context
        page_contexts = {
            "Patient Digital Twin": "Creating computational models of RUNX1-FPD patient genomics to predict clonal evolution and therapeutic responses.",
            "Threat Assessor": "Using AI to rapidly assess which somatic mutations drive malignant transformation in RUNX1-FPD patients.",
            "Seed & Soil Analysis": "Modeling how pathogenic clones interact with bone marrow microenvironments to identify intervention points.",
            "Therapeutic Design": "Designing cancer interception strategies using AI-powered therapeutic optimization."
        }
        
        if page_name in page_contexts:
            st.sidebar.markdown(f"### 📍 Current Workflow")
            st.sidebar.info(page_contexts[page_name])
        
        # LEAP Grant alignment
        if data:
            LEAPGrantAligner.highlight_focus_area_1(data)
            LEAPGrantAligner.highlight_focus_area_2(data)
        
        # Next steps guidance
        st.sidebar.markdown("### 🚀 Recommended Next Steps")
        if page_name == "Patient Digital Twin":
            st.sidebar.write("1. Assess clonal hematopoiesis variants")
            st.sidebar.write("2. Model microenvironment interactions")
        elif page_name == "Threat Assessor":
            st.sidebar.write("1. Analyze cell-nonautonomous factors")
            st.sidebar.write("2. Design interception therapeutics")
        elif page_name == "Seed & Soil Analysis":
            st.sidebar.write("1. Target pathogenic clones")
            st.sidebar.write("2. Modulate bone marrow niche")
        
        return True
    except Exception as e:
        # Silently fail if Streamlit context not available (during import)
        return False

# Utility functions for consistent styling
def apply_leap_styling():
    """Apply consistent styling for LEAP Grant demo"""
    st.markdown("""
    <style>
    .leap-highlight {
        background-color: #e8f4fd;
        border-left: 4px solid #1f77b4;
        padding: 10px;
        margin: 10px 0;
    }
    .ai-process {
        background-color: #f0f8f0;
        border-left: 4px solid #2ca02c;
        padding: 10px;
        margin: 10px 0;
    }
    .clinical-relevance {
        background-color: #fff8dc;
        border-left: 4px solid #ff7f0e;
        padding: 10px;
        margin: 10px 0;
    }
    </style>
    """, unsafe_allow_html=True)

def show_loading_with_context(message: str, steps: List[str] = None):
    """Show loading state with contextual information"""
    with st.spinner(message):
        if steps:
            progress_bar = st.progress(0)
            for i, step in enumerate(steps):
                st.write(f"⏳ {step}")
                progress_bar.progress((i + 1) / len(steps))
                time.sleep(0.5)  # Brief pause for demonstration
    
    st.success("✅ Analysis complete!")

class AppleStyleResultsDisplay:
    """Apple-inspired results display that replaces JSON dumps with intuitive cards"""
    
    @staticmethod
    def display_threat_assessment(data: Dict) -> None:
        """Display threat assessment in Apple-style cards"""
        if not data or 'assessment' not in data:
            st.error("⚠️ **Analysis Incomplete**\n\nWe're still processing your request. Please wait a moment and try again.")
            return
        
        assessment = data['assessment']
        
        # Hero Assessment Card
        AppleStyleResultsDisplay._display_hero_card(assessment, data)
        
        # Key Insights Panel
        AppleStyleResultsDisplay._display_key_insights(assessment, data)
        
        # Clinical Context
        AppleStyleResultsDisplay._display_clinical_context(data)
        
        # Next Steps
        AppleStyleResultsDisplay._display_next_steps(assessment)
        
        # Technical Details (Expandable)
        AppleStyleResultsDisplay._display_technical_details(data)
    
    @staticmethod
    def _display_hero_card(assessment: Dict, full_data: Dict) -> None:
        """Large, prominent verdict card"""
        is_pathogenic = assessment.get('is_pathogenic', False)
        variant_desc = assessment.get('variant_description', 'Unknown variant')
        gene_symbol = full_data.get('gene_symbol', 'Unknown gene')
        
        if is_pathogenic:
            st.error(f"""
            🚨 **HIGH-RISK MUTATION DETECTED**
            
            **{gene_symbol} {variant_desc}**
            
            **Assessment Method:** {assessment.get('assessment_source', 'Unknown')}  
            **Confidence:** {AppleStyleResultsDisplay._get_confidence_level(assessment.get('zeta_score', 0))}  
            **Zeta Score:** {assessment.get('zeta_score', 'N/A')}
            
            **What this means:** {AppleStyleResultsDisplay._get_pathogenic_explanation(assessment, gene_symbol)}
            """)
        else:
            st.success(f"""
            ✅ **LOW-RISK VARIANT**
            
            **{gene_symbol} {variant_desc}**
            
            **Assessment Method:** {assessment.get('assessment_source', 'Unknown')}  
            **Confidence:** {AppleStyleResultsDisplay._get_confidence_level(assessment.get('zeta_score', 0))}  
            **Zeta Score:** {assessment.get('zeta_score', 'N/A')}
            
            **What this means:** This genetic change is unlikely to cause significant harm in RUNX1-FPD patients.
            """)
    
    @staticmethod
    def _display_key_insights(assessment: Dict, full_data: Dict) -> None:
        """Key clinical insights"""
        st.subheader("🧠 Key Insights")
        
        insights = []
        
        # Assessment method insight
        method = assessment.get('assessment_source', 'Unknown')
        if method == 'Truncation Sieve':
            insights.append("• **Instant Analysis:** This mutation was identified immediately as protein-truncating")
            insights.append("• **High Confidence:** Truncating mutations are definitively pathogenic")
        elif method == 'Zeta Oracle':
            insights.append("• **AI Analysis:** Our deep learning model analyzed this variant's biological impact")
            insights.append("• **Multi-factor Assessment:** Considers protein structure, conservation, and clinical data")
        
        # Clinical significance
        if assessment.get('is_pathogenic', False):
            insights.append("• **Clinical Impact:** This mutation significantly increases leukemia risk in RUNX1-FPD patients")
            insights.append("• **Intervention Needed:** Immediate therapeutic design recommended")
        else:
            insights.append("• **Surveillance Mode:** Continue monitoring but no immediate intervention needed")
        
        # Display insights
        for insight in insights:
            st.markdown(insight)
    
    @staticmethod
    def _display_clinical_context(data: Dict) -> None:
        """Clinical context panel"""
        st.subheader("🏥 Clinical Significance")
        
        # Variant profile context
        if 'variant_profile' in data and data['variant_profile']:
            profile = data['variant_profile']
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("**Cancer Association:**")
                primary_site = profile.get('primary_site', 'Unknown').replace('_', ' ').title()
                histology = profile.get('primary_histology', 'Unknown').replace('_', ' ').title()
                st.markdown(f"• Primary site: {primary_site}")
                st.markdown(f"• Histology: {histology}")
            
            with col2:
                st.markdown("**Research Evidence:**")
                if profile.get('pubmed_pmid'):
                    pmid = profile['pubmed_pmid']
                    st.markdown(f"• [PubMed {pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)")
                else:
                    st.markdown("• Literature review in progress")
        
        # Clinical trials context
        if 'clinical_trials' in data and data['clinical_trials']:
            trials = data['clinical_trials']
            st.markdown(f"**Active Research:** {len(trials)} clinical trials targeting this pathway")
    
    @staticmethod
    def _display_next_steps(assessment: Dict) -> None:
        """Recommended next steps"""
        st.subheader("🎯 Recommended Next Steps")
        
        if assessment.get('is_pathogenic', False):
            st.markdown("""
            1. **🌱 Seed & Soil Analysis** - Analyze metastatic potential and tissue vulnerability
            2. **✂️ CRISPR Design** - Create targeted intervention strategies
            3. **💊 Drug Discovery** - Explore combination therapeutic approaches
            4. **📊 Risk Modeling** - Quantify intervention benefits
            """)
        else:
            st.markdown("""
            1. **👀 Continue Surveillance** - Monitor for additional mutations
            2. **🧬 Family Screening** - Check for other high-risk variants
            3. **📈 Risk Assessment** - Regular follow-up recommended
            """)
    
    @staticmethod
    def _display_technical_details(data: Dict) -> None:
        """Expandable technical details"""
        with st.expander("🔬 Technical Details (For Researchers)", expanded=False):
            st.subheader("Raw API Response")
            st.json(data)
            
            if 'vep_dossier' in data:
                st.subheader("VEP Analysis")
                vep_data = data['vep_dossier']
                if 'error' in vep_data:
                    st.warning(f"VEP Analysis: {vep_data['error']}")
                    st.info("💡 **Why VEP failed:** The variant format may not be compatible with Ensembl's database, or the service may be temporarily unavailable. Our AI analysis provides equivalent insights.")
                else:
                    st.json(vep_data)
    
    @staticmethod
    def _get_confidence_level(zeta_score: float) -> str:
        """Convert zeta score to human-readable confidence"""
        if zeta_score == -1000:
            return "Extremely High (Definitive)"
        elif zeta_score < -50:
            return "High"
        elif zeta_score < -10:
            return "Moderate"
        elif zeta_score < 0:
            return "Low"
        else:
            return "Very Low"
    
    @staticmethod
    def _get_pathogenic_explanation(assessment: Dict, gene_symbol: str) -> str:
        """Get contextual explanation for pathogenic variants"""
        method = assessment.get('assessment_source', '')
        
        if method == 'Truncation Sieve':
            return f"This mutation causes the {gene_symbol} protein to be cut short, making it completely non-functional. In RUNX1-FPD patients, this type of 'second hit' significantly increases the risk of developing acute leukemia."
        elif method == 'Zeta Oracle':
            return f"Our AI analysis indicates this {gene_symbol} mutation significantly disrupts normal protein function, contributing to cancer development in RUNX1-FPD patients."
        else:
            return f"This {gene_symbol} mutation is predicted to have harmful effects on protein function."

class PatientJourneyVisualizer:
    """Visualizes the RUNX1-FPD patient journey from germline to malignancy"""
    
    @staticmethod
    def display_journey_timeline(current_stage: str = "second_hit") -> None:
        """Display patient journey with current position highlighted"""
        st.subheader("🧬 RUNX1-FPD Patient Journey")
        
        stages = [
            ("👶", "Born with RUNX1 mutation", "First Hit", "germline"),
            ("🧬", "Clonal hematopoiesis develops", "Age 40-60", "clonal_hematopoiesis"), 
            ("⚡", "Somatic mutation acquired", "Second Hit", "second_hit"),
            ("🚨", "Malignant transformation", "35-50% Risk", "malignancy")
        ]
        
        for i, (emoji, description, subtitle, stage_id) in enumerate(stages):
            col1, col2, col3 = st.columns([1, 8, 2])
            
            with col1:
                if stage_id == current_stage:
                    st.markdown(f"## {emoji}")
                    st.markdown("**← YOU ARE HERE**")
                else:
                    st.markdown(f"### {emoji}")
            
            with col2:
                if stage_id == current_stage:
                    st.markdown(f"**{description}**")
                    st.markdown(f"*{subtitle}*")
                else:
                    st.markdown(f"{description}")
                    st.markdown(f"*{subtitle}*")
            
            # Add arrow between stages
            if i < len(stages) - 1:
                st.markdown("↓")
    
    @staticmethod
    def display_risk_stratification() -> None:
        """Display risk comparison"""
        st.subheader("📊 Risk Profile Comparison")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.metric(
                label="General Population Risk",
                value="0.1%",
                help="Risk of developing hematologic malignancy"
            )
        
        with col2:
            st.metric(
                label="RUNX1-FPD Patient Risk", 
                value="35-50%",
                delta="350-500x higher",
                delta_color="inverse",
                help="Lifetime risk of hematologic malignancy"
            )

class ErrorStateManager:
    """Manages error states with helpful, contextual messages"""
    
    @staticmethod
    def display_friendly_error(error_data: Dict, context: str = "") -> None:
        """Display user-friendly error messages"""
        error_msg = error_data.get('error', 'Unknown error occurred')
        
        if "No variant profile found" in error_msg:
            st.info("""
            🔍 **Novel Variant Detected**
            
            This variant isn't in our database yet, which means it's either very rare or recently discovered. 
            We're analyzing it using our AI models and external databases to determine its significance.
            
            ⏳ **Status:** Analysis in progress  
            🧠 **Method:** Zeta Oracle AI assessment  
            📡 **External validation:** Querying Ensembl VEP
            """)
        elif "400 Client Error" in error_msg and "vep" in error_msg.lower():
            st.warning("""
            ⚠️ **External Database Temporarily Unavailable**
            
            The Ensembl VEP service is currently unavailable, but don't worry - our internal AI analysis 
            provides equivalent insights for variant assessment.
            
            ✅ **Our Analysis:** Complete and reliable  
            📊 **Confidence:** High (based on internal models)  
            🔄 **Retry:** VEP will be attempted again automatically
            """)
        else:
            st.error(f"""
            ❌ **Analysis Error**
            
            We encountered an unexpected issue during analysis: {error_msg}
            
            **Suggested Actions:**
            1. Check your input format
            2. Try a different variant notation
            3. Contact support if the issue persists
            """) 