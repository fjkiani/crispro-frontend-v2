import streamlit as st
import requests
import os
import json
from tools.ui_enhancements import (
    DemoFlowManager, 
    ContextualExplainer, 
    ScientificNarrator, 
    LEAPGrantAligner,
    AppleStyleResultsDisplay,
    PatientJourneyVisualizer,
    ErrorStateManager,
    create_enhanced_sidebar
)
# Import RUNX1 integration tools
from tools.runx1_integration_plan import (
    integrate_runx1_analysis,
    design_runx1_intervention,
    create_runx1_genomic_browser,
    generate_runx1_clinical_report
)

st.set_page_config(
    page_title="Patient Digital Twin",
    page_icon="🧬",
    layout="wide",
)

# Initialize demo manager
demo_manager = DemoFlowManager()

# Enhanced sidebar with context
create_enhanced_sidebar("Patient Digital Twin", {"gene_symbol": "RUNX1", "focus_area": "digital_twin_modeling"})

# Progress indicator for demo mode
demo_manager.show_progress_indicator()

# --- Environment & API Configuration ---
def get_command_center_url():
    """Retrieves the Command Center URL from environment variables."""
    return "https://crispro--command-center-commandcenter-api.modal.run"

COMMAND_CENTER_URL = get_command_center_url()
ASSESS_THREAT_ENDPOINT = "/workflow/assess_threat"

# --- API Communication ---
def assess_threat(gene_symbol, protein_change):
    """Calls the CommandCenter to run the hardened threat assessment."""
    payload = {"gene_symbol": gene_symbol, "protein_change": protein_change}
    try:
        response = requests.post(f"{COMMAND_CENTER_URL}{ASSESS_THREAT_ENDPOINT}", json=payload, timeout=60)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        st.error(f"API Error: Failed to connect to CommandCenter. Details: {e}")
        return None

# --- RUNX1 Analysis Functions ---
def analyze_runx1_variant(gene_symbol, protein_change, variant_type="germline"):
    """Analyze RUNX1 variant using comprehensive AI analysis."""
    try:
        # Create variant dictionary
        variant = {
            "gene": gene_symbol,
            "protein_change": protein_change,
            "variant_type": variant_type,
            "id": f"{gene_symbol}_{protein_change}",
            "position": 36207648,  # Default RUNX1 position
            "consequence": "frameshift_variant" if "fs" in protein_change else "missense_variant"
        }
        
        # Run comprehensive RUNX1 analysis
        if variant_type == "germline":
            result = integrate_runx1_analysis(variant, None, 35.0)
        else:
            # For somatic variants, use as second hit
            germline_variant = st.session_state.get("germline_variant", variant)
            result = integrate_runx1_analysis(germline_variant, variant, 35.0)
        
        return result
    except Exception as e:
        st.error(f"RUNX1 Analysis Error: {e}")
        return None

def display_runx1_analysis_results(result, variant_type="germline"):
    """Display RUNX1 analysis results with enhanced visualization."""
    if not result:
        return
    
    # Get the appropriate analysis section
    if variant_type == "germline":
        analysis = result.get("germline_analysis", {})
        st.markdown("#### 🧬 Germline RUNX1 Analysis")
    else:
        analysis = result.get("somatic_analysis", {})
        st.markdown("#### ⚡ Somatic Mutation Analysis")
    
    # Display variant information
    variant = analysis.get("variant", {})
    st.markdown(f"**Gene:** {variant.get('gene', 'Unknown')}")
    st.markdown(f"**Variant:** {variant.get('protein_change', 'Unknown')}")
    st.markdown(f"**Type:** {variant.get('variant_type', 'Unknown')}")
    
    # Display functional impact
    functional_impact = analysis.get("functional_impact", {})
    impact_level = functional_impact.get("impact", "unknown")
    
    if impact_level == "high":
        st.error(f"🚨 **High Impact:** {functional_impact.get('description', 'Significant functional disruption')}")
    elif impact_level == "moderate":
        st.warning(f"⚠️ **Moderate Impact:** {functional_impact.get('description', 'Moderate functional disruption')}")
    else:
        st.info(f"ℹ️ **Low Impact:** {functional_impact.get('description', 'Minimal functional disruption')}")
    
    # Display functional domains affected
    domains = functional_impact.get("functional_domains", [])
    if domains:
        st.markdown("**Affected Domains:**")
        for domain in domains:
            st.markdown(f"- {domain.get('name', 'Unknown')}: {domain.get('description', 'No description')}")
    
    # Display genomic context
    genomic_context = analysis.get("genomic_context", {})
    if genomic_context:
        st.markdown("**Genomic Context:**")
        st.markdown(f"- Position: {genomic_context.get('position', 'Unknown')}")
        st.markdown(f"- Context: {genomic_context.get('context_description', 'No description')}")

# Demo mode introduction
if demo_manager.demo_mode:
    st.info("""
    🎯 **LEAP Grant Demo Mode Active** 
    
    **What you'll see:** This workflow creates a computational "Digital Twin" of RUNX1-FPD patients, 
    modeling the complete progression from germline predisposition to malignant transformation.
    
    **LEAP Grant Relevance:** This directly addresses Focus Area 1 by quantifying both cell-autonomous 
    mechanisms (germline + somatic mutations) and their temporal interaction in disease progression.
    
    **🚀 NEW: Enhanced with AI-Powered RUNX1 Analysis**
    - Comprehensive variant impact assessment
    - Functional domain mapping
    - Progression modeling with two-hit hypothesis
    - Precision intervention design
    """)

# --- Main UI ---
st.title("🧬 Patient Digital Twin: RUNX1-FPD Precision Medicine")

st.markdown("""
### 🎯 AI-Powered Modeling of RUNX1-FPD Progression

The **Patient Digital Twin** transforms the abstract "Two-Hit Hypothesis" into a precise, quantifiable model 
of how RUNX1-FPD patients progress from germline predisposition to malignant transformation.

**🚀 Enhanced AI Analysis:**
1. **🧬 Germline Analysis** - Comprehensive RUNX1 variant assessment with functional domain mapping
2. **⚡ Somatic Analysis** - Second-hit mutation analysis with progression modeling  
3. **📊 Synergistic Modeling** - Two-hit hypothesis validation with intervention opportunities
4. **🎯 Precision Intervention** - AI-powered therapeutic design and timeline estimation
""")

# RUNX1-FPD Context
st.markdown("""
### 🩸 RUNX1-FPD: From Germline Risk to Malignant Reality

RUNX1 Familial Platelet Disorder patients carry a germline mutation that creates a vulnerable hematopoietic environment. 
The critical question is: **when and how does this vulnerability become malignancy?**

Our AI-powered Digital Twin provides unprecedented insight into this progression.
""")

# Patient Journey Context
if demo_manager.demo_mode:
    PatientJourneyVisualizer.display_journey_timeline("clonal_hematopoiesis")
    
    st.markdown("""
    ### 🔬 LEAP Grant Focus Area 1: Mechanistic Understanding
    
    **Cell-Autonomous Factors We Model:**
    - Germline RUNX1 mutation effects on hematopoietic stem cells
    - Somatic mutation acquisition and clonal expansion
    - Synergistic effects of multiple mutations
    
    **Cell-Nonautonomous Factors:**
    - Bone marrow microenvironment changes
    - Inflammatory signals promoting clonal evolution
    - Age-related selective pressures
    """)

st.info("ℹ️ This analysis uses **live AI models** and may take a moment to complete.")
st.markdown("---")

# --- Enhanced Analysis Workflow ---
col1, col2 = st.columns(2)

with col1:
    st.subheader("🧬 Step 1: Germline RUNX1 Analysis")
    
    if demo_manager.demo_mode:
        st.info("""
        **💡 Enhanced AI Analysis:**
        Our AI system analyzes the germline RUNX1 mutation's impact on:
        - Functional domain disruption (Runt domain, TAD domain)
        - Hematopoietic stem cell function
        - Clonal expansion potential
        - Intervention opportunities
        """)
    
    st.markdown("**Analyze the patient's inherited germline mutation:**")
    
    with st.form("germline_analysis_form"):
        gene_1 = st.text_input("Gene Symbol", "RUNX1", key="gene_1", help="The gene carrying the germline mutation")
        p_change_1 = st.text_input("Protein Change (HGVSp)", "p.Arg135fs", key="p_change_1", help="Specific germline mutation")
        submitted_1 = st.form_submit_button("🔍 Run AI Analysis")

    if submitted_1:
        with st.spinner("🧠 Running comprehensive RUNX1 analysis..."):
            result = analyze_runx1_variant(gene_1, p_change_1, "germline")
            if result:
                st.session_state.runx1_germline_result = result
                st.session_state.germline_variant = {
                    "gene": gene_1,
                    "protein_change": p_change_1,
                    "variant_type": "germline",
                    "id": f"{gene_1}_{p_change_1}",
                    "position": 36207648,
                    "consequence": "frameshift_variant" if "fs" in p_change_1 else "missense_variant"
                }
                st.session_state.somatic_result = None # Reset second hit

if 'runx1_germline_result' in st.session_state and st.session_state.runx1_germline_result:
    with col1:
        result = st.session_state.runx1_germline_result
        st.success("✅ **Germline Analysis Complete**")
        
        # Display enhanced results
        display_runx1_analysis_results(result, "germline")
        
        # Show progression model
        progression_model = result.get("progression_model", {})
        if progression_model:
            st.markdown("#### 📊 Progression Risk Assessment")
            transformation_prob = progression_model.get("transformation_probability", {})
            current_prob = transformation_prob.get("current_probability", 0.0)
            
            if current_prob > 0.3:
                st.error(f"🚨 **High Risk:** {current_prob:.1%} transformation probability")
            elif current_prob > 0.1:
                st.warning(f"⚠️ **Moderate Risk:** {current_prob:.1%} transformation probability")
            else:
                st.info(f"ℹ️ **Low Risk:** {current_prob:.1%} transformation probability")
        
        # Mark progress
        if demo_manager.demo_mode:
            demo_manager.mark_step_complete("germline_analysis")

with col2:
    if 'runx1_germline_result' in st.session_state and st.session_state.runx1_germline_result:
        st.subheader("⚡ Step 2: Somatic Mutation Analysis")
        
        if demo_manager.demo_mode:
            st.info("""
            **💡 Enhanced AI Analysis:**
            Our AI system models the somatic "second hit" by analyzing:
            - Clonal evolution dynamics
            - Synergistic effects with germline mutation
            - Malignant transformation timeline
            - Intervention windows
            """)
        
        st.markdown("**Analyze the somatic mutation found in leukemic cells:**")
        
        with st.form("somatic_analysis_form"):
            gene_2 = st.text_input("Gene Symbol", "ASXL1", key="gene_2", help="Gene carrying the somatic mutation")
            p_change_2 = st.text_input("Protein Change (HGVSp)", "p.Gly646fs", key="p_change_2", help="Specific somatic mutation")
            submitted_2 = st.form_submit_button("⚡ Run AI Analysis")

        if submitted_2:
            with st.spinner("🧠 Modeling malignant transformation..."):
                result = analyze_runx1_variant(gene_2, p_change_2, "somatic")
                if result:
                    st.session_state.runx1_somatic_result = result

if 'runx1_somatic_result' in st.session_state and st.session_state.runx1_somatic_result:
    with col2:
        result = st.session_state.runx1_somatic_result
        st.success("✅ **Somatic Analysis Complete**")
        
        # Display enhanced results
        display_runx1_analysis_results(result, "somatic")
        
        # Show enhanced progression model
        progression_model = result.get("progression_model", {})
        if progression_model:
            st.markdown("#### 📊 Two-Hit Progression Model")
            transformation_prob = progression_model.get("transformation_probability", {})
            current_prob = transformation_prob.get("current_probability", 0.0)
            
            st.markdown(f"**Current Risk:** {current_prob:.1%}")
            st.markdown(f"**Risk Tier:** {result.get('risk_stratification', {}).get('risk_tier', 'Unknown')}")
            st.markdown(f"**Surveillance:** {result.get('risk_stratification', {}).get('surveillance_frequency', 'Unknown')}")
        
        # Mark progress
        if demo_manager.demo_mode:
            demo_manager.mark_step_complete("somatic_analysis")

st.markdown("---")

# Enhanced Digital Twin Summary
if 'runx1_somatic_result' in st.session_state and st.session_state.runx1_somatic_result:
    st.subheader("🎯 Digital Twin Complete: Enhanced Patient Model")
    
    # Get both results
    germline_result = st.session_state.runx1_germline_result
    somatic_result = st.session_state.runx1_somatic_result
    
    # Generate clinical report
    clinical_report = generate_runx1_clinical_report(somatic_result, {"age": 35, "id": "RUNX1_001"})
    
    # Display comprehensive analysis
    col3, col4 = st.columns(2)
    
    with col3:
        st.markdown("#### 🧬 Germline Profile")
        germline_analysis = germline_result.get("germline_analysis", {})
        variant = germline_analysis.get("variant", {})
        functional_impact = germline_analysis.get("functional_impact", {})
        
        st.markdown(f"**Gene:** {variant.get('gene', 'Unknown')}")
        st.markdown(f"**Variant:** {variant.get('protein_change', 'Unknown')}")
        st.markdown(f"**Impact:** {functional_impact.get('impact', 'Unknown')}")
        st.markdown(f"**Domains:** {len(functional_impact.get('functional_domains', []))} affected")
    
    with col4:
        st.markdown("#### ⚡ Somatic Profile")
        somatic_analysis = somatic_result.get("somatic_analysis", {})
        variant = somatic_analysis.get("variant", {})
        functional_impact = somatic_analysis.get("functional_impact", {})
        
        st.markdown(f"**Gene:** {variant.get('gene', 'Unknown')}")
        st.markdown(f"**Variant:** {variant.get('protein_change', 'Unknown')}")
        st.markdown(f"**Impact:** {functional_impact.get('impact', 'Unknown')}")
        st.markdown(f"**Domains:** {len(functional_impact.get('functional_domains', []))} affected")
    
    # Clinical interpretation with AI insights
    st.markdown("#### 🏥 AI-Powered Clinical Interpretation")
    
    progression_model = somatic_result.get("progression_model", {})
    risk_stratification = somatic_result.get("risk_stratification", {})
    
    risk_tier = risk_stratification.get("risk_tier", "unknown")
    current_prob = progression_model.get("transformation_probability", {}).get("current_probability", 0.0)
    
    if risk_tier == "high" or current_prob > 0.3:
        st.error(f"""
        🚨 **High-Risk Profile Confirmed**
        
        **Transformation Risk:** {current_prob:.1%}
        **Risk Tier:** {risk_tier.title()}
        **Surveillance:** {risk_stratification.get('surveillance_frequency', 'Unknown')}
        
        Both mutations create a synergistic high-risk profile for malignant transformation. 
        **Immediate precision intervention is recommended.**
        """)
    elif risk_tier == "moderate" or current_prob > 0.1:
        st.warning(f"""
        ⚠️ **Moderate-Risk Profile**
        
        **Transformation Risk:** {current_prob:.1%}
        **Risk Tier:** {risk_tier.title()}
        **Surveillance:** {risk_stratification.get('surveillance_frequency', 'Unknown')}
        
        Significant risk profile requiring enhanced monitoring and intervention planning.
        """)
    else:
        st.info(f"""
        ℹ️ **Lower-Risk Profile**
        
        **Transformation Risk:** {current_prob:.1%}
        **Risk Tier:** {risk_tier.title()}
        **Surveillance:** {risk_stratification.get('surveillance_frequency', 'Unknown')}
        
        Current mutation profile suggests manageable risk with regular surveillance.
        """)
    
    # Show intervention opportunities
    intervention_opportunities = somatic_result.get("intervention_opportunities", [])
    if intervention_opportunities:
        st.markdown("#### 🎯 Precision Intervention Opportunities")
        for i, opportunity in enumerate(intervention_opportunities, 1):
            st.markdown(f"**{i}. {opportunity.get('intervention_type', 'Unknown')}**")
            st.markdown(f"   - Target: {opportunity.get('target', 'Unknown')}")
            st.markdown(f"   - Timing: {opportunity.get('timing', 'Unknown')}")
            st.markdown(f"   - Success Rate: {opportunity.get('success_probability', 'Unknown')}")
    
    # LEAP Grant context
    if demo_manager.demo_mode:
        st.markdown("""
        ### 🎯 LEAP Grant Impact: AI-Powered Precision Medicine
        
        **Focus Area 1 Achievement:**
        - ✅ **Mechanistic Understanding:** Comprehensive analysis of germline + somatic interactions
        - ✅ **Temporal Modeling:** Quantified progression from predisposition to malignancy
        - ✅ **Functional Mapping:** Identified specific domains and pathways disrupted
        - ✅ **Risk Stratification:** Patient-specific risk assessment with surveillance recommendations
        
        **Focus Area 2 Foundation:**
        - 🎯 **Intervention Design:** AI-powered therapeutic target identification
        - 📊 **Precision Timing:** Optimal intervention window calculation
        - ⏱️ **Success Prediction:** Evidence-based success probability modeling
        """)
    
    # Enhanced next steps with precision intervention
    st.markdown("#### 🚀 AI-Recommended Next Steps")
    
    col5, col6 = st.columns(2)
    
    with col5:
        if st.button("🎯 Design Precision Intervention", type="primary"):
            with st.spinner("🧠 Designing CRISPR therapeutic intervention..."):
                # Get the target variant (somatic for intervention)
                target_variant = somatic_result.get("somatic_analysis", {}).get("variant", {})
                intervention_result = design_runx1_intervention(target_variant, "crispr_correction")
                
                if intervention_result:
                    st.session_state.intervention_result = intervention_result
                    st.success("✅ **Precision Intervention Designed!**")
                    
                    # Display intervention summary
                    guide_design = intervention_result.get("guide_design", {})
                    st.markdown("**Intervention Summary:**")
                    st.markdown(f"- Guides Found: {len(guide_design.get('all_guides', []))}")
                    st.markdown(f"- Success Probability: {intervention_result.get('success_probability', 0.0):.1%}")
                    st.markdown(f"- Timeline: {intervention_result.get('timeline_estimate', {}).get('total_duration_weeks', 'Unknown')} weeks")
                    
                    # Store for next page
                    st.session_state.runx1_intervention_ready = True
    
    with col6:
        if st.button("📊 Generate Clinical Report", type="secondary"):
            with st.spinner("📋 Generating comprehensive clinical report..."):
                # Clinical report already generated above
                st.success("✅ **Clinical Report Generated!**")
                
                # Display key recommendations
                recommendations = somatic_result.get("clinical_recommendations", {})
                immediate_actions = recommendations.get("immediate_actions", [])
                
                if immediate_actions:
                    st.markdown("**Immediate Actions:**")
                    for action in immediate_actions[:3]:  # Show top 3
                        st.markdown(f"- {action}")
                
                # Store for download/sharing
                st.session_state.clinical_report_ready = True
    
    # Navigation to next steps
    st.markdown("---")
    st.markdown("#### 🌟 Continue Your RUNX1-FPD Journey")
    
    nav_col1, nav_col2, nav_col3 = st.columns(3)
    
    with nav_col1:
        if st.button("🧬 Genomic Browser", help="Visualize RUNX1 variants in genomic context"):
            st.session_state.genomic_browser_ready = True
            st.info("🧬 Genomic Browser integration coming in Week 2!")
    
    with nav_col2:
        if st.button("🌱 Seed & Soil Analysis", help="Analyze metastatic potential"):
            if 'runx1_somatic_result' in st.session_state:
                st.session_state.validated_threat = st.session_state.runx1_somatic_result
                st.switch_page("pages/4_🌱_Seed_and_Soil_Analysis.py")
    
    with nav_col3:
        if st.button("⚔️ Threat Assessor", help="Advanced threat analysis"):
            if 'runx1_somatic_result' in st.session_state:
                st.switch_page("pages/2_⚔️_Threat_Assessor.py")

# Footer with enhanced capabilities
if demo_manager.demo_mode:
    st.markdown("---")
    st.markdown("""
    ### 🏆 Enhanced Platform Capabilities Demonstrated
    
    This AI-powered RUNX1-FPD Digital Twin showcases:
    - **🧠 Multi-AI Integration:** Combining multiple AI models for comprehensive analysis
    - **🎯 Precision Medicine:** Patient-specific risk assessment and intervention design
    - **📊 Quantitative Modeling:** Evidence-based progression and success probability calculations
    - **⚡ Real-time Analysis:** Live AI processing with professional-grade results
    - **🏥 Clinical Integration:** Physician-ready reports and actionable recommendations
    
    *This represents the future of precision medicine: AI-powered, patient-specific, intervention-focused.*
    """) 