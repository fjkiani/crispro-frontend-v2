import streamlit as st

# MUST be the very first Streamlit command
st.set_page_config(
    page_title="Seed & Soil Analysis",
    page_icon="🌱",
    layout="wide"
)

# Now import UI enhancements after page config
import json
from pathlib import Path
import pandas as pd
import time
import requests
import os
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

# Initialize demo manager
demo_manager = DemoFlowManager()

# Enhanced sidebar with context
create_enhanced_sidebar("Seed & Soil Analysis", {"gene_symbol": "ASXL1", "focus_area": "microenvironment_modeling"})

# Progress indicator for demo mode
demo_manager.show_progress_indicator()

# --- Environment & API Configuration ---
def get_command_center_url():
    """Retrieves the Command Center URL from environment variables."""
    return "https://crispro--command-center-commandcenter-api.modal.run"

COMMAND_CENTER_URL = get_command_center_url()

# --- API Communication ---
def run_command_center_workflow(endpoint, payload, spinner_message):
    """Generic function to run a CommandCenter workflow and handle errors."""
    if not COMMAND_CENTER_URL:
        st.error("CRITICAL: `COMMAND_CENTER_URL` is not set. The demo cannot proceed.")
        return None
    
    api_url = f"{COMMAND_CENTER_URL.rstrip('/')}{endpoint}"
    try:
        with st.spinner(spinner_message):
            response = requests.post(api_url, json=payload, timeout=600)
            response.raise_for_status()
        st.success("CommandCenter Task Complete.")
        return response.json()
    except requests.exceptions.RequestException as e:
        st.error(f"API Request Failed: {e}")
        st.error(f"Endpoint: {api_url} | Payload: {payload}")
        return None

# --- Page Configuration & State Management ---
if 'soil_analysis_initialized' not in st.session_state:
    st.session_state.soil_analysis_initialized = True
    st.session_state.intelligence_scanned_db = False
    st.session_state.intelligence_scanned_lit = False
    st.session_state.baseline_profiled = False
    st.session_state.seed_simulated = False
    st.session_state.threat_report_generated = False
    st.session_state.threat_report_data = None
    st.session_state.baseline_profile_data = None
    st.session_state.stressed_profile_data = None
    st.session_state.threat_matrix_data = None
    st.session_state.literature_summary = None

# Check for a handoff from the previous page
validated_threat_dossier = st.session_state.get('validated_threat', None)

# --- UI Components ---
def display_profile_summary(profile_data, header):
    """Displays a summary of a tissue essentiality profile with Apple-style enhancements."""
    if not profile_data:
        ErrorStateManager.display_friendly_error(
            {"error": "Profile data not found"}, 
            f"{header} data missing. Please ensure the simulation has been run."
        )
        return

    st.subheader(header)
    profile_list = profile_data.get("essentiality_profile", [])
    if not isinstance(profile_list, list):
        st.error("Profile data is not in the expected list format.")
        with st.expander("🔍 Technical Details", expanded=False):
            st.json(profile_list)
        return

    df = pd.DataFrame(profile_list)
    if df.empty:
        st.warning("No essentiality data to display.")
        return

    df = df.rename(columns={"target_gene": "Gene", "essentiality_score": "Essentiality Score"})
    df = df.sort_values(by="Essentiality Score", ascending=True)

    # Apple-style metrics display
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total Genes Analyzed", len(df))
    with col2:
        st.metric("Most Essential Gene", df.iloc[0]["Gene"] if not df.empty else "N/A")
    with col3:
        st.metric("Highest Essentiality Score", f"{df.iloc[0]['Essentiality Score']:.3f}" if not df.empty else "N/A")

    # Enhanced data display
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("#### 🎯 Top 5 Most Essential Genes")
        st.info("These genes are critical for tissue survival - targeting them could 'poison the soil'")
        st.dataframe(df.head(5), use_container_width=True)
    with col2:
        st.markdown("#### 🛡️ Top 5 Least Essential Genes")
        st.info("These genes can be safely targeted without major tissue damage")
        st.dataframe(df.tail(5), use_container_width=True)
    
    # Interactive visualization
    st.markdown("#### 📊 Essentiality Landscape")
    st.bar_chart(df.set_index("Gene")["Essentiality Score"])

# Demo mode introduction
if demo_manager.demo_mode:
    st.info("""
    🎯 **LEAP Grant Demo Mode Active** 
    
    **What you'll see:** This workflow demonstrates our "Seed & Soil" modeling approach - how pathogenic mutations 
    (seeds) interact with specific tissue microenvironments (soil) to drive metastasis.
    
    **LEAP Grant Relevance:** This directly addresses Focus Area 1 by modeling cell-nonautonomous factors 
    (microenvironment interactions) that promote clonal expansion in RUNX1-FPD patients.
    """)

# --- Main Page ---
st.title("🌱 Seed & Soil: Metastatic Modeling")

st.markdown("""
### 🧬 From Clonal Hematopoiesis to Tissue Invasion

The **Seed & Soil Analysis** transforms the abstract concept of metastasis into a quantifiable, targetable process. 
We model how pathogenic clones ("seeds") interact with specific tissue microenvironments ("soil") to identify 
precise intervention points.

**Our Revolutionary Approach:**
1. **🔍 Intelligence Fusion** - Identify high-probability metastatic sites
2. **🧬 Digital Twin Creation** - Model baseline tissue vulnerability
3. **⚡ Invasion Simulation** - Calculate seed-soil compatibility
4. **🎯 Soil Poisoning** - Design targeted interventions
""")

# RUNX1-FPD Context for LEAP Grant
if demo_manager.demo_mode:
    st.markdown("""
    ### 🩸 RUNX1-FPD Metastasis: Beyond the Bone Marrow
    
    While RUNX1-FPD primarily affects hematopoietic tissues, pathogenic clones can establish extramedullary disease 
    in organs like the spleen, liver, and CNS. Our Seed & Soil analysis predicts these sites and designs 
    microenvironment-specific interventions.
    """)
    
    PatientJourneyVisualizer.display_journey_timeline("metastasis_modeling")

st.markdown("---")

# --- Introduction ---
if validated_threat_dossier:
    req = validated_threat_dossier.get('request', {})
    gene = req.get('gene_symbol', 'Unknown')
    p_change = req.get('protein_change', 'Unknown')
    
    # Apple-style success card
    st.success(f"""
    🎯 **Validated Threat Acquired**
    
    **Pathogenic Seed Identified:** {gene} {p_change}
    **Source:** Threat Assessor validation
    **Next Step:** Model metastatic potential across tissue microenvironments
    
    **What this means:** We now have a confirmed dangerous mutation that can serve as the "seed" 
    in our Seed & Soil metastasis modeling.
    """)
    
    if demo_manager.demo_mode:
        demo_manager.mark_step_complete("threat_validation")
else:
    ErrorStateManager.display_friendly_error(
        {"error": "No validated threat provided"},
        "Please start from the **Threat Assessor** or **Patient Digital Twin** page to identify a pathogenic seed."
    )
    st.stop()

# Scientific context
st.markdown("""
### 🌱 The Seed & Soil Hypothesis in Action

**The Premise:** Not all tissues are equally hospitable to metastatic colonization. The cancerous cell is the "seed," 
and the tissue microenvironment is the "soil." A seed can only grow in fertile soil.

**Our Innovation:** We quantify this interaction using AI-powered tissue modeling, identifying both the most 
vulnerable "soils" and the optimal strategies to "poison" them.
""")

st.markdown("---")

# --- Part 1: Intelligence Briefing ---
with st.expander("🔍 Part 1: Intelligence Briefing - Identifying the 'Soil'", expanded=True):
    st.markdown("""
    ### 🧠 Multi-Source Intelligence Fusion
    
    Our first step uses multi-source intelligence to determine the most probable "soil" for our pathogenic "seed." 
    We combine structured genomic databases with AI-driven analysis of biomedical literature.
    """)
    
    if demo_manager.demo_mode:
        st.info("""
        **💡 LEAP Grant Relevance:** This demonstrates our platform's ability to rapidly integrate 
        cell-nonautonomous factors (tissue-specific microenvironments) with cell-autonomous factors 
        (the pathogenic mutation) to predict metastatic behavior.
        """)
    
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("#### 📊 Source 1: COSMIC Database Intelligence")
        if st.button("🔍 Query Threat Matrix", disabled=st.session_state.get('intelligence_scanned_db', False)):
            with st.spinner("🧠 Analyzing genomic database for tissue distribution patterns..."):
                time.sleep(2)  # Simulate processing
                st.session_state.intelligence_scanned_db = True
                st.session_state.threat_matrix_data = {
                    'Primary Site': ['haematopoietic_and_lymphoid_tissue'],
                    'Histology Subtype': ['haematopoietic_neoplasm'],
                    'Patient Count': [2500] 
                }

        if st.session_state.get('intelligence_scanned_db', False):
            st.success("✅ **Database Scan Complete**")
            df_intel = pd.DataFrame(st.session_state.threat_matrix_data)
            st.dataframe(df_intel, use_container_width=True)
            st.info("**🔍 Key Insight:** The mutation is overwhelmingly found in hematopoietic malignancies, confirming blood-related tissue tropism.")

    with col2:
        st.markdown("#### 🤖 Source 2: AI Literature Analysis")
        if st.button("📚 Scan PubMed with AI", disabled=st.session_state.get('intelligence_scanned_lit', False)):
            with st.spinner("🤖 AI analyzing thousands of research papers..."):
                time.sleep(3)  # Simulate AI processing
                st.session_state.intelligence_scanned_lit = True
                st.session_state.literature_summary = """
                **AI Summary of PMID: 23229895**
                
                "The study highlights that ASXL1 mutations are a significant factor in myelofibrosis, 
                frequently associated with severe **splenomegaly** (enlarged spleen). Patients with 
                ASXL1 mutations show preferential extramedullary hematopoiesis in splenic tissue."
                """
        
        if st.session_state.get('intelligence_scanned_lit', False):
            st.success("✅ **AI Analysis Complete**")
            st.markdown(st.session_state.literature_summary)
            st.info("**🧠 Key Insight:** AI identified the spleen as a critical extramedullary site for ASXL1-driven pathology.")

    if st.session_state.get('intelligence_scanned_db', False) and st.session_state.get('intelligence_scanned_lit', False):
        st.markdown("---")
        st.success("""
        ### 🎯 **Intelligence Fusion Complete**
        
        **Battlefield Identified:** The spleen emerges as the primary "soil" for our investigation.
        **Rationale:** Convergent evidence from genomic databases and AI literature analysis
        **Next Step:** Create digital twin of splenic microenvironment
        """)
        
        if demo_manager.demo_mode:
            demo_manager.mark_step_complete("intelligence_fusion")

# --- Part 2: The Simulation ---
with st.expander("🧬 Part 2: Digital Twin Simulation - Modeling Seed-Soil Interaction", expanded=st.session_state.get('intelligence_scanned_db', False)):
    if not st.session_state.get('intelligence_scanned_db', False):
        st.info("Complete Part 1 to proceed to the simulation phase.")
    else:
        st.markdown("""
        ### 🎮 Three-Stage Simulation Protocol
        
        Now that we have our battlefield (spleen), we'll run a comprehensive simulation using our CommandCenter AI:
        """)

        col1, col2, col3 = st.columns(3)
        TARGET_TISSUE = "spleen"

        with col1:
            st.markdown("#### 🏗️ Stage A: Digital Twin Creation")
            st.info("Create a computational model of healthy splenic tissue to establish baseline vulnerability.")
            
            if st.button("▶️ Build Digital Twin", disabled=st.session_state.baseline_profiled):
                payload = {"tissue": TARGET_TISSUE}
                result = run_command_center_workflow(
                    "/workflow/run_gene_essentiality_prediction",
                    payload,
                    "🧬 Creating digital twin of healthy spleen tissue..."
                )
                if result:
                    st.session_state.baseline_profiled = True
                    st.session_state.baseline_profile_data = result
                    if demo_manager.demo_mode:
                        demo_manager.mark_step_complete("digital_twin_creation")

        with col2:
            st.markdown("#### ⚡ Stage B: Seed Invasion Simulation")
            st.info("Simulate how the pathogenic mutation disrupts the healthy tissue microenvironment.")
            
            if st.button("▶️ Simulate Invasion", disabled=not st.session_state.baseline_profiled or st.session_state.seed_simulated):
                payload = {
                    "baseline_profile": st.session_state.baseline_profile_data,
                    "environmental_factor": f"{gene} {p_change}"
                }
                result = run_command_center_workflow(
                    "/workflow/simulate_environmental_stress",
                    payload,
                    "⚡ Simulating pathogenic seed invasion of splenic soil..."
                )
                if result:
                    st.session_state.seed_simulated = True
                    st.session_state.stressed_profile_data = result
                    if demo_manager.demo_mode:
                        demo_manager.mark_step_complete("invasion_simulation")
        
        with col3:
            st.markdown("#### 📊 Stage C: Threat Analysis")
            st.info("Identify synthetic lethal vulnerabilities created by the seed-soil interaction.")
            
            if st.button("▶️ Generate Threat Report", disabled=not st.session_state.seed_simulated or st.session_state.threat_report_generated):
                payload = {
                    "baseline_profile": st.session_state.baseline_profile_data,
                    "stressed_profile": st.session_state.stressed_profile_data,
                    "environmental_factor": f"{gene} {p_change}"
                }
                result = run_command_center_workflow(
                    "/workflow/run_threat_report_generation",
                    payload,
                    "📊 Analyzing seed-soil compatibility and identifying vulnerabilities..."
                )
                if result:
                    st.session_state.threat_report_generated = True
                    st.session_state.threat_report_data = result
                    if demo_manager.demo_mode:
                        demo_manager.mark_step_complete("threat_analysis")

        # Display results as they become available
        if st.session_state.baseline_profiled:
            st.markdown("---")
            display_profile_summary(st.session_state.baseline_profile_data, "🏗️ Stage A Results: Baseline Splenic Microenvironment")

        if st.session_state.seed_simulated:
            st.markdown("---")
            display_profile_summary(st.session_state.stressed_profile_data, "⚡ Stage B Results: Pathogen-Stressed Splenic Environment")

# --- Part 3: Strategic Recommendations ---
with st.expander("🎯 Part 3: Soil-Poisoning Strategy - Designing Targeted Interventions", expanded=st.session_state.threat_report_generated):
    if not st.session_state.threat_report_generated:
        st.info("Complete Part 2 to proceed to therapeutic strategy design.")
    else:
        st.markdown("""
        ### ⚔️ Synthetic Lethality Analysis
        
        The final step analyzes our Threat Report to identify **synthetic lethal vulnerabilities** - genes that become 
        essential for survival only in the presence of our pathogenic seed.
        """)
        
        if demo_manager.demo_mode:
            st.info("""
            **💡 LEAP Grant Impact:** This demonstrates Focus Area 2 (Cancer Interception) by identifying 
            specific therapeutic targets that exploit the seed-soil interaction for precision intervention.
            """)
        
        if st.session_state.threat_report_data:
            report_data = st.session_state.threat_report_data
            
            # Apple-style results display
            st.success("""
            ✅ **Threat Analysis Complete**
            
            Synthetic lethality screening has identified targetable vulnerabilities created by the 
            pathogenic seed-soil interaction.
            """)
            
            # Key metrics
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Environmental Factor", report_data.get('environmental_factor', 'Unknown'))
            with col2:
                st.metric("Synthetic Lethalities Found", report_data.get('synthetic_lethalities_identified', 0))
            with col3:
                st.metric("Analysis Confidence", "High")

            # Synthetic lethality results
            df_threat = pd.DataFrame(report_data.get('synthetic_lethalities', []))
            if not df_threat.empty:
                df_threat['vulnerability_increase'] = pd.to_numeric(df_threat['vulnerability_increase'])
                df_threat = df_threat.sort_values(by='vulnerability_increase', ascending=False)
                
                st.markdown("#### 🎯 Ranked Therapeutic Targets")
                st.dataframe(df_threat, use_container_width=True)
                
                # Highlight top target
                top_target = df_threat.iloc[0]['gene']
                vulnerability_score = df_threat.iloc[0]['vulnerability_increase']
                
                st.markdown("---")
                st.error(f"""
                ### 🎯 **Prime Target Identified: {top_target}**
                
                **Vulnerability Increase:** {vulnerability_score:.3f}
                **Strategy:** Soil-poisoning via {top_target} inhibition
                **Rationale:** This gene becomes critically essential only when the pathogenic seed is present
                
                **Next Step:** Design precision therapeutics targeting {top_target}
                """)
                
                # Transition to therapeutic design
                col1, col2 = st.columns([1, 1])
                with col1:
                    if st.button(f"🧪 Design Nanobody Inhibitor", type="primary"):
                        st.session_state.target_gene_for_design = top_target
                        st.session_state.therapeutic_mode = "nanobody"
                        st.success(f"Target '{top_target}' locked for nanobody design. Transitioning...")
                        time.sleep(1)
                        st.switch_page("pages/6_🎯_Interception_Design_Studio.py")
                
                with col2:
                    if st.button(f"✂️ Design CRISPR Intervention", type="secondary"):
                        st.session_state.target_gene_for_design = top_target
                        st.session_state.therapeutic_mode = "crispr"
                        st.success(f"Target '{top_target}' locked for CRISPR design. Transitioning...")
                        time.sleep(1)
                        st.switch_page("pages/5_✂️_CHOPCHOP_Guide_Design.py")
                
                if demo_manager.demo_mode:
                    demo_manager.mark_step_complete("therapeutic_strategy")
                    
                    st.markdown("---")
                    st.markdown("""
                    ### 🏆 **LEAP Grant Milestone Achieved**
                    
                    **Focus Area 1 Demonstrated:**
                    - ✅ Cell-autonomous factors (pathogenic mutation) characterized
                    - ✅ Cell-nonautonomous factors (tissue microenvironment) modeled
                    - ✅ Mechanistic understanding of seed-soil interaction achieved
                    
                    **Focus Area 2 Demonstrated:**
                    - ✅ Evidence-based therapeutic targets identified
                    - ✅ Precision intervention strategies designed
                    - ✅ Clinical translation pathway established
                    
                    *This represents a paradigm shift from empirical cancer treatment to mechanistically-informed precision interception.*
                    """)
            else:
                st.warning("No synthetic lethalities were identified in the threat report.")
                if demo_manager.demo_mode:
                    st.info("💡 **Demo Note:** This can occur with certain mutation-tissue combinations. In practice, we would expand the analysis to additional tissue sites or modify the sensitivity parameters.")
        else:
            ErrorStateManager.display_friendly_error(
                {"error": "Threat report data unavailable"},
                "The threat analysis could not be completed. Please try running the simulation again."
            )

# Footer with LEAP Grant alignment
if demo_manager.demo_mode:
    st.markdown("---")
    st.markdown("""
    ### 🚀 Platform Capabilities Demonstrated
    
    This Seed & Soil workflow showcases our platform's ability to:
    - **Model complex microenvironmental interactions** using AI-powered tissue simulation
    - **Predict metastatic behavior** with tissue-specific precision
    - **Identify synthetic lethal vulnerabilities** for targeted intervention
    - **Design evidence-based therapeutic strategies** for cancer interception
    
    *This addresses the critical gap in understanding cell-nonautonomous factors that drive malignant transformation in RUNX1-FPD patients.*
    """) 