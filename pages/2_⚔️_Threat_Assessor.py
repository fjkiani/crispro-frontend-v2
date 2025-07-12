import streamlit as st

# MUST be the very first Streamlit command
st.set_page_config(
    page_title="Threat Assessor",
    page_icon="⚔️",
    layout="wide"
)

# Now import UI enhancements after page config
import requests
import json
import pandas as pd

# Import our new enhancement utilities
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tools.ui_enhancements import (
    DemoFlowManager, ContextualExplainer, ScientificNarrator, 
    ResultsVisualizer, create_enhanced_sidebar, apply_leap_styling,
    show_loading_with_context,
    LEAPGrantAligner,
    AppleStyleResultsDisplay,
    PatientJourneyVisualizer,
    ErrorStateManager
)

# Apply consistent styling
apply_leap_styling()

# Initialize demo flow manager
demo_manager = DemoFlowManager()

# Enhanced sidebar with context
create_enhanced_sidebar("Threat Assessor", {"gene_symbol": "ASXL1", "focus_area": "mechanism_dissection"})

# Progress indicator for demo mode
demo_manager.show_progress_indicator()

# Cell-autonomous vs cell-nonautonomous explanation
if demo_manager.demo_mode:
    LEAPGrantAligner.explain_cell_autonomous_vs_nonautonomous()

# Main content with enhanced explanations
if demo_manager.demo_mode:
    st.info("""
    🎯 **LEAP Grant Demo Mode Active** 
    
    **What you'll see:** This workflow demonstrates our "Triumvirate Protocol" - a three-layer approach to variant assessment 
    that combines instant bioinformatic screening with AI-powered deep analysis.
    
    **LEAP Grant Relevance:** This directly addresses Focus Area 1 (Dissecting Mechanisms) by providing rapid, 
    multi-source analysis of genetic variants and their cell-autonomous effects.
    """)

# Introduction with dynamic context
st.markdown("""
# ⚔️ Threat Assessor

The **Threat Assessor** uses our proprietary **Triumvirate Protocol** to rapidly evaluate genetic variants:

1. **🔍 Truncation Sieve** - Instant detection of protein-breaking mutations
2. **🧠 Zeta Oracle** - AI-powered analysis of subtle variants  
3. **📚 Intelligence Fusion** - Integration with clinical databases and literature

This multi-layered approach ensures no dangerous mutation goes undetected while avoiding false alarms.
""")

# RUNX1-FPD Context
st.markdown("""
### 🧬 RUNX1-FPD Context: The "Two-Hit" Model

In RUNX1 Familial Platelet Disorder, patients are born with a germline RUNX1 mutation (the "first hit"). 
The critical question is: **which somatic mutations act as the dangerous "second hit" that triggers malignant transformation?**

Our Threat Assessor identifies these high-risk second hits with unprecedented speed and accuracy.
""")

# Patient Journey Context
if demo_manager.demo_mode:
    PatientJourneyVisualizer.display_journey_timeline("second_hit")
    PatientJourneyVisualizer.display_risk_stratification()

# Analysis Interface
st.header("🔬 Variant Analysis")

col1, col2 = st.columns([2, 1])

with col1:
    st.subheader("Step 1: Analyze the 'First Hit'")
    
    # First Hit Analysis
    gene_symbol_1 = st.text_input(
        "Gene Symbol (First Hit)", 
        value="RUNX1",
        help="The germline mutation in RUNX1-FPD patients"
    )
    
    protein_change_1 = st.text_input(
        "Protein Change (Optional)", 
        value="p.Arg135fs",
        help="Specific mutation if known, e.g., p.Arg135fs"
    )
    
    if st.button("🔍 Assess First Hit Threat", key="first_hit"):
        if gene_symbol_1:
            with st.spinner("🧠 Analyzing germline mutation..."):
                try:
                    # API call
                    api_url = "https://crispro--command-center-commandcenter-api.modal.run/workflow/assess_threat"
                    payload = {
                        "gene_symbol": gene_symbol_1,
                        "protein_change": protein_change_1 if protein_change_1 else None
                    }
                    
                    response = requests.post(api_url, json=payload, timeout=60)
                    
                    if response.status_code == 200:
                        result = response.json()
                        
                        # Use Apple-style display instead of JSON dump
                        AppleStyleResultsDisplay.display_threat_assessment(result)
                        
                        # Mark progress
                        if 'demo_progress' not in st.session_state:
                            st.session_state.demo_progress = {}
                        st.session_state.demo_progress['first_hit_analysis'] = True
                        
                    else:
                        st.error(f"❌ Analysis failed: {response.status_code} - {response.text}")
                        
                except requests.exceptions.RequestException as e:
                    st.error(f"❌ Connection error: {str(e)}")
                except Exception as e:
                    st.error(f"❌ Unexpected error: {str(e)}")

with col2:
    if demo_manager.demo_mode:
        st.info("""
        **💡 Demo Tip**
        
        The "first hit" is the germline RUNX1 mutation that RUNX1-FPD patients are born with. 
        
        Try analyzing "RUNX1 p.Arg135fs" to see how our system identifies this foundational mutation.
        """)

# Second Hit Analysis
st.subheader("Step 2: Analyze the 'Second Hit'")

col3, col4 = st.columns([2, 1])

with col3:
    gene_symbol_2 = st.text_input(
        "Gene Symbol (Second Hit)", 
        value="ASXL1",
        help="Suspected somatic mutation gene"
    )
    
    protein_change_2 = st.text_input(
        "Protein Change (Optional)", 
        value="p.G646fs*12", 
        help="Specific somatic mutation, e.g., p.G646fs*12"
    )
    
    if st.button("⚡ Assess Second Hit Threat", key="second_hit"):
        if gene_symbol_2:
            with st.spinner("🧠 Analyzing somatic mutation..."):
                try:
                    # API call
                    api_url = "https://crispro--command-center-commandcenter-api.modal.run/workflow/assess_threat"
                    payload = {
                        "gene_symbol": gene_symbol_2,
                        "protein_change": protein_change_2 if protein_change_2 else None
                    }
                    
                    response = requests.post(api_url, json=payload, timeout=60)
                    
                    if response.status_code == 200:
                        result = response.json()
                        
                        # Use Apple-style display instead of JSON dump
                        AppleStyleResultsDisplay.display_threat_assessment(result)
                        
                        # Mark progress
                        if 'demo_progress' not in st.session_state:
                            st.session_state.demo_progress = {}
                        st.session_state.demo_progress['second_hit_analysis'] = True
                        
                    else:
                        st.error(f"❌ Analysis failed: {response.status_code} - {response.text}")
                        
                except requests.exceptions.RequestException as e:
                    st.error(f"❌ Connection error: {str(e)}")
                except Exception as e:
                    st.error(f"❌ Unexpected error: {str(e)}")

with col4:
    if demo_manager.demo_mode:
        st.info("""
        **💡 Demo Tip**
        
        The "second hit" is the somatic mutation that transforms pre-malignant clonal hematopoiesis into acute leukemia.
        
        Try "ASXL1 p.G646fs*12" to see how our Truncation Sieve instantly identifies catastrophic mutations.
        """)

# Comparative Analysis Section
if demo_manager.demo_mode:
    st.header("🔬 Understanding the Analysis")
    
    st.markdown("""
    ### 🧠 Why This Matters for LEAP Grant Objectives
    
    **Focus Area 1 - Mechanistic Understanding:**
    - **Cell-Autonomous Factors:** Each mutation's direct effect on protein function
    - **Temporal Progression:** How germline + somatic mutations interact over time
    - **Molecular Mechanisms:** Precise pathways disrupted by each "hit"
    
    **Focus Area 2 - Cancer Interception:**
    - **Early Detection:** Identify dangerous second hits before malignant transformation
    - **Targeted Intervention:** Design therapies specific to the mutation profile
    - **Prevention Strategy:** Intercept the progression from clonal hematopoiesis to leukemia
    """)

# LEAP Grant alignment
if demo_manager.demo_mode:
    LEAPGrantAligner.highlight_focus_area_1({"gene_symbol": gene_symbol_2})

# Footer with contextual help
st.markdown("---")
st.markdown("""
### 🔗 What's Next?

After identifying threatening mutations, proceed to:
- **🌱 Seed & Soil Analysis** - Model how these mutations spread
- **✂️ CRISPR Design** - Create targeted interventions  
- **💊 Drug Discovery** - Find combination therapies

*This analysis forms the foundation for precision cancer interception in RUNX1-FPD patients.*
""")

# Additional context and help
with st.expander("ℹ️ About the Threat Assessor", expanded=False):
    st.markdown("""
    ### How It Works
    
    The Threat Assessor implements our **Intelligence Fusion** doctrine:
    
    1. **Internal Intelligence:** Queries our curated COSMIC database for known cancer variants
    2. **External Intelligence:** Fetches real-time annotations from Ensembl VEP
    3. **AI Analysis:** Applies our Triumvirate Protocol for comprehensive assessment
    4. **Evidence Synthesis:** Integrates clinical trials, efficacy data, and literature
    
    ### Key Features
    
    - **🔍 Truncation Sieve:** Instant detection of frameshift/nonsense mutations
    - **🧠 Zeta Oracle:** AI-powered impact prediction for complex variants
    - **🌐 Multi-Source Integration:** COSMIC + VEP + Literature + Clinical data
    - **📊 Clinical Context:** Active trials and therapeutic evidence
    
    ### LEAP Grant Alignment
    
    This tool directly addresses **Focus Area 1: Dissecting Mechanisms of Cancer Initiation** by providing:
    - Multi-scale analysis from genomic to clinical levels
    - Integration of environmental and cell-autonomous factors  
    - Rapid mechanistic insights that traditionally require extensive research
    """)

# Footer with technical details
if demo_manager.demo_mode:
    st.markdown("---")
    st.caption("🎯 **LEAP Grant Demo Mode** | Powered by CasPro Intelligence Platform | Real-time API integration with CommandCenter") 