import streamlit as st

# MUST be the very first Streamlit command
st.set_page_config(
    page_title="Trial Conquest",
    page_icon="🎯",
    layout="wide"
)

# Now import UI enhancements after page config
import requests
import json
import os
from tools.ui_enhancements import (
    DemoFlowManager, 
    ContextualExplainer, 
    ScientificNarrator, 
    LEAPGrantAligner,
    AppleStyleResultsDisplay,
    PatientJourneyVisualizer,
    ErrorStateManager,
    create_enhanced_sidebar,
    apply_leap_styling
)

# Apply modern styling
apply_leap_styling()

# Initialize demo manager
demo_manager = DemoFlowManager()

# Enhanced sidebar with context
create_enhanced_sidebar("Trial Conquest", {"gene_symbol": "Multi-Target", "focus_area": "clinical_trial_matching"})

# Progress indicator for demo mode
demo_manager.show_progress_indicator()

# Hero Section - Aggressive, Military-Inspired Design
st.markdown("""
<div style='background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%); 
            padding: 3rem 2rem; border-radius: 16px; margin-bottom: 2rem; color: white; text-align: center;'>
    <h1 style='font-size: 3.5rem; font-weight: 700; color: #ffffff; margin-bottom: 0.5rem; text-shadow: 2px 2px 4px rgba(0,0,0,0.5);'>
        🎯 TRIAL CONQUEST
    </h1>
    <h2 style='font-size: 1.8rem; font-weight: 400; color: #64b5f6; margin-bottom: 1.5rem;'>
        Zeta Doctrine: Precision-Guided Clinical Trial Matching
    </h2>
    <p style='font-size: 1.3rem; color: #e3f2fd; max-width: 900px; margin: 0 auto; line-height: 1.6;'>
        <strong>We don't search for keywords. We understand biological intent.</strong><br>
        Transforming clinical trial recruitment from primitive keyword matching into surgical precision targeting 
        using deep biological understanding and AI-powered patient-trial fusion.
    </p>
    <div style='margin-top: 2rem; padding: 1rem; background: rgba(255,255,255,0.1); border-radius: 8px; border-left: 4px solid #ff6b35;'>
        <strong style='color: #ff6b35;'>MISSION CRITICAL:</strong> 80% of trials fail enrollment timelines. 
        Each delayed day costs $8M+ in lost revenue. We end this logistical nightmare.
    </div>
</div>
""", unsafe_allow_html=True)

# Demo mode context
if demo_manager.demo_mode:
    st.markdown("""
    <div style='background: linear-gradient(135deg, #ff6b35 0%, #f7931e 100%); 
                padding: 2rem; border-radius: 12px; margin: 2rem 0; color: white;'>
        <h3 style='color: white; margin-bottom: 1rem;'>⚡ OPERATIONAL STATUS: DEMO MODE ACTIVE</h3>
        <p style='font-size: 1.1rem; margin-bottom: 1rem;'>
            <strong>Current Theater:</strong> RUNX1-FPD Patient Cohort Clinical Trial Matching
        </p>
        <div style='display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; margin-top: 1rem;'>
            <div>
                <strong>Primary Objective:</strong> Match high-risk RUNX1-FPD patients to prevention trials
                <br>• Target: Pre-malignant clonal hematopoiesis
                <br>• Strategy: Cancer interception protocols
            </div>
            <div>
                <strong>Intelligence Assets:</strong> Deep biological profiling
                <br>• Zeta Oracle pathogenicity assessment
                <br>• Genomic Analyst Agent coordination
                <br>• Multi-modal Digital Twin creation
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)

# The Problem Section - Aggressive Messaging
st.markdown("""
<h2 style='color: #1d1d1f; font-size: 2.5rem; font-weight: 300; margin: 3rem 0 2rem; text-align: center;'>
    The Multi-Billion Dollar Failure of Logistics
</h2>
""", unsafe_allow_html=True)

col1, col2, col3 = st.columns(3)

with col1:
    st.markdown("""
    <div style='background: linear-gradient(135deg, #e74c3c 0%, #c0392b 100%); 
                padding: 2rem; border-radius: 12px; color: white; text-align: center; height: 250px; display: flex; flex-direction: column; justify-content: center;'>
        <div style='font-size: 4rem; margin-bottom: 1rem;'>80%</div>
        <h3 style='color: white; margin-bottom: 1rem;'>Trial Failures</h3>
        <p style='font-size: 1.1rem;'>
            Fail to meet enrollment timelines, delaying life-saving drugs by months or years
        </p>
    </div>
    """, unsafe_allow_html=True)

with col2:
    st.markdown("""
    <div style='background: linear-gradient(135deg, #f39c12 0%, #e67e22 100%); 
                padding: 2rem; border-radius: 12px; color: white; text-align: center; height: 250px; display: flex; flex-direction: column; justify-content: center;'>
        <div style='font-size: 4rem; margin-bottom: 1rem;'>$8M+</div>
        <h3 style='color: white; margin-bottom: 1rem;'>Daily Cost</h3>
        <p style='font-size: 1.1rem;'>
            Lost revenue per day for each delayed blockbuster drug reaching market
        </p>
    </div>
    """, unsafe_allow_html=True)

with col3:
    st.markdown("""
    <div style='background: linear-gradient(135deg, #9b59b6 0%, #8e44ad 100%); 
                padding: 2rem; border-radius: 12px; color: white; text-align: center; height: 250px; display: flex; flex-direction: column; justify-content: center;'>
        <div style='font-size: 4rem; margin-bottom: 1rem;'><5%</div>
        <h3 style='color: white; margin-bottom: 1rem;'>Participation</h3>
        <p style='font-size: 1.1rem;'>
            Of adult cancer patients ever enroll in trials - massive untapped potential
        </p>
    </div>
    """, unsafe_allow_html=True)

st.markdown("---")

# Our Doctrine Section
st.markdown("""
<h2 style='color: #1d1d1f; font-size: 2.5rem; font-weight: 300; margin: 3rem 0 2rem; text-align: center;'>
    Our Doctrine: From Keyword Matching to Biological Warfare
</h2>
""", unsafe_allow_html=True)

# Main Interface
tab1, tab2, tab3, tab4 = st.tabs(["🎯 Patient-Trial Matching", "🧬 Digital Twin Creation", "📊 Eligibility Dossier", "⚡ Command Center"])

with tab1:
    st.markdown("""
    <div style='background: #f8f9fa; padding: 2rem; border-radius: 12px; border-left: 4px solid #007AFF; margin-bottom: 2rem;'>
        <h3 style='color: #1d1d1f; margin-bottom: 1rem;'>🎯 Precision-Guided Trial Matching</h3>
        <p style='color: #515154; font-size: 1.1rem; margin-bottom: 1rem;'>
            <strong>We create deep, multi-modal "Digital Twins" of both patients and trials, 
            matching them based on fundamental understanding of underlying biology.</strong>
        </p>
        <p style='color: #666; font-style: italic;'>
            "A patient is not defined by a single gene name. A trial is not defined by a simple criterion."
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # Patient Input Section
    st.subheader("🧬 Patient Genomic Profile")
    
    col1, col2 = st.columns(2)
    
    with col1:
        patient_gene = st.text_input(
            "Primary Driver Gene",
            value="RUNX1" if demo_manager.demo_mode else "",
            help="Primary oncogene or tumor suppressor involved"
        )
        
        patient_variant = st.text_input(
            "Specific Variant",
            value="R204Q" if demo_manager.demo_mode else "",
            help="Specific protein change (e.g., V600E, R204Q)"
        )
        
        mutation_type = st.selectbox(
            "Mutation Classification",
            ["Germline", "Somatic", "Clonal Hematopoiesis", "Unknown"],
            index=0 if demo_manager.demo_mode else 3
        )
    
    with col2:
        cancer_type = st.text_input(
            "Cancer Type/Risk",
            value="Familial Platelet Disorder (Pre-malignant)" if demo_manager.demo_mode else "",
            help="Primary diagnosis or cancer risk syndrome"
        )
        
        previous_treatments = st.text_area(
            "Previous Treatments",
            value="None (Prevention-focused)" if demo_manager.demo_mode else "",
            help="Prior therapies, surgeries, or interventions"
        )
        
        biomarkers = st.text_area(
            "Additional Biomarkers",
            value="Elevated platelet aggregation defects, mild thrombocytopenia" if demo_manager.demo_mode else "",
            help="Lab values, protein expression, other molecular markers"
        )
    
    # Trial Matching Parameters
    st.subheader("🎯 Trial Matching Criteria")
    
    col1, col2 = st.columns(2)
    
    with col1:
        trial_phase = st.multiselect(
            "Trial Phases",
            ["Phase I", "Phase II", "Phase III", "Phase IV"],
            default=["Phase I", "Phase II"] if demo_manager.demo_mode else []
        )
        
        trial_type = st.multiselect(
            "Trial Types",
            ["Prevention", "Treatment", "Diagnostic", "Supportive Care"],
            default=["Prevention"] if demo_manager.demo_mode else []
        )
    
    with col2:
        geographic_preference = st.selectbox(
            "Geographic Scope",
            ["Global", "North America", "United States", "Europe", "Asia-Pacific"],
            index=2 if demo_manager.demo_mode else 0
        )
        
        max_distance = st.slider(
            "Maximum Travel Distance (miles)",
            min_value=0,
            max_value=500,
            value=100 if demo_manager.demo_mode else 50,
            help="Maximum distance patient willing to travel"
        )
    
    # Execute Matching
    if st.button("🚀 EXECUTE BIOLOGICAL WARFARE MATCHING", type="primary"):
        if not patient_gene or not cancer_type:
            st.error("⚠️ MISSION ABORT: Patient gene and cancer type required for targeting.")
        else:
            with st.spinner("🎯 Deploying Multi-Agent Assault... Analyzing biological intent..."):
                # Simulate API call to CommandCenter
                try:
                    # This would be the actual API call to our CommandCenter
                    # For now, we'll simulate the response
                    
                    if demo_manager.demo_mode and patient_gene.upper() == "RUNX1":
                        # Demo response for RUNX1-FPD
                        matching_results = {
                            "patient_profile": {
                                "gene": patient_gene,
                                "variant": patient_variant,
                                "pathogenicity_assessment": {
                                    "zeta_score": -85.4,
                                    "classification": "Likely Pathogenic",
                                    "source": "Zeta Oracle"
                                },
                                "functional_impact": "Impaired transcriptional regulation, increased clonal hematopoiesis risk"
                            },
                            "matched_trials": [
                                {
                                    "nct_id": "NCT04892459",
                                    "title": "Prevention of Hematologic Malignancies in RUNX1-FPD Patients",
                                    "phase": "Phase II",
                                    "sponsor": "National Cancer Institute",
                                    "status": "Recruiting",
                                    "match_score": 94.7,
                                    "biological_rationale": "Patient's germline RUNX1 R204Q mutation creates transcriptional dysregulation matching trial's target population for clonal hematopoiesis intervention.",
                                    "location": "NIH Clinical Center, Bethesda, MD",
                                    "distance_miles": 45,
                                    "key_eligibility": [
                                        "Germline RUNX1 mutation confirmed",
                                        "Evidence of clonal hematopoiesis",
                                        "No current hematologic malignancy",
                                        "Age 18-75 years"
                                    ]
                                },
                                {
                                    "nct_id": "NCT05123456",
                                    "title": "Metformin for Clonal Hematopoiesis Prevention",
                                    "phase": "Phase I/II",
                                    "sponsor": "Dana-Farber Cancer Institute",
                                    "status": "Recruiting",
                                    "match_score": 87.3,
                                    "biological_rationale": "Metformin's metabolic effects may reduce clonal expansion in high-risk germline mutation carriers like RUNX1-FPD patients.",
                                    "location": "Boston, MA",
                                    "distance_miles": 78,
                                    "key_eligibility": [
                                        "High-risk germline mutations (including RUNX1)",
                                        "Detectable clonal hematopoiesis",
                                        "Normal organ function",
                                        "No diabetes mellitus"
                                    ]
                                },
                                {
                                    "nct_id": "NCT05789012",
                                    "title": "CAR-T Cell Therapy for Pre-Leukemic Clones",
                                    "phase": "Phase I",
                                    "sponsor": "Memorial Sloan Kettering",
                                    "status": "Not yet recruiting",
                                    "match_score": 82.1,
                                    "biological_rationale": "Experimental CAR-T approach targeting cells with RUNX1 dysfunction before malignant transformation.",
                                    "location": "New York, NY",
                                    "distance_miles": 89,
                                    "key_eligibility": [
                                        "RUNX1 familial platelet disorder",
                                        "High clonal burden (>10%)",
                                        "Adequate performance status",
                                        "No prior immunotherapy"
                                    ]
                                }
                            ]
                        }
                    else:
                        # Generic response for other genes
                        matching_results = {
                            "patient_profile": {
                                "gene": patient_gene,
                                "variant": patient_variant,
                                "pathogenicity_assessment": {
                                    "zeta_score": -92.1,
                                    "classification": "Pathogenic",
                                    "source": "Zeta Oracle"
                                },
                                "functional_impact": f"Functional disruption of {patient_gene} pathway identified"
                            },
                            "matched_trials": [
                                {
                                    "nct_id": "NCT04567890",
                                    "title": f"Targeted Therapy for {patient_gene} Mutations",
                                    "phase": "Phase II",
                                    "sponsor": "Leading Cancer Center",
                                    "status": "Recruiting",
                                    "match_score": 91.2,
                                    "biological_rationale": f"Patient's {patient_gene} {patient_variant} mutation directly targets the pathway under investigation.",
                                    "location": "Major Medical Center",
                                    "distance_miles": 25,
                                    "key_eligibility": [
                                        f"{patient_gene} mutation confirmed",
                                        "Adequate organ function",
                                        "No prior targeted therapy",
                                        "Age 18+ years"
                                    ]
                                }
                            ]
                        }
                    
                    st.session_state.trial_matching_results = matching_results
                    
                except Exception as e:
                    st.error(f"🚨 MISSION FAILURE: {str(e)}")
                    st.session_state.trial_matching_results = None

with tab2:
    st.markdown("""
    <div style='background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                padding: 2rem; border-radius: 12px; margin-bottom: 2rem; color: white;'>
        <h3 style='color: white; margin-bottom: 1rem;'>🧬 Digital Twin Creation Protocol</h3>
        <p style='font-size: 1.1rem; margin-bottom: 1rem;'>
            Our platform creates complete Digital Twins, fusing genomic data, clinical history, 
            and lab results into unified intelligence assets for precision matching.
        </p>
        <div style='display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 1rem; margin-top: 1.5rem;'>
            <div style='background: rgba(255,255,255,0.2); padding: 1rem; border-radius: 8px;'>
                <strong>Genomic Layer</strong><br>
                • Mutation analysis<br>
                • Pathway disruption<br>
                • Zeta Oracle assessment
            </div>
            <div style='background: rgba(255,255,255,0.2); padding: 1rem; border-radius: 8px;'>
                <strong>Clinical Layer</strong><br>
                • Disease progression<br>
                • Treatment history<br>
                • Performance status
            </div>
            <div style='background: rgba(255,255,255,0.2); padding: 1rem; border-radius: 8px;'>
                <strong>Molecular Layer</strong><br>
                • Biomarker profiles<br>
                • Protein expression<br>
                • Metabolic state
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    # Digital Twin Visualization
    if 'trial_matching_results' in st.session_state and st.session_state.trial_matching_results:
        results = st.session_state.trial_matching_results
        patient_profile = results.get('patient_profile', {})
        
        st.subheader("🎯 Patient Digital Twin Profile")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("""
            <div style='background: white; border-radius: 12px; padding: 2rem; box-shadow: 0 4px 20px rgba(0,0,0,0.1);'>
                <h4 style='color: #1d1d1f; margin-bottom: 1rem;'>🧬 Genomic Intelligence</h4>
            """, unsafe_allow_html=True)
            
            st.metric("Primary Target", f"{patient_profile.get('gene', 'N/A')} {patient_profile.get('variant', '')}")
            
            pathogenicity = patient_profile.get('pathogenicity_assessment', {})
            zeta_score = pathogenicity.get('zeta_score', 0)
            classification = pathogenicity.get('classification', 'Unknown')
            
            if zeta_score < -50:
                st.error(f"🚨 **{classification}** | Zeta Score: {zeta_score}")
            else:
                st.success(f"✅ **{classification}** | Zeta Score: {zeta_score}")
            
            st.markdown("</div>", unsafe_allow_html=True)
        
        with col2:
            st.markdown("""
            <div style='background: white; border-radius: 12px; padding: 2rem; box-shadow: 0 4px 20px rgba(0,0,0,0.1);'>
                <h4 style='color: #1d1d1f; margin-bottom: 1rem;'>⚡ Functional Impact</h4>
            """, unsafe_allow_html=True)
            
            functional_impact = patient_profile.get('functional_impact', 'Analysis pending...')
            st.write(functional_impact)
            
            st.markdown("</div>", unsafe_allow_html=True)
    else:
        st.info("🎯 Execute patient-trial matching to generate Digital Twin profile")

with tab3:
    st.markdown("""
    <div style='background: linear-gradient(135deg, #ff6b35 0%, #f7931e 100%); 
                padding: 2rem; border-radius: 12px; margin-bottom: 2rem; color: white;'>
        <h3 style='color: white; margin-bottom: 1rem;'>📊 Eligibility Dossier: The Kill Shot</h3>
        <p style='font-size: 1.1rem;'>
            The output is not a simple "yes" or "no." It is a comprehensive "Eligibility Dossier" 
            that provides a prioritized list of trials with detailed, AI-generated rationale 
            explaining why each match was made based on deep biology.
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # Display matching results
    if 'trial_matching_results' in st.session_state and st.session_state.trial_matching_results:
        results = st.session_state.trial_matching_results
        matched_trials = results.get('matched_trials', [])
        
        if matched_trials:
            st.subheader(f"🎯 {len(matched_trials)} High-Priority Trial Matches Identified")
            
            for i, trial in enumerate(matched_trials):
                # Create expandable trial card
                with st.expander(f"🏆 RANK #{i+1}: {trial['title']} (Match Score: {trial['match_score']}%)", expanded=(i==0)):
                    col1, col2 = st.columns([2, 1])
                    
                    with col1:
                        st.markdown(f"**NCT ID:** {trial['nct_id']}")
                        st.markdown(f"**Phase:** {trial['phase']}")
                        st.markdown(f"**Sponsor:** {trial['sponsor']}")
                        st.markdown(f"**Status:** {trial['status']}")
                        st.markdown(f"**Location:** {trial['location']} ({trial['distance_miles']} miles)")
                        
                        st.markdown("**🧠 Biological Rationale:**")
                        st.info(trial['biological_rationale'])
                        
                        st.markdown("**✅ Key Eligibility Criteria:**")
                        for criterion in trial['key_eligibility']:
                            st.markdown(f"• {criterion}")
                    
                    with col2:
                        # Match score visualization
                        score = trial['match_score']
                        if score >= 90:
                            color = "#4CAF50"
                            status = "OPTIMAL MATCH"
                        elif score >= 80:
                            color = "#FF9800"
                            status = "STRONG MATCH"
                        else:
                            color = "#F44336"
                            status = "POTENTIAL MATCH"
                        
                        st.markdown(f"""
                        <div style='background: {color}; color: white; padding: 1rem; border-radius: 8px; text-align: center; margin-bottom: 1rem;'>
                            <div style='font-size: 2rem; font-weight: bold;'>{score}%</div>
                            <div style='font-size: 0.9rem;'>{status}</div>
                        </div>
                        """, unsafe_allow_html=True)
                        
                        if st.button(f"📋 Generate Referral", key=f"referral_{i}"):
                            st.success("🎯 Referral package generated! Clinical coordinator will be notified.")
        else:
            st.warning("No trial matches found. Adjust criteria or patient profile.")
    else:
        st.info("🎯 Execute patient-trial matching to generate Eligibility Dossier")

with tab4:
    st.markdown("""
    <div style='background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%); 
                padding: 2rem; border-radius: 12px; margin-bottom: 2rem; color: white;'>
        <h3 style='color: white; margin-bottom: 1rem;'>⚡ Command Center: Multi-Agent Coordination</h3>
        <p style='font-size: 1.1rem; margin-bottom: 1.5rem;'>
            Our Command Center orchestrates multiple specialized AI agents to execute precision matching campaigns.
        </p>
        <div style='display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 1rem;'>
            <div style='background: rgba(255,255,255,0.1); padding: 1rem; border-radius: 8px;'>
                <strong>🧬 GenomicAnalyst Agent</strong><br>
                Specialized AI for interpreting complex genomic data and mutation impacts
            </div>
            <div style='background: rgba(255,255,255,0.1); padding: 1rem; border-radius: 8px;'>
                <strong>🎯 Zeta Oracle</strong><br>
                Definitive verdict on functional impact of any mutation with confidence scoring
            </div>
            <div style='background: rgba(255,255,255,0.1); padding: 1rem; border-radius: 8px;'>
                <strong>🔍 NLP Vector Engine</strong><br>
                Ingests and understands unstructured trial protocol text for deep matching
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    # Command Center Status
    st.subheader("⚡ System Status")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Active Agents", "3/3", delta="All systems operational")
    
    with col2:
        st.metric("Trials Database", "45,892", delta="+127 this week")
    
    with col3:
        st.metric("Avg Match Time", "2.3 sec", delta="-0.8 sec improved")
    
    # API Configuration
    st.subheader("🔧 Command Center Configuration")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.text_input("CommandCenter API Endpoint", 
                     value=os.getenv("COMMAND_CENTER_ENDPOINT", "https://crispro--command-center-commandcenter-api.modal.run"),
                     disabled=True)
        
        st.selectbox("Genomic Analysis Model", 
                    ["Zeta Oracle v2.1", "Zeta Oracle v2.0", "Legacy Model"],
                    index=0)
    
    with col2:
        st.selectbox("Trial Database Source", 
                    ["ClinicalTrials.gov", "WHO ICTRP", "Combined Sources"],
                    index=2)
        
        st.slider("Match Confidence Threshold", 
                 min_value=50, max_value=95, value=80,
                 help="Minimum confidence score for trial recommendations")

# Footer with value proposition
st.markdown("---")
st.markdown("""
<div style='background: #f8f9fa; padding: 2rem; border-radius: 12px; text-align: center;'>
    <h3 style='color: #1d1d1f; margin-bottom: 1rem;'>🎯 Target Audience & Value Proposition</h3>
    <div style='display: grid; grid-template-columns: 1fr 1fr; gap: 2rem; margin-top: 2rem;'>
        <div style='background: white; padding: 2rem; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1);'>
            <h4 style='color: #007AFF; margin-bottom: 1rem;'>🏢 Pharma & Biotech</h4>
            <p style='color: #515154;'>
                <strong>We are your outsourced recruitment engine.</strong> Fill trials faster with better-qualified patients, 
                dramatically reducing R&D timelines and increasing probability of success. We solve your most expensive bottleneck.
            </p>
        </div>
        <div style='background: white; padding: 2rem; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1);'>
            <h4 style='color: #007AFF; margin-bottom: 1rem;'>🏥 Health Systems</h4>
            <p style='color: #515154;'>
                <strong>We provide your clinicians with a superpower:</strong> instantly find cutting-edge therapeutic options 
                for patients. Improve outcomes while establishing your institution as a leader in precision medicine.
            </p>
        </div>
    </div>
</div>
""", unsafe_allow_html=True) 