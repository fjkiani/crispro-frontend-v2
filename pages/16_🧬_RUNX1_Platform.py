import streamlit as st
import sys
import os

# Add the tools directory to the path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))

try:
    from runx1_integration_plan import RUNX1IntegrationPlan
    from runx1_data_loader import RUNX1DataLoader
    from runx1_progression_modeler import RUNX1ProgressionModeler
    from runx1_genomic_browser import RUNX1GenomicBrowser
    from runx1_demo_scenarios import RUNX1DemoScenarios
    RUNX1_AVAILABLE = True
except ImportError as e:
    st.error(f"RUNX1 components not available: {e}")
    RUNX1_AVAILABLE = False

st.set_page_config(
    page_title="RUNX1-FPD Platform",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 RUNX1-FPD Precision Medicine Platform")
st.markdown("### AI-Powered Clinical Decision Support for RUNX1 Familial Platelet Disorder")

if not RUNX1_AVAILABLE:
    st.error("RUNX1 platform components are not available. Please check the installation.")
    st.stop()

# Initialize session state
if 'runx1_integration' not in st.session_state:
    st.session_state.runx1_integration = RUNX1IntegrationPlan()
if 'runx1_data_loader' not in st.session_state:
    st.session_state.runx1_data_loader = RUNX1DataLoader()
if 'runx1_progression_modeler' not in st.session_state:
    st.session_state.runx1_progression_modeler = RUNX1ProgressionModeler()
if 'runx1_genomic_browser' not in st.session_state:
    st.session_state.runx1_genomic_browser = RUNX1GenomicBrowser()
if 'runx1_demo_scenarios' not in st.session_state:
    st.session_state.runx1_demo_scenarios = RUNX1DemoScenarios()

# Sidebar for navigation
st.sidebar.title("🧬 RUNX1 Platform")
mode = st.sidebar.selectbox(
    "Select Mode",
    ["Patient Analysis", "Demo Scenarios", "Genomic Browser", "Integration Testing"]
)

if mode == "Patient Analysis":
    st.header("🔬 Patient Analysis Workflow")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Patient Information")
        patient_id = st.text_input("Patient ID", "RUNX1_001")
        patient_age = st.number_input("Age", min_value=0, max_value=120, value=45)
        family_history = st.selectbox("Family History", ["Yes", "No", "Unknown"])
        
        st.subheader("Genetic Variants")
        germline_variant = st.text_input("Germline Variant", "RUNX1_R204Q")
        somatic_variants = st.text_area("Somatic Variants (one per line)", "ASXL1_G646fs")
        
        if st.button("🚀 Run Complete Analysis"):
            with st.spinner("Running comprehensive RUNX1 analysis..."):
                try:
                    # Parse somatic variants
                    somatic_list = [v.strip() for v in somatic_variants.split('\n') if v.strip()]
                    
                    # Run progression modeling
                    progression_result = st.session_state.runx1_progression_modeler.model_clonal_evolution(
                        germline_variant, somatic_list, patient_age
                    )
                    
                    # Run integration analysis
                    integration_result = st.session_state.runx1_integration.run_comprehensive_analysis(
                        patient_id, germline_variant, somatic_list, patient_age
                    )
                    
                    st.session_state.analysis_results = {
                        'progression': progression_result,
                        'integration': integration_result
                    }
                    
                    st.success("✅ Analysis complete!")
                    
                except Exception as e:
                    st.error(f"Analysis failed: {str(e)}")
    
    with col2:
        st.subheader("Analysis Results")
        if 'analysis_results' in st.session_state:
            results = st.session_state.analysis_results
            
            # Display progression results
            if 'progression' in results:
                prog = results['progression']
                st.metric("Transformation Risk", f"{prog.get('transformation_probability', 0):.1%}")
                st.metric("Risk Category", prog.get('risk_category', 'Unknown'))
                
                # Display timeline
                if 'timeline' in prog:
                    st.subheader("Disease Timeline")
                    for event in prog['timeline']:
                        st.write(f"**{event['stage']}**: {event['description']}")
            
            # Display integration results
            if 'integration' in results:
                integ = results['integration']
                st.subheader("AI Analysis Summary")
                st.write(integ.get('summary', 'No summary available'))
                
                # Display interventions
                if 'interventions' in integ:
                    st.subheader("Recommended Interventions")
                    for intervention in integ['interventions']:
                        st.write(f"• {intervention}")

elif mode == "Demo Scenarios":
    st.header("🎭 Demo Scenarios")
    
    scenario_names = st.session_state.runx1_demo_scenarios.get_scenario_names()
    selected_scenario = st.selectbox("Select Demo Scenario", scenario_names)
    
    if st.button("🎬 Load Demo Scenario"):
        with st.spinner("Loading demo scenario..."):
            try:
                scenario = st.session_state.runx1_demo_scenarios.get_scenario(selected_scenario)
                st.session_state.current_scenario = scenario
                st.success(f"✅ Loaded scenario: {selected_scenario}")
            except Exception as e:
                st.error(f"Failed to load scenario: {str(e)}")
    
    if 'current_scenario' in st.session_state:
        scenario = st.session_state.current_scenario
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("Patient Profile")
            st.write(f"**Name**: {scenario.get('patient_name', 'Unknown')}")
            st.write(f"**Age**: {scenario.get('age', 'Unknown')}")
            st.write(f"**Clinical Presentation**: {scenario.get('clinical_presentation', 'Unknown')}")
            
            if 'genetic_profile' in scenario:
                st.subheader("Genetic Profile")
                genetic = scenario['genetic_profile']
                st.write(f"**Germline**: {genetic.get('germline_variant', 'Unknown')}")
                st.write(f"**Somatic**: {', '.join(genetic.get('somatic_variants', []))}")
        
        with col2:
            st.subheader("Clinical Analysis")
            if 'analysis_results' in scenario:
                analysis = scenario['analysis_results']
                st.metric("Transformation Risk", f"{analysis.get('transformation_probability', 0):.1%}")
                st.metric("Risk Category", analysis.get('risk_category', 'Unknown'))
                
                if 'clinical_recommendations' in analysis:
                    st.subheader("Clinical Recommendations")
                    for rec in analysis['clinical_recommendations']:
                        st.write(f"• {rec}")

elif mode == "Genomic Browser":
    st.header("🔍 RUNX1 Genomic Browser")
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.subheader("Navigation")
        position = st.number_input("Genomic Position", min_value=1, max_value=261501, value=100000)
        window_size = st.selectbox("Window Size", [1000, 5000, 10000, 50000])
        
        if st.button("🔍 Navigate"):
            with st.spinner("Loading genomic region..."):
                try:
                    region_data = st.session_state.runx1_genomic_browser.get_region_data(
                        position, window_size
                    )
                    st.session_state.current_region = region_data
                    st.success("✅ Region loaded!")
                except Exception as e:
                    st.error(f"Failed to load region: {str(e)}")
    
    with col2:
        st.subheader("Genomic View")
        if 'current_region' in st.session_state:
            region = st.session_state.current_region
            
            # Display sequence info
            st.write(f"**Region**: {region.get('start', 0)}-{region.get('end', 0)}")
            st.write(f"**Length**: {region.get('length', 0)} bp")
            
            # Display features
            if 'features' in region:
                st.subheader("Genomic Features")
                for feature in region['features']:
                    st.write(f"• **{feature['type']}**: {feature['name']} ({feature['start']}-{feature['end']})")
            
            # Display variants
            if 'variants' in region:
                st.subheader("Variants in Region")
                for variant in region['variants']:
                    st.write(f"• **{variant['type']}**: {variant['position']} {variant['change']}")

elif mode == "Integration Testing":
    st.header("🧪 Integration Testing")
    
    if st.button("🚀 Run Integration Tests"):
        with st.spinner("Running integration tests..."):
            try:
                from runx1_integration_test import RUNX1IntegrationTest
                test_runner = RUNX1IntegrationTest()
                
                # Run all tests
                test_results = test_runner.run_all_tests()
                
                st.subheader("Test Results")
                for test_name, result in test_results.items():
                    status = "✅ PASSED" if result['success'] else "❌ FAILED"
                    duration = result.get('duration', 0)
                    st.write(f"**{test_name}**: {status} ({duration:.2f}s)")
                    
                    if not result['success']:
                        st.error(f"Error: {result.get('error', 'Unknown error')}")
                    
                    if 'details' in result:
                        with st.expander(f"Details for {test_name}"):
                            st.json(result['details'])
                
                # Summary
                total_tests = len(test_results)
                passed_tests = sum(1 for r in test_results.values() if r['success'])
                st.metric("Test Success Rate", f"{passed_tests}/{total_tests} ({passed_tests/total_tests*100:.1f}%)")
                
            except Exception as e:
                st.error(f"Integration testing failed: {str(e)}")

# Footer
st.markdown("---")
st.markdown("**🧬 RUNX1-FPD Platform** | Precision Medicine for Familial Platelet Disorder | Built with AI 🚀") 