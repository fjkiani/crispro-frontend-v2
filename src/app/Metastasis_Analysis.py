import streamlit as st
from tools.metastasis_analyzer import MetastasisAnalyzer
import os
import pandas as pd

st.set_page_config(
    page_title="Metastatic Risk Analysis",
    layout="wide"
)

st.title("🔬 AI-Powered Metastatic Risk Analysis")

st.markdown("""
Welcome to the AI-Powered Metastasis Prevention Platform. 

This tool leverages the Evo2 biological foundation model to analyze a patient's genomic data (from a VCF file) and assess metastatic risk across the 8-step cascade.

**How to use:**
1.  Upload a VCF file containing somatic variants.
2.  Click the "Run Analysis" button.
3.  Review the generated "Metastatic Potential Report".
""")

# File Uploader
uploaded_file = st.file_uploader("Upload your VCF file", type=["vcf"])

if uploaded_file is not None:
    # Save the uploaded file to a temporary location to be used by the analyzer
    temp_file_path = os.path.join("/tmp", uploaded_file.name)
    with open(temp_file_path, "wb") as f:
        f.write(uploaded_file.getbuffer())

    st.success(f"File '{uploaded_file.name}' uploaded successfully.")

    if st.button("Run Analysis", key="run_analysis"):
        with st.spinner("Analyzing patient data... This may take a moment."):
            try:
                analyzer = MetastasisAnalyzer()
                analysis_results = analyzer.analyze_patient_data(temp_file_path)

                st.session_state['analysis_results'] = analysis_results
            except Exception as e:
                st.error(f"An error occurred during analysis: {e}")
                st.stop()

# Display results if they exist in the session state
if 'analysis_results' in st.session_state:
    results = st.session_state['analysis_results']
    st.header("Metastatic Potential Report")

    if results:
        # In a real app, we would calculate a meaningful risk score.
        # For now, we generate a placeholder score.
        st.metric("Overall Metastatic Risk Score", "7.8 / 10.0", "High Risk")
        
        st.subheader("Stage-by-Stage Risk Breakdown")
        
        for step, result in results.items():
            step_name = step.replace('_', ' ').title()
            with st.expander(f"**{step_name}**"):
                if result:
                    if result.get("status") == "placeholder":
                        st.info("Analysis for this step is not yet implemented.")
                    else:
                        # This assumes the result is a JSON object.
                        # In a real app, we would format this nicely.
                        st.json(result)
                else:
                    st.warning("No result was obtained for this step. Please check the API connection and logs.")

    else:
        st.warning("Analysis did not produce any results.") 