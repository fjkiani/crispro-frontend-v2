import streamlit as st
import requests
import os
import sys
import json
import time
import pandas as pd
import streamlit_gosling as st_gos
from typing import Dict

# Add the project root to the Python path to allow for tool imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from tools.ui_enhancements import create_enhanced_sidebar
from tools.runx1_integration_plan import run_v4_intervention_design

# --- Page Config ---
st.set_page_config(
    page_title="Patient Digital Twin v4",
    page_icon="🧬",
    layout="wide",
)

create_enhanced_sidebar("Patient Digital Twin", {"gene_symbol": "RUNX1", "focus_area": "digital_twin_modeling"})

def main_digital_twin_page_v4():
    st.title("Patient Digital Twin v4: RUNX1-FPD")
    st.markdown("### V4 Kill Chain Integration")
    st.info("This interface is now directly connected to the AI General Command Center's V4 workflow.")

    # Example hardcoded patient data for the demo
    patient_id = "Patient-RUNX1-001"
    germline_variant_data = {
        "id": "var-runx1-1",
        "gene": "RUNX1",
        "protein_change": "p.Arg135fs",
        "sequence_for_perplexity": "CCTGCCTGGGCTCCCCAGCCCCAGCAGCCCTGCAGCCCAGCTCTGGCCTGCCGGAGGCCACGCGC"
    }

    st.markdown("---")
    st.markdown("#### Patient Profile")
    col1, col2 = st.columns(2)
    with col1:
        st.write(f"**Patient ID:**")
        st.code(patient_id, language='text')
    with col2:
        st.write(f"**Germline Variant:**")
        st.code(f"{germline_variant_data['gene']} {germline_variant_data['protein_change']}", language='text')

    st.markdown("---")

    if st.button("🚀 Design V4 Precision Intervention", type="primary", use_container_width=True):
        st.session_state.v4_intervention_results = run_v4_intervention_design(germline_variant_data, patient_id)
        st.rerun()

    if 'v4_intervention_results' in st.session_state:
        results = st.session_state.v4_intervention_results
        st.markdown("### V4 Intervention Design Results")
        
        if not results or "error" in results:
            st.error(f"Intervention Design Failed: {results.get('error', 'Unknown error')}")
        else:
            guides = results.get("guides", [])
            if not guides:
                st.warning("The campaign succeeded but returned no guide candidates.")
            else:
                st.success(f"V4 Kill Chain Complete! {len(guides)} interventions designed and ranked.")
                df = pd.DataFrame(guides)
                df_display = df[['guide_sequence', 'assassin_score']].rename(columns={
                    "guide_sequence": "Guide Sequence",
                    "assassin_score": "Assassin Score"
                }).sort_values(by="Assassin Score", ascending=False)
                
                st.dataframe(df_display, use_container_width=True)
                st.info("Higher 'Assassin Score' indicates a better-predicted balance of efficacy and safety.")
                
                with st.expander("🔬 View Detailed Guide Dossiers"):
                    st.json(guides)

if __name__ == "__main__":
    main_digital_twin_page_v4() 