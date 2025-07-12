import streamlit as st
import pandas as pd
import json
import sys
import os

# Add the project root to the Python path to allow importing from 'tools'
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
from tools.command_center_client import get_all_patients, get_latest_battle_plan_for_patient

st.set_page_config(
    page_title="AI General Command Deck",
    page_icon="⚔️",
    layout="wide",
)

st.title("⚔️ AI General Command Deck")
st.markdown("A High-Fidelity, Dynamic Battlefield Simulation of Patient Disease")

# --- DATA LOADING ---
@st.cache_data
def load_patients():
    """Load patient list from the backend."""
    patients = get_all_patients()
    if patients:
        return [p['patient_identifier'] for p in patients]
    return []

patient_list = load_patients()

if not patient_list:
    st.warning("Could not retrieve patient list from the Command Center. Displaying mock data.")
    st.info("Ensure the CommandCenter Modal service is deployed and the URL is set in your environment.")
    patient_list = ["Patient_A (TP53)", "Patient_B (BRAF)", "Patient_C (KRAS)"]


# --- UI COMPONENTS ---

# Sidebar for patient selection
st.sidebar.title("Patient Selection")
selected_patient_id = st.sidebar.selectbox("Select a Patient Case", patient_list)

# Main content area
st.header(f"Digital Twin: {selected_patient_id}")

# Load the battle plan for the selected patient
if selected_patient_id:
    battle_plan = get_latest_battle_plan_for_patient(selected_patient_id)

    if battle_plan and "error" not in battle_plan:
        # Confidence Score
        st.subheader("Confidence Score")
        confidence = battle_plan.get("perplexity_score", 0.0) # Using perplexity as confidence
        st.metric(label="AI General Confidence (1 - Perplexity)", value=f"{1 - confidence:.2%}", delta=None)
        st.progress(1 - confidence)

        # Main Intelligence Tabs
        tab1, tab2, tab3 = st.tabs(["🧬 Genomic Integrity", "🎯 Ensemble Assessment", "🗺️ Metastasis Forecast"])

        with tab1:
            st.subheader("Structural Variant Scan")
            sv_data = battle_plan.get("structural_variants", {})
            if not sv_data or not any(sv_data.values()):
                st.success("No significant structural variants detected.")
            else:
                for sv_type, variants in sv_data.items():
                    if variants:
                        st.write(f"**{sv_type.replace('_', ' ').title()}:**")
                        for variant in variants:
                            st.code(json.dumps(variant, indent=2), language="json")

        with tab2:
            st.subheader("Ensemble Variant Assessment")
            ensemble_data = battle_plan.get("ensemble_scores", {})
            if ensemble_data:
                df = pd.DataFrame.from_dict(ensemble_data, orient='index')
                st.dataframe(df, use_container_width=True)
            else:
                st.info("No ensemble assessment data available for this plan.")

        with tab3:
            st.subheader("8-Step Metastasis Forecast")
            # This part is a placeholder as the battle_plan model doesn't have this field yet.
            # We will use mock data for now.
            forecast_data = {
                "Primary Tumor Growth": {"score": 0.8, "evidence": "High proliferation rate"},
                "Local Invasion": {"score": 0.6, "evidence": "EMT markers detected"},
            }
            st.info("Metastasis forecast data is currently mocked.")
            for step, details in forecast_data.items():
                st.write(f"**{step}**")
                st.progress(details['score'])
                with st.expander("Evidence"):
                    st.info(details['evidence'])
    else:
        st.error(f"Could not load battle plan for {selected_patient_id}. Error: {battle_plan.get('error', 'Unknown error')}")
else:
    st.info("Select a patient from the sidebar to view their Digital Twin.") 