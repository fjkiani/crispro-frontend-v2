import streamlit as st
import pandas as pd
from PIL import Image
import sys
import os
import requests

# MUST be the very first Streamlit command
st.set_page_config(
    page_title="RadOnc Co-Pilot",
    page_icon="☢️",
    layout="wide"
)

# Add project root to path to allow import of cousin modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tools.ui_enhancements import (
    apply_leap_styling,
    create_enhanced_sidebar,
    show_loading_with_context,
    AppleStyleResultsDisplay
)

# Apply consistent styling
apply_leap_styling()

# --- Page Content ---

st.markdown("# ☢️ RadOnc Co-Pilot")
st.markdown("### Predicting Radiotherapy Response in Lung Cancer")

st.markdown("""
This page presents the findings from our deep-dive analysis into the TCGA dataset, exploring how a patient's **TP53 mutation status** correlates with their survival outcome after receiving **radiation therapy**.
""")

# --- Display Key Results ---

st.subheader("Key Finding: TP53 Pathogenicity Predicts Poor Radiotherapy Response")

try:
    # Load the survival plot
    survival_plot_path = "analyses/rad_onc_tp53/results/survival_plot.png"
    image = Image.open(survival_plot_path)
    st.image(image, caption="Kaplan-Meier survival curves for four patient cohorts.", use_column_width=True)

    st.markdown("""
    **Interpretation:** The orange line represents patients with a pathogenic (harmful) TP53 mutation who received radiation. Their dramatically lower survival probability demonstrates that for this patient group, radiation therapy is not only ineffective but potentially detrimental. This is a powerful, data-driven biomarker that could help personalize cancer treatment.
    """)

except FileNotFoundError:
    st.error(f"❌ **Error:** The survival plot was not found. Please ensure the analysis has been run and the file exists at `{survival_plot_path}`.")


# --- Interactive Data Exploration ---

st.subheader("Explore the Full Patient Dataset")

try:
    # Load the full analysis data
    master_data_path = "analyses/rad_onc_tp53/results/pathogenicity_analysis.tsv"
    df = pd.read_csv(master_data_path, sep='\\t')

    st.markdown("Use the table below to search, sort, and filter the full dataset from the study.")
    st.dataframe(df)

except FileNotFoundError:
    st.error(f"❌ **Error:** The analysis data was not found. Please ensure the analysis has been run and the file exists at `{master_data_path}`.")
except Exception as e:
    st.error(f"An error occurred while loading the data: {e}")

# --- Live Assessment Widget ---

st.subheader("Live Assessment: Analyze a New Variant")

st.markdown("""
This demonstrates the co-pilot's ability to act as a live decision support tool. Enter a TP53 variant to get a real-time pathogenicity assessment from our platform.
""")

col1, col2 = st.columns([1, 1])

with col1:
    gene_input = st.text_input("Gene Symbol", value="TP53", key="live_gene")

with col2:
    variant_input = st.text_input("Protein Change", value="p.R248Q", key="live_variant")

if st.button("⚡ Assess Variant in Real-Time"):
    if not gene_input or not variant_input:
        st.warning("Please provide both a gene symbol and a protein change.")
    else:
        with st.spinner(f"🧠 Assessing {gene_input} {variant_input}... This may take a moment."):
            try:
                # Construct the API call
                api_url = "https://crispro--command-center-v2-commandcenter-api.modal.run/workflow/assess_threat"
                payload = {
                    "gene_symbol": gene_input,
                    "protein_change": variant_input
                }

                response = requests.post(api_url, json=payload, timeout=120)

                if response.status_code == 200:
                    result = response.json()
                    # Use the Apple-style display to show the result
                    AppleStyleResultsDisplay.display_threat_assessment(result)
                else:
                    st.error(f"❌ Analysis failed: {response.status_code} - {response.text}")

            except requests.exceptions.RequestException as e:
                st.error(f"❌ Connection error: {str(e)}")
            except Exception as e:
                st.error(f"❌ An unexpected error occurred: {str(e)}") 