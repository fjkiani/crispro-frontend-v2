"""
The primary user interface for designing and executing CRISPR-based
interventions against identified genetic threats.
"""
import streamlit as st
import os
import requests
from dotenv import load_dotenv
import time

# --- Page Configuration ---
st.set_page_config(
    page_title="Interception Design Studio",
    page_icon="🎯",
    layout="wide"
)

# --- Environment and API Setup ---
load_dotenv()
COMMAND_CENTER_URL = os.getenv("COMMAND_CENTER_URL")
if not COMMAND_CENTER_URL:
    st.error("COMMAND_CENTER_URL not set in .env file. Please configure it to connect to the backend.")
    st.stop()

NANOBODY_ENDPOINT = f"{COMMAND_CENTER_URL}/workflow/design_nanobody"
DRUG_SCREEN_ENDPOINT = f"{COMMAND_CENTER_URL}/workflow/run_drug_repurposing"

# --- Session State ---
if 'nanobody_results' not in st.session_state:
    st.session_state.nanobody_results = None
if 'drug_screen_results' not in st.session_state:
    st.session_state.drug_screen_results = None

# --- Main UI ---
st.title("🎯 Interception Design Studio")
st.markdown("""
A validated target is only the beginning. This studio provides a suite of tools to design multiple, distinct therapeutic interventions against that target. We can attack the problem from multiple angles, from gene-level editing to protein-level inhibition.
""")
st.markdown("---")

# Check for a handoff from another page
target_gene = st.session_state.get('target_gene_for_design', "ASXL1")
if 'target_gene_for_design' in st.session_state:
    st.success(f"🎯 **Target Acquired:** Ready to design interventions for **{target_gene}**.", icon="✅")
    # We can delete the key now as its value is stored in target_gene
    del st.session_state.target_gene_for_design


tab1, tab2, tab3 = st.tabs(["**Biologic Inhibitor (Nanobody)**", "**Drug Repurposing Screen**", "**Gene Editing (CRISPR)**"])

# --- Tab 1: Nanobody Design ---
with tab1:
    st.header("Generate Novel Biologic Inhibitor")
    st.markdown("This module uses our **Zeta Oracle (Generative)** to design a novel nanobody—a small, highly specific antibody fragment—to bind to and inhibit the protein product of our target gene.")

    with st.form("nanobody_form"):
        protein_name = st.text_input("Target Protein Name", value=target_gene)
        protein_domain = st.text_input("Target Protein Domain", value="Active Site", help="e.g., 'Active Site', 'Binding Domain', 'SH2 Domain'")
        submitted = st.form_submit_button("Design Nanobody")

    if submitted:
        with st.spinner("Invoking Zeta Oracle... Forging novel nanobody sequence..."):
            payload = {"target_protein_name": protein_name, "target_domain": protein_domain}
            try:
                response = requests.post(NANOBODY_ENDPOINT, json=payload, timeout=300)
                response.raise_for_status()
                st.session_state.nanobody_results = response.json()
                st.success("Nanobody design complete.")
            except Exception as e:
                st.error(f"An error occurred: {e}")

    if st.session_state.nanobody_results:
        st.subheader("Forged Nanobody Blueprint")
        results = st.session_state.nanobody_results
        st.text(f"Designed for: {results.get('design_parameters', {}).get('target_protein')}")

        dna_sequence = results.get('generated_sequence', {}).get('dna_sequence', '')
        st.code(dna_sequence, language="fasta")
        st.info(results.get('generated_sequence', {}).get('comment'))


# --- Tab 2: Drug Repurposing ---
with tab2:
    st.header("Screen for Repurposable Drugs")
    st.markdown("This module simulates a high-throughput virtual screen, using the predicted structure of the target protein to find existing, FDA-approved drugs that could be repurposed to inhibit it.")

    with st.form("drug_screen_form"):
        protein_structure = st.text_area("Target Protein Structure (PDB format)", value=f"PLACEHOLDER_PDB_FOR_{target_gene}", help="In a real scenario, this would be the PDB data from an AlphaFold prediction.")
        drug_library = st.selectbox("Drug Library to Screen", ["DrugBank_approved", "ZINC", "ChEMBL"])
        submitted = st.form_submit_button("Run Repurposing Screen")

    if submitted:
        with st.spinner("Invoking Virtual Screening service..."):
            payload = {"mutant_protein_structure": protein_structure, "drug_library": drug_library}
            try:
                response = requests.post(DRUG_SCREEN_ENDPOINT, json=payload, timeout=300)
                response.raise_for_status()
                st.session_state.drug_screen_results = response.json()
                st.success("Virtual screen complete.")
            except Exception as e:
                st.error(f"An error occurred: {e}")

    if st.session_state.drug_screen_results:
        st.subheader("Top Drug Hits")
        results = st.session_state.drug_screen_results
        st.dataframe(results.get('top_hits'), use_container_width=True)


# --- Tab 3: Gene Editing ---
with tab3:
    st.header("Design Gene Knockout Therapy (CRISPR)")
    st.markdown("For a direct, permanent solution, we can design a gene editing therapy to knock out the target gene at the DNA level. This leverages our high-precision `CHOPCHOP` guide designer.")

    st.warning("This capability has been integrated into a seamless workflow.")
    st.markdown("To design a gene knockout therapy, please use the following workflow:")
    st.markdown("1. Start at the **Future Threat Simulation** page to identify a threat.")
    st.markdown("2. Proceed to the **Seed & Soil Analysis** page to validate the target in a specific tissue.")
    st.markdown("3. The final step of the Seed & Soil analysis will automatically hand you off to the **CHOPCHOP Guide Designer** with the correct target pre-loaded.")

    if st.button("Go to CHOPCHOP Guide Designer"):
        st.switch_page("pages/5_✂️_CHOPCHOP_Guide_Design.py") 