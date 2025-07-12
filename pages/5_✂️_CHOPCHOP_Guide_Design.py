import streamlit as st
import os
import requests
from dotenv import load_dotenv
import pandas as pd

# --- Page Configuration ---
st.set_page_config(
    page_title="CHOPCHOP Guide Designer",
    page_icon="✂️",
    layout="wide"
)

# --- Environment and API Setup ---
load_dotenv()
COMMAND_CENTER_URL = os.getenv("COMMAND_CENTER_URL")
if not COMMAND_CENTER_URL:
    st.error("COMMAND_CENTER_URL not set in .env file. Please configure it to connect to the backend.")
    st.stop()
    
# This endpoint is versatile and can be used for simple guide design by providing a genomic context.
DESIGN_ENDPOINT = f"{COMMAND_CENTER_URL}/workflow/design_germline_correction" 
OFF_TARGET_ENDPOINT = f"{COMMAND_CENTER_URL}/workflow/off_target_analysis"
NANOBODY_DESIGN_ENDPOINT = f"{COMMAND_CENTER_URL}/workflow/design_nanobody"

# --- Session State Initialization ---
if 'chopchop_results' not in st.session_state:
    st.session_state.chopchop_results = None
if 'off_target_reports' not in st.session_state:
    st.session_state.off_target_reports = {}
if 'nanobody_design_result' not in st.session_state:
    st.session_state.nanobody_design_result = None

# This is the key for the handoff from the Seed & Soil page.
# We do not initialize it here, we only check for its existence.

# --- API Communication ---
def run_off_target_analysis(guide_sequence):
    """Calls the CommandCenter to run off-target analysis for a given guide."""
    payload = {"guide_sequence": guide_sequence}
    try:
        with st.spinner(f"Running off-target analysis for {guide_sequence[:10]}..."):
            response = requests.post(OFF_TARGET_ENDPOINT, json=payload, timeout=300)
            response.raise_for_status()
        st.success("Off-target analysis complete!")
        return response.json()
    except requests.exceptions.RequestException as e:
        st.error(f"Off-target analysis API Request Failed: {e}")
        return None

def design_nanobody_inhibitor(target_protein_name, target_domain):
    """Calls the CommandCenter to design a nanobody inhibitor."""
    payload = {
        "target_protein_name": target_protein_name,
        "target_domain": target_domain
    }
    try:
        with st.spinner(f"Designing nanobody for {target_protein_name}..."):
            response = requests.post(NANOBODY_DESIGN_ENDPOINT, json=payload, timeout=600) # Longer timeout for generative model
            response.raise_for_status()
        st.success("Nanobody design complete!")
        return response.json()
    except requests.exceptions.RequestException as e:
        st.error(f"Nanobody Design API Request Failed: {e}")
        st.error(f"Please ensure the CommandCenter is deployed and the `COMMAND_CENTER_URL` is correct.")
        return None

# --- Main UI ---
st.title("✂️ CHOPCHOP Guide RNA Designer")
st.markdown("""
This tool designs and validates guide RNAs for gene editing. It serves two primary functions:
1.  **Germline Correction:** Designing therapies to *correct* a faulty gene.
2.  **Gene Knockout:** Designing guides to *disable* a target gene (e.g., a "soil-poisoning" target identified in metastasis modeling).
""")
st.markdown("---")

# Using tabs for different functionalities
guide_designer_tab, nanobody_designer_tab = st.tabs(["Guide RNA Designer", "Nanobody Designer"])

with guide_designer_tab:
    # --- Input Form ---
    # Check if a target gene was passed from another page via session_state
    initial_target = st.session_state.get('target_gene_for_design', "")

    # Display a prominent message if a handoff has occurred.
    if initial_target:
        st.success(f"🎯 **Target Acquired from Seed & Soil Analysis:** A therapy will be designed to knock out the **{initial_target}** gene, poisoning the 'soil' to make it inhospitable for cancer cells.", icon="🎯")

    with st.form("guide_design_form"):
        st.subheader("Design Parameters")
        
        # Use the passed-in gene name as the default, otherwise use a placeholder
        target_gene = st.text_input(
            "Target Gene Name", 
            value=initial_target or "ASXL1", 
            help="Enter the official symbol for the gene you want to target."
        )
        
        # The genomic context is the most critical piece of information.
        # We provide a default that is known to work for the demo.
        target_context = st.text_area(
            "Genomic Context (Coordinates or Sequence)", 
            "chr20:31022420-31022480", # A default valid coordinate for ASXL1
            height=100,
            help="Provide the genomic coordinates (e.g., chr1:1-1000) or a raw DNA sequence for the target region."
        )
        
        submitted = st.form_submit_button("Design & Validate Guides")

    # --- Logic and Display ---
    if submitted:
        if not target_gene or not target_context:
            st.error("Please provide both a target gene and a genomic context.")
        else:
            with st.spinner("Invoking CommandCenter... Designing and validating guides..."):
                # For a knockout, the 'corrected_sequence' is irrelevant but required by the endpoint.
                # We send a clear signal of our intent. The backend handles the logic.
                # The most important part is the target_site_context, which the guides are designed against.
                payload = {
                    "corrected_sequence": "KNOCKOUT_TARGET", 
                    "target_site_context": target_context
                }
                try:
                    response = requests.post(DESIGN_ENDPOINT, json=payload, timeout=300)
                    response.raise_for_status()
                    st.session_state.chopchop_results = response.json()
                    st.success("Guide design and validation complete!")
                    
                    # Critical step: Clear the handoff state variable after it has been used.
                    # This prevents the message from showing up on subsequent visits to the page.
                    if 'target_gene_for_design' in st.session_state:
                        del st.session_state.target_gene_for_design

                except requests.exceptions.HTTPError as http_err:
                    st.error(f"HTTP error occurred: {http_err}")
                    st.error(f"Response Body: {response.text}")
                except requests.exceptions.RequestException as req_err:
                    st.error(f"Request error occurred: {req_err}")

    # --- Results Display ---
    if st.session_state.chopchop_results:
        st.markdown("---")
        st.subheader("Design & Validation Report")
        
        results = st.session_state.chopchop_results
        
        if results.get("error"):
            st.error(f"An error occurred in the backend: {results.get('error')}")
        else:
            col1, col2 = st.columns(2)

            with col1:
                # Display Guide RNA Design Results
                guide_design = results.get("guide_rna_design", {})
                st.markdown("#### ⚔️ Guide RNA Candidates")
                if guide_design.get("error"):
                     st.error(f"Could not design guides: {guide_design.get('error')}")
                else:
                    st.metric("Total Candidates Found", guide_design.get("guides_found", 0))
                    best_guide = guide_design.get("best_guide")
                    if best_guide:
                        st.code(best_guide, language=None)
                    
                    candidates = guide_design.get("candidates", [])
                    if candidates:
                        st.markdown("##### All Candidate Guides:")
                        for i, guide_seq in enumerate(candidates):
                            guide_col, button_col = st.columns([0.7, 0.3])
                            with guide_col:
                                st.markdown(f"`{guide_seq}`")
                            with button_col:
                                if st.button("Validate Safety", key=f"validate_btn_{guide_seq}"):
                                    off_target_report = run_off_target_analysis(guide_seq)
                                    st.session_state.off_target_reports[guide_seq] = off_target_report
                                    st.rerun()
                            
                            # Display off-target results for this specific guide if available
                            if guide_seq in st.session_state.off_target_reports and st.session_state.off_target_reports[guide_seq] is not None:
                                report = st.session_state.off_target_reports[guide_seq]
                                if report.get("error"):
                                    st.error(f"Off-target report error: {report.get('error')}")
                                else:
                                    num_hits = report.get("total_off_targets_found", 0)
                                    if num_hits > 0:
                                        st.warning(f"⚠️ Found {num_hits} potential off-target sites.")
                                        with st.expander(f"View Details for {guide_seq}"):
                                            st.dataframe(pd.DataFrame(report["off_target_hits"]), use_container_width=True)
                                    else:
                                        st.success(f"✅ No significant off-target sites found for {guide_seq}.")
                                    st.markdown("---") # Separator


            with col2:
                # Display Off-Target Analysis (from initial submission, which currently only validates the best_guide)
                # This section remains for the overview, but individual buttons now provide detailed reports.
                off_target_analysis = results.get("off_target_analysis", {})
                st.markdown("#### 🛡️ Initial Off-Target Overview")
                if off_target_analysis.get("error"):
                    st.warning(f"Could not perform initial off-target analysis for best guide: {off_target_analysis.get('error')}")
                else:
                    hits = off_target_analysis.get("off_target_hits", [])
                    if not hits:
                        st.success("✅ No significant off-target sites found for the best guide.")
                    else:
                        st.metric("Potential Off-Target Sites (Best Guide)", len(hits))
                        with st.expander("View Initial Off-Target Details"):
                            st.json(hits)

with nanobody_designer_tab:
    st.subheader("🔬 Novel Biologic Inhibitor (Nanobody) Designer")
    st.markdown("""
    This module allows you to design novel nanobody sequences targeted to specific proteins and their functional domains. These nanobodies can act as highly specific therapeutic inhibitors.
    """)

    with st.form("nanobody_design_form"):
        target_protein_name = st.text_input(
            "Target Protein Name", 
            value="ASXL1", 
            help="Enter the name of the target protein (e.g., 'KRAS', 'PD-L1')."
        )
        target_domain = st.text_input(
            "Target Domain", 
            value="PHD Finger Domain", 
            help="Specify the functional domain to inhibit (e.g., 'Active Site', 'Binding Pocket')."
        )
        design_submitted = st.form_submit_button("Design Nanobody")

    if design_submitted:
        if not target_protein_name or not target_domain:
            st.error("Please provide both a target protein name and a target domain.")
        else:
            result = design_nanobody_inhibitor(target_protein_name, target_domain)
            st.session_state.nanobody_design_result = result
            st.rerun()

    if st.session_state.nanobody_design_result:
        st.markdown("---")
        st.subheader("🧬 Generated Nanobody Sequence")
        design_res = st.session_state.nanobody_design_result

        if design_res.get("error"):
            st.error(f"An error occurred during nanobody design: {design_res.get('error')}")
        else:
            generated_sequence_data = design_res.get("generated_sequence", {})
            nanobody_dna_sequence = generated_sequence_data.get("dna_sequence", "")
            comment = generated_sequence_data.get("comment", "")

            if nanobody_dna_sequence:
                st.code(nanobody_dna_sequence, language="dna")
                if comment:
                    st.info(comment)
            else:
                st.warning("No nanobody sequence was generated.") 