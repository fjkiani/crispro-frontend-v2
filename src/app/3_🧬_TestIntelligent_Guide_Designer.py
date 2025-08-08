import streamlit as st
import pandas as pd
import re
import json
import time
import httpx
import os
from tools.st_utils import (
    init_session_state,
    create_educational_sidebar,
    apply_custom_css,
    next_steps_module
)

# --- Page Config ---
st.set_page_config(
    page_title="Intelligent Guide Designer v4",
    page_icon="⚔️",
    layout="wide",
)

apply_custom_css()
init_session_state()
create_educational_sidebar()

# --- Client for CommandCenter ---
COMMAND_CENTER_URL = os.environ.get("COMMAND_CENTER_URL", "https://crispro--crispr-assistant-command-center-v3-commandcenter-api.modal.run")

def get_command_center_client():
    """Initializes a httpx client for the command center."""
    return httpx.Client(base_url=COMMAND_CENTER_URL, timeout=30.0)

def intelligent_guide_designer_page_v4():
    """
    UI for the Intelligent Guide Designer page, now powered by the v4 Kill Chain.
    """
    st.markdown('<div class="main-header">Intelligent Guide Designer v4</div>', unsafe_allow_html=True)
    st.write("This tool leverages the AI General Command Center to design and rank CRISPR guide RNAs using the v4 Kill Chain.")

    # Handoff from other tools
    hugo, prot_change, locus_from_handoff, sequence_from_handoff = "", "", "", ""
    if 'active_mutation' in st.session_state:
        mutation = st.session_state.active_mutation
        coord_full = mutation.get("genomic_coordinate_hg38", "")
        match = re.match(r'^(chr[A-Za-z0-9]+:\d+)', coord_full)
        if match:
            locus_from_handoff = match.group(1)
        hugo = mutation.get("hugo_gene_symbol", "N/A")
        prot_change = mutation.get("protein_change", "N/A")
        sequence_from_handoff = mutation.get("sequence_for_perplexity", "") # Assuming this field exists from handoff
            st.info(f"🧬 Received target from Digital Twin: Designing guides for {hugo} {prot_change} at locus {locus_from_handoff}.")

    # Input Form
    with st.form("guide_design_form_v4"):
        st.subheader("Design Inputs")
        patient_id_input = st.text_input("Patient Identifier", value="Patient-001")
        locus_input = st.text_input("Genomic Locus", value=locus_from_handoff, help="e.g., `chr7:140753336`")
        hgvs_input = st.text_input("Mutation (HGVS-p)", value=prot_change, help="e.g., `p.Val600Glu`")
        sequence_input = st.text_area("Target Sequence (for guide generation)", value=sequence_from_handoff, height=150)
        
        submitted = st.form_submit_button("⚔️ Execute V4 Design Campaign")

    if submitted:
        if not all([patient_id_input, locus_input, hgvs_input, sequence_input]):
            st.error("Please fill all input fields.")
        else:
            with st.spinner("Initiating V4 Kill Chain... Contacting Command Center..."):
                try:
                    client = get_command_center_client()
                    # Step 1: Formulate Battle Plan
                    formulate_payload = {
                        "patient_identifier": patient_id_input,
                        "gene": hugo,
                        "mutation_hgvs_p": hgvs_input,
                        "sequence_for_perplexity": sequence_input
                    }
                    response = client.post("/v3/workflow/formulate_battle_plan", json=formulate_payload)
                    response.raise_for_status()
                    plan_data = response.json()
                    battle_plan_id = plan_data.get("id")

                    if not battle_plan_id:
                        st.error(f"Failed to formulate battle plan: {plan_data.get('detail', 'Unknown error')}")
                        return
                    
                    st.session_state.battle_plan_id = battle_plan_id
                    st.success(f"Battle Plan #{battle_plan_id} formulated. Executing kill chain...")

                    # Step 2: Execute the kill chain
                    execute_response = client.post(f"/v4/workflow/execute_guide_design_campaign/{battle_plan_id}")
                    execute_response.raise_for_status()

                    st.session_state.design_job_running = True
                    st.rerun()

                except httpx.HTTPStatusError as e:
                    st.error(f"Command Center Error: {e.response.status_code} - {e.response.text}")
                except Exception as e:
                    st.error(f"An unexpected error occurred: {e}")

    # Polling and Results Display
    if st.session_state.get("design_job_running"):
        poll_and_display_results()

def poll_and_display_results():
    """Polls the command center for job status and displays results when complete."""
    battle_plan_id = st.session_state.battle_plan_id
    st.info(f"Monitoring Battle Plan #{battle_plan_id}. Please wait, this may take several minutes...")
    
    progress_bar = st.progress(0)
    status_text = st.empty()

    for i in range(100): # Poll for a max of ~5 minutes
        try:
            client = get_command_center_client()
            response = client.get(f"/v3/battle_plan/{battle_plan_id}")
            response.raise_for_status()
            data = response.json()
            status = data.get("status")

            status_text.text(f"Current Status: {status}")
            progress = (i + 1) / 100
            progress_bar.progress(progress)

            if status == "interventions_designed":
                st.success("V4 Kill Chain Complete! Interventions designed.")
                st.session_state.design_job_running = False
                st.session_state.final_results = data.get("results", {}).get("guides", [])
                st.rerun()
            elif status == "design_failed":
                st.error("Guide design campaign failed. Check Command Center logs for details.")
                st.session_state.design_job_running = False
                st.rerun()

            time.sleep(3)
        except Exception as e:
            st.error(f"Error while polling for results: {e}")
            st.session_state.design_job_running = False
            st.rerun()
            
    if st.session_state.design_job_running:
         st.warning("Operation timed out. The job may still be running in the background. Please check back later.")
         st.session_state.design_job_running = False
    
# Display final results if they exist
if 'final_results' in st.session_state:
    results = st.session_state.final_results
    st.subheader("Ranked Guide RNA Candidates (V4)")
    if not results:
        st.warning("No guide candidates could be generated.")
    else:
    df = pd.DataFrame(results)
        df_display = df[['guide_sequence', 'assassin_score']].rename(columns={
        "guide_sequence": "Guide Sequence",
            "assassin_score": "Assassin Score"
    })
    st.dataframe(df_display, use_container_width=True)
        st.info("Higher 'Assassin Score' indicates a better-predicted balance of efficacy and safety.")

        with st.expander("🔬 View Detailed Dossiers"):
            st.json(results)

if __name__ == "__main__":
    intelligent_guide_designer_page_v4() 