import streamlit as st
import pandas as pd
import re
import json
from tools.st_utils import (
    init_session_state,
    create_educational_sidebar,
    apply_custom_css,
    guide_finder_module,
    next_steps_module
)

# --- Page Config ---
st.set_page_config(
    page_title="Intelligent Guide Designer",
    page_icon="🧬",
    layout="wide",
)

apply_custom_css()
init_session_state()
create_educational_sidebar()


# --- Main Page UI ---

def intelligent_guide_designer_page():
    """
    UI for the Intelligent Guide Designer page.
    This function orchestrates the user interaction for designing guides with the AI backend.
    """
    st.markdown('<div class="main-header">Intelligent Guide Designer</div>', unsafe_allow_html=True)
    st.write("This tool leverages AI to design and rank CRISPR guide RNAs based on predicted on-target efficacy and potential off-target effects.")

    # Handoff from other tools
    locus_from_handoff = ""
    if 'active_mutation' in st.session_state and st.session_state.active_mutation.get("genomic_coordinate_hg38"):
        coord_full = st.session_state.active_mutation["genomic_coordinate_hg38"]
        match = re.match(r'^(chr[A-Za-z0-9]+:\d+)', coord_full)
        if match:
            locus_from_handoff = match.group(1)
            hugo = st.session_state.active_mutation.get("hugo_gene_symbol", "N/A")
            prot_change = st.session_state.active_mutation.get("protein_change", "N/A")
            st.info(f"🧬 Received target from Digital Twin: Designing guides for {hugo} {prot_change} at locus {locus_from_handoff}.")

    # Input Form
    with st.form("guide_design_form"):
        st.subheader("Design Inputs")
        locus_input = st.text_input(
            "Genomic Locus",
            value=locus_from_handoff,
            help="Enter a genomic coordinate (e.g., `chr7:140753336`) or a region (e.g., `chr7:140,753,000-140,754,000`)."
        )
        genome_input = st.selectbox("Reference Genome", ["hg38", "hg19"], index=0)
        submitted = st.form_submit_button("Design Guides")

    if submitted:
        if not locus_input:
            st.error("Please provide a genomic locus.")
        else:
            st.session_state.intelligent_guides_running = True
            st.session_state.intelligent_guides_results = None
            with st.spinner("Designing intelligent guides... This may take several minutes (UCSC, BLAST, AI calls)."):
                results = guide_finder_module.find_intelligent_guides(locus=locus_input, genome=genome_input)
                st.session_state.intelligent_guides_results = results
                st.session_state.intelligent_guides_running = False
                st.rerun()

    # Results Display
    if st.session_state.get('intelligent_guides_results') is not None:
        render_results(locus_input)

def render_results(locus_input):
    """Renders the guide design results and next steps advisor."""
    st.subheader("Ranked Guide RNA Candidates")
    results = st.session_state.intelligent_guides_results

    if not results:
        st.warning("No guide candidates could be generated for the given input. Please check the locus and genome, then try again.")
        return

    df = pd.DataFrame(results)
    df_display = df.rename(columns={
        "guide_sequence": "Guide Sequence",
        "on_target_score": "On-Target Score (AI)",
        "off_target_hits": "Off-Target Hits"
    })
    st.dataframe(df_display, use_container_width=True)
    st.info("Higher 'On-Target Score' suggests better efficacy. Lower 'Off-Target Hits' suggests better safety.")

    with st.expander("🔬 Experimental Plan Advisor"):
        if st.button("Generate Next Steps"):
            with st.spinner("Asking the AI advisor for next steps..."):
                recommendation_str = next_steps_module.get_next_steps_recommendation(
                    editing_type="Gene Knockout",
                    result_data={"status": "Guide Design Complete"},
                    issue_list=[],
                    therapeutic_context={"goal": f"Targeting locus {locus_input}"},
                    ranked_guides=results
                )
                try:
                    recommendation = json.loads(recommendation_str)
                    st.success(f"**{recommendation.get('recommendation_title', 'Recommendation')}**")
                    st.markdown(f"**Summary:** {recommendation.get('summary', '')}")
                    st.markdown("**Next Steps:**")
                    for step in recommendation.get('next_steps', []):
                        st.markdown(f"- {step}")
                except (json.JSONDecodeError, TypeError):
                    st.markdown(recommendation_str)


if __name__ == "__main__":
    intelligent_guide_designer_page() 