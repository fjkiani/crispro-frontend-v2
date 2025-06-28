import streamlit as st
from tools.intelligent_guide_finder import find_intelligent_guides
from tools.next_steps import get_next_steps_recommendation
import pandas as pd
import json

st.set_page_config(
    page_title="Intelligent Guide Designer",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 Intelligent Guide Designer")
st.write("This tool leverages AI to design and rank CRISPR guide RNAs based on predicted on-target efficacy and potential off-target effects.")

# --- Session State Initialization ---
if 'intelligent_guides_results' not in st.session_state:
    st.session_state.intelligent_guides_results = None
if 'intelligent_guides_running' not in st.session_state:
    st.session_state.intelligent_guides_running = False

# --- Handoff from Digital Twin ---
locus_from_twin = ""
if "targeted_variant_for_crispr" in st.session_state and st.session_state.targeted_variant_for_crispr:
    variant_data = st.session_state.targeted_variant_for_crispr
    locus_from_twin = variant_data.get("target_locus", "")
    st.info(f"🧬 Received target from Digital Twin: Designing guides for {variant_data.get('variant', 'N/A')}.")
    # Clear the state variable so it doesn't persist across reloads
    del st.session_state.targeted_variant_for_crispr

# --- Input Form ---
with st.form("guide_design_form"):
    st.subheader("Design Inputs")
    locus_input = st.text_input(
        "Genomic Locus", 
        value=locus_from_twin, 
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
        
        with st.spinner("Designing intelligent guides... This may take several minutes as it involves multiple API calls (UCSC, BLAST, and AI models)."):
            results = find_intelligent_guides(locus=locus_input, genome=genome_input)
            st.session_state.intelligent_guides_results = results
            st.session_state.intelligent_guides_running = False

# --- Results Display ---
if st.session_state.intelligent_guides_results is not None:
    st.subheader("Ranked Guide RNA Candidates")
    results = st.session_state.intelligent_guides_results
    
    if results:
        df = pd.DataFrame(results)
        df.rename(columns={
            "guide_sequence": "Guide Sequence",
            "on_target_score": "On-Target Score (Efficacy)",
            "off_target_hits": "Off-Target Hits (Safety)"
        }, inplace=True)
        
        st.dataframe(df, use_container_width=True)
        st.info("Higher 'On-Target Score' suggests better efficacy. Lower 'Off-Target Hits' suggests better safety.")

        # --- Next Steps Advisor ---
        st.subheader("🔬 Experimental Plan Advisor")
        if st.button("Generate Next Steps"):
            with st.spinner("Asking the AI advisor for next steps..."):
                # We can create a simplified context for the advisor
                recommendation_str = get_next_steps_recommendation(
                    editing_type="Gene Knockout",
                    result_data={"status": "Guide Design Complete"},
                    issue_list=[],
                    therapeutic_context={"goal": f"Targeting locus {locus_input}"},
                    ranked_guides=results
                )
                
                try:
                    # The new function should return a JSON string
                    recommendation = json.loads(recommendation_str)
                    st.success(f"**{recommendation.get('recommendation_title', 'Recommendation')}**")
                    st.markdown(f"**Summary:** {recommendation.get('summary', '')}")
                    
                    st.markdown("**Next Steps:**")
                    for step in recommendation.get('next_steps', []):
                        st.markdown(f"- {step}")
                except (json.JSONDecodeError, TypeError):
                    # Fallback for non-JSON responses
                    st.markdown(recommendation_str)

    else:
        st.warning("No guide candidates could be generated for the given input. Please check the locus and try again.")
elif st.session_state.intelligent_guides_running:
    st.info("The guide design process is running in the background.") 