import streamlit as st
import json
from tools.myeloma_digital_twin import predict_myeloma_drug_response, RAS_MAPK_PATHWAY_GENES, TP53_GENE

st.set_page_config(
    page_title="Myeloma Digital Twin",
    page_icon="🧬",
    layout="wide"
)

# --- Page Title and Introduction ---
st.title("🔬 Myeloma Digital Twin")
st.markdown("""
This module uses the **Evo 2 AI model** to predict a Multiple Myeloma patient's likely response to therapies targeting the **MAPK pathway**. 
Enter a patient's somatic mutations below to generate a real-time analysis. The model assesses the functional impact of each variant to determine resistance potential.
""")

# --- Helper Functions ---
def get_impact_color(impact_level):
    if impact_level >= 3:
        return "red"
    elif impact_level >= 2:
        return "orange"
    elif impact_level >= 1:
        return "gold"
    else:
        return "grey"

# --- Main UI ---
st.subheader("Patient Mutation Input")

# Initialize session state for mutations
if 'mutations' not in st.session_state:
    st.session_state.mutations = [
        {"id": 1, "gene": "BRAF", "hgvs_p": "p.Val600Glu", "variant_info": "chr7:140753336 A>T", "build": "hg38"},
    ]
if 'next_mutation_id' not in st.session_state:
    st.session_state.next_mutation_id = 2
    
# Display input fields for each mutation
for i, mutation in enumerate(st.session_state.mutations):
    st.markdown(f"**Mutation {i+1}**")
    cols = st.columns([2, 2, 4, 1, 1])
    mutation['gene'] = cols[0].text_input("Gene Symbol", value=mutation['gene'], key=f"gene_{mutation['id']}")
    mutation['hgvs_p'] = cols[1].text_input("Protein Change (HGVS_p)", value=mutation['hgvs_p'], key=f"hgvs_p_{mutation['id']}")
    mutation['variant_info'] = cols[2].text_input("Genomic Coordinate (e.g., chr:pos REF>ALT)", value=mutation['variant_info'], key=f"variant_info_{mutation['id']}")
    mutation['build'] = cols[3].selectbox("Build", ["hg38", "hg19"], key=f"build_{mutation['id']}", index=["hg38", "hg19"].index(mutation['build']))
    
    if cols[4].button("Remove", key=f"remove_{mutation['id']}"):
        st.session_state.mutations.pop(i)
        st.rerun()

# Add/Remove buttons
cols = st.columns(2)
if cols[0].button("＋ Add Mutation"):
    new_id = st.session_state.next_mutation_id
    st.session_state.mutations.append({"id": new_id, "gene": "", "hgvs_p": "", "variant_info": "", "build": "hg38"})
    st.session_state.next_mutation_id += 1
    st.rerun()


# --- Analysis Execution ---
if st.button("🚀 Analyze Drug Response", type="primary", use_container_width=True):
    with st.spinner("Analyzing variants with Evo 2 model... This may take a moment."):
        # Filter out empty mutations before sending to the backend
        mutations_to_analyze = [m for m in st.session_state.mutations if m['gene'] and m['variant_info']]
        if not mutations_to_analyze:
            st.error("Please enter at least one valid mutation.")
        else:
            prediction_output = predict_myeloma_drug_response(mutations_to_analyze)
            st.session_state.prediction_output = prediction_output
            
# --- Results Display ---
if 'prediction_output' in st.session_state:
    output = st.session_state.prediction_output
    
    st.subheader("Analysis Complete")
    
    # Display final prediction
    st.info(f"**Final Prediction:** {output['prediction']}")

    # Display summary metrics
    cols = st.columns(2)
    cols[0].metric("RAS/MAPK Pathway Impact Score", output['summed_impact_ras_pathway'])
    cols[1].metric("TP53 Pathway Impact Score", output['summed_impact_tp53'])
    
    with st.expander("Show Detailed Variant-Level Analysis", expanded=True):
        for detail in output['detailed_analysis']:
            impact_level = detail['calculated_impact_level']
            impact_color = get_impact_color(impact_level)
            variant_info = detail['evo2_result']

            st.markdown(f"""
            <div style="border-left: 6px solid {impact_color}; padding-left: 10px; margin-bottom: 10px; padding-bottom: 10px;">
                <h4>Variant: {detail['variant']}</h4>
                <p>
                    <b>Calculated Impact Level: {impact_level}</b><br>
                    Evo2 Prediction: <i>'{variant_info.get('classification', 'N/A')}'</i><br>
                    Confidence: {variant_info.get('confidence', 0):.2%}<br>
                    Delta Score: {variant_info.get('delta_score', 0):.6f}
                </p>
            </div>
            """, unsafe_allow_html=True)

            # --- Task 1.1: Add "Call to Action" Button ---
            if impact_level >= 2: # Condition for showing the button
                # This check ensures we have the data structure we need
                if "original_variant_data" in detail:
                    cols = st.columns([1, 4]) 
                    with cols[0]:
                        if st.button("🎯 Design CRISPR Therapy", key=f"design_{detail['variant']}"):
                            # --- Task 1.2: Store Context in Session State ---
                            st.session_state['targeted_variant_for_crispr'] = {
                                "gene": detail['gene'],
                                "variant": detail['variant'],
                                # We construct a locus string that our new tool can use
                                "target_locus": f"{detail['chrom']}:{detail['pos']}"
                            }
                            
                            st.success(f"Preparing to design therapy for {detail['variant']}. Navigating to Guide Designer...")
                            # --- Task 1.3: Navigation --_
                            st.switch_page("pages/2_Intelligent_Guide_Designer.py")

# --- Sidebar with Explanations ---
with st.sidebar:
    st.header("About the Myeloma Digital Twin")
    st.markdown("""
    **Pathway Analysis:**
    The model specifically checks for mutations in genes known to drive resistance in Multiple Myeloma.
    - **RAS/MAPK Pathway:** `{', '.join(RAS_MAPK_PATHWAY_GENES)}`
    - **TP53 Pathway:** `{', '.join(TP53_GENE)}`

    **Impact Scoring:**
    The "Impact Level" (0-3) is calculated based on the Evo 2 model's confidence and the severity of the predicted effect (delta score). A higher score indicates a greater likelihood that the variant contributes to drug resistance.

    **Prediction Logic:**
    - **Likely Resistant:** High total impact score in the RAS/MAPK pathway.
    - **Potential Resistance:** Significant impact in the TP53 pathway or moderate RAS/MAPK impact.
    - **Likely Sensitive:** Low or no impact in key resistance pathways.
    """)
    st.markdown("---")
    st.info("This is a research tool. All outputs should be interpreted by qualified experts.") 