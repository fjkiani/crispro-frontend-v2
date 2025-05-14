import streamlit as st
import sys
import os
import time # Import time for potential delays if needed for visual effect
import random # Import random if needed for frontend mock data generation

# Add the project root to the Python path to allow importing ai_research_assistant
# Assuming streamlit_app.py is in the project root, and ai_research_assistant.py is in a subdirectory
# Adjust the path based on your actual project structure if necessary
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if PROJECT_ROOT not in sys.path:
    sys.path.append(PROJECT_ROOT)

try:
    # Import the specific handler function from the ai_research_assistant module
    # Note: The ai_research_assistant module needs to be importable,
    # and handle_design_therapy_interaction needs to be defined within it.
    # We also assume handle_design_therapy_interaction is modified to accept a display_callback.
    from ai_research_assistant import (
        handle_design_therapy_interaction, 
        simulate_component_sequence_refinement,
        simulate_guide_rna_optimization,
        simulate_off_target_analysis,
        simulate_repair_template_optimization,
        simulate_in_vitro_characterization,
        simulate_cell_based_assay,
        simulate_off_target_validation,
        simulate_delivery_system_optimization,
        simulate_therapeutic_development_planning
    )
    # Placeholder for the actual get_llm_chat_response if it's not in ai_research_assistant and needs to be imported separately
    # from tools.llm_api import get_llm_chat_response
except ImportError as e:
    st.error(f"Failed to import `ai_research_assistant` module or `handle_design_therapy_interaction`. Make sure the path is correct and the function exists: {e}")
    # Fallback for critical function if import fails, to prevent app crash
    def handle_design_therapy_interaction(*args, **kwargs):
        st.error("`handle_design_therapy_interaction` is not available due to import error.")
        # Return an empty list of recommendations in case of failure
        return []

# --- Streamlit Page Function ---
def display_mock_therapy_design_page():
    st.set_page_config(page_title="Mock Therapeutic System Design", layout="wide")
    st.title("🧬 Mock Multi-Gene Therapeutic System Design")
    st.markdown("""
    This page allows you to simulate the design of a multi-gene CRISPR-based therapeutic system
    for a target cancer mutation. The process uses mock Evo 2, CHOPCHOP, and AlphaFold 3 data for demonstration.
    **All outputs, scores, and rationales are SIMULATED.**
    """)

    # --- Explanation Expander ---
    with st.expander("ℹ️ Understanding These Mock Results", expanded=False):
        st.markdown("""
        This simulation demonstrates a conceptual workflow for designing multi-gene therapeutic systems
        by integrating simulated data from different tools: Evo 2 (generation & target impact),
        CHOPCHOP (guide validation), and AlphaFold 3 (structural prediction).

        **Workflow Stages (Simulated):**
        1.  **Target Identification & Prioritization:** Analyzing the original mutation's potential impact (simulated Evo 2).
        2.  **Design Goal Specification:** Defining the therapeutic aim (e.g., disrupt, correct).
        3.  **System Component Generation:** Proposing sequences for system parts (simulated Evo 2).
        4.  **Initial Guide Evaluation:** Validating guide RNA components (simulated CHOPCHOP).
        5.  **Structural Evaluation:** Predicting 3D structures of components/complexes (simulated AlphaFold 3).
        6.  **Integrated Scoring:** Combining all simulated data to rank candidates.
        7.  **Recommendation:** Presenting the top-scoring system designs.

        **Supporting Evidence Metrics (Simulated):**
        These are mock values but represent real-world concepts. They provide the data points used for the Confidence Score.

        * **`Mock Original Variant Evo2 Confidence`**: (Simulated Evo2 VEP score) Represents the predicted confidence that the *original* target mutation has a significant biological impact related to cancer (e.g., it's pathogenic, oncogenic, or confers drug resistance). A higher value suggests the target is more likely to be relevant for therapeutic intervention.
        * **`Mock Original Variant Delta Likelihood`**: (Simulated Evo2 VEP score) Predicts how the mutation changes the likelihood of the sequence appearing naturally. Large negative values can indicate variants that are strongly selected against, often including disease-causing or functionally disruptive mutations.
        * **`Mock Overall System AF Score`**: (Simulated AlphaFold3 score) An aggregate structural assessment score for the entire multi-gene system's predicted 3D configuration. A higher score suggests that the system's components are likely to fold correctly and interact as intended, which is crucial for function.
        * **`Comp '<component_name>' Avg AF`**: (Simulated AlphaFold3 metric) An average of mock structural quality metrics (like pLDDT, ranking, and structural similarity) for an individual protein or nucleic acid component within the system. A higher value indicates better predicted structural quality and potential stability for that specific part.
        * **`Mock CHOPCHOP On-target Score`**: (Simulated CHOPCHOP score) Predicts the likelihood or efficiency of the designed guide RNA leading to a cut at the intended genomic target site. This often considers sequence features like GC content and position relative to the PAM. A higher score suggests better on-target editing potential.
        * **`Mock CHOPCHOP Off-target Count`**: (Simulated CHOPCHOP count) The predicted number of sites in the genome (other than the intended target) where the designed guide RNA could potentially bind and cause unintended edits. Minimizing off-target activity is critical for safety. A lower count is better for specificity.
        * **`Mock CHOPCHOP MFE`**: (Simulated CHOPCHOP score) The Minimum Free Energy (MFE) of the designed guide RNA's predicted secondary structure. MFE is a standard measure of RNA folding stability; more negative values indicate a more stable predicted structure. Stable guide RNA structure is important for efficient binding to the Cas protein and the target DNA.

        **Designed System Components (Mock):** Placeholder sequences representing the genetic parts of the designed therapeutic system (e.g., the sequence for a guide RNA, a Cas protein, or a DNA repair template).

        **Mock AlphaFold System Metrics per Component:** Raw simulated metric values from AlphaFold 3 predictions for each individual component's structure.

        """)


    # --- 1. Initialize Session State ---
    # Initialize conversation history if not already present
    if 'conversation_history' not in st.session_state:
        st.session_state.conversation_history = []
    # Initialize active mutation if not already present
    if 'active_mutation' not in st.session_state:
        # Initialize with a default or provide input fields
        st.session_state.active_mutation = {
            "hugo_gene_symbol": "BRAF", "protein_change": "V600E", "variant_type": "Missense_Mutation",
            "genomic_coordinate_hg38": "chr7:140753336A>T", "allele_frequency": 0.45, "mutation_id": "MU100001"
        }
    # Initialize therapy recommendations list
    if 'therapy_recommendations' not in st.session_state:
        st.session_state.therapy_recommendations = []
    # Initialize workflow log list
    if 'workflow_log' not in st.session_state:
        st.session_state.workflow_log = []


    # --- 2. Target Mutation Input/Display ---
    st.sidebar.subheader("🎯 Target Mutation")
    # Allow user to modify or confirm the target mutation (simplified for now)
    gene_symbol = st.sidebar.text_input("Gene Symbol", value=st.session_state.active_mutation.get("hugo_gene_symbol", "BRAF")) # Corrected key
    protein_change = st.sidebar.text_input("Protein Change", value=st.session_state.active_mutation.get("protein_change", "V600E"))
    # Ideally, more fields or a selection mechanism from a list of mutations

    # Update active_mutation if changed by user
    # This is a basic way; a more robust app might have a dedicated mutation selection/management component
    current_target_mutation = {
        "hugo_gene_symbol": gene_symbol,
        "protein_change": protein_change,
        "variant_type": st.session_state.active_mutation.get("variant_type", "Missense_Mutation"), # Keep others for now
        "genomic_coordinate_hg38": st.session_state.active_mutation.get("genomic_coordinate_hg38", ""),
        "allele_frequency": st.session_state.active_mutation.get("allele_frequency", 0.0),
        "mutation_id": st.session_state.active_mutation.get("mutation_id", "DEFAULT_MUT")
    }
    # Update session state if the user changed the input fields
    if (st.session_state.active_mutation.get("hugo_gene_symbol") != gene_symbol or
        st.session_state.active_mutation.get("protein_change") != protein_change):
        st.session_state.active_mutation = current_target_mutation
        # Clear previous recommendations and log if the target changes
        st.session_state.therapy_recommendations = []
        st.session_state.workflow_log = []
        st.rerun() # Rerun to clear the display based on the new target

    st.sidebar.info(f"Designing for: **{current_target_mutation['hugo_gene_symbol']} {current_target_mutation['protein_change']}**")


    # --- 3. Trigger Button & Workflow Execution ---
    # Use a form to prevent rerunning the script on every input change within the form
    with st.form("design_form"):
        st.subheader("⚙️ Run Design Workflow")
        # Add any other design parameters here if needed (e.g., design goal selection)
        # For now, the design goal is hardcoded in the backend mock handler

        submit_button = st.form_submit_button("🚀 Design Mock Therapy System", type="primary", help="Click to start the mock design process")

    if submit_button:
        st.session_state.therapy_recommendations = [] # Clear previous recommendations
        st.session_state.workflow_log = [] # Clear previous log

        # Create a placeholder for the workflow log
        log_placeholder = st.empty()

        # Define the callback function to update the log in Streamlit
        def streamlit_display_callback(message):
            # Append message to the log list in session state
            st.session_state.workflow_log.append(message)
            # Update the display using the placeholder
            with log_placeholder.container(): # Use a container within the placeholder
                 st.subheader("⚙️ Workflow Log")
                 for log_message in st.session_state.workflow_log:
                     st.text(log_message) # Display each log message as text

        # Run the backend design function with the callback
        with st.spinner("Simulating Evo 2, CHOPCHOP, AlphaFold 3, and generating recommendations... This may take a moment."):
            try:
                # Pass the current target mutation and the callback function
                # Note: handle_design_therapy_interaction needs to accept these arguments
                returned_recommendations = handle_design_therapy_interaction(
                    current_target_mutation, # Pass the current target mutation details
                    st.session_state.conversation_history, # Pass the conversation history list
                    display_callback=streamlit_display_callback # Pass the callback function
                )
                st.session_state.therapy_recommendations = returned_recommendations
            except Exception as e:
                st.error(f"An error occurred during the design process: {e}")
                # Append error to the log as well
                st.session_state.workflow_log.append(f"ERROR: {e}")
                st.session_state.therapy_recommendations = []

        # Display final status message
        if st.session_state.therapy_recommendations:
            st.success("✅ Mock therapy system design process completed! See recommendations below.")
        else:
            st.info("ℹ️ No mock recommendations were generated or the process encountered an issue.")

    # --- Display the Workflow Log (will be updated by the callback during execution) ---
    # This section is primarily for showing the log after a rerun or if the user expands it later
    if st.session_state.workflow_log and not submit_button: # Don't show this area if the button was just clicked, let the placeholder handle it
         with st.expander("⚙️ Workflow Log", expanded=False):
              for log_message in st.session_state.workflow_log:
                  st.text(log_message)


    # --- 4. Display Results ---
    if st.session_state.therapy_recommendations:
        st.subheader("📊 Simulated System Recommendations")

        for i, rec in enumerate(st.session_state.therapy_recommendations):
            approach_text = rec.get('recommended_approach', f'Recommendation {i+1}')
            confidence_score = rec.get('confidence_score', 'N/A')

            # Ensure expander_key is not used if it was causing issues, or define it if needed for other purposes
            # For now, removing the key argument from st.expander as it was the source of a TypeError.
            with st.expander(f"Recommendation {i+1}: {approach_text} (Simulated Confidence: {confidence_score:.2f})", expanded=(i == 0)):
                # Display Detailed Rationale
                st.markdown(f"**Detailed Rationale (Simulated):** {rec.get('detailed_rationale', 'Detailed rationale not available.')}")
                st.markdown("---") # Separator after rationale

                tab_evidence, tab_components, tab_metrics, tab_next_steps = st.tabs(["Supporting Evidence", "🧬 System Components", "🔬 Evaluation Metrics", "➡️ Suggested Next Steps"])

                with tab_evidence:
                    st.markdown("**Supporting Evidence (Mock):**")
                    st.caption("Metrics below are simulated. See 'Understanding These Mock Results' expander above for detailed explanations if needed.")
                    supporting_evidence = rec.get('supporting_evidence', [])
                    if supporting_evidence:
                        for evidence_item in supporting_evidence:
                            # Default explanation if no specific match
                            explanation = ""

                            if "Mock Original Variant Evo2 Confidence" in evidence_item:
                                explanation = "_(Simulated Evo2 VEP score. Represents predicted confidence that the original target mutation has a significant biological impact related to cancer. Higher suggests a good therapeutic target.)_"
                            elif "Mock Original Variant Delta Likelihood" in evidence_item:
                                explanation = "_(Simulated Evo2 VEP score. Predicts how the mutation changes the likelihood of the sequence appearing naturally. Large negative values can indicate impactful variants.)_"
                            elif "Mock Overall System AF Score" in evidence_item:
                                explanation = "_(Simulated AlphaFold3 score. Overall structural assessment of the entire multi-gene system. Higher suggests components are likely to fold/interact correctly.)_"
                            elif "Comp " in evidence_item and "Avg AF" in evidence_item:
                                explanation = "_(Simulated average AlphaFold metrics (pLDDT, ranking, similarity) for an individual component. Higher indicates better predicted structural quality for that component.)_"
                            elif "Mock CHOPCHOP On-target Score" in evidence_item:
                                explanation = "_(Simulated CHOPCHOP score. Predicted on-target efficacy of the guide RNA. Higher is better.)_ "
                            elif "Mock CHOPCHOP Off-target Count" in evidence_item:
                                explanation = "_(Simulated CHOPCHOP count. Predicted number of off-target sites for the guide RNA. Lower is better for specificity.)_ "
                            elif "Mock CHOPCHOP MFE" in evidence_item:
                                explanation = "_(Simulated CHOPCHOP score. Minimum Free Energy for the guide RNA's predicted secondary structure. More negative indicates a more stable structure.)_ "
                            # Fallback for pLDDT, Ranking, Struct Similarity if they appear outside "Comp Avg AF" context, though less likely now
                            elif "pLDDT" in evidence_item: # Should be caught by "Comp Avg AF" if it's an average
                                explanation = "_(Simulated protein structure confidence; 0-100, higher >70 is generally better)_ "
                            elif "Ranking" in evidence_item: # Should be caught by "Comp Avg AF"
                                explanation = "_(Simulated quality/binding score; higher is generally better)_ "
                            elif "Struct Similarity" in evidence_item: # Should be caught by "Comp Avg AF"
                                explanation = "_(Simulated structural similarity to a reference; interpretation is context-dependent)_ "

                            st.markdown(f"- {evidence_item} {explanation}")
                    else:
                        st.markdown("  _No supporting evidence provided._")

                with tab_components:
                    st.markdown("**Designed System Components (Mock):**")
                    st.caption("Placeholder sequences representing the genetic parts of the designed therapeutic system. In a real-world scenario, these would be actual DNA or RNA sequences.")
                    components = rec.get('designed_system_components', {})
                    if components:
                        if "system_id" in components:
                            st.markdown(f"- **System ID**: `{components['system_id']}`")
                            st.caption("  _A unique identifier for this mock designed therapeutic system candidate._")
                        if "description" in components:
                            st.markdown(f"- **Description**: {components['description']}")
                            st.caption("  _A brief summary of this mock system candidate and its intended target._")

                        # Specific explanations for common components and display mock sequences
                        st.markdown("**Mock Designed Sequences:**")
                        component_sequences_found = False
                        if "guide_rna_sequence" in components:
                            st.markdown(f"- **Guide RNA Sequence**: `{components['guide_rna_sequence'][:60]}...`")
                            st.caption("""
                                _Mock placeholder sequence for the guide RNA (gRNA). In CRISPR systems, the gRNA is a short synthetic RNA molecule that consists of a scaffold sequence necessary for binding to the Cas protein and a spacer sequence that is complementary to the target DNA sequence. The spacer sequence directs the Cas protein to the specific genomic location where the edit is intended._
                            """)
                            component_sequences_found = True
                        if "cas_protein_sequence" in components:
                            st.markdown(f"- **Cas Protein Sequence**: `{components['cas_protein_sequence'][:60]}...`")
                            st.caption("""
                                _Mock placeholder sequence for the Cas enzyme (CRISPR-associated protein), such as Cas9. The Cas protein is the 'molecular scissors' that, guided by the gRNA, makes a precise cut in the DNA at the target site. Different Cas proteins exist with varying properties and PAM requirements._
                            """)
                            component_sequences_found = True
                        if "repair_template_dna" in components:
                            st.markdown(f"- **Repair Template DNA**: `{components['repair_template_dna'][:60]}...`")
                            st.caption("""
                                _Mock placeholder sequence for a DNA repair template. This template is used in homology-directed repair (HDR), a cellular DNA repair pathway. When a double-strand break is made by the Cas enzyme, the cell can use this provided template as a blueprint to accurately repair the break, allowing for precise gene correction or insertion (knock-in) of new genetic material._
                            """)
                            component_sequences_found = True
                        if "modified_cas_protein_sequence" in components:
                            st.markdown(f"- **Modified Cas Protein Sequence**: `{components['modified_cas_protein_sequence'][:60]}...`")
                            st.caption("""
                                _Mock placeholder sequence for a modified Cas enzyme. Cas proteins can be engineered to alter their properties, such as increasing specificity (reducing off-target edits), changing their size for better delivery, or modifying their enzymatic activity (e.g., creating a 'nickase' that cuts only one DNA strand, or a catalytically inactive dCas9 for gene activation/repression)._
                            """)
                            component_sequences_found = True
                        if "viral_packaging_signal_dna" in components:
                            st.markdown(f"- **Viral Packaging Signal DNA**: `{components['viral_packaging_signal_dna'][:60]}...`")
                            st.caption("""
                                _Mock placeholder sequence for a viral packaging signal. This is a specific DNA sequence required for packaging the therapeutic genetic construct (containing the Cas gene, gRNA gene, repair template, etc.) into a viral vector, such as an Adeno-Associated Virus (AAV), for efficient delivery into target cells._
                            """)
                            component_sequences_found = True
                        if "therapeutic_gene_sequence" in components:
                             st.markdown(f"- **Therapeutic Gene Sequence**: `{components['therapeutic_gene_sequence'][:60]}...`")
                             st.caption("""
                                _Mock placeholder for the actual gene sequence intended to have a therapeutic effect. This could be a gene to be inserted (knock-in), a gene to replace a mutated one (gene correction), or a gene designed to express a therapeutic protein._
                            """)
                             component_sequences_found = True
                        if "promoter_sequence" in components:
                             st.markdown(f"- **Promoter Sequence**: `{components['promoter_sequence'][:60]}...`")
                             st.caption("""
                                _Mock placeholder for a promoter sequence. Promoters are DNA sequences located near genes that act as binding sites for proteins that initiate transcription. They control when and where a gene is expressed, ensuring that the therapeutic genes in the system are turned on in the target cells._
                            """)
                             component_sequences_found = True

                        # Generic fallback for other components if any
                        st.markdown("**Other Mock Sequences:**")
                        other_comp_count = 0
                        for comp_name, comp_val in components.items():
                            if comp_name not in ["system_id", "description", "guide_rna_sequence",
                                                "cas_protein_sequence", "repair_template_dna",
                                                "modified_cas_protein_sequence", "viral_packaging_signal_dna",
                                                "therapeutic_gene_sequence", "promoter_sequence"]:
                                other_comp_count +=1
                                if isinstance(comp_val, str) and comp_val.startswith("MOCK_"):
                                    st.markdown(f"- **{comp_name.replace('_', ' ').title()}**: `{comp_val[:60]}...`")
                                elif isinstance(comp_val, str):
                                    st.markdown(f"- **{comp_name.replace('_', ' ').title()}**: {comp_val}")
                        if other_comp_count == 0 and not component_sequences_found:
                             st.markdown("  _No component details available._")
                        elif other_comp_count > 0:
                             st.caption("  _Above are other mock sequences for necessary parts of a therapeutic system._")


                    else:
                        st.markdown("  _No component details available._")


                with tab_metrics:
                    st.markdown("**Evaluation Metrics (Mock):**")
                    st.caption("Simulated metrics from AlphaFold 3 and CHOPCHOP evaluations.")

                    # Display Overall AlphaFold Score
                    # Corrected: Access 'mock_overall_af_score' from the 'rec' dictionary
                    overall_af_score = rec.get('mock_overall_af_score', 'N/A')
                    st.markdown(f"- **Overall Mock AlphaFold System Score:** {overall_af_score:.2f}")
                    st.caption("  _Simulated AlphaFold3 aggregate score for the entire system's structural integrity. Higher suggests better folding/interaction._")
                    st.markdown("---")

                    # Display CHOPCHOP Metrics (if available)
                    chopchop_metrics_system = rec.get('mock_chopchop_metrics', {})
                    if chopchop_metrics_system:
                        st.markdown("**Mock CHOPCHOP Metrics (Guide RNA):**")
                        # Assuming guide RNA metrics are stored under 'guide_rna_sequence' key
                        guide_chopchop_metrics = chopchop_metrics_system.get('guide_rna_sequence', {})
                        if guide_chopchop_metrics:
                            st.markdown(f"  - Mock On-target Score: {guide_chopchop_metrics.get('chopchop_on_target_score', 'N/A'):.2f}")
                            st.caption("    _Simulated CHOPCHOP score. Predicted on-target efficacy of the guide RNA. This score often considers sequence features like GC content and position relative to the PAM sequence. A higher score suggests a greater likelihood of the guide RNA leading to a successful cut at the intended genomic location._")
                            st.markdown(f"  - Mock Off-target Count: {guide_chopchop_metrics.get('chopchop_off_target_count', 'N/A')}")
                            st.caption("    _Simulated CHOPCHOP count. This is the predicted number of sites in the genome (other than the intended target) where the designed guide RNA could potentially bind and cause unintended edits. Minimizing off-target activity is crucial for the safety and specificity of CRISPR therapies. A lower count indicates higher specificity and reduced risk of unwanted edits._")
                            st.markdown(f"  - Mock MFE: {guide_chopchop_metrics.get('chopchop_mfe', 'N/A'):.2f}")
                            st.caption("    _Simulated Minimum Free Energy (MFE) for the guide RNA's predicted secondary structure. MFE is a standard measure in RNA folding; more negative values indicate a more stable predicted RNA structure. A stable guide RNA structure is important for its proper interaction with the Cas protein and efficient binding to the target DNA._")
                        else:
                             st.markdown("  _No mock CHOPCHOP metrics available for the guide RNA component._")
                        st.markdown("---")


                    # Display Component-Specific AlphaFold Metrics
                    component_metrics_af = rec.get('mock_alphafold_system_metrics', {})
                    if component_metrics_af:
                        st.markdown("**Mock AlphaFold Metrics per Component:**")
                        for comp_name, metrics_dict in component_metrics_af.items():
                            # Skip system_id and description as they don't have these bio-metrics
                            if comp_name.lower() in ["system_id", "description"]:
                                continue

                            st.markdown(f"**{comp_name.replace('_', ' ').title()} Metrics:**")
                            if metrics_dict:
                                if "mock_plddt_score" in metrics_dict:
                                    st.markdown(f"  - Mock pLDDT Score: {metrics_dict['mock_plddt_score']:.2f}")
                                    st.caption("    _Simulated AlphaFold metric (ranging from 0 to 100). Predicted Local Distance Difference Test (pLDDT) indicates the confidence in the local accuracy of the predicted 3D structure for a protein component. Scores above 70 generally indicate a good prediction, while scores above 90 suggest very high confidence in the predicted structure._")
                                if "mock_ranking_score" in metrics_dict:
                                    st.markdown(f"  - Mock Ranking Score: {metrics_dict['mock_ranking_score']:.2f}")
                                    st.caption("    _Simulated AlphaFold metric. This score often reflects the overall quality of a predicted structure or the confidence in a predicted interaction (like protein-DNA binding). A higher score generally indicates a more reliable prediction._")
                                if "mock_structural_similarity_score" in metrics_dict:
                                    st.markdown(f"  - Mock Structural Similarity Score: {metrics_dict['mock_structural_similarity_score']:.2f}")
                                    st.caption("    _Simulated AlphaFold metric. Measures the structural similarity of the designed component's predicted structure to a reference structure (e.g., the wild-type protein). High similarity is desirable for designs aiming to restore function, while low similarity might be expected for designs intended to disrupt protein structure._")
                                if "mock_mfe" in metrics_dict: # Specific to RNA components
                                    st.markdown(f"  - Mock MFE: {metrics_dict['mock_mfe']:.2f}")
                                    st.caption("    _Simulated Minimum Free Energy (MFE) for RNA structure. More negative values indicate a more stable predicted RNA secondary structure. This is relevant for RNA components like guide RNAs, where proper folding is essential for function._")
                                # Check if any expected metrics were found for this component, otherwise print a generic message
                                if not any(k in metrics_dict for k in ["mock_plddt_score", "mock_ranking_score", "mock_structural_similarity_score", "mock_mfe"]):
                                    st.caption("    _No specific AlphaFold-simulated metrics (pLDDT, Ranking, MFE, Similarity) for this component._")
                            else:
                                st.markdown("  _No metrics available for this component._")
                            st.markdown(" ")
                    else:
                        st.markdown("  _No component-specific AlphaFold metrics available._")

                with tab_next_steps:
                    st.markdown("**Suggested Next Steps Based on Simulation:**")
                    st.caption("These are conceptual suggestions to guide further refinement and planning based on the simulated design and evaluation metrics. Experimental validation is always crucial.")

                    # Display In Silico Refinement Suggestions
                    st.markdown("### 🛠️ In Silico Refinement Suggestions")
                    st.markdown("""
                    Based on the simulated evaluation metrics, the following *in silico* (computational) steps could be considered to potentially improve the design:
                    * **Refine Component Sequences:** For components with lower simulated structural confidence (e.g., low pLDDT for a protein like the Cas enzyme) or other unfavorable metrics, explore generating or optimizing alternative sequences using tools like Evo 2. This could involve:
                        * Investigating regions of predicted disorder or potential misfolding.
                        * Considering sequence variants that are predicted to improve stability or desired interactions.
                        * For protein components, exploring codon optimization for the intended host expression system to improve protein yield and folding.
                    * **Optimize Guide RNA:** If the simulated CHOPCHOP metrics for the guide RNA are not ideal (e.g., high off-target count, suboptimal on-target score, unstable MFE), explore designing alternative guide sequences targeting the same region or nearby sites.
                    * **Advanced Off-Target Analysis:** The simulated CHOPCHOP off-target count provides an initial estimate. For therapeutic applications, a more thorough *in silico* off-target analysis using tools that search against the entire genome with more sophisticated algorithms is crucial. This involves identifying all potential off-target sites and assessing their likelihood of being edited.
                    * **Refine Repair Template (if applicable):** If the design involves a repair template for gene correction or knock-in, optimize the template sequence, including adjusting the length and sequence of homology arms, to improve the efficiency and accuracy of homology-directed repair (HDR).
                    """)

                    # Display Experimental Validation Plan
                    st.markdown("### 🧪 Experimental Validation Plan (Conceptual)")
                    st.markdown("""
                    Promising *in silico* designs must be validated experimentally. The following steps outline a typical conceptual plan:
                    * **In Vitro Characterization:**
                        * Synthesize or transcribe the designed guide RNA(s).
                        * Express and purify the designed (or selected) Cas protein.
                        * Perform *in vitro* cleavage assays using PCR amplicons of the target genomic region to confirm that the Cas-gRNA complex can efficiently cut the target DNA at the predicted site.
                    * **Cell-Based Assays:**
                        * Introduce the therapeutic system components (e.g., as plasmids, mRNA, or pre-assembled ribonucleoproteins - RNPs) into relevant cancer cell lines (e.g., cell lines derived from the specific cancer type with the target mutation).
                        * Assess on-target editing efficiency using methods like Sanger sequencing followed by TIDE/ICE analysis, or more quantitative next-generation sequencing (NGS)-based methods.
                        * If the goal is gene knockout, confirm reduced or absent protein expression using techniques like Western blot or flow cytometry (FACS).
                        * If the goal is gene correction or knock-in, assess the frequency of precise edits and potentially confirm expression of the corrected or inserted gene.
                    * **Experimental Off-Target Validation:**
                        * Once on-target activity is confirmed, rigorously assess off-target effects. Use unbiased genome-wide methods like GUIDE-seq, CIRCLE-seq, or SITE-seq to identify actual off-target cleavage sites in cells.
                        * Validate the most likely off-target sites identified by these methods (and potentially high-scoring sites from *in silico* analysis) using targeted deep sequencing to quantify editing frequency.
                    * **Functional Validation:**
                        * Assess the functional consequences of the editing in relevant cell models. This depends on the design goal:
                            * For oncogene disruption: Measure cell viability, proliferation (e.g., MTT, CellTiter-Glo, colony formation assays), and induction of apoptosis (e.g., Annexin V staining) to see if editing inhibits cancer cell growth or promotes cell death.
                            * For tumor suppressor restoration: Assess markers of restored protein function or downstream pathway activity.
                    """)

                    # Display Broader Therapeutic Considerations
                    st.markdown("### 🌍 Broader Therapeutic Considerations")
                    st.markdown("""
                    Translating a successful design from cell lines to a therapy involves additional critical considerations:
                    * **Delivery Strategy:** Define and develop a suitable method to deliver the multi-component system into the target cells/tissues in a patient. The choice of delivery vector (e.g., viral vectors like Adeno-Associated Virus (AAV) or lentivirus, non-viral methods like lipid nanoparticles (LNPs) or electroporation, or even bacterial delivery) is critical and has implications for packaging capacity, cell/tissue targeting (tropism), manufacturing, and potential immune responses.
                    * **Pre-clinical Studies:** If *in vitro* and cell line experiments are promising, the next major step is testing in preclinical animal models that accurately mimic the human disease (e.g., patient-derived xenografts (PDX), genetically engineered mouse models (GEMMs)). These studies assess efficacy in a living system, toxicity, biodistribution of the therapeutic, and optimal dosing.
                    * **Safety & Immunogenicity:** Evaluate the potential for the therapeutic system components (especially the bacterial-derived Cas enzyme) to trigger an undesirable immune response in patients. Develop strategies to mitigate immunogenicity if necessary (e.g., using humanized or engineered Cas variants, transient delivery methods).
                    * **Long-term Effects:** Evaluate the long-term durability of the therapeutic effect and monitor for potential late-onset off-target effects or other unintended consequences over time.
                    """)
                    
                    # Add Interactive AI Agents section
                    st.markdown("### 🤖 Interactive AI Agents")
                    st.markdown("""
                    The following interactive agents can help you simulate how AI-powered research assistants might 
                    optimize components of your design or validate its effectiveness. Select an agent to see a 
                    simulation of the process and results.
                    """)
                    
                    # Create columns for agent buttons
                    col1, col2, col3 = st.columns(3)
                    
                    # Function to create a buffer for agent output
                    def create_agent_output_area():
                        return st.empty()
                    
                    # Streamlit callback for agent display
                    def streamlit_display_callback(message):
                        if 'agent_output_area' in st.session_state:
                            current_content = st.session_state.agent_output_content
                            st.session_state.agent_output_content = current_content + "\n" + message
                            st.session_state.agent_output_area.markdown(f"```\n{st.session_state.agent_output_content}\n```")
                    
                    # Display area for agent output
                    agent_output_container = st.container()
                    
                    # Generate a mock recommendation if needed
                    # Use the current recommendation as input for the agents
                    current_recommendation = rec
                    
                    # Guide RNA Optimization Agent
                    with col1:
                        if st.button("Guide RNA Optimizer", key=f"guide_optimizer_{i}"):
                            # Create output area
                            with agent_output_container:
                                st.subheader("Guide RNA Optimization Agent")
                                st.session_state.agent_output_area = create_agent_output_area()
                                st.session_state.agent_output_content = "Starting Guide RNA optimization process..."
                                st.session_state.agent_output_area.markdown(f"```\n{st.session_state.agent_output_content}\n```")
                                
                                # Run the agent
                                result = run_interactive_agent(
                                    agent_type="guide_optimization",
                                    recommendation=current_recommendation,
                                    editing_type="knockout",  # Default to knockout, could be different based on design_goal
                                    display_callback=streamlit_display_callback
                                )
                                
                                # Show comparison of results
                                if result:
                                    st.write("### Results Comparison")
                                    
                                    col1, col2 = st.columns(2)
                                    with col1:
                                        st.write("**Original Metrics**")
                                        st.write(f"On-target Score: {result['changes']['original_metrics']['chopchop_on_target_score']:.2f}")
                                        st.write(f"Off-target Count: {result['changes']['original_metrics']['chopchop_off_target_count']}")
                                        st.write(f"MFE: {result['changes']['original_metrics']['chopchop_mfe']:.2f}")
                                    
                                    with col2:
                                        st.write("**Optimized Metrics**")
                                        st.write(f"On-target Score: {result['changes']['new_metrics']['chopchop_on_target_score']:.2f}")
                                        st.write(f"Off-target Count: {result['changes']['new_metrics']['chopchop_off_target_count']}")
                                        st.write(f"MFE: {result['changes']['new_metrics']['chopchop_mfe']:.2f}")
                                    
                                    st.write("### Explanation")
                                    st.write(result['explanation'])
                    
                    # Protein Structure Optimization Agent
                    with col2:
                        if st.button("Protein Component Optimizer", key=f"protein_optimizer_{i}"):
                            # Create output area
                            with agent_output_container:
                                st.subheader("Protein Component Optimization Agent")
                                st.session_state.agent_output_area = create_agent_output_area()
                                st.session_state.agent_output_content = "Starting Protein Component optimization process..."
                                st.session_state.agent_output_area.markdown(f"```\n{st.session_state.agent_output_content}\n```")
                                
                                # Run the agent
                                result = run_interactive_agent(
                                    agent_type="component_optimization",
                                    recommendation=current_recommendation,
                                    component_name="cas_protein_sequence",
                                    editing_type="knockout",  # Default to knockout
                                    display_callback=streamlit_display_callback
                                )
                                
                                # Show comparison of results
                                if result:
                                    st.write("### Results Comparison")
                                    
                                    col1, col2 = st.columns(2)
                                    with col1:
                                        st.write("**Original Metrics**")
                                        st.write(f"pLDDT Score: {result['changes']['original_metrics']['mock_plddt_score']:.2f}")
                                        st.write(f"Ranking Score: {result['changes']['original_metrics']['mock_ranking_score']:.2f}")
                                    
                                    with col2:
                                        st.write("**Optimized Metrics**")
                                        st.write(f"pLDDT Score: {result['changes']['new_metrics']['mock_plddt_score']:.2f}")
                                        st.write(f"Ranking Score: {result['changes']['new_metrics']['mock_ranking_score']:.2f}")
                                    
                                    st.write("### Explanation")
                                    st.write(result['explanation'])
                    
                    # In Vitro Validation Agent
                    with col3:
                        if st.button("In Vitro Validation", key=f"in_vitro_validation_{i}"):
                            # Create output area
                            with agent_output_container:
                                st.subheader("In Vitro Validation Agent")
                                st.session_state.agent_output_area = create_agent_output_area()
                                st.session_state.agent_output_content = "Starting In Vitro Validation process..."
                                st.session_state.agent_output_area.markdown(f"```\n{st.session_state.agent_output_content}\n```")
                                
                                # Run the agent
                                result = run_interactive_agent(
                                    agent_type="in_vitro_validation",
                                    recommendation=current_recommendation,
                                    editing_type="knockout",  # Default to knockout
                                    display_callback=streamlit_display_callback
                                )
                                
                                # Show results
                                if result:
                                    st.write("### Validation Results")
                                    
                                    # Create a progress bar for cleavage efficiency
                                    cleavage_efficiency = result['results']['cleavage_efficiency']
                                    st.write(f"**Cleavage Efficiency**: {cleavage_efficiency:.1%}")
                                    st.progress(cleavage_efficiency)
                                    
                                    # Show interpretation
                                    st.write("### Interpretation")
                                    st.write(result['results']['interpretation'])
                                    
                                    st.write("### Explanation")
                                    st.write(result['explanation'])


    # --- 5. Display AI Chat Assistant Interaction ---
    # This section can display the conversational output from the LLM if you integrate it
    # For this specific workflow log request, we are focusing on the structured log above
    # If you want to display the conversational history here, uncomment and use your history state
    # st.subheader("💬 AI Assistant Interaction Log")
    # chat_container = st.container()
    # with chat_container:
    #     for message in st.session_state.conversation_history:
    #         role = message.get("role", "user")
    #         content = message.get("content", "")
    #         with st.chat_message(role):
    #             st.markdown(content)


# --- Run the Streamlit Page ---
if __name__ == "__main__":
    display_mock_therapy_design_page()
