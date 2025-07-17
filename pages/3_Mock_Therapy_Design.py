import streamlit as st
import sys
import os
import time # Import time for potential delays if needed for visual effect
import random # Import random if needed for frontend mock data generation
import re

# Moved st.set_page_config to be the first Streamlit command
st.set_page_config(page_title="CasPro 🧬 - Mock Therapeutic System Design", layout="wide")

# Add the project root to the Python path to allow importing ai_research_assistant
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if PROJECT_ROOT not in sys.path:
    sys.path.append(PROJECT_ROOT)

try:
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
        simulate_therapeutic_development_planning,
        run_interactive_agent
    )
    # Ensure get_llm_chat_response is imported
    from tools.llm_api import get_llm_chat_response
except ImportError as e:
    st.error(f"Failed to import modules. Ensure `ai_research_assistant.py` and `tools/llm_api.py` are accessible and correct: {e}")
    # Fallback for critical functions if import fails
    def handle_design_therapy_interaction(*args, **kwargs):
        st.error("`handle_design_therapy_interaction` is not available due to import error.")
        return []
    def run_interactive_agent(*args, **kwargs):
        st.error("`run_interactive_agent` is not available due to import error.")
        return {"error": "Agent not available"}
    def get_llm_chat_response(*args, **kwargs):
        st.error("`get_llm_chat_response` is not available due to import error.")
        return "LLM service is currently unavailable."
    # Add fallbacks for simulate_... functions as needed

# --- Streamlit Page Function ---
def display_mock_therapy_design_page():
    # st.set_page_config was here, MOVED TO TOP

    # --- Set Page-Specific Context for LLM ---
    page_name = "Mock Therapeutic System Design"
    page_description = ("""
    This page allows you to simulate the design of a multi-gene CRISPR-based therapeutic system
    for a target cancer mutation. The process uses mock Evo 2, CHOPCHOP, and AlphaFold 3 data for demonstration.
    Key activities include defining a target mutation, running a simulated design workflow, and reviewing detailed recommendations.
    Recommendations include system components, supporting evidence (simulated Evo2, AlphaFold, CHOPCHOP metrics),
    evaluation metrics, and suggested next steps with interactive AI agent simulations.
    All outputs, scores, and rationales are SIMULATED.
    """)
    active_target_display = "Not set"
    if 'active_mutation' in st.session_state and st.session_state.active_mutation:
        active_target_display = f"{st.session_state.active_mutation.get('hugo_gene_symbol', '?')} {st.session_state.active_mutation.get('protein_change', '?')}"

    st.session_state.current_page_info = {
        "name": page_name,
        "description": page_description,
        "current_target": active_target_display
    }
    # --- End Set Page-Specific Context ---

    st.title("CasPro 🧬 Mock Multi-Gene Therapeutic System Design")
    st.markdown("""
    This page allows you to simulate the design of a multi-gene CRISPR-based therapeutic system
    for a target cancer mutation. The process uses mock Evo 2, CHOPCHOP, and AlphaFold 3 data for demonstration.
    **All outputs, scores, and rationales are SIMULATED.**
    """)

    # --- New Expander for Workflow Explanation ---
    with st.expander("⚙️ What Happens When You Click 'Design Mock Therapy System'?", expanded=False):
        st.markdown("""
        When you click the **'Design Mock Therapy System'** button, the system initiates our mock multi-gene therapy system design workflow. This is a simulation of how our platform might leverage powerful AI models (conceptually similar to Evo2 and AlphaFold 3) to propose and evaluate complex genetic constructs for potential treatment.

        Here's a breakdown of what's happening behind the scenes in this mock scenario:

        1.  **Identifying the Target:**
            *   The system first focuses on your selected target mutation (e.g., BRAF V600E).
            *   It retrieves some *simulated Evo2 data* for this original variant, like a mock confidence score and delta likelihood. This mock data conceptually helps confirm if the selected mutation is a relevant and impactful therapeutic target.

        2.  **Simulating AI-Guided Generation (Mock Evo2):**
            *   Next, our `MockEvo2Generator` component (a placeholder for a more sophisticated AI) kicks in.
            *   This simulates the capability of a real generative model like Evo2, which could be trained to design multi-gene systems.
            *   Based on your target mutation and a mock *design goal* (e.g., 'disrupt oncogene', 'correct mutation'), this generator proposes several *candidate multi-gene systems*.
            *   Each system isn't just a single sequence but a collection of *mock sequences* representing different components potentially needed for a therapy – such as a guide RNA sequence, a modified Cas protein sequence, a repair template, or even elements for a delivery vector.

        3.  **Simulating Structural & Biological Evaluation (Mock AlphaFold 3 & Mock CHOPCHOP):**
            *   For each of these generated system candidates, the platform then simulates running evaluations using tools conceptually similar to AlphaFold 3 (for structure) and CHOPCHOP (for guide RNA validation):
                *   **Mock AlphaFold 3 Analysis (`MockAlphaFoldPredictor`):** This simulates AlphaFold 3 predicting the 3D structure and interactions of the components within the designed system. It generates mock metrics such as:
                    *   `pLDDT` (simulating confidence in local protein/RNA structure).
                    *   `Ranking Scores` (simulating overall confidence in the predicted structure or interactions).
                    *   `Structural Similarity` (simulating how the designed components compare to known reference structures).
                    This step is conceptually crucial because the 3D shape and how components fit together directly impacts whether a real therapeutic system would function correctly.
                *   **Mock CHOPCHOP Analysis (for guide RNAs):** For components like guide RNAs, we also simulate validation using metrics similar to those from CHOPCHOP. This includes:
                    *   Predicted *on-target efficacy scores*.
                    *   Predicted *off-target counts*.
                    *   *Minimum Free Energy (MFE)* for guide RNA stability.
                    This is vital for ensuring the specificity (hitting only the intended target) and safety of a CRISPR-based therapy.

        4.  **Simulating Integrated Scoring:**
            *   All these varied simulated metrics – the mock Evo2 data for the original variant, the mock AlphaFold 3 structural scores for the designed system's components and overall structure, and the mock CHOPCHOP guide validation scores – are fed into our `calculate_mock_system_design_score` function.
            *   This function simulates a comprehensive scoring logic, combining these different types of evidence into a single *mock confidence score* for each entire system candidate.
            *   This score represents the simulated likelihood of that designed system achieving the therapeutic goal, based on all the available mock data.

        5.  **Presenting Mock Recommendations:**
            *   Finally, the system uses an LLM (Large Language Model) to present the top-scoring mock system designs to you.
            *   It shows the recommended approach (the description of the mock system), the simulated confidence score, a mock rationale explaining why that score was given, and the detailed simulated supporting evidence from the mock Evo2, AlphaFold 3, and CHOPCHOP evaluations.

        **Key Takeaway:** This entire process is a *simulation* designed to mimic how an advanced AI-driven platform might approach therapeutic design. The data is mock, but the concepts represent real-world considerations in genetic medicine development.
        """)
    # --- End New Expander ---

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
    gene_symbol = st.sidebar.text_input("Gene Symbol", value=st.session_state.active_mutation.get("hugo_gene_symbol", "BRAF"))
    protein_change = st.sidebar.text_input("Protein Change", value=st.session_state.active_mutation.get("protein_change", "V600E"))

    # --- CRISPR Education Center (placeholder, can be expanded or moved) ---
    # st.sidebar.subheader("📚 CRISPR Education Center")
    # with st.sidebar.expander("🧬 CRISPR Terminology"):
    #     st.markdown("Explore key CRISPR terms and definitions.")
    # with st.sidebar.expander("💡 Core Concepts"):
    #     st.markdown("Understand the fundamental principles of CRISPR technology.")
    # with st.sidebar.expander("🤖 Ask AI Chat Assistant (Edu Center)"):
    #     st.markdown("Have questions? Ask our AI assistant for help.")
    #     if st.sidebar.button("Open Chat (Edu)"):
    #         st.info("Chat functionality to be implemented here.")
    # --- End CRISPR Education Center ---

    # Update active_mutation if changed by user
    current_target_mutation = {
        "hugo_gene_symbol": gene_symbol,
        "protein_change": protein_change,
        "variant_type": st.session_state.active_mutation.get("variant_type", "Missense_Mutation"),
        "genomic_coordinate_hg38": st.session_state.active_mutation.get("genomic_coordinate_hg38", ""),
        "allele_frequency": st.session_state.active_mutation.get("allele_frequency", 0.0),
        "mutation_id": st.session_state.active_mutation.get("mutation_id", "DEFAULT_MUT")
    }
    if (st.session_state.active_mutation.get("hugo_gene_symbol") != gene_symbol or
        st.session_state.active_mutation.get("protein_change") != protein_change):
        st.session_state.active_mutation = current_target_mutation
        st.session_state.therapy_recommendations = []
        st.session_state.workflow_log = []
        st.rerun()

    if current_target_mutation['hugo_gene_symbol'] == "BRAF" and current_target_mutation['protein_change'] == "V600E":
        with st.sidebar.expander("Understanding BRAF V600E Context", expanded=False):
            st.markdown("""
            **Target: BRAF V600E**
            *   **BRAF Gene:** A key gene involved in normal cell growth signaling pathways (MAPK pathway).
            *   **V600E Mutation:** A specific alteration where Valine (V) at amino acid position 600 is replaced by Glutamic acid (E). This is an **activating mutation**, meaning it causes the BRAF protein to be constantly "on," leading to uncontrolled cell division and tumor growth. It's a well-known **oncogenic driver** found in various cancers, most notably melanoma (around 50% of cases), but also in colorectal cancer, thyroid cancer, and others.
            *   **Therapeutic Goal (Simulated on this page):** The simulation aims to design a CRISPR-based system to counteract this cancer-causing mutation. Common strategies simulated could include:
                *   **Disrupting/Knocking out** the mutated BRAF V600E allele.
                *   **Correcting** the V600E mutation back to the wild-type sequence via Homology Directed Repair (HDR) using a repair template.
                *   Developing a delivery system for these components.

            **Why is this context important for the "Mock Therapeutic System Design" page?**
            Understanding the specific target (BRAF V600E) helps interpret the simulated design choices and evaluation metrics. For example:
            *   High `Mock Original Variant Evo2 Confidence` for BRAF V600E would reinforce its validity as a therapeutic target.
            *   If the design goal is 'disrupt_oncogene', the system components (guide RNA, Cas) must effectively target and inactivate BRAF V600E.
            *   If 'correct_mutation' is the goal, the system would also include a repair template, and its structural properties (e.g., from mock AlphaFold) become important.
            *   Off-target scores for the guide RNA are critical, as editing the normal BRAF allele or other genes could have significant side effects.

            This page simulates the *in silico* (computational) design and initial evaluation phase for such a complex therapeutic.
            """)

    # --- Ask About CRISPR section (MOVED TO SIDEBAR) ---
    with st.sidebar:
        st.subheader("💬 Ask About CRISPR")

        if 'chat_question' not in st.session_state:
            st.session_state.chat_question = ""
        if 'chat_answer' not in st.session_state:
            st.session_state.chat_answer = ""

        user_question_sidebar = st.text_input("Ask a question:", key="crispr_chat_input_sidebar")

        if st.button("Submit Question", key="crispr_chat_submit_sidebar"):
            if user_question_sidebar:
                st.session_state.chat_question = user_question_sidebar
                try:
                    # --- Construct detailed system prompt ---
                    page_info = st.session_state.get('current_page_info', {})
                    page_name = page_info.get('name', 'Mock Therapeutic System Design')

                    active_mut = st.session_state.get('active_mutation', {})
                    target_gene = active_mut.get('hugo_gene_symbol', 'N/A')
                    protein_change = active_mut.get('protein_change', 'N/A')
                    active_target_display = f"{target_gene} {protein_change}"

                    workflow_status = "No therapy design has been run yet."
                    if st.session_state.get('therapy_recommendations'):
                        workflow_status = f"Mock therapy system recommendations are currently displayed for {active_target_display}."
                    
                    therapeutic_context_toggle_status = "Enabled" if st.session_state.get('provide_therapeutic_context', False) else "Disabled"

                    # Corrected f-string definition for system_prompt
                    system_prompt = f"""
System: You are an AI assistant for CasPro, a CRISPR therapeutic design simulation tool.
User is currently on the '{page_name}' page.

**Current Page Context & Workflow Overview:**

This page simulates the initial steps of designing a CRISPR-based therapy for a target cancer mutation, using a series of mock tools and analyses. The goal is to make the complex process of early-stage therapeutic design more tangible, even if all the underlying science is a simulation.

**Current Active Target Mutation:**
- Gene: {target_gene}
- Protein Change: {protein_change}

**Simulated Workflow Stages for Designing a Therapy for the Active Target:**

1.  **Target Confirmation (Simulated Evo 2):**
    -   Input: The active target mutation.
    -   Process: The system (mock) analyzes this mutation to confirm its known oncogenic potential. Metrics like 'Mock Original Variant Evo2 Confidence' and 'Delta Likelihood' are simulated outputs.
    -   Interpretation: High confidence reinforces that the target is a strong candidate for therapeutic intervention.

2.  **Therapeutic Strategy & System Component Generation (Simulated Evo 2):**
    -   Process: Based on the target and a chosen therapeutic strategy (e.g., 'disrupt oncogene'), the system (mock) generates candidate components for the CRISPR therapy (e.g., mock sequences for guide RNA, Cas protein, potentially repair templates).

3.  **Guide RNA Validation (Simulated CHOPCHOP):**
    -   Process: The generated mock guide RNA sequence is (mock) evaluated for potential effectiveness and safety.
    -   Output Metrics: Mock CHOPCHOP On-target Score, Off-target Count, MFE (Minimum Free Energy).

4.  **Structural Evaluation of Components (Simulated AlphaFold 3):**
    -   Process: The 3D structures of protein/RNA components are (mock) predicted and evaluated.
    -   Output Metrics: Mock Overall System AF Score, individual component scores (pLDDT, Rank, Similarity).

5.  **Integrated Scoring & Recommendation:**
    -   Process: All simulated metrics are combined to give an overall 'Confidence Score' for each complete therapeutic system candidate.
    -   Output: The page presents these mock recommendations, ranked by simulated confidence, with explanations of what the mock numbers conceptually represent.

**Current Workflow Status on Page:**
- {workflow_status}

**Therapeutic Context Mode for Explanations:**
- {therapeutic_context_toggle_status}

Based on this detailed context and your general CRISPR knowledge, please answer the user's question.
"""
                    # --- End detailed system prompt ---

                    messages = [
                        {"role": "system", "content": system_prompt},
                        {"role": "user", "content": user_question_sidebar}
                    ]
                    llm_response = get_llm_chat_response(messages, provider="gemini") 
                    st.session_state.chat_answer = llm_response
                except Exception as e:
                    st.session_state.chat_answer = f"Error getting response from LLM: {e}"
                    st.error(f"LLM Error: {e}")
            else:
                st.session_state.chat_answer = "Please ask a question."

        if st.session_state.chat_question:
            st.markdown(f"**You asked:** {st.session_state.chat_question}")
        if st.session_state.chat_answer:
            with st.chat_message("assistant"):
                st.markdown(st.session_state.chat_answer)
    # --- END: Ask About CRISPR section in sidebar ---

    # --- 3. Trigger Button & Workflow Execution ---
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
            st.session_state.workflow_log.append(message)
            with log_placeholder.container():
                st.subheader("⚙️ Workflow Log (Live Updates)")
                # Use a new container with a fixed height and make it scrollable
                with st.container(height=350): # Adjust height as needed
                    # Iterate in reverse to show newest messages at the top in the live view
                    for log_message in reversed(st.session_state.workflow_log):
                        if "ERROR:" in log_message:
                            st.error(log_message)
                        elif "WARNING:" in log_message or "Warning:" in log_message:
                            st.warning(log_message)
                        elif "SUCCESS:" in log_message or "completed!" in log_message.lower() or "finished generating" in log_message.lower() or "prediction complete" in log_message.lower() or "validation complete" in log_message.lower() or "scoring complete" in log_message.lower():
                            st.success(log_message)
                        elif "AI Research Assistant (Scoring)" in log_message:
                            st.markdown(f"📊 _{log_message}_")
                        elif "AI Research Assistant (ChopChopSim)" in log_message:
                            st.markdown(f"✂️ _{log_message}_")
                        elif "AI Research Assistant (AlphaFoldSim)" in log_message:
                            st.markdown(f"🧬 _{log_message}_")
                        elif "AI Research Assistant (Evo2Sim)" in log_message:
                            st.markdown(f"🌱 _{log_message}_")
                        elif "Simulating" in log_message or "Processing" in log_message or "Retrieving" in log_message or "Focusing" in log_message or "Initiating" in log_message:
                            st.markdown(f"⏳ {log_message}")
                        elif "Generated mock" in log_message or "details:" in log_message or "Inputs for scoring:" in log_message:
                            st.markdown(f"📄 _{log_message}_")
                        else:
                            st.text(log_message) # Default for messages not caught by specific styling

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

    # --- Display the Full Workflow Log in an Expander ---
    # This section is primarily for showing the log after a rerun or if the user expands it later
    if st.session_state.workflow_log: # Show if there's any log content
         with st.expander("📜 View Full Workflow Log", expanded=False): # Default to collapsed
              # Display messages in chronological order for the full log
              for log_message in st.session_state.workflow_log:
                  if "ERROR:" in log_message:
                      st.error(log_message)
                  elif "WARNING:" in log_message or "Warning:" in log_message :
                      st.warning(log_message)
                  elif "SUCCESS:" in log_message or "completed!" in log_message.lower() or "finished generating" in log_message.lower() or "prediction complete" in log_message.lower() or "validation complete" in log_message.lower() or "scoring complete" in log_message.lower():
                      st.success(log_message)
                  elif "AI Research Assistant (Scoring)" in log_message:
                      st.markdown(f"📊 _{log_message}_")
                  elif "AI Research Assistant (ChopChopSim)" in log_message:
                      st.markdown(f"✂️ _{log_message}_")
                  elif "AI Research Assistant (AlphaFoldSim)" in log_message:
                      st.markdown(f"🧬 _{log_message}_")
                  elif "AI Research Assistant (Evo2Sim)" in log_message:
                      st.markdown(f"🌱 _{log_message}_")
                  elif "Simulating" in log_message or "Processing" in log_message or "Retrieving" in log_message or "Focusing" in log_message or "Initiating" in log_message:
                      st.markdown(f"⏳ {log_message}")
                  elif "Generated mock" in log_message or "details:" in log_message or "Inputs for scoring:" in log_message:
                      st.markdown(f"📄 _{log_message}_")
                  else:
                      st.text(log_message) # Default


    # --- 4. Display Results ---
    if st.session_state.therapy_recommendations:
        st.subheader("📊 Simulated System Recommendations")

        for i, rec in enumerate(st.session_state.therapy_recommendations):
            approach_text = rec.get('recommended_approach', f'Recommendation {i+1}')
            confidence_score = rec.get('confidence_score', 'N/A')

            with st.expander(f"Recommendation {i+1}: {approach_text} (Simulated Confidence: {confidence_score:.2f})", expanded=(i == 0)):
                st.markdown(f"**Detailed Rationale (Simulated):** {rec.get('detailed_rationale', 'Detailed rationale not available.')}")
                st.markdown("---") # Separator after rationale

                tab_evidence, tab_components, tab_metrics, tab_next_steps = st.tabs(["Supporting Evidence", "🧬 System Components", "🔬 Evaluation Metrics", "➡️ Suggested Next Steps"])

                with tab_evidence:
                    st.markdown("### 📊 Supporting Evidence (Simulated Data)")
                    st.caption("Below are key simulated metrics that contribute to this system candidate\'s evaluation. Higher scores are generally better, unless specified (e.g., off-target counts). Refer to the \'Understanding These Mock Results\' expander for general metric definitions.")
                    st.divider()

                    # --- Contextual Explanation for BRAF V600E ---
                    # active_gene = st.session_state.active_mutation.get("hugo_gene_symbol", "BRAF")
                    # active_protein_change = st.session_state.active_mutation.get("protein_change", "V600E")
                    # if active_gene == "BRAF" and active_protein_change == "V600E":
                    #     st.markdown(\"""
                    #         #### Interpreting Evidence for BRAF V600E Therapeutic Design:
                    #
                    #         You are currently evaluating a mock therapeutic system designed to target the **BRAF V600E** mutation, a common driver in various cancers, including melanoma. The goal of such a therapy is often to **disrupt or correct the oncogenic activity** of this mutated protein.
                    #
                    #         As you review the metrics below:
                    #         - **High `Mock Original Variant Evo2 Confidence` and significant `Mock Original Variant Delta Likelihood`** would confirm that BRAF V600E itself is a strong, impactful target, justifying intervention.
                    #         - **`Mock Overall System AF Score` and individual component pLDDT/Rank scores** (e.g., for the Cas protein and guide RNA) are crucial. For a *disruption* strategy (like knockout), the guide RNA and Cas protein must form a functional complex to edit the gene. For a *correction* strategy, the repair template also needs to be structurally sound and interact correctly.
                    #         - **`Mock CHOPCHOP On-target Score`** should be high, indicating the guide RNA is likely to cut efficiently at the BRAF V600E locus.
                    #         - **`Mock CHOPCHOP Off-target Count`** should be ideally zero, or very low. Off-target edits are a major safety concern, especially for a gene like BRAF involved in normal cellular processes.
                    #         - **`Mock CHOPCHOP MFE`** (Minimum Free Energy) for the guide RNA suggests structural stability, which is important for its function with the Cas protein.
                    #
                    #         The mock component similarities (e.g., for the Cas protein) are less critical for a knockout strategy but would be important if aiming to restore a wild-type function or if comparing to a known functional structure. For BRAF V600E, if a repair template is used for correction, its similarity to the wild-type sequence in the corrected region would be paramount.
                    #     \"\"\")
                    #     st.divider()
                    # --- End Contextual Explanation ---
                    supporting_evidence = rec.get('supporting_evidence', [])
                    if supporting_evidence:
                        # Create two columns for the grid layout
                        col1, col2 = st.columns(2)
                        
                        # Distribute items between columns
                        for index, evidence_item_full_string in enumerate(supporting_evidence):
                            target_column = col1 if index % 2 == 0 else col2 # Alternate columns
                            with target_column:
                                parts = evidence_item_full_string.split(':')
                                metric_name = parts[0].strip()
                                metric_value_str = parts[1].strip() if len(parts) > 1 else "N/A"
                                
                                explanation = "_No specific explanation available for this exact metric line-item, but it contributes to the overall assessment._"
                                emoji = "ℹ️"
                                value_to_display = metric_value_str # Default to string
                                numeric_val_for_assessment = None

                                try:
                                    numeric_val = float(metric_value_str)
                                    numeric_val_for_assessment = numeric_val # Store for assessment
                                    value_to_display = f"{numeric_val:.2f}"
                                    if metric_name == "Mock Original Variant Delta Likelihood":
                                        value_to_display = f"{numeric_val:.7f}"
                                    elif "Count" in metric_name or "Off-target" in metric_name:
                                        value_to_display = int(numeric_val)
                                except ValueError:
                                    pass

                                st.markdown(f"##### {metric_name}")

                                if "Mock Original Variant Evo2 Confidence" == metric_name:
                                    emoji = "🌱"
                                    explanation = ("_Simulated Evo 2 score indicating confidence in the **original target mutation\\'s** biological impact (e.g., pathogenicity). "
                                                 "Higher scores suggest the mutation is a relevant therapeutic target._")
                                elif "Mock Original Variant Delta Likelihood" == metric_name:
                                    emoji = "📉"
                                    explanation = ("_Simulated Evo 2 score predicting how the **original mutation** changes sequence likelihood. Large negative values can indicate impactful variants "
                                                 "(e.g., disease-causing) that are selected against._")
                                elif "Mock Overall System AF Score" == metric_name:
                                    emoji = "🏗️"
                                    explanation = ("_Simulated AlphaFold 3 aggregate score for the **entire system\\'s** structural integrity. Higher scores suggest components are likely to fold and interact correctly, which is crucial for function._")
                                elif metric_name.startswith("Comp ") and " pLDDT" in metric_name:
                                    emoji = "🧬"
                                    comp_actual_name = metric_name.split(" pLDDT")[0].replace("Comp ","").strip("'")
                                    explanation = (f"_Simulated AlphaFold pLDDT score (0-100) for the **{comp_actual_name}** component. Indicates confidence in the local accuracy of its predicted 3D structure. "
                                                 f"Scores >70 are generally good; >90 suggest high confidence in the fold._")
                                elif metric_name.startswith("Comp ") and " Rank" in metric_name:
                                    emoji = "🏅"
                                    comp_actual_name = metric_name.split(" Rank")[0].replace("Comp ","").strip("'")
                                    explanation = (f"_Simulated AlphaFold ranking score for the **{comp_actual_name}** component. This can reflect overall structural quality or confidence in predicted interactions. Higher is generally better._")
                                elif metric_name.startswith("Comp ") and " Similarity" in metric_name:
                                    emoji = "🔗"
                                    comp_actual_name = metric_name.split(" Similarity")[0].replace("Comp ","").strip("'")
                                    explanation = (f"_Simulated AlphaFold structural similarity for the **{comp_actual_name}** to a reference. High similarity is key if restoring function (e.g., to wild-type). For disruption, lower similarity to an active form might be expected._")
                                elif "Mock CHOPCHOP On-target Score" == metric_name:
                                    emoji = "🎯"
                                    explanation = ("_Simulated CHOPCHOP score for the **guide RNA\\'s** predicted on-target efficacy. Considers sequence features like GC content and PAM proximity. Higher scores suggest a greater likelihood of successful DNA cleavage at the intended site._")
                                elif "Mock CHOPCHOP Off-target Count" == metric_name:
                                    emoji = "⚠️"
                                    explanation = ("_Simulated CHOPCHOP count for the **guide RNA\\'s** predicted off-target sites. Lower counts are crucial for specificity and safety, indicating fewer unintended edits._")
                                elif "Mock CHOPCHOP MFE" == metric_name:
                                    emoji = "➰"
                                    explanation = ("_Simulated CHOPCHOP Minimum Free Energy (MFE) for the **guide RNA\\'s** predicted secondary structure. More negative values indicate a more stable structure, important for Cas protein interaction and DNA binding._")
                                
                                st.metric(label=emoji, value=value_to_display)
                                st.caption(explanation)

                                # --- Add Conditional Assessment ---
                                assessment_text = ""
                                assessment_emoji = ""
                                design_goal = rec.get('editing_type', '').lower() # Example: "Mock System Design (disrupt_oncogene)"
                                
                                if numeric_val_for_assessment is not None:
                                    # Default thresholds
                                    good_threshold = 0.7
                                    bad_threshold = 0.4
                                    is_lower_better = False

                                    if "Mock Original Variant Evo2 Confidence" == metric_name:
                                        good_threshold = 0.7; bad_threshold = 0.3
                                    elif "Mock Original Variant Delta Likelihood" == metric_name:
                                        # Larger negative is more impactful, so "good" means more negative
                                        # For simplicity, we'll treat it as higher absolute value is good if negative
                                        good_threshold = -0.001; bad_threshold = -0.0001 
                                        # This logic is tricky with "good" for more negative. Let's simplify for now.
                                        # Or we can check if it's very negative.
                                        if numeric_val_for_assessment < -0.0005: assessment_emoji = "👍"; assessment_text = "Potentially impactful variant signal."
                                        else: assessment_emoji = "🤔"; assessment_text = "Less clear impact signal."
                                    elif "Mock Overall System AF Score" == metric_name:
                                        good_threshold = 0.75; bad_threshold = 0.5
                                    elif "pLDDT" in metric_name:
                                        good_threshold = 80; bad_threshold = 60 # pLDDT is 0-100
                                    elif "Rank" in metric_name: # Assuming higher rank is better
                                        good_threshold = 0.8; bad_threshold = 0.5 
                                    elif "Similarity" in metric_name:
                                        if "disrupt" in design_goal: # For disruption, lower similarity to active form might be okay/good
                                            good_threshold = 0.3; bad_threshold = 0.6 # Inverted logic example
                                            is_lower_better = True 
                                        else: # For correction, higher similarity to WT is good
                                            good_threshold = 0.7; bad_threshold = 0.4
                                    elif "Mock CHOPCHOP On-target Score" == metric_name:
                                        good_threshold = 0.8; bad_threshold = 0.5
                                    elif "Mock CHOPCHOP Off-target Count" == metric_name:
                                        good_threshold = 1; bad_threshold = 5 # Lower is better
                                        is_lower_better = True
                                    elif "Mock CHOPCHOP MFE" == metric_name:
                                        good_threshold = -20; bad_threshold = -10 # More negative is better
                                        # This needs careful handling: lower MFE (more negative) is better.
                                        if numeric_val_for_assessment <= good_threshold: assessment_emoji = "👍"; assessment_text = "Favorable (stable structure)."
                                        elif numeric_val_for_assessment > bad_threshold: assessment_emoji = "👎"; assessment_text = "Less favorable (potentially unstable)."
                                        else: assessment_emoji = "🤔"; assessment_text = "Neutral/acceptable stability."
                                    
                                    # Apply general threshold logic if not handled by specific MFE/DeltaLikelihood
                                    if not assessment_text: # If not already set by specific logic
                                        if is_lower_better:
                                            if numeric_val_for_assessment <= good_threshold: assessment_emoji = "👍"; assessment_text = "Good (meets/exceeds target for lower-is-better)."
                                            elif numeric_val_for_assessment > bad_threshold: assessment_emoji = "👎"; assessment_text = "Potentially concerning (does not meet target for lower-is-better)."
                                            else: assessment_emoji = "🤔"; assessment_text = "Acceptable/Neutral."
                                        else: # This corresponds to the 'if is_lower_better:'
                                            if numeric_val_for_assessment >= good_threshold: assessment_emoji = "👍"; assessment_text = "Good (meets/exceeds target)."
                                            elif numeric_val_for_assessment < bad_threshold: assessment_emoji = "👎"; assessment_text = "Potentially concerning (below target)."
                                            else: assessment_emoji = "🤔"; assessment_text = "Acceptable/Neutral."
                                # This 'else' corresponds to 'if not assessment_text:'
                                # It also handles the case where metric_value is not numeric initially
                                else: # Non-numeric metric OR assessment_text was already set by specific logic
                                    # If assessment_text is already set, we don't overwrite it.
                                    # If it's not set AND it's non-numeric, then it's a qualitative point.
                                    if not assessment_text:
                                        assessment_emoji = "ℹ️"
                                        assessment_text = "Qualitative data point."
                                
                                if assessment_text: # Ensure assessment_text has a value before displaying
                                    st.markdown(f"**Assessment:** {assessment_emoji} _{assessment_text}_")
                                # --- End Conditional Assessment ---
                                st.divider() 
                    else:
                        st.markdown("  _No supporting evidence provided for this candidate._")

                with tab_components:
                    st.markdown("### 🧬 System Components (Mock)")
                    st.markdown(
                        "_These are the core biological parts that make up the proposed therapeutic system. "
                        "Each plays a distinct role in achieving the desired genetic modification for "
                        f"targeting **{current_target_mutation['hugo_gene_symbol']} "
                        f"{current_target_mutation['protein_change']}**. "
                        "The explanations below detail each mock component and its significance in this context._"
                    )
                    st.divider()

                    components = rec.get('designed_system_components', {})
                    design_goal_raw = rec.get('editing_type', 'unknown strategy') # e.g., "Mock System Design (disrupt_oncogene)"
                    strategy_match = re.search(r'\\((.*?)\\)', design_goal_raw)
                    design_goal = strategy_match.group(1).strip() if strategy_match else "disrupt_oncogene" # Default or extracted

                    if components:
                        # --- System ID and Description ---
                        if "system_id" in components:
                            st.markdown(f"#### System Identifier: `{components['system_id']}`")
                            st.caption(f"_disrupt_oncogene: This is the chosen therapeutic strategy. In this case, the simulated goal is to disrupt the function of the BRAF oncogene (which is driving cancer growth due to the V600E mutation). (e.g., {design_goal}) for **{current_target_mutation['hugo_gene_symbol']} {current_target_mutation['protein_change']}**._")
                            st.markdown("---")
                        if "description" in components:
                            st.markdown(f"#### System Description: {components['description']}")
                            st.caption(f"_A concise summary of this system's intended purpose and target. For example, this system is designed to '{design_goal}' related to the **{current_target_mutation['hugo_gene_symbol']} {current_target_mutation['protein_change']}** mutation. Clear descriptions aid in collaborative efforts and maintaining focus on the therapeutic objective._")
                            st.markdown("---")

                        # --- Core CRISPR Components ---
                        st.markdown("#### Key Therapeutic Sequences (Mock):")
                        component_sequences_found = False

                        if "guide_rna_sequence" in components:
                            st.markdown(f"##### Guide RNA (gRNA) Sequence: `{components['guide_rna_sequence'][:60]}...`")
                            st.caption(f"_This mock sequence represents the gRNA, the component responsible for **targeting specificity**. It contains a crucial 'spacer' region (typically ~20 nucleotides) that is designed to be complementary to a specific DNA sequence within or near the **{current_target_mutation['hugo_gene_symbol']}** gene, specifically at the site relevant to the **{current_target_mutation['protein_change']}** mutation. Its precise sequence dictates where the Cas enzyme will be directed to make a cut. For a '{design_goal}' strategy, this targeting is paramount for efficacy and safety._")
                            component_sequences_found = True
                            st.markdown(" ") # Add a little space

                        if "cas_protein_sequence" in components:
                            st.markdown(f"##### Cas Protein Sequence: `{components['cas_protein_sequence'][:60]}...`")
                            st.caption(f"_This mock sequence represents the Cas enzyme (e.g., Cas9, Cas12). This protein acts as the 'molecular scissors' that performs the DNA cut once guided by the gRNA. The choice of Cas protein can influence editing efficiency, the type of cut made (e.g., double-strand break for NHEJ-mediated '{design_goal}', or nick for base editing), and PAM sequence requirements. The sequence itself can be optimized for better activity, stability, or reduced immunogenicity._")
                            component_sequences_found = True
                            st.markdown(" ")

                        if "repair_template_dna" in components:
                            st.markdown(f"##### Repair Template DNA Sequence: `{components['repair_template_dna'][:60]}...`")
                            st.caption(f"_This mock sequence is for a DNA repair template, essential if the therapeutic goal is precise gene editing via Homology Directed Repair (HDR), such as correcting the **{current_target_mutation['protein_change']}** mutation in **{current_target_mutation['hugo_gene_symbol']}**. It contains the desired corrected sequence flanked by 'homology arms' that match the DNA sequence around the target site, guiding the cell's repair machinery to incorporate the change accurately. This component is critical for strategies like 'correct_mutation' or 'insert_gene'._")
                            component_sequences_found = True
                            st.markdown(" ")

                        if "modified_cas_protein_sequence" in components:
                            st.markdown(f"##### Modified Cas Protein Sequence: `{components['modified_cas_protein_sequence'][:60]}...`")
                            st.caption(f"_This mock sequence represents an engineered Cas protein. Modifications can include: base editors (to make C>T or A>G changes without a full double-strand break), prime editors (for more complex edits like small insertions/deletions or all base conversions), nickases (to cut only one DNA strand, often used with paired guides for increased specificity), or dCas9 (catalytically inactive Cas9 used for gene activation/repression by fusing it to effector domains). The specific modification would be chosen based on the '{design_goal}' for **{current_target_mutation['hugo_gene_symbol']}**._")
                            component_sequences_found = True
                            st.markdown(" ")
                        
                        # --- Delivery & Expression Components (if applicable) ---
                        delivery_components_present = any(k in components for k in ["viral_packaging_signal_dna", "therapeutic_gene_sequence", "promoter_sequence"])
                        if delivery_components_present:
                            st.markdown("#### Delivery & Expression Elements (Mock):")

                        if "viral_packaging_signal_dna" in components:
                            st.markdown(f"##### Viral Packaging Signal DNA: `{components['viral_packaging_signal_dna'][:60]}...`")
                            st.caption("_This mock DNA sequence is essential for viral vector-based delivery systems (e.g., AAV, Lentivirus). It enables the therapeutic payload (containing the CRISPR components) to be encapsulated within viral particles for efficient delivery to target cells/tissues. The choice of signal depends on the viral vector used._")
                            component_sequences_found = True
                            st.markdown(" ")

                        if "therapeutic_gene_sequence" in components:
                             st.markdown(f"##### Therapeutic Gene Sequence (Payload): `{components['therapeutic_gene_sequence'][:60]}...`")
                             st.caption(f"_This mock sequence represents the actual genetic material intended to exert the therapeutic effect, distinct from the CRISPR machinery itself. This could be a corrected version of **{current_target_mutation['hugo_gene_symbol']}**, a cDNA for a missing protein, or a shRNA. This is relevant for 'gene_addition' or 'gene_replacement' strategies._")
                             component_sequences_found = True
                             st.markdown(" ")

                        if "promoter_sequence" in components:
                             st.markdown(f"##### Promoter Sequence: `{components['promoter_sequence'][:60]}...`")
                             st.caption(f"_This mock DNA sequence controls the expression of the linked therapeutic genes (e.g., Cas enzyme, gRNA, or a therapeutic payload gene). The choice of promoter (e.g., constitutive, tissue-specific, inducible) is critical for ensuring the CRISPR system or therapeutic gene is active in the correct cells and at appropriate levels, which is vital for both efficacy and safety when targeting **{current_target_mutation['hugo_gene_symbol']}**._")
                             component_sequences_found = True
                             st.markdown(" ")

                        # --- Other/Generic Components ---
                        st.markdown("#### Other System Components (Mock):")
                        other_comp_count = 0
                        for comp_name, comp_val in components.items():
                            # Filter out already handled components
                            if comp_name not in ["system_id", "description", "guide_rna_sequence",
                                                "cas_protein_sequence", "repair_template_dna",
                                                "modified_cas_protein_sequence", "viral_packaging_signal_dna",
                                                "therapeutic_gene_sequence", "promoter_sequence"]:
                                other_comp_count += 1
                                st.markdown(f"##### {comp_name.replace('_', ' ').title()}: `{str(comp_val)[:60]}...`")
                                st.caption(f"_This mock sequence represents an additional component, '{comp_name.replace('_', ' ')}'. In a real system, this could be a linker, a tag for purification/detection, a regulatory element, or another part of a multi-component delivery vehicle. Its specific role would be defined by the overall system architecture and therapeutic strategy for **{current_target_mutation['hugo_gene_symbol']} {current_target_mutation['protein_change']}**._")
                                st.markdown(" ")
                        
                        if other_comp_count == 0 and not component_sequences_found and not delivery_components_present:
                             st.markdown("  _No specific component sequences listed for this candidate beyond ID and description._")
                        elif other_comp_count == 0 and (component_sequences_found or delivery_components_present):
                            st.caption("  _All identified mock components listed above._")
                        elif other_comp_count > 0:
                             st.caption("  _The components above are mock placeholders for other necessary parts of a complex therapeutic system, each with a defined role contributing to the overall function._")

                    else:
                        st.markdown("  _No component details available for this candidate._")
                    st.markdown("---") # Final divider for the tab

                with tab_metrics:
                    st.markdown("### 🔬 Evaluation Metrics (Mock)")
                    st.caption("Below are the raw simulated metric values for each component's structure.")
                    st.divider()

                    # Display AlphaFold metrics
                    alpha_fold_metrics = rec.get('alpha_fold_metrics', {})
                    for component_name, component_metrics in alpha_fold_metrics.items():
                        st.markdown(f"**{component_name}:**")
                        for metric_name, metric_value in component_metrics.items():
                            st.markdown(f"- {metric_name}: {metric_value}")
                        st.divider()

                with tab_next_steps:
                    st.markdown("### ➡️ Suggested Next Steps")
                    st.caption("Below are potential next steps for refining this system candidate.")
                    st.divider()

                    # --- In Silico Refinement Opportunities ---
                    st.subheader("🛠️ In Silico Refinement Opportunities")
                    st.markdown("""
                    <div style="background-color: #f0f0f0; padding: 10px; border-radius: 5px;">
                        <p>🧬 <strong>Optimize Guide RNA</strong></p>
                        <ul>
                            <li>📊 <strong>Common Approaches:</strong></li>
                            <ul>
                                <li>📊 <em>Analyze guide RNA sequence features</em> (e.g., GC content, PAM proximity) to optimize on-target efficiency.</li>
                                <li>💻 <em>Use machine learning models</em> (like CHOPCHOP) to predict off-target sites and minimize their impact.</li>
                                <li>🗺️ <em>Explore guide RNA sequence space</em> to find alternative sequences with similar on-target efficacy but lower off-target potential.</li>
                                <li>🔬 <em>Evaluate guide RNA stability</em> using tools like Mfold or RNAfold to ensure it remains functional in the cellular environment.</li>
                                <li>➰ <em>Optimize guide RNA secondary structure</em> to improve its binding to the Cas protein and DNA target.</li>
                                <li>📈 <em>Compare multiple guide RNA designs</em> to find the most effective one for the target mutation.</li>
                                <li>🖥️ <em>Utilize in silico tools</em> like CRISPRscan or CCTop to assess guide RNA performance.</li>
                                <li>📋 <em>Consider using multiple guide RNAs</em> in combination to increase editing efficiency and specificity.</li>
                            </ul>
                        </ul>
                    </div>
                    """, unsafe_allow_html=True)
                    st.divider()

                    st.subheader("🧪 Experimental Validation Plan")
                    st.markdown("""
                    <div style="background-color: #f0f0f0; padding: 10px; border-radius: 5px;">
                        <p>🧪 <strong>Perform In Vitro Characterization</strong></p>
                        <ul>
                            <li>📊 <strong>Common Approaches:</strong></li>
                            <ul>
                                <li>🧪 <em>Test guide RNA activity</em> in vitro using assays like Cas9 cleavage or T7 endonuclease assays.</li>
                                <li>🔬 <em>Evaluate guide RNA specificity</em> using techniques like Sanger sequencing or deep sequencing.</li>
                                <li>📈 <em>Assess guide RNA efficiency</em> by measuring the percentage of edited cells or alleles.</li>
                                <li>🌡️ <em>Test guide RNA stability</em> under various conditions (e.g., temperature, pH) to ensure its performance in vivo.</li>
                                <li>🧬 <em>Analyze off-target effects</em> using techniques like genome-wide sequencing or targeted amplicon sequencing.</li>
                                <li>📋 <em>Compare multiple guide RNA designs</em> in vitro to select the most promising ones for further testing.</li>
                            </ul>
                        </ul>
                    </div>
                    """, unsafe_allow_html=True)
                    st.divider()

                    st.subheader("🌍 Broader Therapeutic Considerations")
                    st.markdown("""
                    <div style="background-color: #f0f0f0; padding: 10px; border-radius: 5px;">
                        <p>🌎 <strong>Consider Off-Target Effects</strong></p>
                        <ul>
                            <li>📊 <strong>Common Approaches:</strong></li>
                            <ul>
                                <li>🧪 <em>Evaluate off-target effects</em> in relevant cell lines or animal models.</li>
                                <li>🌡️ <em>Test guide RNA specificity</em> under various conditions to ensure minimal off-target activity.</li>
                                <li>📋 <em>Compare multiple guide RNA designs</em> for their off-target profiles to select the safest ones.</li>
                                <li>📈 <em>Assess the potential for on-target resistance</em> and develop strategies to mitigate it.</li>
                                <li>🌐 <em>Consider the broader impact</em> of the therapy on the patient's health and well-being.</li>
                            </ul>
                        </ul>
                    </div>
                    """, unsafe_allow_html=True)
                    st.divider()

                    st.subheader("🤖 Interactive AI Agents")
                    st.markdown("""
                    <div style="background-color: #f0f0f0; padding: 10px; border-radius: 5px;">
                        <p>🤖 <strong>Utilize AI Agents for Guided Optimization</strong></p>
                        <ul>
                            <li>📊 <strong>Common Approaches:</strong></li>
                            <ul>
                                <li>🤖 <em>Leverage AI agents</em> to automate the design and optimization of guide RNAs.</li>
                                <li>💻 <em>Use machine learning models</em> to predict guide RNA performance and suggest improvements.</li>
                                <li>📈 <em>Iteratively refine guide RNA designs</em> based on AI agent feedback and experimental data.</li>
                                <li>🌐 <em>Consider the broader impact</em> of the therapy on the patient's health and well-being.</li>
                            </ul>
                        </ul>
                    </div>
                    """, unsafe_allow_html=True)
                    st.divider()

                    # --- Call to Action ---
                    st.subheader("📝 Call to Action")
                    st.markdown("""
                    <div style="background-color: #f0f0f0; padding: 10px; border-radius: 5px;">
                        <p>📝 <strong>Next Steps:</strong></p>
                        <ul>
                            <li>📊 <em>Refine the guide RNA sequence</em> based on the in silico analysis and experimental validation.</li>
                            <li>🧪 <em>Perform additional in vitro characterization</em> to ensure the guide RNA's performance and safety.</li>
                            <li>🌍 <em>Consider the broader therapeutic implications</em> and potential off-target effects.</li>
                            <li>🤖 <em>Utilize AI agents</em> to guide the optimization process and suggest improvements.</li>
                        </ul>
                    </div>
                    """, unsafe_allow_html=True)
                    st.divider()

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
                    st.divider()

                    st.markdown("### 🛠️ In Silico Refinement Opportunities")
                    st.markdown("Explore these computational steps to potentially improve your design before wet-lab experiments. Review each area and consider simulating improvements with our Interactive AI Agents where applicable.")
                    st.divider()

                    # --- Guide RNA Optimization Section ---
                    st.subheader("🧬 Optimize Guide RNA")
                    st.info("**Importance:** Maximizing on-target editing efficacy and minimizing off-target effects are critical for therapeutic safety and success.")
                    st.markdown("""
                        **🛠️ Common Approaches:**
                        - 📊 **Analyze Metrics:** Scrutinize *on-target scores* (e.g., from CHOPCHOP), off-target predictions, and structural stability (MFE).
                        - 💻 **Predict Impact:** Utilize *in silico tools* to forecast effects of sequence modifications.
                        - 🗺️ **Explore Alternatives:** Investigate *different guide sequences* targeting the same or nearby genomic loci.
                        - 🔬 **Key Factors:** Consider GC content, protospacer adjacent motif (PAM) compatibility, and potential RNA secondary structures that could impede function.
                    """)
                    guide_metrics = rec.get('mock_chopchop_metrics', {}).get('guide_rna_sequence', {})
                    if guide_metrics:
                        on_target_score = guide_metrics.get('chopchop_on_target_score', 0.0)
                        off_target_count = guide_metrics.get('chopchop_off_target_count', 0)
                        mfe = guide_metrics.get('chopchop_mfe', 0.0)
                        
                        score_col1, score_col2, score_col3 = st.columns(3)
                        score_col1.metric("Current On-Target", f"{on_target_score:.2f}")
                        score_col2.metric("Current Off-Targets", off_target_count)
                        score_col3.metric("Current MFE", f"{mfe:.2f}")

                        if on_target_score < 0.75 or off_target_count > 3:
                            st.warning(f"Simulated metrics suggest potential for improvement.")
                        else:
                            st.success("Current simulated guide metrics appear strong!")
                    else:
                        st.markdown("_No specific guide metrics found in the current recommendation to display._")
                    st.markdown("👉 **Simulate this:** Use the `Guide RNA Optimizer` agent in the 'Interactive AI Agents' section below.")
                    st.divider()

                    # --- Component Sequence Refinement Section ---
                    st.subheader("🔄 Refine Component Sequences (e.g., Cas Protein, Repair Template)")
                    st.info("**Importance:** Ensuring the stability, correct folding (for proteins), and functionality of all system components (Cas enzymes, repair templates, etc.) is vital for overall efficacy and can impact immunogenicity.")
                    st.markdown("""
                        **🛠️ Common Approaches:**
                        - 🧬 **Proteins (e.g., Cas):** Analyze structural predictions (e.g., pLDDT from AlphaFold). Identify regions of low confidence or potential misfolding. Explore *sequence variants* (mutations, truncations) predicted to improve stability, activity, or specificity. Consider *codon optimization* for the target expression system.
                        - 📄 **DNA (e.g., Repair Templates):** Optimize *homology arm lengths* and sequences for HDR. Ensure correct edit incorporation and minimize undesired integration or template-induced mutations.
                        - ➰ **RNA (e.g., viral genome components):** Ensure structural integrity and proper function of elements like ITRs or packaging signals if applicable.
                    """)
                    component_metrics_af = rec.get('mock_alphafold_system_metrics', {})
                    low_confidence_proteins = {
                        name: metrics.get('mock_plddt_score') 
                        for name, metrics in component_metrics_af.items() 
                        if isinstance(metrics, dict) and "protein" in name.lower() and metrics.get('mock_plddt_score', 1.0) < 0.7
                    }
                    if low_confidence_proteins:
                        st.warning("Found protein components with lower simulated structural confidence (pLDDT < 0.7):")
                        for name, plddt in low_confidence_proteins.items():
                            st.markdown(f"  - `{name}`: pLDDT {plddt:.2f}")
                    else:
                        st.success("Simulated structural confidence for protein components appears good (pLDDT ≥ 0.7 where applicable).")
                    st.markdown("👉 **Simulate this:** Use the `Protein Component Optimizer` agent (for Cas/proteins) or conceptually refine template designs based on literature for DNA components.")
                    st.divider()

                    # --- Advanced Off-Target Analysis Section ---
                    st.subheader("🎯 Perform Advanced Off-Target Analysis")
                    st.info("**Importance:** Comprehensive off-target analysis is a cornerstone of safety for CRISPR-based therapies, aiming to identify and mitigate risks of unintended genomic modifications.")
                    st.markdown("""
                        **🛠️ Common Approaches:**
                        - 📈 **Initial Scan:** Start with predictions from tools like CHOPCHOP.
                        - 🖥️ **Deep Dive Tools:** Employ advanced *in silico* tools (e.g., CasOFFinder, DeepHF, CFD-Score based utilities) to search the entire genome with sophisticated algorithms. These consider mismatch tolerance, DNA accessibility (sometimes), and specific Cas enzyme profiles.
                        - 📋 **Prioritize & Validate:** Generate a *prioritized list* of potential off-target sites for subsequent experimental validation.
                    """)
                    # current_off_targets is already defined if guide_metrics exists
                    if guide_metrics and 'chopchop_off_target_count' in guide_metrics:
                        current_off_targets = guide_metrics.get('chopchop_off_target_count', -1)
                        st.markdown(f"Current CHOPCHOP predicted off-targets: **{current_off_targets}**.")
                        if current_off_targets > 5:
                            st.warning("An advanced off-target scan is highly recommended due to the number of initial predictions.")
                    else:
                        st.markdown("_Guide metrics for off-target count not available._")
                    st.markdown("👉 **Simulate/Consider this:** While a direct agent for *advanced* genome-wide scans isn't simulated here, the *concept* is vital. Plan for using specialized bioinformatics tools for this analysis.")
                    st.divider()

                    # Display Experimental Validation Plan (original markdown)
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

                # --- EXECUTIVE SUMMARY SECTION ---
                st.divider()
                st.subheader("📝 Executive Summary of Simulated Findings")
                
                confidence_score_value = rec.get('confidence_score', 0.0)
                simulated_findings = rec.get('simulated_findings', [])
                
                # Safely extract strategy from editing_type, defaulting if not found
                editing_type_raw = rec.get('editing_type', 'unknown strategy')
                strategy_match = re.search(r'\\((.*?)\\)', editing_type_raw) # Regex to find content in parentheses
                strategy_text = strategy_match.group(1).strip() if strategy_match else editing_type_raw

                summary_intro = (
                    f"This mock therapeutic system, **{rec.get('recommended_approach', f'Candidate {i+1}')}**, designed for "
                    f"targeting **{current_target_mutation['hugo_gene_symbol']} {current_target_mutation['protein_change']}** via a "
                    f"'{strategy_text}' strategy, achieved a simulated confidence score of **{confidence_score_value:.2f}**. "
                    f"This score reflects an aggregation of various simulated metrics."
                )
                st.markdown(summary_intro)
                
                if simulated_findings:
                    st.markdown("---")
                    st.markdown("##### Key Simulated Findings:")
                    strengths = [f for f in simulated_findings if f.get('type') == 'strength']
                    weaknesses = [f for f in simulated_findings if f.get('type') == 'weakness']
                    
                    if strengths:
                        st.markdown("<div style='margin-left: 20px;'><strong>Positive Indicators (Strengths):</strong></div>", unsafe_allow_html=True)
                        for finding in strengths:
                            st.markdown(f"<div style='margin-left: 40px;'>✅ <strong>{finding.get('metric', 'N/A').replace('_', ' ').title()} ({finding.get('component', 'N/A').replace('_', ' ').title()})</strong>: {finding.get('description', 'N/A')}</div>", unsafe_allow_html=True)
                    else:
                        st.markdown("<div style='margin-left: 20px;'><em>No specific positive indicators highlighted in simulated findings.</em></div>", unsafe_allow_html=True)
                        
                    if weaknesses:
                        st.markdown("<div style='margin-left: 20px;'><strong>Areas for Attention (Weaknesses):</strong></div>", unsafe_allow_html=True)
                        for finding in weaknesses:
                            st.markdown(f"<div style='margin-left: 40px;'>⚠️ <strong>{finding.get('metric', 'N/A').replace('_', ' ').title()} ({finding.get('component', 'N/A').replace('_', ' ').title()})</strong>: {finding.get('description', 'N/A')}</div>", unsafe_allow_html=True)
                    else:
                        st.markdown("<div style='margin-left: 20px;'><em>No specific weaknesses highlighted in simulated findings, suggesting overall moderate to good simulated performance based on these criteria.</em></div>", unsafe_allow_html=True)
                else:
                    st.markdown("<em>No structured simulated findings were available to detail strengths or weaknesses. The 'Detailed Rationale' above provides an overall narrative.</em>")
                
                st.markdown("---")
                st.markdown("##### Overall Simulated Assessment:")
                if confidence_score_value >= 0.8:
                    st.success(f"This candidate demonstrates **strong promise** based on simulated data for targeting {current_target_mutation['hugo_gene_symbol']} {current_target_mutation['protein_change']}. Key strengths support its potential, though experimental validation is paramount.")
                elif confidence_score_value >= 0.5:
                    st.warning(f"This candidate shows **moderate promise** for {current_target_mutation['hugo_gene_symbol']} {current_target_mutation['protein_change']}. While some aspects are positive, certain simulated metrics indicate areas requiring further investigation or optimization before proceeding.")
                else:
                    st.error(f"This candidate has a **lower simulated confidence** for {current_target_mutation['hugo_gene_symbol']} {current_target_mutation['protein_change']}. Significant refinement and review of foundational metrics (e.g., target validity, component scores) are likely necessary based on these mock results.")
                st.caption("_Reminder: All findings and scores are part of a simulation and require rigorous experimental validation._")
                # --- END EXECUTIVE SUMMARY SECTION ---


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

    # --- 7. Display Workflow Log ---
    if st.session_state.workflow_log:
        with st.expander("🔍 View Full Workflow Log", expanded=False):
            for log_entry in st.session_state.workflow_log:
                if isinstance(log_entry, tuple) and len(log_entry) == 2: # Check if it's a tuple of (timestamp, message)
                    timestamp, message = log_entry
                    st.text(f"{timestamp}: {message}")
                elif isinstance(log_entry, str): # Handle plain string log entries
                    st.text(log_entry)
                else:
                    st.text(str(log_entry)) # Fallback for other types
    else:
        st.info("Workflow log is empty. Run a design process to populate it.")

    # --- 8. Simulated System Recommendations ---
    st.header("🏆 Simulated System Recommendations")
    if not st.session_state.therapy_recommendations:
        st.info("No mock therapy system recommendations generated yet. Click the 'Design Mock Therapy System' button above.")
    else:
        # Sort recommendations by confidence score, highest first
        sorted_recommendations = sorted(st.session_state.therapy_recommendations, key=lambda x: x.get('confidence_score', 0.0), reverse=True)
        
        # Use st.columns to create a layout for recommendations
        # Determine number of columns (e.g., 2 or 3 based on preference)
        # num_cols = 2 # Or 3
        # cols = st.columns(num_cols)
        
        active_recommendation = None
        if 'active_recommendation_id' in st.session_state and st.session_state.active_recommendation_id is not None:
            for rec_item in sorted_recommendations: # Use a different variable name here
                if rec_item.get("recommended_approach") == st.session_state.active_recommendation_id: 
                    active_recommendation = rec_item # Assign to active_recommendation
                    break
        
        if active_recommendation:
            main_col, details_col = st.columns([0.6, 0.4])
        else:
            main_col = st.container() # Use a container for the main column if no details view

        with main_col:
            st.subheader("Select a System to Explore Further:")
            for i, rec in enumerate(sorted_recommendations):
                rec_id = rec.get("recommended_approach", f"System_{i}") # Use a unique ID
                
                with st.container(): # Removed border=True for cleaner look, can be re-added
                    st.markdown(f"##### **{i+1}. {rec.get('recommended_approach', 'N/A')}**")
                    st.metric(label="Mock Confidence Score", value=f"{rec.get('confidence_score', 0.0):.2f}")
                    
                    # Truncate or summarize components for brevity in this main list
                    components = rec.get('designed_system_components', {})
                    component_names = list(components.keys())
                    if component_names:
                        st.markdown(f"**Key Components (Simulated):** {', '.join(component_names[:3])}{'...' if len(component_names) > 3 else ''}")
                    
                    # Button to select this recommendation for detailed view
                    if st.button(f"Explore System: {rec.get('recommended_approach', 'N/A')}", key=f"explore_{rec_id}"):
                        st.session_state.active_recommendation_id = rec.get("recommended_approach")
                        st.session_state.agent_results = {} # Clear previous agent results
                        st.rerun() # Rerun to update the details column

        if active_recommendation: # Only populate details_col if it exists
            with details_col:
                st.subheader(f"🔬 Details for: {active_recommendation.get('recommended_approach', 'N/A')}")
                
                # Use a more structured display for details
                # Confidence Score (already shown in metric, but can repeat or add context)
                st.markdown(f"**Mock Confidence Score:** {active_recommendation.get('confidence_score', 0.0):.2f}")

                # Editing Type
                st.markdown(f"**Editing Type (Simulated):** {active_recommendation.get('editing_type', 'N/A')}")

                # Designed System Components with Sequences (using expander for long sequences)
                with st.expander("🧬 Designed System Components (Simulated Sequences)", expanded=False):
                    components = active_recommendation.get('designed_system_components', {})
                    if components:
                        for comp_name, comp_seq in components.items():
                            st.code(f"{comp_name}:\n{comp_seq}", language="text")
                    else:
                        st.caption("No detailed components provided.")
                
                # Detailed Rationale
                with st.expander("Detailed Rationale (Simulated)", expanded=True):
                    st.markdown(active_recommendation.get('detailed_rationale', 'Not available.'))

                # Simulated Findings (Strengths/Weaknesses)
                with st.expander("🔍 Simulated Findings (Strengths & Weaknesses)", expanded=False):
                    findings = active_recommendation.get('simulated_findings', [])
                    if findings:
                        for finding in findings:
                            emoji = "👍" if finding.get('type') == 'strength' else "👎"
                            st.markdown(f"{emoji} **{finding.get('type', '').capitalize()} in {finding.get('metric', '')} for '{finding.get('component', '')}' (Value: {finding.get('value', 'N/A')}):** {finding.get('description', '')}")
                    else:
                        st.caption("No specific findings provided.")

                # Supporting Evidence (Raw Metrics)
                with st.expander("📊 Supporting Evidence (Simulated Raw Metrics)", expanded=False):
                    evidence = active_recommendation.get('supporting_evidence', [])
                    if evidence:
                        for item in evidence:
                            st.markdown(f"- {item}")
                    else:
                        st.caption("No supporting evidence provided.")

                # Mock AlphaFold System Metrics (Per Component)
                with st.expander("🧱 Mock AlphaFold Component Metrics (Simulated)", expanded=False):
                    af_metrics = active_recommendation.get('mock_alphafold_system_metrics', {})
                    if af_metrics:
                        for comp, metrics_dict in af_metrics.items():
                            st.markdown(f"**Component: {comp}**")
                            for key, val in metrics_dict.items():
                                st.markdown(f"  - `{key}`: {val}")
                    else:
                        st.caption("No AlphaFold component metrics provided.")
                
                # Mock CHOPCHOP Metrics (for guide RNAs)
                with st.expander("✂️ Mock CHOPCHOP Guide Metrics (Simulated)", expanded=False):
                    cc_metrics_all = active_recommendation.get('mock_chopchop_metrics', {})
                    if cc_metrics_all and 'guide_rna_sequence' in cc_metrics_all : # Assuming structure { 'guide_rna_sequence': {metrics} }
                        guide_metrics = cc_metrics_all['guide_rna_sequence']
                        st.markdown(f"**For Guide RNA:**") # Assuming one primary guide per system for this mock
                        for key, val in guide_metrics.items():
                             st.markdown(f"  - `{key}`: {val}")
                    elif cc_metrics_all: # If structure is just {metrics} directly
                        for key, val in cc_metrics_all.items():
                             st.markdown(f"  - `{key}`: {val}")
                    else:
                        st.caption("No CHOPCHOP metrics provided.")

                # --- Suggested Next Steps ---
                st.markdown("<h4 style='color: black;'>➡️ Suggested Next Steps</h4>", unsafe_allow_html=True)
                st.caption("Interactive agents for these steps have been moved to the 'Interactive Agents' page in the navigation panel.")

                # The entire "Interactive Agents / Next Steps" section that was here has been removed.
                # This section previously handled displaying potential agent suggestions and their execution.
                # It started with a check for active_recommendation and 'potential_agent_suggestions'
                # and included loops for buttons, agent calls, and results display.

    # --- 9. Educational Content / LLM Chat - Already in Sidebar ---
    # The main content of the display_mock_therapy_design_page function might end here or before fallbacks.

# --- Main execution ---
if __name__ == "__main__":
    # This is where the Streamlit page content is actually rendered
    display_mock_therapy_design_page()

