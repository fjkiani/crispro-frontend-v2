import streamlit as st
import sys
import os
import time
import random

# Page Configuration
st.set_page_config(page_title="CasPro 🧬 - Interactive Agents", layout="wide")

# Add the project root to the Python path
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if PROJECT_ROOT not in sys.path:
    sys.path.append(PROJECT_ROOT)

# Attempt to import the run_interactive_agent function and other necessary components
try:
    from ai_research_assistant import run_interactive_agent
    # If other specific simulation functions were directly called by agents, import them too,
    # though the design pattern is for run_interactive_agent to handle this.
except ImportError as e:
    st.error(f"Failed to import `run_interactive_agent` from `ai_research_assistant`. Ensure the module is accessible: {e}")
    def run_interactive_agent(*args, **kwargs): # Fallback
        st.error("`run_interactive_agent` is not available due to import error.")
        # Ensure the display_callback is handled if it's a kwarg
        display_callback = kwargs.get('display_callback', st.text)
        display_callback("Error: Agent simulation unavailable.")
        return {"error": "Agent not available", "explanation": "Import failed."}

# --- Initialize Session State Keys specific to this page or if not already globally managed ---
if 'agent_results' not in st.session_state: # This should ideally be managed where recommendations are selected
    st.session_state.agent_results = {}
if 'active_recommendation_id' not in st.session_state: # This ID links to the selected recommendation
    st.session_state.active_recommendation_id = None
if 'recommendations' not in st.session_state: # Full list of recommendations
    st.session_state.recommendations = []

# --- Helper for display callback ---
def agent_display_callback(message, container=None):
    if container:
        container.info(message) # Display within a specific container if provided
    else:
        st.info(message) # Default display


# --- Main Page Display ---
st.title("🤖 Interactive Therapeutic Agent Simulations")
st.markdown("""
This page allows you to run simulated follow-up actions ("agents") based on the detailed system recommendation
you explored on the "Mock Therapeutic System Design" page. Select an agent to simulate its process and see mock results.
""")

active_recommendation = None
if st.session_state.active_recommendation_id and st.session_state.recommendations:
    for rec in st.session_state.recommendations:
        # Assuming 'recommended_approach' was used as a unique ID
        if rec.get("recommended_approach") == st.session_state.active_recommendation_id:
            active_recommendation = rec
            break

if not active_recommendation:
    st.warning("Please first select a 'Mock System Recommendation' to explore on the 'Mock Therapeutic System Design' page. "
               "The interactive agents are specific to a selected system.")
else:
    st.subheader(f"🔬 Agents for System: {active_recommendation.get('recommended_approach', 'N/A')}")
    st.markdown(f"**Original Mock Confidence Score:** {active_recommendation.get('confidence_score', 0.0):.2f}")
    
    potential_agents = active_recommendation.get('potential_agent_suggestions', [])

    if not potential_agents:
        st.info("No interactive agent suggestions available for this selected system recommendation.")
    else:
        st.markdown("Select an agent below to simulate the next step:")

        # Display area for agent's operational messages
        agent_op_container = st.empty() 

        num_cols = 2 # Adjust as needed
        cols = st.columns(num_cols)
        col_idx = 0

        for agent_suggestion in potential_agents:
            agent_type = agent_suggestion.get('type') # e.g., "refinement_suggestion", "experimental_suggestion"
            agent_text = agent_suggestion.get('text', 'Unnamed Agent')
            agent_details = agent_suggestion.get('details', 'No details provided.')
            # Construct a unique key for the button based on agent text and type to avoid collisions
            button_key = f"agent_button_{agent_type}_{agent_text.replace(' ','_').lower()}"
            
            # Use the agent_type from the suggestion, map to what run_interactive_agent expects
            # This mapping needs to be robust. Example: 'refinement_suggestion' might map to 'guide_optimization' or 'component_optimization'
            # For now, we'll assume the 'text' or a new field in 'potential_agent_suggestions' more directly maps to agent_type for run_interactive_agent
            
            # A simple heuristic for mapping:
            run_agent_type_mapping = {
                "Optimize guide RNA on-target score": "guide_optimization",
                "General in silico refinement": "component_optimization", # Or a more generic one
                "Investigate Cas9 protein sequence stability": "component_optimization", # Assuming Cas9 is the component_name
                "Address high off-target count": "off_target_prediction", # Or could be guide_optimization
                "Prioritize off-target validation": "experimental_off_target_validation",
                "Plan in vitro characterization": "in_vitro_validation",
                "Design cell-based assays": "cell_based_assay_simulation",
                "Outline off-target analysis strategy": "experimental_off_target_validation", # Or in_silico prediction
                "Assess HDR efficiency (if applicable)": "hdr_template_optimization", # Needs more context
                "Discuss delivery strategy": "delivery_optimization_simulation",
                "Plan pre-clinical studies": "full_development_planning", # This is broad
                "Consider immunogenicity": "component_optimization", # (e.g., deimmunize protein)
                "Evaluate long-term effects": "full_development_planning" # Broad
            }
            # Default to a generic component_optimization if text doesn't match a specific agent, or use a new field if available
            # For example, if agent_suggestion had a 'agent_tool_id' field:
            # current_agent_tool_id = agent_suggestion.get('agent_tool_id', 'component_optimization')

            # Fallback: try to infer from text, or use a generic name if agent_suggestion has agent_key
            current_agent_tool_id = agent_suggestion.get('agent_key') # Prefer this if present
            if not current_agent_tool_id:
                 current_agent_tool_id = run_agent_type_mapping.get(agent_text, "component_optimization")


            with cols[col_idx % num_cols]:
                if st.button(f"🚀 Simulate: {agent_text}", key=button_key, help=agent_details):
                    st.session_state.agent_results[button_key] = None # Clear previous result for this button
                    agent_op_container.info(f"Simulating: {agent_text}...")
                    
                    # Prepare context for run_interactive_agent
                    # This might need specific component names or other data from `active_recommendation`
                    # For "component_optimization", we need a component_name
                    component_name_for_agent = None
                    if "guide" in agent_text.lower() or "rna" in agent_text.lower():
                         # Try to get the actual guide RNA sequence name if available, otherwise default
                        component_name_for_agent = next((k for k in active_recommendation.get('designed_system_components', {}) if "guide" in k.lower()), None)
                        if not component_name_for_agent and 'guide_rna_sequence' in active_recommendation.get('designed_system_components', {}):
                            component_name_for_agent = 'guide_rna_sequence'

                    elif "protein" in agent_text.lower() or "cas" in agent_text.lower():
                        component_name_for_agent = next((k for k in active_recommendation.get('designed_system_components', {}) if "protein" in k.lower() or "cas" in k.lower()), None)
                        if not component_name_for_agent and 'cas_protein_sequence' in active_recommendation.get('designed_system_components', {}):
                            component_name_for_agent = 'cas_protein_sequence'
                    # Add more specific component name extractions if needed for other agent types
                    
                    with st.spinner(f"Running agent: {agent_text}..."):
                        agent_result = run_interactive_agent(
                            agent_type=current_agent_tool_id, # This is the crucial mapping
                            recommendation=active_recommendation, # Pass the full recommendation
                            editing_type=active_recommendation.get('editing_type'),
                            component_name=component_name_for_agent, # Pass identified component name
                            display_callback=lambda msg: agent_display_callback(msg, agent_op_container)
                        )
                    st.session_state.agent_results[button_key] = agent_result
                    agent_op_container.empty() # Clear the "Simulating..." message
                    st.rerun() # Rerun to display results immediately
            col_idx += 1

        # Display Agent Results
        st.markdown("---")
        st.subheader("Agent Simulation Results:")
        
        if not st.session_state.agent_results:
            st.caption("No agent simulations run yet for this system recommendation.")
        
        for button_key, result in st.session_state.agent_results.items():
            if result: # If a result exists for this button
                agent_text_from_key = button_key.replace("agent_button_", "").replace("_", " ").title() # Reconstruct approx agent name
                # Try to find the original agent_text for a better title
                original_agent_text = agent_text_from_key # Default
                for ag_sugg in potential_agents:
                    reconstructed_key = f"agent_button_{ag_sugg.get('type')}_{ag_sugg.get('text', '').replace(' ','_').lower()}"
                    if reconstructed_key == button_key:
                        original_agent_text = ag_sugg.get('text', agent_text_from_key)
                        break
                
                with st.expander(f"📜 Result for: Simulate: {original_agent_text}", expanded=True):
                    if isinstance(result, dict) and result.get("error"):
                        st.error(f"Error: {result.get('error')}")
                        if result.get("explanation"):
                            st.caption(result.get("explanation"))
                    elif isinstance(result, dict):
                        st.success("Simulation Complete!")
                        # Generic display of dictionary results
                        for res_key, res_value in result.items():
                            if res_key.lower() == "explanation":
                                st.markdown(f"**Explanation:** {res_value}")
                            elif isinstance(res_value, dict):
                                st.markdown(f"**{res_key.replace('_',' ').title()}:**")
                                for sub_key, sub_val in res_value.items():
                                    st.markdown(f"  - `{sub_key}`: {sub_val}")
                            elif isinstance(res_value, list):
                                st.markdown(f"**{res_key.replace('_',' ').title()}:**")
                                for item in res_value:
                                    if isinstance(item, dict):
                                        for sub_k, sub_v in item.items():
                                             st.markdown(f"    - `{sub_k}`: {sub_v}")
                                    else:
                                        st.markdown(f"  - {item}")
                            else:
                                st.markdown(f"**{res_key.replace('_',' ').title()}:** {res_value}")
                    else:
                        st.write(result) # Fallback for non-dict results

# Add some spacing at the bottom
st.markdown("<br><br>", unsafe_allow_html=True)

# To run this page directly if needed for testing (streamlit run pages/4_Interactive_Agents.py)
if __name__ == "__main__":
    pass
    # Minimal setup for direct run if needed for isolated testing
    # For example, you might call a main function for this page if it were structured that way:
    # main_interactive_agents_page_logic() 