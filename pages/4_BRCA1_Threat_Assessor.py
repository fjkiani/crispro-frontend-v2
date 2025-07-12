import streamlit as st
import requests
import os
import sys

# Ensure the app can find the tools directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tools.ui_enhancements import (
    apply_leap_styling, 
    create_enhanced_sidebar, 
    AppleStyleResultsDisplay
)

# --- Page Configuration ---
st.set_page_config(
    page_title="BRCA1 Threat Assessor",
    page_icon="🛡️",
    layout="wide"
)

# Apply consistent styling
apply_leap_styling()

# --- Sidebar ---
# Using a simplified version of the sidebar for this focused tool
create_enhanced_sidebar("Guardian Digital Twin", {"gene_symbol": "BRCA1", "focus_area": "risk_assessment"})

# --- Main Content ---
st.markdown("# 🛡️ BRCA1 Threat Assessor")
st.markdown("---")
st.markdown("""
This tool provides a real-time risk assessment for any `BRCA1` variant using our validated **`CommandCenter`** API. 
Enter a protein change below to see the **Triumvirate Protocol** in action.
""")

# --- Analysis Interface ---
st.header("🔬 Analyze a BRCA1 Variant")

# Use a form to group the input and button
with st.form(key='brca1_assessment_form'):
    protein_change = st.text_input(
        label="Enter BRCA1 Protein Change",
        value="p.C61G",
        help="e.g., p.C61G (a missense mutation) or p.I1108* (a truncating mutation)"
    )
    submit_button = st.form_submit_button(label="⚡ Assess Threat")

# --- Backend Logic & Results Display ---
if submit_button and protein_change:
    st.markdown("---")
    # Display a spinner while the API call is in progress
    with st.spinner("Calling `CommandCenter`... Simulating mutation and running Truncation Sieve..."):
        try:
            # The validated endpoint for our service
            api_url = "https://crispro--command-center-v2-commandcenter-api.modal.run/workflow/assess_threat"
            payload = {
                "gene_symbol": "BRCA1",
                "protein_change": protein_change
            }
            
            response = requests.post(api_url, json=payload, timeout=60)
            
            if response.status_code == 200:
                result = response.json()
                
                # Use the same high-quality display component from the other page
                AppleStyleResultsDisplay.display_threat_assessment(result)
                
                # --- INJECTING THE COUNTER-OFFENSIVE PROTOCOL UI ---
                # Check if the assessment is pathogenic from the nested dictionary
                if result.get("assessment", {}).get("is_pathogenic"):
                    st.markdown("---")
                    st.markdown("## 🔥 Metastatic Counter-Offensive Protocol 🔥")
                    st.warning("🚨 **ACTION REQUIRED:** Pathogenic variant detected. The following protocols from the 8-Step Metastasis Kill Chain are now available. Select a counter-measure to dispatch a design order to the Zeta Forge.")

                    # --- Kill Chain Step 4: Invasion ---
                    with st.expander("💣 STEP 4: DESIGN INVASION BLOCKER"):
                        st.markdown("Prevent the enemy from breaking containment by destroying their tissue-dissolving enzymes (MMPs).")
                        if st.button("Task Zeta Forge: Design MMP Neutralizer", key="mmp_neutralizer"):
                            st.success("✅ COMMAND ACKNOWLEDGED. Zeta Forge is designing a gRNA weapon to neutralize the MMP gene family.")

                    # --- Kill Chain Step 6: Homing ---
                    with st.expander("📡 STEP 6: JAM HOMING BEACONS"):
                        st.markdown("Blind the cancer's guidance systems (e.g., CXCR4 receptors) so it can't find a new organ to invade.")
                        if st.button("Task Zeta Forge: Design Receptor Jammer", key="receptor_jammer"):
                            st.success("✅ COMMAND ACKNOWLEDGED. Zeta Forge is designing a nanobody to physically block the CXCR4 receptor.")
                    
                    # --- Kill Chain Step 5: Immune Evasion ---
                    with st.expander("🛡️ STEP 5: DISABLE ENEMY CLOAKING"):
                        st.markdown("Disable the cancer's 'cloaking device' (e.g. PD-L1) to make it visible to our forces (the immune system).")
                        if st.button("Task Zeta Forge: Design PD-L1 Neutralizer", key="pd_l1_neutralizer"):
                            st.success("✅ COMMAND ACKNOWLEDGED. Zeta Forge is designing a gRNA weapon to disable the PD-L1 gene.")
                
            else:
                st.error(f"❌ Analysis failed with status code {response.status_code}:")
                st.json(response.json())
                
        except requests.exceptions.RequestException as e:
            st.error(f"❌ A connection error occurred: {str(e)}")
        except Exception as e:
            st.error(f"❌ An unexpected error occurred: {str(e)}")

# --- Explainer ---
with st.expander("ℹ️ How this works", expanded=False):
    st.markdown("""
    This tool demonstrates the first layer of our **Triumvirate Protocol**:
    
    1.  **Input:** You provide a `BRCA1` protein change.
    2.  **API Call:** The page calls our deployed `CommandCenter` service.
    3.  **Truncation Sieve:** The service fetches the `BRCA1` DNA sequence and simulates the mutation's effect. It translates the mutated DNA to a protein and checks for a premature STOP codon.
    4.  **Verdict:** If a premature stop is found, it's flagged as **PATHOGENIC**. If not, it's passed to the next stage (the Zeta-Oracle AI, which is part of Phase 3 of our plan).
    """) 