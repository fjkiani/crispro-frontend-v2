"""
Platinum Response Data Hunt - Main Streamlit App

Beautiful, modular dashboard showcasing extracted platinum response data.
"""

import sys
from pathlib import Path

# Add dashboard directory to path for imports
dashboard_dir = Path(__file__).parent
if str(dashboard_dir) not in sys.path:
    sys.path.insert(0, str(dashboard_dir))

import streamlit as st
from config import PAGE_CONFIG
from pages import overview, patients, overlap, validation

# Page configuration
st.set_page_config(**PAGE_CONFIG)

# Sidebar navigation
st.sidebar.title("⚔️ Platinum Response Dashboard")
st.sidebar.markdown("---")

# Page selection
page = st.sidebar.radio(
    "Navigate",
    ["Overview", "Patients", "Overlap Analysis", "Validation Readiness"],
    index=0,
)

st.sidebar.markdown("---")
st.sidebar.markdown("### 📊 Quick Stats")
st.sidebar.info("""
**Data Sources:**
- Jr2's Dataset: 469 patients
- Zo's Dataset: 200 patients
- Overlap: 161 patients

**Status:** ✅ Validation Ready
""")

# Render selected page
if page == "Overview":
    overview.render_overview_page()
elif page == "Patients":
    patients.render_patients_page()
elif page == "Overlap Analysis":
    overlap.render_overlap_page()
elif page == "Validation Readiness":
    validation.render_validation_page()

# Footer
st.markdown("---")
st.markdown(
    """
    <div style='text-align: center; color: #86868b; padding: 2rem;'>
        <p>⚔️ Platinum Response Data Hunt Dashboard</p>
        <p style='font-size: 0.9rem;'>Powered by CrisPRO.ai | Research Use Only</p>
    </div>
    """,
    unsafe_allow_html=True,
)

