"""
Overlap Page - Detailed overlap analysis
"""

import sys
from pathlib import Path

# Add dashboard directory to path
dashboard_dir = Path(__file__).parent.parent
if str(dashboard_dir) not in sys.path:
    sys.path.insert(0, str(dashboard_dir))

import streamlit as st
from data.loader import load_datasets
from data.processor import compute_overlap
from components.overlap_analysis import render_overlap_analysis


def render_overlap_page():
    """Render the overlap analysis page"""
    st.title("🔗 Overlap Analysis")
    st.markdown("### Detailed analysis of dataset overlap and matching")
    
    # Load data
    with st.spinner("Computing overlap..."):
        jr2_data, zo_data = load_datasets()
        
        if not jr2_data or not zo_data:
            st.error("❌ Failed to load datasets. Please check data files.")
            return
        
        overlap_data = compute_overlap(jr2_data, zo_data)
    
    # Overlap analysis
    render_overlap_analysis(overlap_data)
    
    # Additional metrics
    st.markdown("---")
    st.subheader("📊 Additional Overlap Metrics")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric(
            "Jr2 Only",
            f"{overlap_data.get('jr2_only', 0):,}",
            "Patients with response but no mutations",
        )
    
    with col2:
        st.metric(
            "Zo Only",
            f"{overlap_data.get('zo_only', 0):,}",
            "Patients with mutations but no response",
        )
    
    with col3:
        st.metric(
            "Both",
            f"{overlap_data.get('overlap_count', 0):,}",
            "Patients with both response and mutations",
        )

