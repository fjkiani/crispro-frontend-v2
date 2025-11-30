"""
Validation Page - Validation readiness and statistical power
"""

import sys
from pathlib import Path

# Add dashboard directory to path
dashboard_dir = Path(__file__).parent.parent
if str(dashboard_dir) not in sys.path:
    sys.path.insert(0, str(dashboard_dir))

import streamlit as st
from data.loader import load_datasets, merge_datasets
from data.processor import compute_statistics, compute_overlap
from components.validation_status import render_validation_status


def render_validation_page():
    """Render the validation readiness page"""
    st.title("⚔️ Validation Readiness")
    st.markdown("### Statistical power analysis and validation checklist")
    
    # Load data
    with st.spinner("Loading validation data..."):
        jr2_data, zo_data = load_datasets()
        
        if not jr2_data or not zo_data:
            st.error("❌ Failed to load datasets. Please check data files.")
            return
        
        merged_data = merge_datasets(jr2_data, zo_data)
        stats = compute_statistics(jr2_data, zo_data, merged_data)
        overlap_data = compute_overlap(jr2_data, zo_data)
    
    # Validation status
    render_validation_status(stats, overlap_data)
    
    # Response distribution for validation
    st.markdown("---")
    st.subheader("📊 Response Distribution for Validation Cohort")
    
    if merged_data:
        response_counts = {}
        for patient in merged_data:
            response = patient.get("platinum_response", "Unknown")
            response_counts[response] = response_counts.get(response, 0) + 1
        
        col1, col2, col3 = st.columns(3)
        for i, (response, count) in enumerate(response_counts.items()):
            with [col1, col2, col3][i % 3]:
                st.metric(
                    label=response.title(),
                    value=f"{count:,}",
                    delta=f"{count / len(merged_data) * 100:.1f}%",
                )

