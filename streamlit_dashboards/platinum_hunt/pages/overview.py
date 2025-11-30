"""
Overview Page - Main dashboard with hero metrics and summary
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
from components.hero_metrics import render_hero_metrics
from components.response_charts import render_response_distribution
from components.validation_status import render_validation_status


def render_overview_page():
    """Render the overview/dashboard page"""
    st.title("⚔️ Platinum Response Data Hunt Dashboard")
    st.markdown("### Comprehensive visualization of extracted platinum response data")
    
    # Load data
    with st.spinner("Loading datasets..."):
        jr2_data, zo_data = load_datasets()
        
        if not jr2_data or not zo_data:
            st.error("❌ Failed to load datasets. Please check data files.")
            return
        
        merged_data = merge_datasets(jr2_data, zo_data)
        stats = compute_statistics(jr2_data, zo_data, merged_data)
        overlap_data = compute_overlap(jr2_data, zo_data)
    
    # Hero metrics
    st.markdown("---")
    render_hero_metrics(stats)
    
    # Response distribution
    st.markdown("---")
    render_response_distribution(stats, jr2_data)
    
    # Validation status
    st.markdown("---")
    render_validation_status(stats, overlap_data)
    
    # Quick stats
    st.markdown("---")
    st.subheader("📊 Quick Statistics")
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Jr2's Total", f"{stats.get('jr2_total', 0):,}")
    
    with col2:
        st.metric("Zo's Total", f"{stats.get('zo_total', 0):,}")
    
    with col3:
        st.metric("Overlap", f"{stats.get('overlap_count', 0):,}")
    
    with col4:
        st.metric("Match Rate", f"{stats.get('match_rate', 0):.1f}%")
    
    # Source breakdown
    if stats.get("source_breakdown"):
        st.subheader("📦 Data Source Breakdown")
        source_breakdown = stats.get("source_breakdown", {})
        
        col1, col2 = st.columns(2)
        
        with col1:
            for source, count in source_breakdown.items():
                st.write(f"**{source}:** {count:,} patients")
        
        with col2:
            # Source pie chart
            import plotly.express as px
            fig = px.pie(
                values=list(source_breakdown.values()),
                names=list(source_breakdown.keys()),
                title="Patients by Source",
            )
            st.plotly_chart(fig, use_container_width=True)

