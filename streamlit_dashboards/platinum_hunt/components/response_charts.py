"""
Response Distribution Charts Component
"""

import sys
from pathlib import Path

# Add dashboard directory to path
dashboard_dir = Path(__file__).parent.parent
if str(dashboard_dir) not in sys.path:
    sys.path.insert(0, str(dashboard_dir))

import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
from collections import Counter
from config import RESPONSE_COLORS, RESPONSE_LABELS


def render_response_distribution(stats: dict, jr2_data: dict = None):
    """
    Render response distribution charts
    
    Args:
        stats: Dictionary with response_distribution
        jr2_data: Optional full Jr2 dataset for detailed breakdown
    """
    response_dist = stats.get("response_distribution", {})
    
    if not response_dist:
        st.warning("⚠️ No response distribution data available")
        return
    
    # Create two columns for side-by-side charts
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("📊 Response Distribution (Pie Chart)")
        
        # Prepare data for pie chart
        labels = [RESPONSE_LABELS.get(k, k.title()) for k in response_dist.keys()]
        values = list(response_dist.values())
        colors_list = [RESPONSE_COLORS.get(k, "#808080") for k in response_dist.keys()]
        
        fig_pie = px.pie(
            values=values,
            names=labels,
            color=labels,
            color_discrete_map={RESPONSE_LABELS.get(k, k): RESPONSE_COLORS.get(k, "#808080") 
                               for k in response_dist.keys()},
            hole=0.4,  # Donut chart
        )
        fig_pie.update_traces(
            textposition='inside',
            textinfo='percent+label',
            hovertemplate='<b>%{label}</b><br>Count: %{value}<br>Percentage: %{percent}<extra></extra>',
        )
        fig_pie.update_layout(
            showlegend=True,
            height=400,
            margin=dict(l=20, r=20, t=20, b=20),
        )
        st.plotly_chart(fig_pie, use_container_width=True)
    
    with col2:
        st.subheader("📊 Response Distribution (Bar Chart)")
        
        # Prepare data for bar chart
        labels = [RESPONSE_LABELS.get(k, k.title()) for k in response_dist.keys()]
        values = list(response_dist.values())
        colors_list = [RESPONSE_COLORS.get(k, "#808080") for k in response_dist.keys()]
        
        fig_bar = px.bar(
            x=labels,
            y=values,
            color=labels,
            color_discrete_map={RESPONSE_LABELS.get(k, k): RESPONSE_COLORS.get(k, "#808080") 
                               for k in response_dist.keys()},
            labels={'x': 'Response Type', 'y': 'Patient Count'},
        )
        fig_bar.update_traces(
            hovertemplate='<b>%{x}</b><br>Count: %{y}<extra></extra>',
            text=values,
            textposition='outside',
        )
        fig_bar.update_layout(
            showlegend=False,
            height=400,
            margin=dict(l=20, r=20, t=20, b=20),
            xaxis_title="Response Type",
            yaxis_title="Patient Count",
        )
        st.plotly_chart(fig_bar, use_container_width=True)
    
    # Response statistics table
    st.subheader("📋 Response Statistics")
    col1, col2, col3 = st.columns(3)
    
    total = sum(response_dist.values())
    
    for i, (response_type, count) in enumerate(response_dist.items()):
        percentage = (count / total * 100) if total > 0 else 0
        with [col1, col2, col3][i % 3]:
            st.metric(
                label=RESPONSE_LABELS.get(response_type, response_type.title()),
                value=f"{count:,}",
                delta=f"{percentage:.1f}%",
            )

