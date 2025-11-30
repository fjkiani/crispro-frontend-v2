"""
Hero Metrics Component - Displays key statistics in beautiful cards
"""

import sys
from pathlib import Path

# Add dashboard directory to path
dashboard_dir = Path(__file__).parent.parent
if str(dashboard_dir) not in sys.path:
    sys.path.insert(0, str(dashboard_dir))

import streamlit as st
from config import COLORS


def render_hero_metrics(stats: dict):
    """
    Render hero metrics cards
    
    Args:
        stats: Dictionary with statistics (jr2_total, overlap_count, match_rate, validation_ready)
    """
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric(
            label="📊 Total Patients Extracted",
            value=f"{stats.get('jr2_total', 0):,}",
            delta=f"From {stats.get('source_breakdown', {}).get('gdc_xml', 0)} GDC files",
            delta_color="normal",
        )
    
    with col2:
        st.metric(
            label="🔗 Overlap with Zo's Dataset",
            value=f"{stats.get('overlap_count', 0):,}",
            delta=f"{stats.get('match_rate', 0):.1f}% match rate",
            delta_color="normal" if stats.get('overlap_count', 0) >= 40 else "off",
        )
    
    with col3:
        match_rate = stats.get('match_rate', 0)
        st.metric(
            label="📈 Match Rate",
            value=f"{match_rate:.1f}%",
            delta="34.3% overall" if match_rate > 30 else "Below target",
            delta_color="normal" if match_rate > 30 else "off",
        )
    
    with col4:
        validation_status = "✅ Ready" if stats.get('validation_ready', False) else "⏸️ Pending"
        validation_color = "normal" if stats.get('validation_ready', False) else "off"
        overlap_count = stats.get('overlap_count', 0)
        threshold = 40
        delta_value = f"{overlap_count - threshold} over threshold" if overlap_count >= threshold else f"{threshold - overlap_count} needed"
        
        st.metric(
            label="⚔️ Validation Status",
            value=validation_status,
            delta=delta_value,
            delta_color=validation_color,
        )

