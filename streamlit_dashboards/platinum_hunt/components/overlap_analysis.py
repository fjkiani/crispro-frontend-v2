"""
Overlap Analysis Component - Venn diagrams and match rate analysis
"""

import sys
from pathlib import Path

# Add dashboard directory to path
dashboard_dir = Path(__file__).parent.parent
if str(dashboard_dir) not in sys.path:
    sys.path.insert(0, str(dashboard_dir))

import streamlit as st
import plotly.graph_objects as go
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import io
import base64
from config import COLORS


def render_overlap_analysis(overlap_data: dict):
    """
    Render overlap analysis with Venn diagram and metrics
    
    Args:
        overlap_data: Dictionary with overlap metrics from compute_overlap()
    """
    if not overlap_data:
        st.warning("⚠️ No overlap data available")
        return
    
    # Key metrics
    jr2_total = overlap_data.get("jr2_total", 0)
    zo_total = overlap_data.get("zo_total", 0)
    overlap_count = overlap_data.get("overlap_count", 0)
    match_rate_jr2 = overlap_data.get("match_rate_jr2", 0)
    match_rate_zo = overlap_data.get("match_rate_zo", 0)
    
    # Display metrics
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric(
            label="Jr2's Dataset",
            value=f"{jr2_total:,} patients",
            delta=f"{match_rate_jr2:.1f}% matched",
        )
    
    with col2:
        st.metric(
            label="Overlap",
            value=f"{overlap_count:,} patients",
            delta=f"{overlap_count - 40} over threshold" if overlap_count >= 40 else f"{40 - overlap_count} needed",
            delta_color="normal" if overlap_count >= 40 else "off",
        )
    
    with col3:
        st.metric(
            label="Zo's Dataset",
            value=f"{zo_total:,} patients",
            delta=f"{match_rate_zo:.1f}% matched",
        )
    
    # Venn diagram
    st.subheader("🔗 Dataset Overlap Visualization")
    
    try:
        # Create Venn diagram using matplotlib
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Calculate sets
        jr2_only = overlap_data.get("jr2_only", 0)
        zo_only = overlap_data.get("zo_only", 0)
        
        # Create Venn diagram
        venn = venn2(
            subsets=(jr2_only, zo_only, overlap_count),
            set_labels=("Jr2's Dataset\n(Platinum Response)", "Zo's Dataset\n(Mutations)"),
            ax=ax,
            set_colors=(COLORS["primary"], COLORS["secondary"]),
            alpha=0.7,
        )
        
        # Style the diagram
        for text in venn.set_labels:
            text.set_fontsize(14)
            text.set_fontweight('bold')
        
        for text in venn.subset_labels:
            if text:
                text.set_fontsize(12)
                text.set_fontweight('bold')
        
        ax.set_title("Dataset Overlap Analysis", fontsize=16, fontweight='bold', pad=20)
        
        # Convert to image for Streamlit
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=150, bbox_inches='tight')
        buf.seek(0)
        img_data = base64.b64encode(buf.read()).decode()
        plt.close()
        
        st.image(f"data:image/png;base64,{img_data}", use_container_width=True)
    except ImportError:
        st.warning("⚠️ matplotlib-venn not installed. Install with: `pip install matplotlib-venn`")
        # Fallback: Show text-based overlap
        st.info(f"""
        **Dataset Overlap:**
        - Jr2 Only: {jr2_only:,} patients
        - Overlap: {overlap_count:,} patients
        - Zo Only: {zo_only:,} patients
        """)
    except Exception as e:
        st.error(f"Error rendering Venn diagram: {e}")
        # Fallback: Show text-based overlap
        jr2_only = overlap_data.get("jr2_only", 0)
        zo_only = overlap_data.get("zo_only", 0)
        st.info(f"""
        **Dataset Overlap:**
        - Jr2 Only: {jr2_only:,} patients
        - Overlap: {overlap_count:,} patients
        - Zo Only: {zo_only:,} patients
        """)
    
    # Match rate breakdown
    st.subheader("📊 Match Rate Analysis")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.write("**From Jr2's Perspective:**")
        st.progress(match_rate_jr2 / 100)
        st.caption(f"{overlap_count:,} of {jr2_total:,} patients matched ({match_rate_jr2:.1f}%)")
    
    with col2:
        st.write("**From Zo's Perspective:**")
        st.progress(match_rate_zo / 100)
        st.caption(f"{overlap_count:,} of {zo_total:,} patients matched ({match_rate_zo:.1f}%)")
    
    # Sample ID breakdown
    if overlap_data.get("overlap_samples"):
        with st.expander("🔍 View Overlapping Sample IDs", expanded=False):
            overlap_samples = overlap_data.get("overlap_samples", [])
            st.write(f"**{len(overlap_samples)} overlapping samples:**")
            # Display in columns for better readability
            cols = st.columns(5)
            for i, sample_id in enumerate(overlap_samples[:50]):  # Show first 50
                with cols[i % 5]:
                    st.code(sample_id, language=None)
            
            if len(overlap_samples) > 50:
                st.caption(f"... and {len(overlap_samples) - 50} more samples")

