"""
Validation Status Component - Statistical power and sample size analysis
"""

import sys
from pathlib import Path

# Add dashboard directory to path
dashboard_dir = Path(__file__).parent.parent
if str(dashboard_dir) not in sys.path:
    sys.path.insert(0, str(dashboard_dir))

import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from config import VALIDATION_THRESHOLDS, COLORS


def render_validation_status(stats: dict, overlap_data: dict):
    """
    Render validation readiness status with statistical power analysis
    
    Args:
        stats: Dictionary with statistics
        overlap_data: Dictionary with overlap metrics
    """
    overlap_count = stats.get("overlap_count", 0)
    min_n = VALIDATION_THRESHOLDS["min_n"]
    target_n = VALIDATION_THRESHOLDS["target_n"]
    
    # Validation status
    st.subheader("⚔️ Validation Readiness")
    
    # Status card
    if overlap_count >= min_n:
        st.success(f"✅ **VALIDATION READY** - N={overlap_count:,} exceeds minimum threshold of {min_n}")
    else:
        st.warning(f"⏸️ **PENDING** - N={overlap_count:,} below minimum threshold of {min_n}")
    
    # Sample size comparison
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric(
            label="Current Sample Size",
            value=f"{overlap_count:,}",
            delta=f"{overlap_count - min_n} over minimum" if overlap_count >= min_n else f"{min_n - overlap_count} needed",
            delta_color="normal" if overlap_count >= min_n else "off",
        )
    
    with col2:
        st.metric(
            label="Minimum Required",
            value=f"{min_n:,}",
            delta="For Chi-square test",
        )
    
    with col3:
        st.metric(
            label="Target Sample Size",
            value=f"{target_n:,}",
            delta=f"{overlap_count - target_n} difference" if overlap_count >= target_n else f"{target_n - overlap_count} to target",
            delta_color="normal" if overlap_count >= target_n else "off",
        )
    
    # Progress visualization
    st.subheader("📊 Sample Size Progress")
    
    # Create progress bar
    progress = min(overlap_count / target_n, 1.0) if target_n > 0 else 0
    st.progress(progress)
    st.caption(f"{overlap_count:,} / {target_n:,} patients ({progress*100:.1f}% of target)")
    
    # Statistical power analysis
    st.subheader("🔬 Statistical Power Analysis")
    
    power_info = {
        "Chi-square Test": {
            "min_n": 40,
            "current_n": overlap_count,
            "status": "✅ Sufficient" if overlap_count >= 40 else "❌ Insufficient",
        },
        "Fisher's Exact Test": {
            "min_n": 20,
            "current_n": overlap_count,
            "status": "✅ Sufficient" if overlap_count >= 20 else "❌ Insufficient",
        },
        "Logistic Regression": {
            "min_n": 50,
            "current_n": overlap_count,
            "status": "✅ Sufficient" if overlap_count >= 50 else "❌ Insufficient",
        },
    }
    
    # Display power analysis table
    power_df = pd.DataFrame(power_info).T
    power_df.columns = ["Minimum N", "Current N", "Status"]
    st.dataframe(power_df, use_container_width=True)
    
    # Next steps checklist
    st.subheader("📋 Validation Checklist")
    
    checklist_items = [
        ("✅" if overlap_count >= min_n else "⏸️", f"Sample size ≥{min_n} (Current: {overlap_count})"),
        ("✅" if stats.get("response_distribution") else "⏸️", "Response distribution available"),
        ("✅" if overlap_data else "⏸️", "Overlap analysis complete"),
        ("⏸️", "SAE features computed for matched patients"),
        ("⏸️", "Statistical tests executed"),
        ("⏸️", "Validation report generated"),
    ]
    
    for status, item in checklist_items:
        st.write(f"{status} {item}")

