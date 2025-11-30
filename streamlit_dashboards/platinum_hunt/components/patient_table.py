"""
Patient Table Component - Searchable, filterable patient explorer
"""

import sys
from pathlib import Path

# Add dashboard directory to path
dashboard_dir = Path(__file__).parent.parent
if str(dashboard_dir) not in sys.path:
    sys.path.insert(0, str(dashboard_dir))

import streamlit as st
import pandas as pd
from config import RESPONSE_LABELS, RESPONSE_COLORS


def render_patient_table(df: pd.DataFrame, merged_data: list = None):
    """
    Render searchable, filterable patient table
    
    Args:
        df: DataFrame with patient data
        merged_data: Optional full merged data for detail view
    """
    if df.empty:
        st.warning("⚠️ No patient data available")
        return
    
    st.subheader("🔍 Patient Explorer")
    
    # Filters
    col1, col2, col3 = st.columns(3)
    
    with col1:
        response_filter = st.multiselect(
            "Filter by Response Type",
            options=list(RESPONSE_LABELS.values()),
            default=[],
        )
    
    with col2:
        source_filter = st.multiselect(
            "Filter by Source",
            options=df["source"].unique().tolist() if "source" in df.columns else [],
            default=[],
        )
    
    with col3:
        has_mutations_filter = st.selectbox(
            "Has Mutations",
            options=["All", "Yes", "No"],
            index=0,
        )
    
    # Apply filters
    filtered_df = df.copy()
    
    if response_filter:
        # Convert response labels back to keys for filtering
        response_keys = [k for k, v in RESPONSE_LABELS.items() if v in response_filter]
        filtered_df = filtered_df[filtered_df["platinum_response"].isin(response_keys)]
    
    if source_filter:
        filtered_df = filtered_df[filtered_df["source"].isin(source_filter)]
    
    if has_mutations_filter == "Yes":
        filtered_df = filtered_df[filtered_df["has_mutations"] == True]
    elif has_mutations_filter == "No":
        filtered_df = filtered_df[filtered_df["has_mutations"] == False]
    
    # Search
    search_term = st.text_input("🔍 Search by Sample ID or Patient ID", "")
    if search_term:
        mask = (
            filtered_df["sample_id"].str.contains(search_term, case=False, na=False) |
            filtered_df["patient_id"].str.contains(search_term, case=False, na=False)
        )
        filtered_df = filtered_df[mask]
    
    # Display results count
    st.caption(f"Showing {len(filtered_df)} of {len(df)} patients")
    
    # Display table
    display_columns = ["sample_id", "patient_id", "platinum_response", "source", "has_mutations", "mutation_count"]
    available_columns = [col for col in display_columns if col in filtered_df.columns]
    
    # Format response labels
    display_df = filtered_df[available_columns].copy()
    if "platinum_response" in display_df.columns:
        display_df["platinum_response"] = display_df["platinum_response"].map(
            lambda x: RESPONSE_LABELS.get(x, x.title())
        )
    
    # Styled table
    st.dataframe(
        display_df,
        use_container_width=True,
        height=400,
        hide_index=True,
    )
    
    # Export functionality
    col1, col2 = st.columns(2)
    
    with col1:
        csv = filtered_df.to_csv(index=False)
        st.download_button(
            label="📥 Download Filtered Data (CSV)",
            data=csv,
            file_name="platinum_response_patients.csv",
            mime="text/csv",
        )
    
    with col2:
        if merged_data:
            import json
            json_data = json.dumps([d for d in merged_data if d.get("sample_id") in filtered_df["sample_id"].tolist()], indent=2)
            st.download_button(
                label="📥 Download Filtered Data (JSON)",
                data=json_data,
                file_name="platinum_response_patients.json",
                mime="application/json",
            )

