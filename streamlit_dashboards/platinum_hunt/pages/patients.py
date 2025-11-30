"""
Patients Page - Patient explorer with search and filters
"""

import sys
from pathlib import Path

# Add dashboard directory to path
dashboard_dir = Path(__file__).parent.parent
if str(dashboard_dir) not in sys.path:
    sys.path.insert(0, str(dashboard_dir))

import streamlit as st
import pandas as pd
from data.loader import load_datasets, merge_datasets
from data.processor import create_patients_dataframe
from components.patient_table import render_patient_table


def render_patients_page():
    """Render the patient explorer page"""
    st.title("🔍 Patient Explorer")
    st.markdown("### Search and explore individual patient records")
    
    # Load data
    with st.spinner("Loading patient data..."):
        jr2_data, zo_data = load_datasets()
        
        if not jr2_data or not zo_data:
            st.error("❌ Failed to load datasets. Please check data files.")
            return
        
        merged_data = merge_datasets(jr2_data, zo_data)
        df = create_patients_dataframe(merged_data)
    
    # Patient table
    render_patient_table(df, merged_data)
    
    # Patient detail view (if row selected)
    if not df.empty:
        st.markdown("---")
        st.subheader("📋 Patient Details")
        
        # Sample ID selector
        sample_ids = df["sample_id"].tolist()
        selected_sample = st.selectbox("Select a patient to view details:", [""] + sample_ids)
        
        if selected_sample:
            patient = next((p for p in merged_data if p.get("sample_id") == selected_sample), None)
            
            if patient:
                col1, col2 = st.columns(2)
                
                with col1:
                    st.write("**Basic Information:**")
                    st.json({
                        "sample_id": patient.get("sample_id"),
                        "patient_id": patient.get("patient_id"),
                        "platinum_response": patient.get("platinum_response"),
                        "source": patient.get("source"),
                    })
                
                with col2:
                    st.write("**Mutation Data:**")
                    mutation_count = len(patient.get("mutations", [])) if patient.get("mutations") else 0
                    st.json({
                        "has_mutations": patient.get("has_mutations", False),
                        "mutation_count": mutation_count,
                        "os_months": patient.get("os_months"),
                        "stage": patient.get("stage"),
                    })
                
                # Show mutations if available
                if patient.get("mutations"):
                    st.write("**Mutations:**")
                    mutations_df = pd.DataFrame(patient["mutations"][:10])  # Show first 10
                    st.dataframe(mutations_df, use_container_width=True)
                    if len(patient["mutations"]) > 10:
                        st.caption(f"... and {len(patient['mutations']) - 10} more mutations")

