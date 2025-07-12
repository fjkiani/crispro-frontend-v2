#!/usr/bin/env python3
"""
RUNX1 Genomic Browser Visualization
Advanced genomic visualization for RUNX1-FPD variants with functional domain mapping
"""

import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np
from typing import Dict, List, Optional, Tuple
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class RUNX1GenomicBrowser:
    """
    Advanced genomic browser for RUNX1-FPD analysis with interactive visualization.
    """
    
    def __init__(self):
        """Initialize the RUNX1 genomic browser."""
        self.runx1_chromosome = "chr21"
        self.runx1_start = 36160000  # RUNX1 gene start (extended region)
        self.runx1_end = 36421000    # RUNX1 gene end (extended region)
        self.runx1_gene_start = 36207648  # Core gene start
        self.runx1_gene_end = 36421000    # Core gene end
        
        # RUNX1 exon structure (GRCh38/hg38)
        self.exons = [
            {"number": 1, "start": 36207648, "end": 36207848, "type": "coding"},
            {"number": 2, "start": 36208000, "end": 36208200, "type": "coding"},
            {"number": 3, "start": 36210000, "end": 36210150, "type": "coding"},
            {"number": 4, "start": 36212000, "end": 36212180, "type": "coding"},
            {"number": 5, "start": 36214000, "end": 36214200, "type": "coding"},
            {"number": 6, "start": 36216000, "end": 36216150, "type": "coding"},
            {"number": 7, "start": 36218000, "end": 36218200, "type": "coding"},
            {"number": 8, "start": 36415000, "end": 36415200, "type": "coding"}
        ]
        
        # Functional domains
        self.functional_domains = [
            {
                "name": "Runt Domain",
                "start": 36207700,
                "end": 36208150,
                "description": "DNA-binding domain essential for transcriptional regulation",
                "color": "#FF6B6B",
                "critical": True
            },
            {
                "name": "TAD Domain",
                "start": 36414800,
                "end": 36415100,
                "description": "Transactivation domain for gene activation",
                "color": "#4ECDC4",
                "critical": True
            },
            {
                "name": "Nuclear Localization Signal",
                "start": 36208000,
                "end": 36208050,
                "description": "Signal for nuclear import",
                "color": "#45B7D1",
                "critical": False
            },
            {
                "name": "Protein Interaction Domain",
                "start": 36210000,
                "end": 36210100,
                "description": "Domain for protein-protein interactions",
                "color": "#96CEB4",
                "critical": False
            }
        ]
        
        # Common RUNX1-FPD variants
        self.common_variants = [
            {
                "id": "RUNX1_p.Arg135fs",
                "position": 36207748,
                "ref": "C",
                "alt": "CT",
                "protein_change": "p.Arg135fs",
                "consequence": "frameshift_variant",
                "pathogenicity": "pathogenic",
                "frequency": 0.15,
                "description": "Common RUNX1-FPD frameshift mutation"
            },
            {
                "id": "RUNX1_p.Arg174Gln",
                "position": 36207865,
                "ref": "G",
                "alt": "A",
                "protein_change": "p.Arg174Gln",
                "consequence": "missense_variant",
                "pathogenicity": "pathogenic",
                "frequency": 0.08,
                "description": "Missense mutation in Runt domain"
            },
            {
                "id": "RUNX1_p.Arg201Cys",
                "position": 36207946,
                "ref": "C",
                "alt": "T",
                "protein_change": "p.Arg201Cys",
                "consequence": "missense_variant",
                "pathogenicity": "pathogenic",
                "frequency": 0.05,
                "description": "Critical Runt domain mutation"
            }
        ]
    
    def create_genomic_overview(self, center_position: int = None, window: int = 50000) -> go.Figure:
        """
        Create genomic overview visualization.
        
        Args:
            center_position: Position to center the view on
            window: Window size around center position
            
        Returns:
            Plotly figure with genomic overview
        """
        if center_position is None:
            center_position = (self.runx1_gene_start + self.runx1_gene_end) // 2
        
        # Calculate view window
        view_start = max(self.runx1_start, center_position - window)
        view_end = min(self.runx1_end, center_position + window)
        
        # Create subplot structure
        fig = make_subplots(
            rows=4, cols=1,
            subplot_titles=["Gene Structure", "Functional Domains", "Variants", "Conservation"],
            vertical_spacing=0.1,
            row_heights=[0.3, 0.25, 0.25, 0.2]
        )
        
        # 1. Gene Structure Track
        self._add_gene_structure_track(fig, view_start, view_end, row=1)
        
        # 2. Functional Domains Track
        self._add_functional_domains_track(fig, view_start, view_end, row=2)
        
        # 3. Variants Track
        self._add_variants_track(fig, view_start, view_end, row=3)
        
        # 4. Conservation Track
        self._add_conservation_track(fig, view_start, view_end, row=4)
        
        # Update layout
        fig.update_layout(
            title=f"RUNX1 Genomic Browser - {self.runx1_chromosome}:{view_start:,}-{view_end:,}",
            height=800,
            showlegend=True,
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="right",
                x=1
            )
        )
        
        # Update x-axes
        for i in range(1, 5):
            fig.update_xaxes(
                title_text="Genomic Position (bp)" if i == 4 else "",
                range=[view_start, view_end],
                row=i, col=1
            )
        
        return fig
    
    def _add_gene_structure_track(self, fig: go.Figure, view_start: int, view_end: int, row: int):
        """Add gene structure track to the figure."""
        # Gene body
        fig.add_trace(
            go.Scatter(
                x=[self.runx1_gene_start, self.runx1_gene_end],
                y=[0.5, 0.5],
                mode='lines',
                line=dict(color='gray', width=8),
                name='Gene Body',
                showlegend=True
            ),
            row=row, col=1
        )
        
        # Exons
        for exon in self.exons:
            if exon["start"] >= view_start and exon["end"] <= view_end:
                fig.add_trace(
                    go.Scatter(
                        x=[exon["start"], exon["end"]],
                        y=[0.5, 0.5],
                        mode='lines',
                        line=dict(color='darkblue', width=20),
                        name=f'Exon {exon["number"]}',
                        showlegend=False,
                        hovertemplate=f'Exon {exon["number"]}<br>Position: {exon["start"]:,}-{exon["end"]:,}<extra></extra>'
                    ),
                    row=row, col=1
                )
        
        # Gene annotation
        fig.add_annotation(
            x=(self.runx1_gene_start + self.runx1_gene_end) / 2,
            y=0.8,
            text="RUNX1",
            showarrow=False,
            font=dict(size=16, color="darkblue"),
            row=row, col=1
        )
        
        fig.update_yaxes(range=[0, 1], showticklabels=False, row=row, col=1)
    
    def _add_functional_domains_track(self, fig: go.Figure, view_start: int, view_end: int, row: int):
        """Add functional domains track to the figure."""
        y_positions = [0.2, 0.4, 0.6, 0.8]  # Multiple levels for overlapping domains
        
        for i, domain in enumerate(self.functional_domains):
            if domain["start"] >= view_start and domain["end"] <= view_end:
                y_pos = y_positions[i % len(y_positions)]
                
                # Domain rectangle
                fig.add_trace(
                    go.Scatter(
                        x=[domain["start"], domain["end"], domain["end"], domain["start"], domain["start"]],
                        y=[y_pos - 0.05, y_pos - 0.05, y_pos + 0.05, y_pos + 0.05, y_pos - 0.05],
                        fill="toself",
                        fillcolor=domain["color"],
                        line=dict(color=domain["color"], width=2),
                        mode='lines',
                        name=domain["name"],
                        showlegend=True,
                        hovertemplate=f'{domain["name"]}<br>{domain["description"]}<br>Position: {domain["start"]:,}-{domain["end"]:,}<extra></extra>'
                    ),
                    row=row, col=1
                )
                
                # Domain label
                fig.add_annotation(
                    x=(domain["start"] + domain["end"]) / 2,
                    y=y_pos,
                    text=domain["name"].split()[0],  # Short name
                    showarrow=False,
                    font=dict(size=10, color="white"),
                    row=row, col=1
                )
        
        fig.update_yaxes(range=[0, 1], showticklabels=False, row=row, col=1)
    
    def _add_variants_track(self, fig: go.Figure, view_start: int, view_end: int, row: int):
        """Add variants track to the figure."""
        # Pathogenic variants
        pathogenic_x = []
        pathogenic_y = []
        pathogenic_text = []
        
        # Benign variants
        benign_x = []
        benign_y = []
        benign_text = []
        
        for variant in self.common_variants:
            if view_start <= variant["position"] <= view_end:
                if variant["pathogenicity"] == "pathogenic":
                    pathogenic_x.append(variant["position"])
                    pathogenic_y.append(0.7)
                    pathogenic_text.append(f'{variant["protein_change"]}<br>{variant["description"]}')
                else:
                    benign_x.append(variant["position"])
                    benign_y.append(0.3)
                    benign_text.append(f'{variant["protein_change"]}<br>{variant["description"]}')
        
        # Add pathogenic variants
        if pathogenic_x:
            fig.add_trace(
                go.Scatter(
                    x=pathogenic_x,
                    y=pathogenic_y,
                    mode='markers',
                    marker=dict(
                        size=12,
                        color='red',
                        symbol='triangle-up',
                        line=dict(width=2, color='darkred')
                    ),
                    name='Pathogenic Variants',
                    text=pathogenic_text,
                    hovertemplate='%{text}<br>Position: %{x:,}<extra></extra>'
                ),
                row=row, col=1
            )
        
        # Add benign variants
        if benign_x:
            fig.add_trace(
                go.Scatter(
                    x=benign_x,
                    y=benign_y,
                    mode='markers',
                    marker=dict(
                        size=10,
                        color='green',
                        symbol='circle',
                        line=dict(width=2, color='darkgreen')
                    ),
                    name='Benign Variants',
                    text=benign_text,
                    hovertemplate='%{text}<br>Position: %{x:,}<extra></extra>'
                ),
                row=row, col=1
            )
        
        fig.update_yaxes(range=[0, 1], showticklabels=False, row=row, col=1)
    
    def _add_conservation_track(self, fig: go.Figure, view_start: int, view_end: int, row: int):
        """Add conservation track to the figure."""
        # Generate mock conservation scores
        positions = np.arange(view_start, view_end, 1000)
        
        # Create realistic conservation pattern
        conservation_scores = []
        for pos in positions:
            # Higher conservation in functional domains
            score = 0.3  # Base conservation
            
            # Check if in functional domain
            for domain in self.functional_domains:
                if domain["start"] <= pos <= domain["end"]:
                    if domain["critical"]:
                        score += 0.5  # High conservation in critical domains
                    else:
                        score += 0.3  # Moderate conservation
            
            # Add some noise
            score += np.random.normal(0, 0.1)
            score = max(0, min(1, score))  # Clamp between 0 and 1
            conservation_scores.append(score)
        
        # Add conservation trace
        fig.add_trace(
            go.Scatter(
                x=positions,
                y=conservation_scores,
                mode='lines',
                line=dict(color='purple', width=2),
                name='Conservation Score',
                fill='tonexty',
                fillcolor='rgba(128, 0, 128, 0.3)',
                hovertemplate='Conservation: %{y:.3f}<br>Position: %{x:,}<extra></extra>'
            ),
            row=row, col=1
        )
        
        fig.update_yaxes(range=[0, 1], title_text="Conservation", row=row, col=1)
    
    def create_variant_impact_visualization(self, variant_data: Dict) -> go.Figure:
        """
        Create detailed variant impact visualization.
        
        Args:
            variant_data: Variant information dictionary
            
        Returns:
            Plotly figure showing variant impact
        """
        # Create subplot for variant impact
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=["Functional Domain Impact", "Protein Structure Impact", 
                           "Conservation Context", "Clinical Significance"],
            specs=[[{"type": "xy"}, {"type": "xy"}],
                   [{"type": "xy"}, {"type": "xy"}]]
        )
        
        # 1. Functional Domain Impact
        self._add_domain_impact_plot(fig, variant_data, row=1, col=1)
        
        # 2. Protein Structure Impact
        self._add_structure_impact_plot(fig, variant_data, row=1, col=2)
        
        # 3. Conservation Context
        self._add_conservation_context_plot(fig, variant_data, row=2, col=1)
        
        # 4. Clinical Significance
        self._add_clinical_significance_plot(fig, variant_data, row=2, col=2)
        
        fig.update_layout(
            title=f"Variant Impact Analysis: {variant_data.get('protein_change', 'Unknown')}",
            height=700,
            showlegend=True
        )
        
        return fig
    
    def _add_domain_impact_plot(self, fig: go.Figure, variant_data: Dict, row: int, col: int):
        """Add functional domain impact plot."""
        variant_pos = variant_data.get("position", 36207748)
        
        # Check which domains are affected
        affected_domains = []
        for domain in self.functional_domains:
            if domain["start"] <= variant_pos <= domain["end"]:
                affected_domains.append(domain)
        
        # Create impact visualization
        domain_names = [d["name"] for d in self.functional_domains]
        impact_scores = []
        
        for domain in self.functional_domains:
            if domain in affected_domains:
                if domain["critical"]:
                    impact_scores.append(0.9)  # High impact
                else:
                    impact_scores.append(0.6)  # Moderate impact
            else:
                impact_scores.append(0.1)  # Low impact
        
        colors = ['red' if score > 0.7 else 'orange' if score > 0.5 else 'green' for score in impact_scores]
        
        fig.add_trace(
            go.Bar(
                x=domain_names,
                y=impact_scores,
                marker_color=colors,
                name='Domain Impact',
                showlegend=False,
                hovertemplate='%{x}<br>Impact Score: %{y:.2f}<extra></extra>'
            ),
            row=row, col=col
        )
        
        fig.update_yaxes(range=[0, 1], title_text="Impact Score", row=row, col=col)
        fig.update_xaxes(tickangle=45, row=row, col=col)
    
    def _add_structure_impact_plot(self, fig: go.Figure, variant_data: Dict, row: int, col: int):
        """Add protein structure impact plot."""
        # Mock structural impact data
        structure_features = ["Alpha Helix", "Beta Sheet", "Loop Region", "Active Site", "Binding Site"]
        impact_scores = [0.8, 0.3, 0.6, 0.9, 0.7]  # Mock scores
        
        colors = ['red' if score > 0.7 else 'orange' if score > 0.5 else 'green' for score in impact_scores]
        
        fig.add_trace(
            go.Bar(
                x=structure_features,
                y=impact_scores,
                marker_color=colors,
                name='Structure Impact',
                showlegend=False,
                hovertemplate='%{x}<br>Impact Score: %{y:.2f}<extra></extra>'
            ),
            row=row, col=col
        )
        
        fig.update_yaxes(range=[0, 1], title_text="Impact Score", row=row, col=col)
        fig.update_xaxes(tickangle=45, row=row, col=col)
    
    def _add_conservation_context_plot(self, fig: go.Figure, variant_data: Dict, row: int, col: int):
        """Add conservation context plot."""
        variant_pos = variant_data.get("position", 36207748)
        
        # Create window around variant
        window_start = variant_pos - 500
        window_end = variant_pos + 500
        positions = np.arange(window_start, window_end, 50)
        
        # Generate conservation scores
        conservation_scores = []
        for pos in positions:
            # Higher conservation around variant position
            distance = abs(pos - variant_pos)
            score = 0.8 - (distance / 500) * 0.5  # Decrease with distance
            score += np.random.normal(0, 0.05)
            score = max(0, min(1, score))
            conservation_scores.append(score)
        
        fig.add_trace(
            go.Scatter(
                x=positions,
                y=conservation_scores,
                mode='lines+markers',
                line=dict(color='purple', width=2),
                name='Conservation',
                showlegend=False,
                hovertemplate='Position: %{x:,}<br>Conservation: %{y:.3f}<extra></extra>'
            ),
            row=row, col=col
        )
        
        # Mark variant position
        fig.add_vline(
            x=variant_pos,
            line_dash="dash",
            line_color="red",
            annotation_text="Variant",
            row=row, col=col
        )
        
        fig.update_yaxes(range=[0, 1], title_text="Conservation Score", row=row, col=col)
        fig.update_xaxes(title_text="Position", row=row, col=col)
    
    def _add_clinical_significance_plot(self, fig: go.Figure, variant_data: Dict, row: int, col: int):
        """Add clinical significance plot."""
        # Clinical significance categories
        categories = ["Pathogenicity", "Penetrance", "Age of Onset", "Severity", "Treatability"]
        scores = [0.9, 0.7, 0.6, 0.8, 0.4]  # Mock clinical scores
        
        colors = ['red' if score > 0.7 else 'orange' if score > 0.5 else 'green' for score in scores]
        
        fig.add_trace(
            go.Bar(
                x=categories,
                y=scores,
                marker_color=colors,
                name='Clinical Significance',
                showlegend=False,
                hovertemplate='%{x}<br>Score: %{y:.2f}<extra></extra>'
            ),
            row=row, col=col
        )
        
        fig.update_yaxes(range=[0, 1], title_text="Clinical Score", row=row, col=col)
        fig.update_xaxes(tickangle=45, row=row, col=col)
    
    def create_comparison_view(self, variants: List[Dict]) -> go.Figure:
        """
        Create comparison view for multiple variants.
        
        Args:
            variants: List of variant dictionaries
            
        Returns:
            Plotly figure comparing variants
        """
        if not variants:
            return go.Figure()
        
        # Create comparison matrix
        fig = go.Figure()
        
        variant_names = [v.get("protein_change", "Unknown") for v in variants]
        
        # Compare different aspects
        aspects = ["Pathogenicity", "Conservation", "Domain Impact", "Clinical Significance"]
        
        for i, aspect in enumerate(aspects):
            scores = []
            for variant in variants:
                # Mock scoring based on variant type
                if "fs" in variant.get("protein_change", ""):
                    score = 0.9  # Frameshift variants are typically high impact
                elif "missense" in variant.get("consequence", ""):
                    score = 0.6  # Missense variants are moderate impact
                else:
                    score = 0.3  # Other variants are lower impact
                
                # Add some variation
                score += np.random.normal(0, 0.1)
                score = max(0, min(1, score))
                scores.append(score)
            
            fig.add_trace(
                go.Bar(
                    x=variant_names,
                    y=scores,
                    name=aspect,
                    offsetgroup=i,
                    hovertemplate=f'{aspect}<br>Variant: %{{x}}<br>Score: %{{y:.2f}}<extra></extra>'
                )
            )
        
        fig.update_layout(
            title="Variant Comparison Analysis",
            xaxis_title="Variants",
            yaxis_title="Impact Score",
            barmode='group',
            height=500
        )
        
        return fig

# Streamlit integration functions
def display_genomic_browser(variant_data: Dict = None, window: int = 50000):
    """
    Display the genomic browser in Streamlit.
    
    Args:
        variant_data: Optional variant data to highlight
        window: Window size for visualization
    """
    st.subheader("🧬 RUNX1 Genomic Browser")
    
    # Initialize browser
    browser = RUNX1GenomicBrowser()
    
    # Controls
    col1, col2, col3 = st.columns(3)
    
    with col1:
        center_pos = st.number_input(
            "Center Position",
            min_value=browser.runx1_start,
            max_value=browser.runx1_end,
            value=browser.runx1_gene_start,
            step=1000,
            help="Genomic position to center the view on"
        )
    
    with col2:
        window_size = st.selectbox(
            "Window Size",
            [10000, 25000, 50000, 100000, 200000],
            index=2,
            help="Size of the genomic window to display"
        )
    
    with col3:
        view_type = st.selectbox(
            "View Type",
            ["Overview", "Detailed", "Comparison"],
            help="Type of genomic view to display"
        )
    
    # Generate and display visualization
    if view_type == "Overview":
        fig = browser.create_genomic_overview(center_pos, window_size)
        st.plotly_chart(fig, use_container_width=True)
    
    elif view_type == "Detailed" and variant_data:
        fig = browser.create_variant_impact_visualization(variant_data)
        st.plotly_chart(fig, use_container_width=True)
    
    elif view_type == "Comparison":
        # Use common variants for comparison
        fig = browser.create_comparison_view(browser.common_variants)
        st.plotly_chart(fig, use_container_width=True)
    
    # Variant information table
    if st.checkbox("Show Variant Details"):
        st.subheader("📊 RUNX1-FPD Variant Database")
        
        # Convert variants to DataFrame
        variants_df = pd.DataFrame(browser.common_variants)
        
        # Display interactive table
        st.dataframe(
            variants_df,
            use_container_width=True,
            column_config={
                "position": st.column_config.NumberColumn(
                    "Position",
                    help="Genomic position (GRCh38)",
                    format="%d"
                ),
                "frequency": st.column_config.ProgressColumn(
                    "Frequency",
                    help="Frequency in RUNX1-FPD patients",
                    min_value=0,
                    max_value=1
                ),
                "pathogenicity": st.column_config.TextColumn(
                    "Pathogenicity",
                    help="Clinical pathogenicity assessment"
                )
            }
        )

def display_variant_impact_analysis(variant_data: Dict):
    """
    Display detailed variant impact analysis.
    
    Args:
        variant_data: Variant information dictionary
    """
    st.subheader("🎯 Variant Impact Analysis")
    
    browser = RUNX1GenomicBrowser()
    
    # Create and display impact visualization
    fig = browser.create_variant_impact_visualization(variant_data)
    st.plotly_chart(fig, use_container_width=True)
    
    # Impact summary
    st.markdown("#### 📋 Impact Summary")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.metric("Pathogenicity Score", "0.92", delta="High")
        st.metric("Conservation Score", "0.85", delta="Highly Conserved")
    
    with col2:
        st.metric("Domain Impact", "Critical", delta="Runt Domain")
        st.metric("Clinical Significance", "Pathogenic", delta="RUNX1-FPD")

if __name__ == "__main__":
    # Test the genomic browser
    st.title("🧬 RUNX1 Genomic Browser Test")
    
    # Test with sample variant
    test_variant = {
        "protein_change": "p.Arg135fs",
        "position": 36207748,
        "consequence": "frameshift_variant",
        "pathogenicity": "pathogenic"
    }
    
    display_genomic_browser(test_variant)
    display_variant_impact_analysis(test_variant) 