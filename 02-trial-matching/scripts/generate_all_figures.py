#!/usr/bin/env python3
"""
Master script to generate all publication figures and tables.

Usage:
    python generate_all_figures.py [--output-dir OUTPUT_DIR] [--format FORMAT]

Options:
    --output-dir: Directory to save figures (default: ../figures)
    --format: Output format - pdf, png, or both (default: both)
"""

import argparse
import sys
import os
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent.parent / "oncology-coPilot" / "oncology-backend-minimal"))

# Import figure generation scripts
from figure1_system_architecture import generate_figure1
from figure2_mechanism_fit_performance import generate_figure2
from figure3_clinical_example import generate_figure3
from figure4_ranking_accuracy import generate_figure4
from figure5_shortlist_compression import generate_figure5
from supplementary_figures import generate_all_supplementary
from generate_tables import generate_all_tables

def main():
    parser = argparse.ArgumentParser(description="Generate all publication figures and tables")
    parser.add_argument("--output-dir", type=str, default="../figures", 
                       help="Directory to save figures (default: ../figures)")
    parser.add_argument("--format", type=str, default="both", 
                       choices=["pdf", "png", "both"],
                       help="Output format (default: both)")
    parser.add_argument("--skip-supplementary", action="store_true",
                       help="Skip supplementary figures")
    parser.add_argument("--skip-tables", action="store_true",
                       help="Skip table generation")
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 80)
    print("PUBLICATION FIGURE GENERATION")
    print("=" * 80)
    print(f"Output directory: {output_dir}")
    print(f"Format: {args.format}")
    print()
    
    # Generate main figures
    print("Generating main figures...")
    print("-" * 80)
    
    try:
        print("Figure 1: System Architecture...")
        generate_figure1(output_dir, args.format)
        print("  ✅ Complete")
    except Exception as e:
        print(f"  ❌ Error: {e}")
    
    try:
        print("Figure 2: Mechanism Fit Performance...")
        generate_figure2(output_dir, args.format)
        print("  ✅ Complete")
    except Exception as e:
        print(f"  ❌ Error: {e}")
    
    try:
        print("Figure 3: Clinical Example (MBD4+TP53)...")
        generate_figure3(output_dir, args.format)
        print("  ✅ Complete")
    except Exception as e:
        print(f"  ❌ Error: {e}")
    
    try:
        print("Figure 4: Ranking Accuracy Comparison...")
        generate_figure4(output_dir, args.format)
        print("  ✅ Complete")
    except Exception as e:
        print(f"  ❌ Error: {e}")
    
    try:
        print("Figure 5: Shortlist Compression...")
        generate_figure5(output_dir, args.format)
        print("  ✅ Complete")
    except Exception as e:
        print(f"  ❌ Error: {e}")
    
    # Generate supplementary figures
    if not args.skip_supplementary:
        print()
        print("Generating supplementary figures...")
        print("-" * 80)
        try:
            generate_all_supplementary(output_dir, args.format)
            print("  ✅ Complete")
        except Exception as e:
            print(f"  ❌ Error: {e}")
    
    # Generate tables
    if not args.skip_tables:
        print()
        print("Generating tables...")
        print("-" * 80)
        try:
            table_dir = output_dir.parent / "tables"
            table_dir.mkdir(parents=True, exist_ok=True)
            generate_all_tables(table_dir)
            print("  ✅ Complete")
        except Exception as e:
            print(f"  ❌ Error: {e}")
    
    print()
    print("=" * 80)
    print("GENERATION COMPLETE")
    print("=" * 80)
    print(f"Figures saved to: {output_dir}")
    if not args.skip_tables:
        print(f"Tables saved to: {table_dir}")

if __name__ == "__main__":
    main()


