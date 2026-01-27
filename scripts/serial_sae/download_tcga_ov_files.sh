#!/bin/bash
# Download TCGA-OV files using GDC Data Transfer Tool

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DATA_DIR="$PROJECT_ROOT/data/serial_sae/tcga_ov"
MANIFEST="$DATA_DIR/tcga_ov_rnaseq_manifest.csv"

echo "================================================================================"
echo "Download TCGA-OV Files from GDC Data Portal"
echo "================================================================================"
echo "Date: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Check if gdc-client is installed
if ! command -v gdc-client &> /dev/null; then
    echo "❌ Error: gdc-client not found"
    echo ""
    echo "Install with:"
    echo "  pip install gdc-client"
    echo ""
    echo "Or download from: https://gdc.cancer.gov/access-data/gdc-data-transfer-tool"
    exit 1
fi

echo "✅ GDC client found: $(gdc-client --version | head -1)"
echo ""

# Check if manifest exists
if [ ! -f "$MANIFEST" ]; then
    echo "❌ Error: Manifest file not found: $MANIFEST"
    echo "Run: python3 scripts/serial_sae/download_tcga_ov_gdc.py first"
    exit 1
fi

echo "📋 Manifest file: $MANIFEST"
FILE_COUNT=$(tail -n +2 "$MANIFEST" | wc -l | tr -d ' ')
echo "📊 Files to download: $FILE_COUNT"
echo ""

# Create download directory
DOWNLOAD_DIR="$DATA_DIR/downloads"
mkdir -p "$DOWNLOAD_DIR"

echo "📥 Starting download..."
echo "   Output directory: $DOWNLOAD_DIR"
echo ""

# Download files
cd "$DOWNLOAD_DIR"
gdc-client download -m "$MANIFEST" --no-verify

echo ""
echo "================================================================================"
echo "✅ Download complete"
echo "================================================================================"
echo ""
echo "Files downloaded to: $DOWNLOAD_DIR"
echo ""
echo "Next steps:"
echo "1. Process RNA-seq files for SAE computation"
echo "2. Extract expression data"
echo "3. Compute pathway scores"
