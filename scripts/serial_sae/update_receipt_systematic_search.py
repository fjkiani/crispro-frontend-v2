#!/usr/bin/env python3
"""Update mission receipt with systematic search results."""

import json
from pathlib import Path
from datetime import datetime

# Update mission receipt
receipt_path = Path("receipts/sae_serial_framework.json")

with open(receipt_path) as f:
    receipt = json.load(f)

# Update Phase 2 with systematic search results
receipt["phases_completed"]["phase2"]["status"] = "complete"
receipt["phases_completed"]["phase2"]["systematic_search"] = {
    "databases_searched": ["GENIE-BPC", "GEO/SRA", "dbGAP", "PubMed/PMC", "JCI"],
    "total_datasets_evaluated": 5,
    "validation_studies_reviewed": 2,
    "top_datasets": {
        "msk_spectrum": {
            "score": 9.5,
            "access": "dbGAP phs002857.v3.p1",
            "status": "requires_application",
            "priority": "highest"
        },
        "gse165897": {
            "score": 8.7,
            "access": "GEO (public)",
            "status": "ready_for_download",
            "priority": "high"
        }
    }
}

receipt["phases_completed"]["phase2"]["findings"]["systematic_search_complete"] = True
receipt["phases_completed"]["phase2"]["findings"]["recommended_datasets"] = [
    "MSK_SPECTRUM (dbGAP) - 9.5/10 - Requires application",
    "GSE165897 (GEO) - 8.7/10 - Public, ready for immediate download"
]
receipt["phases_completed"]["phase2"]["findings"]["tcga_ov_correlation"] = {
    "status": "complete",
    "patients_analyzed": 426,
    "significant_pathways": ["HER2", "PI3K", "Efflux", "RAS_MAPK"]
}

receipt["mission_status"] = "Phase 2 complete - Systematic search finished. Top datasets: MSK_SPECTRUM (dbGAP), GSE165897 (GEO). TCGA-OV correlation complete (426 patients)."

with open(receipt_path, "w") as f:
    json.dump(receipt, f, indent=2)

print("✅ Mission receipt updated with systematic search results")
print("\n📊 Phase 2 Summary:")
print("  - Systematic search: COMPLETE")
print("  - Top datasets identified:")
print("    1. MSK_SPECTRUM (dbGAP) - 9.5/10 - Requires application")
print("    2. GSE165897 (GEO) - 8.7/10 - Public, ready for download")
print("  - TCGA-OV correlation: COMPLETE (426 patients)")
print("\n📋 Next Steps:")
print("  1. Download and process GSE165897 (immediate, public access)")
print("  2. Submit dbGAP application for MSK_SPECTRUM")
print("  3. Integrate findings with existing TCGA-OV analysis")
