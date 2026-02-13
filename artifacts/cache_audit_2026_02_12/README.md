# Cache Audit: 2026-02-12

**Source File**: `~/.cache/oncology_ayesha_fit/therapy_fit_previews_v1.json`
**Extraction Date**: 2026-02-12
**Source Timestamp**: 2026-02-11T03:00:41Z (Most likely corresponds to the user's "2/12" working session data)

## Contents

This folder contains the separated payloads extracted from the Ayesha Therapy Fit cache:

- **`L1_payload.json`**: The Level 1 (Baseline) analysis payload.
  - Contains top drugs (Olaparib, etc.), mechanism panel, and confidence scores.

- **`L2_payloads.json`**: The Level 2 (Scenario) analysis payloads.
  - Contains variations like `L2A_HRDhi_TMBhi`, `L2B_HRDhighTMBlow`, etc.

- **`L3_combinations_payloads.json`**: The Level 3 (Deep Simulation) payloads.
  - Contains the combinatorial expansion of L2 scenarios into L3 states (e.g. `L2A` x `L3C`).

- **`Metadata.json`**: Top-level cache metadata (hashes, timestamps, ttl).

## Usage
These files represent the specific "agent-generated" output stored in the backend cache, allowing for offline inspection of the therapy recommendations and mechanism scoring logic.
