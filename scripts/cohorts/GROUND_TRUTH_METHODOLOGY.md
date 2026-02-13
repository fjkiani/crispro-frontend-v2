# Ground Truth Methodology: Timing Engine Validation

**Date:** January 31, 2026
**Subject:** Independence of Ground Truth PFI vs Computed PFI

## 1. Goal
To verify that the calculated Platinum-Free Interval (PFI) from the Timing Engine is accurate and NOT circular (i.e., not just reading the answer key).

## 2. Data Source
*   **Dataset**: Villalobos 2018 (TCGA-OV Analysis)
*   **File**: `oncology-coPilot/oncology-backend-minimal/data/TCGA-OV/ds_cci.17.00096-1.xlsx`
*   **Sheet**: 'Master clinical dataset'

## 3. Independence Check (No Leakage)

### A. The Engine's Input (Calculated)
The Timing Engine (`timing_chemo_features.py`) calculates PFI using **only** these derived dates:
1.  **Last Platinum End Date**: Derived from `last_platinum_dose_date` OR inferred from `regimen_start + TTF`.
2.  **Progression Date**: Derived from `progression_date` field.
3.  **Calculation**: `PFI = Progression Date - Last Platinum End Date`

### B. The Ground Truth (The "Answer Key")
The validation script compares against a pre-calculated column in the dataset:
*   **Column**: `Days off platinum prior to recurrence 1st line`
*   **Origin**: Manually curated by Villalobos et al. authors.
*   **Independence**: This specific column is **NEVER** passed to the Timing Engine. It is loaded separately in `validate_timing_engine_REAL.py` and used ONLY for the final `merge` comparison step.

## 4. Methodology Proof
The validation script `validate_timing_engine_REAL.py` demonstrates this separation:
1.  Loads Excel file.
2.  Extracts `ground_truth_pfi` into a separate dataframe `ground_truth`.
3.  Constructs `RegimenTable` using *other* columns (TTF, dates) but explicitly **excludes** the PFI column.
4.  Calls `build_timing_chemo_features(...)` with the clean `RegimenTable`.
5.  Merges the Engine Output with the `ground_truth` dataframe for scoring.

## 5. Circularity Check
*   Does `build_timing_chemo_features` read the Excel file? **NO**.
*   Does the input `regimen_table` contain the answer PFI? **NO**.
*   Is the calculation trivial? **NO**, it requires correctly identifying the "Platinum" regimen among multiple lines and matching it to the correct "Progression" event date.
