# Project Plan: Operation Zeta Shield

## 1.0 Project Title

Operation Zeta Shield: AI-Powered Stratification of Cancer Patients for Radiation Therapy

## 2.0 Objective

To validate the CrisPRO platform's ability to meaningfully stratify cancer patients for radiation therapy. We will achieve this by using the proprietary **Zeta Oracle** AI to score the functional impact of TP53 mutations and correlating these scores with patient survival outcomes from real-world clinical data.

The successful completion of this project will provide the first concrete evidence that our AI can generate clinically relevant, predictive biomarkers, forming a cornerstone for the LEAP Grant proposal.

## 3.0 Hypothesis

**Patients with TP53 mutations that receive a high-damage "Zeta Score" from our AI oracle will have significantly different survival outcomes when treated with radiation therapy compared to patients with low-damage score mutations or wild-type TP53.**

Specifically, we hypothesize that high-damage mutations impair the DNA damage response pathway, sensitizing tumors to radiation and leading to a different survival curve than tumors with functional TP53.

## 4.0 Data Sources

All data required for this analysis is located within the project repository.

1.  **Genomic Data:**
    *   **File:** `data/data_mutations.txt`
    *   **Source:** The Cancer Genome Atlas (TCGA) Pan-Cancer Atlas cohort.
    *   **Use:** To identify all non-synonymous TP53 mutations within the lung adenocarcinoma (LUAD) patient cohort.
    *   **Key Columns:** `Hugo_Symbol`, `HGVSp_Short`, `Tumor_Sample_Barcode`.

2.  **Clinical Data:**
    *   **File:** `data/luad_tcga_pan_can_atlas_2018_clinical_data.tsv`
    *   **Source:** TCGA Pan-Cancer Atlas clinical data for the LUAD cohort.
    *   **Use:** To provide patient survival outcomes and, crucially, to identify the cohort of patients who received radiation therapy.
    *   **Key Columns:** `Patient ID`, `Overall Survival (Months)`, `Overall Survival Status`, `Radiation Therapy`.

3.  **Reference Sequence Data:**
    *   **File:** `data/reference/tp53_sequence.fasta`
    *   **Source:** Curated reference coding sequence (CDS) for human TP53.
    *   **Use:** To serve as the "ground truth" wild-type sequence against which all mutated sequences are compared by the Zeta Oracle.

## 5.0 Methodology (Execution Plan)

This operation will be conducted with a rigorous, verification-first approach to prevent the failures of the previous "Operation Zeta Strike".

**Phase 1: System Validation (Complete)**

*   **1.1: Oracle Logic Repair:** The bug within `services/command_center/main.py` that passed incorrect data to the Zeta Oracle has been fixed. The system now correctly generates a mutated DNA sequence based on the specific protein change requested.
*   **1.2: Analysis Script Reinforcement:** The analysis script, `tools/zeta_striker.py`, has been reinforced with a `_perform_sanity_check` function. This function will automatically verify the scores of known pathogenic TP53 mutations (e.g., R175H, R249M) and will halt the entire analysis if the scores are not plausible.

**Phase 2: Data Processing & AI Analysis**

*   **2.1: Full Execution:** The reinforced `zeta_striker.py` script will be executed.
*   **2.2: Mutation Scoring:** The script will iterate through all unique TP53 mutations found in the genomic data file. For each mutation, it will make a request to the `CommandCenter`'s `/workflow/assess_threat` endpoint.
*   **2.3: Intelligence Fusion:** The returned Zeta Scores will be merged with the clinical data, creating a master analysis file: `results/master_clinical_zeta_scores_corrected.tsv`. Each patient will have an associated Zeta Score (0 for Wild-Type, or the score of their specific mutation).

**Phase 3: Survival Analysis & Hypothesis Testing**

*   **3.1: Cohort Isolation:** The master analysis file will be filtered to isolate only those patients who received radiation therapy (`Radiation Therapy == 'Yes'`).
*   **3.2: Damage Stratification:** This radiation cohort will be stratified into two groups:
    *   **High Damage:** Patients with a `zeta_score <= -50`.
    *   **Low Damage:** Patients with a `zeta_score > -50` (this includes Wild-Type patients).
*   **3.3: Statistical Analysis:** A Kaplan-Meier survival analysis will be performed to visualize the survival curves for the High Damage vs. Low Damage groups. A log-rank test will be used to calculate a p-value, determining if the difference between the curves is statistically significant.

## 6.0 Success Criteria

The operation will be deemed a success if the following criteria are met:

1.  **Technical Success:** The `zeta_striker.py` script runs to completion without being halted by the sanity checks.
2.  **Qualitative Validation:** The output `TP53_zeta_scores_corrected.json` shows meaningful, non-zero scores for known pathogenic missense mutations.
3.  **Primary Objective Met:** The final `zeta_shield_analysis.png` Kaplan-Meier plot shows a clear separation between the survival curves of the "High Damage" and "Low Damage" cohorts.
4.  **Hypothesis Validated (Stretch Goal):** The log-rank test returns a p-value less than 0.05, providing statistical validation of our core hypothesis. 