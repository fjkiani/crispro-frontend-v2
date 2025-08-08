# After-Action Report: Operation Guardian

## 1. Mission Objective

To validate the `CommandCenter` API and the Triumvirate assessment protocol by conducting an end-to-end survival analysis on the TCGA Ovarian Cancer cohort, specifically testing the hypothesis that pathogenic `BRCA1` mutations correlate with patient survival outcomes.

## 2. Executive Summary

Operation Guardian was a resounding technical success, but was fraught with a cascading series of complex, intertwined bugs that tested the limits of our debugging capabilities. The mission ultimately succeeded, proving the viability of our automated analysis engine. However, the path to success revealed critical weaknesses in our client-server contracts, error handling visibility, and initial server-side logic. The core achievement was the transformation of a multi-day manual bioinformatics task into a sub-minute, fully automated pipeline. The key lesson is that client-side parsing of server responses is as critical as the server logic itself, and failures in this contract can create deeply misleading "ghost" bugs.

## 3. Timeline of Failures & Resolutions

The operation can be characterized by a "bug cascade," where fixing one issue immediately revealed a deeper, underlying problem.

*   **Failure 1: `500 Internal Server Error` (Initial API Calls)**
    *   **Symptom:** The first analysis script failed with generic 500 errors.
    *   **Root Cause:** The deployed `CommandCenter` service was missing the `BRCA1.fasta` reference sequence file in its container.
    *   **Resolution:** Added the `data/reference` directory to the Modal app's image build.

*   **Failure 2: `KeyError: 'Patient ID'` (Data Merge)**
    *   **Symptom:** The unified analysis script failed immediately during the data merging step.
    *   **Root Cause:** A mismatch between patient identifier columns. The mutation file used `Tumor_Sample_Barcode`, while the clinical file used `Patient ID`.
    *   **Resolution:** Implemented logic in the client script to create a common `Patient ID` key by slicing the `Tumor_Sample_Barcode`.

*   **Failure 3: The "Ghost" Bug (`Pathogenic Cohort Size: 0`)**
    *   **Symptom:** The script ran to completion but produced a meaningless result, with zero patients classified as having pathogenic mutations. This was the most persistent and misleading failure.
    *   **Initial Diagnosis (Incorrect):** Assumed file I/O errors or data type mismatches (string vs. boolean) between script executions. This led to multiple failed attempts at fixing the wrong problem.
    *   **Second Diagnosis (Partially Correct):** Discovered a placeholder `_apply_mutation_to_sequence` function on the server that returned un-mutated DNA, short-circuiting the `Truncation Sieve`.
    *   **Third Diagnosis (Partially Correct):** After fixing the placeholder, a `not enough values to unpack` error occurred due to brittle regex, which was hidden by a generic `try...except` block on the server.
    *   **Fourth Diagnosis (Partially Correct):** After fixing the regex, the logic for comparing stop codon positions in the `Truncation Sieve` was found to be flawed.
    *   **FINAL ROOT CAUSE:** After fixing all server-side logic, the problem persisted. Server logs finally proved the `CommandCenter` was correctly identifying mutations as pathogenic and returning `True`. The final bug was in the **client script**. The client was not parsing the nested JSON response from the server (`{"assessment": {"is_pathogenic": true}}`) and was instead looking for `is_pathogenic` at the top level, causing it to default to `False` every time.

## 4. Lessons Learned

1.  **Define and Enforce API Contracts:** The client-server contract is paramount. The final, most elusive bug was a simple discrepancy between the structure of the JSON the server sent and what the client expected. This contract must be documented and validated.
2.  **Fail Loudly and Specifically:** The generic `try...except Exception` block in the `CommandCenter`'s main assessment function was disastrous. It swallowed critical errors (like the regex `unpack` error) and returned a default "non-pathogenic" result, making debugging nearly impossible. Errors should be allowed to propagate or be caught and re-raised with more context.
3.  **The Client is Part of the System:** We spent days debugging the server when it was, in the end, working perfectly. The assumption that the server was the sole point of failure was incorrect. A full-system view, including the client, is essential.
4.  **Log-based Debugging is Invaluable:** The breakthrough only occurred when we added detailed, step-by-step logging to the server-side functions. This, combined with streaming the service logs (`modal app logs`), was the only way to prove the server was working and pivot our attention back to the client. 