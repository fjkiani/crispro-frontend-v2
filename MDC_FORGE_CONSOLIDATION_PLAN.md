# Forge/Boltz/Gauntlet Doctrine Consolidation Plan

## Objective
To consolidate the `forge_boltz_master.mdc` and `master_forge_boltz_gauntlet_doctrine.mdc` into a single, authoritative source of truth for the Forge -> Boltz -> Gauntlet protocol, and to eliminate redundant individual doctrine files.

## Current State Analysis

### Directories
- `.cursor/rules/forge_boltz/`: Contains `forge_boltz_master.mdc` and an `archive/` sub-directory with individual doctrine files.
- `.cursor/rules/forge-doctrine/`: Contains `master_forge_boltz_gauntlet_doctrine.mdc` and individual doctrine files.

### Key Files and Findings

1.  **`forge_boltz_master.mdc`** (`.cursor/rules/forge_boltz/forge_boltz_master.mdc`)
    -   **Last Updated:** 2025-01-XX
    -   **Status:** ✅ CONSOLIDATED SOURCE OF TRUTH
    -   Lists specific files it was consolidated from in its metadata.
    -   Contains sections like "Structural Integrity Protocol," "Truth or Treason Protocol," "Boltz Integration Doctrine," "Forge Service Architecture," "Gauntlet Appropriation Strategy," "Failure Analysis & Remediation," and "Advanced Therapeutic Design."

2.  **`master_forge_boltz_gauntlet_doctrine.mdc`** (`.cursor/rules/forge-doctrine/master_forge_boltz_gauntlet_doctrine.mdc`)
    -   **Last Updated:** 2024
    -   **Status:** OPERATIONAL DIRECTIVE
    -   Organized into "CHAPTERS" (e.g., "CHAPTER 1: THE STRUCTURAL INTEGRITY PROTOCOL").
    -   Contains "Strategic Implications" section (Competitive Advantage, Market Position, Future Roadmap) not present in `forge_boltz_master.mdc`.

3.  **Individual Doctrine Files (Duplicates)**
    -   The following files are identical in content between `.cursor/rules/forge_boltz/archive/` and `.cursor/rules/forge-doctrine/`:
        -   `forge_boltz_failure_analysis_doctrine.mdc`
        -   `gauntlet_doctrine.mdc`
        -   `zeta_forge_doctrine.mdc`
        -   `zeta_forge_reclamation_doctrine.mdc`

## Consolidation Strategy

### Phase 1: Identify Master Candidate
-   The `.cursor/rules/forge-doctrine/master_forge_boltz_gauntlet_doctrine.mdc` seems to be more current ("OPERATIONAL DIRECTIVE" vs "CONSOLIDATED SOURCE OF TRUTH" which might imply it's an older consolidation) and includes "Strategic Implications" which is valuable. It also has a more structured chapter format.
-   We will designate `master_forge_boltz_gauntlet_doctrine.mdc` as the primary master file to be enhanced.

### Phase 2: Merge Content (Blocked)
Due to `.cursorignore` rules, I cannot directly read the content of `.mdc` files using the `read_file` tool. Therefore, I cannot programmatically integrate the unique content from `forge_boltz_master.mdc` into `master_forge_boltz_gauntlet_doctrine.mdc` at this time.

**Manual Action Required:**
A manual review and merge of content from `forge_boltz_master.mdc` into `master_forge_boltz_gauntlet_doctrine.mdc` is required to ensure all relevant information is consolidated. This includes:
-   Reviewing the "consolidated from" list in `forge_boltz_master.mdc` to ensure all relevant information from those archived files is present in `master_forge_boltz_gauntlet_doctrine.mdc`.
-   Reviewing section content for any details or phrasing unique to `forge_boltz_master.mdc` that should be preserved.
-   Adding the explicit "Consolidated From" list from `forge_boltz_master.mdc` to the `master_forge_boltz_gauntlet_doctrine.mdc` metadata.

### Phase 3: File Management
1.  **Archive `forge_boltz_master.mdc`:** Once its unique content is merged (manually), move `.cursor/rules/forge_boltz/forge_boltz_master.mdc` to an `archive/` sub-directory within `.cursor/rules/forge_boltz/`.
2.  **Delete Duplicate Individual Doctrines (Blocked):**
    -   Due to `.cursorignore` rules, I cannot directly delete the duplicate individual doctrine files from `.cursor/rules/forge_boltz/archive/` or `.cursor/rules/forge-doctrine/` programmatically.
    -   **Manual Action Required:** Manually delete the duplicate individual doctrine files from `.cursor/rules/forge_boltz/archive/` (e.g., `forge_boltz_failure_analysis_doctrine.mdc`, `gauntlet_doctrine.mdc`, etc.) since identical copies exist in `.cursor/rules/forge-doctrine/`. This will reduce redundancy and improve clarity.
3.  **Retain `.cursor/rules/forge-doctrine/` individual doctrines:** These will serve as the source documents for historical reference, as per the `README.md` in that directory.

## Action Plan (Revised)

1.  **Review `forge_boltz_master.mdc` for unique content and explicit "Consolidated From" list.** (Already partially done in previous steps).
2.  **MANUALLY Edit `.cursor/rules/forge-doctrine/master_forge_boltz_gauntlet_doctrine.mdc` to incorporate unique content from `forge_boltz_master.mdc`.**
3.  **MANUALLY Move `forge_boltz_master.mdc` to archive.**
4.  **MANUALLY Delete duplicate individual doctrine files from `.cursor/rules/forge_boltz/archive/`.**
5.  **Update `MDC_FORGE_CONSOLIDATION_PLAN.md` with completion status (This will be done by me once the manual steps are confirmed).

### Phase 3: File Management
1.  **Archive `forge_boltz_master.mdc`:** Once its unique content is merged, move `.cursor/rules/forge_boltz/forge_boltz_master.mdc` to an `archive/` sub-directory within `.cursor/rules/forge_boltz/`.
2.  **Delete Duplicate Individual Doctrines:**
    -   Delete the duplicate individual doctrine files from `.cursor/rules/forge_boltz/archive/` (e.g., `forge_boltz_failure_analysis_doctrine.mdc`, `gauntlet_doctrine.mdc`, etc.) since identical copies exist in `.cursor/rules/forge-doctrine/`. This will reduce redundancy and improve clarity.
3.  **Retain `.cursor/rules/forge-doctrine/` individual doctrines:** These will serve as the source documents for historical reference, as per the `README.md` in that directory.

## Action Plan

1.  **Review `forge_boltz_master.mdc` for unique content and explicit "Consolidated From" list.** (Already partially done in previous steps).
2.  **Edit `.cursor/rules/forge-doctrine/master_forge_boltz_gauntlet_doctrine.mdc` to incorporate unique content from `forge_boltz_master.mdc`.**
3.  **Move `forge_boltz_master.mdc` to archive.**
4.  **Delete duplicate individual doctrine files from `.cursor/rules/forge_boltz/archive/`.**
5.  **Update `MDC_FORGE_CONSOLIDATION_PLAN.md` with completion status.**
