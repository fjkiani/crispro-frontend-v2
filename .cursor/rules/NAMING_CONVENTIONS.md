# Naming Conventions Document

## Objective
To establish clear and consistent naming conventions for files and directories within the `.cursor/rules/` directory to improve readability, searchability, and overall organization.

## General Principles

1.  **Clarity and Conciseness**: Names should clearly convey the content or purpose of the file/directory without being excessively long.
2.  **Consistency**: Adhere to a single naming style across the entire directory.
3.  **Readability**: Avoid abbreviations where possible, or use widely recognized ones.
4.  **Searchability**: Use keywords that are likely to be searched for.

## File Naming Conventions

### Style: `SNAKE_CASE_WITH_CAPS.mdc` or `snake_case_with_no_caps.md`

-   **Default for Doctrines/Plans/Masters**: Use `SNAKE_CASE_WITH_CAPS` for critical doctrine, master, or plan files (e.g., `MASTER_FORGE_BOLTZ_GAUNTLET_DOCTRINE.mdc`, `PLATFORM_EVOLUTION_MASTER_PLAN.mdc`).
-   **Default for General Documentation/Blogs**: Use `snake_case_with_no_caps` for general documentation, blogs, reports, or less critical files (e.g., `operation_guardian_AAR.md`, `tox101-blog.mdc`).
-   **Consistency in Extensions**: Use `.mdc` for Cursor-specific markdown files with metadata, and `.md` for standard markdown files.

### Prefixes and Suffixes

-   **`MASTER_`**: For files that consolidate information from multiple sources or serve as a primary reference (e.g., `CLINICAL_TRIALS_MASTER_DOCUMENT.md`).
-   **`_DOCTRINE`**: For files detailing core principles, protocols, or frameworks (e.g., `ZETA_FORGE_DOCTRINE.mdc`).
-   **`_PLAN`**: For files outlining execution steps, strategies, or roadmaps (e.g., `DOCTRINE_MODULARIZATION_PLAN.md`).
-   **`_REPORT` / `_AAR`**: For after-action reports or summary documents (e.g., `operation_guardian_AAR.md`).
-   **`00_MASTER_INDEX`**: For primary index files within a directory (e.g., `00_MASTER_INDEX.mdc`).

### Dates and Versions

-   Avoid including dates or version numbers in file names unless absolutely necessary for historical tracking, and prefer to include them in the file's metadata or content.

## Directory Naming Conventions

### Style: `kebab-case` or `snake_case`

-   **Default**: Use `kebab-case` for most directories (e.g., `forge-doctrine`, `clinical-trials-agents`).
-   **Exceptions**: `snake_case` can be used for directories that contain a specific logical grouping (e.g., `forge_boltz`). Strive for consistency within sub-directories.

### Structure

-   **Plural Nouns**: Use plural nouns for directories containing multiple related files (e.g., `doctrines/`, `agents/`, `concepts/`).
-   **Singular Nouns**: Use singular nouns for directories representing a single entity or a distinct, focused topic (e.g., `ayesha/`, `MOAT/`).

## Examples of Good Naming (Existing)

-   `.cursor/rules/forge-doctrine/master_forge_boltz_gauntlet_doctrine.mdc`
-   `.cursor/rules/MM/INDEX.md`
-   `.cursor/rules/after_action_reports/operation_guardian_AAR.md`
-   `.cursor/rules/clinical_trials_agents/CLINICAL_TRIALS_MASTER_DOCUMENT.md`
-   `.cursor/rules/research/00_MASTER_INDEX.mdc`

## Initial Files to Consider for Renaming (based on a quick scan)

-   `.cursor/rules/CrisPRO_Command_Center/00_MASTER_INDEX.mdc` (previously `README.md`, renamed for consistency)
-   `.cursor/rules/MOAT/README.md` (already addressed, distinction documented)
-   `.cursor/rules/MODULARIZATION_COMPLETE_REPORT.mdc` (previously `MODULARIZATION_COMPLETE_SUMMARY.md`; requires manual content consolidation from `MODULARIZATION_PROJECT_SUMMARY.md` due to `.cursorignore` restrictions)
-   `.cursor/rules/MODULARIZATION_PROJECT_SUMMARY.md` (previously `MODULARIZATION_COMPLETE.md`; requires manual content consolidation into `MODULARIZATION_COMPLETE_REPORT.mdc` due to `.cursorignore` restrictions)
-   `.cursor/rules/P0_INTEGRATION_SUMMARY.md` (previously `P0_INTEGRATION_COMPLETE.md`, renamed for clarity)

## Status

-   **Draft**: Initial version for review and iteration.
-   **Last Updated**: January 2025
