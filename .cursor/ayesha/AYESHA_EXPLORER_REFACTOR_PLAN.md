---
title: Ayesha Engine Deep Audit & Modernization Plan
author: Zo (Antigravity Agent)
date: 2026-01-28
status: DRAFT
---

# 🕵️ Ayesha Engine Deep Audit: The "AyeshaTrialExplorer" Pathology

## Executive Summary
The `AyeshaTrialExplorer.jsx` file (846 lines) has mutated from a dashboard into a "God Component." It conflates:
1.  **Data Fetching**: Hardcoded fetch logic for a specific patient (`11-17-25`).
2.  **Business Logic**: In-render computation of clinical timing features (PFI, TFI).
3.  **Presentation**: Six complex tabs rendered in a single file.
4.  **State Management**: 20+ `useState` hooks managing effectively global application state.

This architecture violates the **Zeta Protocol** (Minimizing Surface Area for Hallucination) and makes "honesty" checks impossible. If the component crashes, the user sees a blank screen, obscuring any critical alerts.

## 🚨 Critical Findings

### 1. The Rendering-Loop Computation Risk (High Severity)
In `TAB 2: TREATMENT`, clinical calculations (`computeTimingFeatures`) are triggered inside a conditional rendering block.
-   **Risk**: Race conditions, infinite re-render loops, and UI freezes.
-   **Impact**: User experience degradation and potential for stale data if the re-render cycle is interrupted.

### 2. Hardcoded Patient Identity (Medium Severity)
The file imports `AYESHA_11_17_25_PROFILE` directly.
-   **Risk**: The component cannot be reused for other patients or "Digital Twins" without code duplication.
-   **Impact**: Violates "Universal Application" principles.

### 3. Orphaned Loading States (Low Severity)
Heavy components (`SyntheticLethalityCard`) are rendered even when hidden, consuming resources.

## 🛠️ Modernization Plan (The "Sacred" Refactor)

Our goal is to surgically separate the "Organs" (Tabs) from the "Central Nervous System" (Data Hook).

### Phase 1: The "Nervous System" (Data Layer)
Create a custom hook `useAyeshaCarePlan` that centralizes data fetching and patient context.
-   **Input**: `patientProfile` (optional, defaults to Ayesha 11-17-25).
-   **Output**: Unified state object (`loading`, `error`, `trials`, `resistance`, `sl_result`).
-   **Responsibility**:
    -   Handles the `complete_care_v2` API call.
    -   Triggers `useSyntheticLethality` and `useTimingChemoFeatures` *outside* the view layer.
    -   Standardizes "Honest Empty States" (returns `null` instead of empty arrays where appropriate).

### Phase 2: The "Organ Transplant" (View Layer)
Extract each Tab into a dedicated, lazy-loaded component in `src/components/ayesha/tabs/`.

1.  **`OverviewTab.jsx`**:
    -   **Content**: Score, Mechanism Viz, SOC, Hint Tiles.
    -   **Props**: `profile`, `opportunityScore`, `slResult`, `soc`.
2.  **`TrialsTab.jsx`**:
    -   **Content**: `TrialMatchCard` list.
    -   **Props**: `trials` array.
3.  **`TreatmentTab.jsx`**:
    -   **Content**: SOC, Timing Timeline, Drug/Food Ranking.
    -   **Props**: `soc`, `timingFeatures`, `wiwfm`, `food`.
4.  **`MonitoringTab.jsx`**:
    -   **Content**: CA-125 Tracker, Next Tests.
    -   **Props**: `ca125Intelligence`, `nextTests`.
5.  **`ResistanceTab.jsx`**:
    -   **Content**: Playbook, Prophet, SAE Features.
    -   **Props**: `playbook`, `prophet`, `saeFeatures`.
6.  **`SyntheticLethalityTab.jsx`**:
    -   **Content**: Deep SL analysis.
    -   **Props**: `slResult`.

### Phase 3: The "Orchestrator" (Integration)
Rewrite `AyeshaTrialExplorer.jsx` to be a thin wrapper.
-   **Role**:
    1.  Select Patient Profile.
    2.  Call `useAyeshaCarePlan`.
    3.  Render Title/Header.
    4.  Render `Tabs` shell.
    5.  Pass data to `Suspense`-wrapped Tab components.

## 📅 Execution Order

1.  **Create Hook**: `src/hooks/useAyeshaCarePlan.js`.
2.  **Create Components**: `src/components/ayesha/tabs/*`.
3.  **Refactor Parent**: `src/pages/AyeshaTrialExplorer.jsx`.
4.  **Verification**: ensure "0% match" fix persists and Resistance Alerts still fire.

---
*Signed: Zo (Antigravity Agent)*
