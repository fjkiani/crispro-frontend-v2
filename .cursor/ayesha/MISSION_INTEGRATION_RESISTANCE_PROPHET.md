# Mission: Resistance Prophet Integration (The Bridge)

**Status:** 🚦 **PROPOSED**
**Priority:** High
**owner:** Resistance Lab Team

## 🎯 Objective
Eliminate the "heuristic orphan" status of the Resistance Prophet Service by wiring it to consume the rigorously verified Timing & Chemosensitivity Engine.

## 🧱 Current State (The Gap)
*   **Timing Engine**: ✅ **Online**. Produces rigorous PFI, PTPI, TFI, and KELIM metrics.
*   **Ayesha Dashboard**: ✅ **Wired**. Displays these metrics using `useTimingChemoFeatures`.
*   **Resistance Prophet**: ❌ **Orphaned**. Uses simplified heuristics in `resistance_prophet_service.py` (`mode="heuristic"`) and ignores the Timing Engine.

## 🛠️ Scope of Work

### 1. API Contract Definition
    - [x] Define `PatientRegimenFeatureTable` schema (Unified Input).
    - [x] Update Prophet Service to accept this table as an optional input `provenance`.

### 2. Prophet Service Refactor
    - [x] **Deprecate Heuristics**: Remove duplicate/simplified PFI logic in Prophet.
    - [x] **Consume Real Features**:
        - If `TimingFeatures` provided -> Use them for Logic Rules.
        - Rule: `PFI < 6m` (from Engine) -> Resistance Risk.
        - Rule: `KELIM < 1.0` (from Engine) -> Chemo-Refractory Signal.

### 3. Frontend "Resistance Lab" Integration
    - [x] Update `ResistanceLab.jsx` to fetch Timing Features *first*.
    - [x] Transform and inject `timing_context` into `/api/ayesha/resistance/simulate` payload.

### 4. Verification
    - [x] **Backend Bridge Test**: Verify `predict_resistance` accepts `PatientRegimenFeatureTable` and triggers Engine Logic.
    - [x] **Frontend Simulation**: Verify Simulation uses Real Engine data (e.g. Verified KELIM).
    - [x] **E2E Test**: Verify that Dashboard (Tab 2) and Prophet (Resistance Lab) produce **identical** conclusions for the same patient.
    - [x] **Audit**: Prove no "data drift" between the two components.

## ⚠️ Constraints & Guardrails
*   **RUO Only**: The Timing Engine is RUO. Prophet predictions must clearly state this limitation.
*   **Fallback**: If Timing Engine fails or data is missing, Prophet must fail gracefully (or explicitly fallback to heuristics with a warning).
