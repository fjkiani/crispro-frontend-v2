# 🔮 Resistance Prophet: Master Documentation

**Status:** ✅ Production Ready (Validated)
**Version:** 2.0 (Consolidated)
**Last Updated:** February 2025

**Supersedes:** `.cursor/ayesha/RESISTANCE_PREDICTION_VALIDATED.md`, `RESISTANCE_PROPHET_AUDIT.md`, `RESISTANCE_PROPHET_INTEGRATION_AUDIT.md`, `RESISTANCE_PROPHET_INTEGRATION_COMPLETE.md`, `RESISTANCE_PROPHET_REFACTOR_PLAN.md`, `RESISTANCE_PROPHET_REFACTOR_PLAN_FINAL.md` (archived)

---

## 1. Executive Summary
Resistance Prophet is the core predictive engine of MOAT, responsible for identifying potential drug resistance mechanisms before treatment begins. It integrates sequence analysis, pathway topology, and clinical evidence.

**Key Capabilities:**
*   **Gene-Level Prediction:** Validated predictions for specific variants (e.g., DIS3 in MM).
*   **Pathway-Level Prediction:** Validated pathway aggregation (e.g., MAPK pathway in Ovarian).
*   **Evidence Tiering:** Automatic classification of evidence strength (Supported, Consider, Insufficient).

---

## 2. Validated Metrics (The Ledger)

The following claims are strictly validated with reproducibility receipts.

### Multiple Myeloma (MM)
| Biomarker | Cohort | Relative Risk (RR) | p-value | Status |
|-----------|--------|-------------------|---------|--------|
| **DIS3** | MMRF CoMMpass (n=995) | **2.08** | **0.0145** | ✅ Significant |
| **TP53** | MMRF CoMMpass (n=995) | **1.90** | 0.11 | ⚠️ Clinical Trend |

### Ovarian Cancer (OV)
| Pathway | Cohort | Relative Risk (RR) | p-value | Finding |
|---------|--------|-------------------|---------|---------|
| **NF1** | TCGA-OV (n=469) | **2.10** | <0.05 | Predicts Platinum Resistance |
| **MAPK** | TCGA-OV (n=469) | **1.97** | <0.05 | Predicts Platinum Resistance |
| **PI3K** | TCGA-OV (n=469) | **1.39** | 0.02 | Predicts Platinum Resistance |

**DDR (DNA repair) for outcome prediction:** Not validated as a standalone prognostic marker in our cohorts. Do not claim DDR for outcome prediction. MAPK and NF1 are the validated OV resistance signals.

### Production Reliability
*   **Resistance E2E Fixtures:** 100% Pass (5/5). Validated handling of L0/L1/L2 input completeness and confidence caps.
*   **Gating Logic:** Validated deterministic behavior (penalty vs rescue).

---

## 3. Implementation Status

### Core Components
*   **Service:** `ResistanceProphetService` (wrapper) + modular package `api/services/resistance_prophet/`:
    *   `engine.py` – main prediction logic; `aggregation.py` – risk aggregation
    *   `signals/ovarian.py` – OV signals (MAPK, NF1, PI3K); `signals/mm.py` – MM signals (DIS3, TP53)
    *   `actions.py`, `shim.py`, `constants.py`, `schemas.py`, `spine.py`, `timing_engine.py`, `prognosis.py`
*   **Playbook:** `ResistancePlaybookService`. Implemented.
*   **Endpoints:** `/api/resistance/predict`, `/api/care/resistance_playbook`.

### Integration
*   **Orchestrator:** Fully wired into `AGENT_03`.
*   **Frontend:** `ResistancePlaybook.jsx` (Placeholder/Partial) - *See Remaining Work*.

### Mars Protocol Refactor (Complete)
Modular package layout implemented. Logic moved from monolithic `resistance_prophet_service.py` into `api/services/resistance_prophet/` with engine, signals (ovarian, mm), aggregation, actions, shim, spine, timing_engine, and prognosis.

### Baseline Artifact Fix (Ayesha HIGH-Risk Overcall)
**Root cause:** When baseline SAE features are missing, the system used population average instead, which inflated resistance probability and produced spurious HIGH risk.
**Fix:**
*   Apply a **30% probability penalty** when baseline is missing (population average used).
*   Require **≥3 signals** for HIGH risk when baseline is missing (prevents single-signal HIGH overcall).

---

## 4. Mission: Multiple Myeloma (MM)

**Goal:** Complete resistance prediction ecosystem for Multiple Myeloma.

### Status
*   **Achieved:** Validated signal for DIS3 and TP53.
*   **Gap:** PSMB5 (Proteasome Inhibitor) and CRBN (IMiD) specific resistance mutations are not yet fully implemented/validated due to data sparsity in standard cohorts.

**Consolidation Strategy:**
*   This document replaces `MISSION_MM_RESISTANCE...` and `MM_RESISTANCE_IMPLEMENTATION...`.
*   Future work should focus on extracting PSMB5/CRBN data from focused lit search or specialized cohorts.

---

## 5. Remaining Work (Backlog)

### Critical
- [ ] **Frontend:** Implement full `ResistancePlaybook.jsx` to visualize resistance pathways (currently a placeholder).
- [ ] **Validation:** Expand MM validation to include PSMB5 and CRBN if data permits.
- [ ] **Evo2 Integration:** Ensure Evo2 service is called as a secondary signal for variants of uncertain significance (VUS).

### Enhancements
- [ ] **Mechanism Explanation:** Improve "Why" text generation for predicted resistance.
