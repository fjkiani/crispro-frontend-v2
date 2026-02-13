# SAE Intelligence System: The "Moat" Architecture

**Status**: PRODUCTION (Refactored Jan 2026)
**Core Logic**: Validated Biology (MFAP4) + Heuristic Proxy (Mutations)
**Deprecated**: "True SAE" Features (Failed Validation)

---

## 1. System Overview

The **SAE Intelligence System** is the decision engine for the Resistance Prophet. It aggregates signals from 4 layers of resistance biology to predict patient outcomes and recommend next-line therapies.

### The Core Thesis
We do not rely on "black box" neural features. We rely on **Validated Biological Signals**.

| Layer | Component | Status | Signal Strength |
| :--- | :--- | :--- | :--- |
| **L1: Intrinsic** | **MFAP4 Biomarker** | 🏆 **Gold** | AUROC 0.763 (External) |
| **L2: Genetic** | **Proxy SAE** | ✅ **Silver** | AUROC 0.628 (Baseline) |
| **L3: Adaptive** | **CA-125 Kinetics** | ⚠️ **Alert** | >25% Risk |
| **L4: Clearance** | **ABCB1 Efflux** | 🧪 *Exp* | Future |

---

## 2. The Logic Engine

### Step 1: Layer 1 Check (Cell State)
*   **Input**: RNA Expression (MFAP4).
*   **Logic**: If MFAP4 > 1.5 SD -> **EMT Phenotype Detected**.
*   **Result**: High Probability of Platinum Resistance.

### Step 2: Layer 2 Check (Mutations)
*   **Input**: Mutation Profile (VCF).
*   **Logic (Proxy)**:
    *   `Map Mutations -> Pathways` (e.g. NF1 -> MAPK).
    *   `Sum Pathway Burden`.
    *   `Compare to Drug Target`.
*   **Result**: "Mechanism Match" or "Escape Detected".

### Step 3: Layer 3 Check (Dynamics)
*   **Input**: CA-125 Series.
*   **Logic**: Calculate K-elimination and Delta from Nadir.
*   **Result**: "Relapse Imminent" vs "Stable".

---

## 3. Deprecated Components (Do Not Use)

*   **Diamond SAE Features**: The 9 features (27607, etc.) previously claimed to predict resistance have been **debunked** (Audit Jan 2026). They are random noise.
*   **True SAE Vectors**: Direct use of SAE latent vectors is suspended until N > 500.

---

## 4. Validated Drug Combinations (Smart Logic)
The system identifies patients who benefit from multi-angle attacks based on **Proxy Pathway Vectors**:

1.  **PARP + Bevacizumab (Avastin)**
    *   **Logic**: `High DDR Mutation Load AND High VEGF Expression`
    *   **Rationale**: Attack DNA repair while starving tumor (PAOLA-1 Mechanism).
    
2.  **PARP + ATR Inhibitor**
    *   **Logic**: `High DDR` + `Recurrent Disease`
    *   **Rationale**: Synthetic lethality overlap.

3.  **IO + PARP (The "Topacio" Signal)**
    *   **Logic**: `High TMB` OR `High DDR`
    *   **Rationale**: DNA damage primes immune recognition. 
