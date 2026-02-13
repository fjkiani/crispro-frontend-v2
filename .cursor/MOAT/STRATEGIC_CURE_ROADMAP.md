# 🎯 Strategic Roadmap: High Value Targets to Cure Cancer

**Mission:** Move beyond "predicting outcomes" to **changing outcomes** by intercepting resistance before it happens.
**Status:** Strategic Plan (Draft)
**Source:** Synthesis of SAE Intelligence Audit & Resistance Prophet Metrics.

---

## 💎 The "Big 3" High Value Targets

Research confirms three specific areas where MOAT can deliver disproportionate value. These are not just "features"; they are potential cure-enablers.

### 1. The "Intrinsic Lock": MFAP4 & EMT State
*   **The Problem:** Platinum resistance is often intrinsic (pre-existing) in High-Grade Serous Ovarian Cancer (HGSOC), leading to futile 1st-line toxicity.
*   **The Solution:** Deploy **MFAP4** as a definitive pre-treatment biomarker.
*   **Value:** Avoid futile platinum rounds. Switch immediately to non-platinum alternatives (e.g., trials, PARP-frontline) for resistant patients.
*   **Metric:** Validated AUROC **0.763** (beats standard-of-care prediction).

### 2. The "Adaptive Radar": Serial CA-125 Kinetics
*   **The Problem:** Resistance emerges adaptively during treatment. Static biomarkers (DNA) miss this dynamic evolution.
*   **The Solution:** Implement **Sequence-agnostic Serial Monitoring**. Use kinetic modeling of CA-125 decay/rise to flag resistance *weeks* before radiological progression.
*   **Value:** "Switch Time" reduction. Catching resistance 4 weeks earlier can mean the difference between widespread metastasis and treatable disease.

### 3. The "Myeloma Shield": Complete Resistance Ecosystem
*   **The Problem:** Multiple Myeloma has clear resistance drivers (DIS3, BRAF, MEK) but they are often ignored in community settings.
*   **The Solution:** Close the gap on **PSMB5 (Proteasome Inhibitors)** and **CRBN (IMiDs)**.
*   **Value:** 100% of MM patients receive these drugs. Predicting resistance to them is the "Holy Grail" of MM precision medicine.
*   **Status:** DIS3/TP53 validated (**RR=2.08**). PSMB5/CRBN validation is the final mile.

---

## 🗺️ Operation "Cure" Deliverables

| Target | Deliverable | Type | Dependencies |
| :--- | :--- | :--- | :--- |
| **1. Intrinsic Lock** | **MFAP4 Module Integration** | Code | `SAE_INTELLIGENCE` |
| | *Action*: Wire MFAP4 score into `BiomarkerAgent` output. | | |
| **2. Adaptive Radar** | **Serial Kinetics Engine** | Code | `PatientState` History |
| | *Action*: Build `def analyze_kinetics(history)` in `MonitoringAgent`. | | |
| **3. Myeloma Shield** | **PSMB5/CRBN Validation** | Research | External Data / Lit |
| | *Action*: Targeted extraction of PSMB5/CRBN variants outcomes. | | |

---

## 🚫 Strategic "Stop" List (Low Value / Distractions)

To focus on the cure, we must cut the noise:

1.  **STOP** Diamond/True SAE Features: Audit proved them to be noise (AUROC 0.55). Do not waste compute on them.
2.  **STOP** Broad-Spectrum "AI" Markers: Without biological grounding (like MFAP4), they are hallucinations.
3.  **STOP** General Chatbots: Focus on **Agentic Actions** (matching, predicting, alerting), not conversation.

---

## 🚀 Execution Plan

1.  **Immediate (Week 1):** Integrate MFAP4 into the production pipeline (`BiomarkerAgent`).
2.  **Near Term (Week 2):** Prototype Serial Kinetics for Ovarian Cancer.
3.  **Strategic (Month 1):** Launch the "Myeloma Shield" data campaign for PSMB5/CRBN.
