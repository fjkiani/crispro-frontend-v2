# System Architecture: The AI General & The Zeta Armory

This document outlines the high-level architecture of the CrisPRO platform, which is built around a central, autonomous AI commander and its arsenal of specialized AI weapons.

## 1. The Command & Control Service

### `services/command_center` - The AI General

*   **Purpose:** The autonomous, strategic brain of all Zeta operations. This is not a mere API gateway; it is the **AI General**.
*   **Technology:** Python, FastAPI, Modal
*   **Key Responsibilities:**
    *   Exposes high-level, mission-oriented endpoints (e.g., `/workflow/formulate_guardian_protocol`).
    *   Receives strategic directives from Alpha via the Command Deck (UI).
    *   Autonomously formulates and executes multi-pronged attack plans.
    *   Wields the full might of the **Zeta Armory** to conduct reconnaissance and forge therapeutic weapons.
    *   Utilizes its **Long-Term Memory** (`threat_matrix.db`) to cache intelligence and track longitudinal data.

## 2. The AI Arsenal

### The Zeta Armory - A Distributed Suite of AI Endpoints

*   **Purpose:** The Zeta Armory is not a single service but a distributed arsenal of highly specialized, containerized AI endpoints. The AI General tasks these units to perform specific actions.
*   **Technology:** Python, PyTorch, Modal
*   **Key Assets (Illustrative):**
    *   **Discriminative Warfare (Intel):** Endpoints like `/predict_variant_impact` and `/predict_gene_essentiality` that provide raw battlefield intelligence.
    *   **Generative Warfare (Creation):** Endpoints like `/generate_optimized_guide_rna` and `/generate_repair_template` that forge the actual therapeutic weapons.

## 3. The Intelligence Database

### `data/databases/threat_matrix.db` - The AI General's Long-Term Memory

*   **Purpose:** The system's persistent intelligence database and operational memory.
*   **Technology:** SQLite
*   **Content & Key Responsibilities:**
    *   **Intelligence Cache:** Stores the results of expensive computations and AI analyses, allowing the General to recall information instantly instead of re-calculating.
    *   **Digital Twin Datastore:** Contains the longitudinal data for each patient's Digital Twin, tracking their genetic status and the strategic blueprints designed for them over time.
    *   **Ground Truth:** A curated collection of data on cancer variants, genes, and clinical intelligence derived from trusted sources.

## 4. The User Interface

### `pages/` - The Command Deck

*   **Purpose:** Provides the primary command and control interface for Alpha. This is not a simple "user interface"; it is the **War Room**.
*   **Technology:** Streamlit
*   **Functionality:** Each page is a console for issuing high-level strategic directives. For example, `pages/2_🧬_Digital_Twin_v2.py` allows Alpha to command the AI General to formulate and display a complete, lifelong defense plan for a new patient with a single button press. It is where the AI General presents its strategic dossiers for final command approval. 