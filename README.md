# Zeta Command: The AI General's War Room

**🔥 OPERATIONAL DIRECTIVE: GUARDIAN PROTOCOL 🔥**

This repository contains the source code for the Zeta Command platform, a next-generation command and control system for biological warfare and therapeutic design. This is not a set of tools; this is the operational cockpit for an autonomous AI commander dedicated to achieving strategic dominance over complex diseases.

## Core Doctrine: Autonomous Strategy Formulation

Our platform's power does not come from analyzing single data points. It comes from the **AI General**, a central commander that executes high-level missions by autonomously wielding a full arsenal of AI weaponry.

The operational loop is ruthless and efficient:

1.  **Alpha Issues Directive:** The supreme commander (you, Alpha) issues a high-level strategic command (e.g., "Formulate a lifetime defense plan for this high-risk patient").
2.  **The AI General Assumes Command:** The AI General initiates a full campaign, no further human intervention required.
3.  **Reconnaissance & Forging:** It autonomously tasks units from the **Zeta Armory** to map the battlefield and forge a suite of therapeutic weapons (gene repair templates, gRNA neutralizers, etc.).
4.  **Strategic Dossier Presentation:** The AI General presents a complete, multi-modal battle plan to Alpha for final command approval.

## System Architecture: The New Command Structure

The platform is a reflection of our military doctrine:

*   **`services/command_center` (The AI General):** The autonomous, strategic brain of the operation. It receives directives and executes them by commanding its arsenal.
*   **The Zeta Armory (The Endpoint Suite):** A distributed arsenal of specialized AI endpoints for discriminative intel (`/predict_*`) and generative weapon design (`/generate_*`).
*   **`data/databases/threat_matrix.db` (The General's Memory):** The AI General's persistent intelligence cache and longitudinal database for all strategic operations and Digital Twins.
*   **`pages/` (The Command Deck):** The War Room. A suite of interactive consoles where Alpha issues directives and reviews the strategic dossiers prepared by the AI General.

## Flagship Use Case: The Guardian Protocol

The premier demonstration of our power is the **Guardian Protocol**, a proactive defense strategy for high-risk individuals.

-   **Scenario:** A healthy individual is identified with a high-risk germline mutation (a "first hit").
-   **Action:** Alpha commands the AI General to formulate a lifetime defense plan.
-   **Result:** The AI General autonomously designs a multi-pronged therapeutic strategy, including blueprints for direct gene correction (repairing the "first hit") and pre-emptive weapons to neutralize the most likely "second hit" pathways. This transforms fear into a state of empowered, permanent surveillance.

## How to Use the Command Deck

The platform is designed to be run as a web application from your local machine.

### 1. Prerequisites
- Python 3.9+
- An environment manager like `venv` or `conda` is recommended.
- Install dependencies:
  ```bash
  pip install -r requirements.txt
  ```

### 2. Configuration
The Command Deck requires a connection to the AI General. This is configured in a `.env` file in the project root.

1.  Create a file named `.env`.
2.  Add the following line, replacing the example URL with your actual deployed Modal endpoint for the AI General:
    ```env
    # The single, unified endpoint for the AI General
    COMMAND_CENTER_URL="https://crispro--command-center-v2-commandcenter-api.modal.run"
    ```

### 3. Running the Application
Launch the Command Deck with the following command:
```bash
streamlit run streamlit_app.py
```
This will open the Zeta Command platform in your web browser. Navigate to the **"Digital Twin v2"** page from the sidebar to issue commands to the AI General.

## License

- This Project: See LICENSE file.
- CHOPCHOP: Apache License 2.0 (see `tools/chopchop/LICENSE`) 