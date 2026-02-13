# ☁️ SaaS Platform: Master Documentation

**Status:** ✅ Live / 🚧 Security Hardening
**Version:** 2.0 (Consolidated)
**Last Updated:** January 2025

---

## 1. Architecture Overview

### Vision
A unified "Manager" agent (Archon) orchestrating specialized sub-agents for biomarkers, resistance, efficacy, and nutrition.

### Orchestration Pipeline
The system follows a strict **7-Phase Pipeline**:
1.  **INITIALIZED:** Request received.
2.  **EXTRACTING:** Parsing files (VCF, PDF).
3.  **ANALYZING:** Parallel agent execution (Biomarker, Resistance).
4.  **RANKING:** Drug efficacy scoring (S/P/E).
5.  **MATCHING:** Clinical trial matching.
6.  **PLANNING:** Nutrition and Care Plan synthesis.
7.  **MONITORING:** Continuous surveillance setup.

### Integration
*   **PubMearch:** Keyword hotspot analysis (Audited & Ready).
*   **PubMed Parser:** Deep XML parsing for quantitative data (Audited & Ready).

---

## 2. Security Capabilities & Audit

**Status:** Gaps Identified (Roadmap in place).

### Validated Controls
*   **Orchestrator Provenance:** Every decision is tagged with the agent ID and run signature.

### Gap Analysis (From Audit)
*   **Auth:** No robust RBAC implemented yet.
*   **Data:** Patient data encryption at rest needs verification.
*   **Audit Logs:** Centralized audit logging is partial (Orchestrator logs exist, but structured security logging is needed).

---

## 3. Tech Stack & Pattern

### The Golden Rule
**✅ DO:** Import services directly (e.g., `from ..service import X`).
**❌ DON'T:** Make HTTP calls from within the orchestrator to its own API (avoids race conditions/overhead).

### Codebase Structure
*   `api/services/orchestrator/`: Core logic.
*   `api/services/resistance/`: Resistance Prophet.
*   `api/services/drug_efficacy/`: Therapy Fit.
*   `oncology-frontend/`: React-based dashboard.

---

## 4. Remaining Work (Backlog)

### Platform
- [ ] **Orchestrator Wrapper:** Create `api/services/research_intelligence/orchestrator.py` to unify the research tools.
- [ ] **Frontend Mapper:** Create `orchestratorMapper.js` to translate new API responses for legacy components.

### Security (P1)
- [ ] **RBAC:** Implement basic Role-Based Access Control.
- [ ] **Encryption:** Verify encryption at rest for `PatientState`.
- [ ] **Logging:** Implement structured security audit logs.
