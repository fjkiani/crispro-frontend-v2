# CAPABILITIES REGISTRY (ZETA DOCTRINE)

**Status Definitions:**
*   🟢 **OPERATIONAL**: Production-ready, tested, currently in use.
*   🟡 **PARTIAL**: Implemented but requires warnings/caveats (e.g., heuristic only).
*   🔴 **PLACEHOLDER**: Stubbed code, returns static/fake data. MUST BE FLAGGED.
*   💀 **REMOVED**: Purged from codebase.

| Subsystem | Component | Status | Verification/Receipts | Notes |
| :--- | :--- | :--- | :--- | :--- |
| **Resistance Engine** | `ClinicalHeuristicEngine` (Proxy) | 🟢 **OPERATIONAL** | `api/services/resistance_prophet_service.py` | Renamed from "Proxy SAE". production default. |
| **Resistance Engine** | `Evo2SAEExtractor_RUO` (True) | 🔴 **PLACEHOLDER** | *Pending Location* | **BLOCKED.** Research Use Only. Do not call in production. |
| **Ayesha Fit** | `ResistanceLab` (Glass Box) | 🟡 **PARTIAL** | `ResistanceLab.jsx` | Connected but lacks specific visualizations (`EMTRiskGauge`). |
| **Biomarker** | `MFAP4` (Transcriptomic) | 🟢 **OPERATIONAL** | `api/services/ayesha_fit/builder.py` | Provenance: GSE63885 (AUROC 0.763). |
| **Biomarker** | `Diamond SAE` (Features) | 💀 **REMOVED** | N/A | Logic removed from active paths. |
| **Architecture** | `compute_serial_sae` (Kinetics) | ❓ **UNVERIFIED** | TBD | Suspected vaporware. Needs audit. |

**Enforcement Rule:**
Any module marked 🔴 **PLACEHOLDER** must return explicit structure:
```json
{
  "status": "not_implemented",
  "reason": "Research Use Only - Pending Validation",
  "provenance": null
}
```
Silent failures or static "mock" data are strictly prohibited.
