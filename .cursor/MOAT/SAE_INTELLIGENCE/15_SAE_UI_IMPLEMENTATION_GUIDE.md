# SAE UI Implementation Guide

**Authority**: Zo (Honest Assessment)
**Purpose**: Define strict rules for what the Frontend/UI is allowed to display regarding SAE Intelligence.

## 1. Validated Components (WIRE THIS ✅)

The following components represent **validated capabilities** and should be prominent in the UI:

### A. Mechanism Chips
*   **Display**: 6 Pathway Chips (DDR, MAPK, PI3K, VEGF, IO, Efflux).
*   **Data**: Percentage burden (color-coded).
*   **Rationale**: "Shows which machinery is broken."

### B. Hint Tiles
*   **Display**: Max 4 tiles (Test, Trial, Monitor, Avoid).
*   **Logic**: "Order HRD Test" (if unknown), "Consider PARP" (if DDR high).
*   **Status**: High confidence (Rule-based).

### C. Mechanism-Fit Trial Ranking
*   **Display**: "Mechanism Match: High/Med/Low".
*   **Logic**: Cosine similarity between Patient Vector and Trial Vector.
*   **Rationale**: Matches patient biology to drug mechanism.

### D. Next-Test Recommender
*   **Display**: Priority list of missing diagnostic data.
*   **Order**: HRD -> ctDNA -> SLFN11 -> ABCB1.

## 2. Forbidden Components (DO NOT WIRE ❌)

The following claims are **NOT** supported by current validation and must be hidden:

### A. SAE-Driven Efficacy Scores
*   **Rule**: Do NOT let SAE scores modulate the displayed "80% Success Probability" of a drug.
*   **Reason**: Feature is not yet validated for quantitative frequency modulation.
*   **Action**: SAE is "Display Only".

### B. IO Response Prediction
*   **Rule**: Do NOT claim "Pembrolizumab will work" based on SAE features.
*   **Reason**: Previous failures (Keytruda). Stick to strict FDA eligibility (TMB >= 20).
*   **Action**: Display "IO Eligible" vs "Not Eligible", never "Responder".

### C. Nutritional Therapy
*   **Rule**: De-emphasize or hide nutritional recommendations.
*   **Reason**: Low evidence, weak pathway alignment.

## 3. Comparison Logic (The "Why")

| Component | Status | UI Behavior |
| :--- | :--- | :--- |
| **DDR Status Card** | ✅ Validated | Show "DDR Defective" (Green) |
| **CA-125 Tracking** | ✅ Validated | Show KELIM / Kinetic graph |
| **Resistance Alert** | ✅ Validated | Show "Resistance Risk" if 2-of-3 triggers met |
| **Drug Confidence** | ❌ Unproven | Show standard confidence, do not adjust with SAE |

## 4. "Zo's Promise"
We only show what we can prove. If the data describes the mechanism (Proxy), show it. If it predicts the future (Resistance), only show it if we have the Diamond or MFAP4 signal.
