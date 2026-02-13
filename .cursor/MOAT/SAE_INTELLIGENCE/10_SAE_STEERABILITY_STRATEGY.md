# SAE Steerability Strategy

## 1. The Concept: "Steerability"
Steerability is the capability to **intervene** on specific biological features in silico to simulate outcomes. Unlike "Black Box" predictions, a "Glass Box" SAE model allows us to clamp features (e.g., "Zero out the MAPK path") and observe the predicted change in resistance probability.

## 2. The 3-Layer Moat
1.  **De-Noising (Layer 1)**: Converting dense, messy embeddings into sparse, interpretable features. (Status: Built via SAE).
2.  **Audit Trail (Layer 2)**: Providing provenance for every prediction. "Why did we predict resistance? Because Features X, Y, Z are active." (Status: Built via Proxy).
3.  **Steerability (Layer 3 - The Crown Jewel)**: Simulating counterfactuals. "What if we inhibit PARP?" (Status: Research/Future).

## 3. Implementation Tiers

### V0: Proxy Steerability (Current)
*   **Method**: Intervene on the 7D Mechanism Vector.
*   **Example**: Manually set `DDR` score to 0.0 (simulating perfect repair restoration) -> Watch Resistance Probability spike.
*   **Use Case**: Clinician "What-If" scenarios in the dashboard.

### V1: Proxy Steerability (THE STANDARD)
*   **Method**: Intervene on the 7D Mechanism Vector or Specific Genes.
*   **Example**: "Mock `NF1` Mutation" -> Observe MAPK Pathway Spike -> Resistance Risk Increase.
*   **Evidence**: Validated by RR=2.0 (NF1/DIS3).
*   **Status**: ✅ Production Ready.

### V2: True SAE Steerability (DEPRECATED)
*   **Status**: ❌ Paused (Validation Failed).
*   **Reason**: The "Diamond Bin" features (27607) proved to be noise.
*   **Future**: Requires N > 500 to revisit high-dimensional steering.

## 4. Strategic Value
*   **Regulatory**: Approximates "White Box" AI, essential for FDA SAMD (Software as a Medical Device) approval.
*   **Clinical**: Empowers oncologists to test hypotheses before prescribing toxic therapies.
*   **Commercial**: Differentiates from generic "Black Box" predictors (e.g., Tempus/Foundation) which offer raw probabilities without interactive simulation.
