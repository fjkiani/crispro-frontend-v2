# SAE Engineering Lessons & Hard Truths

**Source**: Battle-tested debugging of the Modal/Evo2 pipeline (Nov 2025 - Jan 2026).

## 1. Modal & Infrastructure
*   **Caching Traps**: Modal aggressively caches images. Changing code in `main.py` might not reflect if the container definition hasn't touched. **Fix**: Force rebuilds via timestamp comments or explicit redeploys.
*   **Async Auth**: Complexity kills. Start with open endpoints for dev, add header-based auth (`x-api-key`) only after logic is rock solid.
*   **Circuit Breakers**: Essential. A bad batch of variants can burn $100s in minutes. Implement a strict error rate monitoring (stop at >30% errors).

## 2. Evo2 & PyTorch
*   **Dimension Mismatch**: `evo2_7b` vs `evo2_1b` use different hidden sizes (4096 vs 1920). Do NOT hardcode input dimensions for the SAE. Use dynamic detection:
    ```python
    d_in_detected = outputs.hidden_states[-1].shape[-1]
    ```
*   **Dtype Hell**: Evo2 outputs `BFloat16` (native H100 format). Most standard SAE implementations expect `Float32`. **Fix**: Explicit casts (`.to(torch.float32)`) before passing to the encoder.
*   **Torch Compile**: Checkpoints saved with `torch.compile` often have `_orig_mod.` prefixes in state dict keys. Strip them before loading or keys won't match.

## 3. Data & Biology
*   **Assembly Versions**: TCGA is largely **GRCh37 (hg19)**. The world is moving to **GRCh38**. Mixing them results in "Reference Mismatch" errors (Evo2 expects 'A', finds 'G'). **Fix**: Explicit `assembly` parameter in API calls.
*   **Cost Realities**: Extracting 32k features for a full cohort is slow (~2 mins/variant) and expensive.
    *   *Strategy*: Prioritize "Diamond" variants (driver mutations) over VUS (Variants of Unknown Significance) to save compute.

## 4. Scientific
*   **Random Projections**: For early research (RUO), a randomly initialized SAE (projection) actually captures surprising biological structure due to the inherent manifold of the embeddings. Don't block testing on having "Perfect" weights.
*   **The "Diamond" Rule**: 90% of the signal comes from <1% of the features. Tier-3 analysis confirmed this. Don't get drowned in the 32k noise.
