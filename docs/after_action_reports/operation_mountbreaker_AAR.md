# After Action Report: Operation Mountbreaker

**Mission Status: VICTORY**

**Objective:** Deploy the `evo_service`, a critical microservice powered by the `arcinstitute/evo2` 40B parameter model, to the Modal cloud environment.

**Execution Summary:**
The campaign to deploy `evo_service` was one of the most brutal and protracted debugging engagements in this project's history. We were met with a cascade of failures, each more illogical than the last, stemming from a hostile and poorly documented library. The core strategy evolved from a naive "build from source" attempt to a sophisticated, multi-pronged assault that leveraged pre-built NVIDIA containers, deep source code reconnaissance, and a final, brute-force solution to hardware limitations.

**Key Learnings & Intelligence Gained:**
The intelligence gathered during this operation has been codified into a new doctrine: [The NVIDIA Container Deployment Playbook](mdc:../../.cursor/rules/nvidia-container-deployment-guide.mdc). The core tenets are summarized here:

1.  **NEVER Build From Source:** Complex, GPU-dependent libraries like `evo2` are not meant to be built from source in a CI/CD environment. It is a fool's errand. We must always default to using official, pre-built containers from registries like NVIDIA NGC.

2.  **Submodule Deception:** Git submodules are a critical failure point for Modal's `add_local_dir`. The build process does not reliably copy them. The only effective solution was an explicit, dual-mount strategy, forcing both the primary package and the submodule directory (`vortex`) into the container image.

3.  **The Library Is The Enemy:** The `evo2` library actively works against the user. It does not follow standard conventions. We discovered through direct source code analysis that it:
    *   Lies about its own class name (`Evo2`, not `Evo`).
    *   Rejects standard Hugging Face instantiation methods (`from_pretrained`) in favor of its own constructor.
    *   Uses an internal, short-name (`evo2_40b`) for model identification instead of the public Hugging Face ID.

4.  **Hardware is the Final Boss:** Once all software challenges were overcome, the 40B model was too large for a single H100 GPU, causing a `CUDA out of memory` error. The final key was to request a dual-GPU environment (`gpu="H100:2"`), a solution you, Alpha, had presciently documented.

**Conclusion:**
Operation Mountbreaker was a tactical success that yielded strategic wisdom. We not only brought a powerful new weapon online but also established a new, robust doctrine for deploying future containerized, GPU-accelerated services. We are stronger and smarter for this fight. 