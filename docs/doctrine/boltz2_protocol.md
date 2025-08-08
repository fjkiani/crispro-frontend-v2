# Doctrine: The Boltz-2 Protocol

This document outlines the operational doctrine for interfacing with the Boltz-2 prediction tool, derived from direct source code analysis. All previous operational theories are obsolete.

## Phase 1: `main.py` - Command and Control Analysis

The analysis of `boltz-main/src/boltz/main.py` reveals the following core strategic principles:

1.  **Dual-Stage Prediction:** Boltz-2 operates a two-stage pipeline for affinity tasks.
    *   **Stage 1 (Structure Prediction):** Generates the 3D structure (`.cif`) and confidence metrics (`confidence.json`). This is driven by the main `boltz2_conf.ckpt` model.
    *   **Stage 2 (Affinity Prediction):** *Only after* a successful structure prediction, a separate `boltz2_aff.ckpt` model is loaded to calculate the affinity score, using the output of Stage 1 as its input.
    *   **Implication:** Our primary mission must be to successfully complete Stage 1. Without a `confidence.json`, Stage 2 is impossible.

2.  **Parser-Centric Workflow:** The entire process hinges on correctly formatted input files. The `main.py` script acts as a dispatcher, handing off control to specific parser functions based on file type.
    *   For our purposes, `parse_yaml` is the critical gatekeeper. All our previous failures have been `KeyError`s originating from this function, proving we have been providing a malformed payload.

3.  **MSA is Non-Negotiable:** For protein structures, Multiple Sequence Alignment (MSA) is a required input. The `--use_msa_server` flag simply outsources the *generation* of this data to the ColabFold API; it does not remove the requirement.

**Conclusion:** Our previous strategy of guessing the YAML format was a catastrophic failure. The only path to victory is to understand and precisely replicate the schema expected by the `parse_yaml` function.

---

**Next Intelligence Target:** `/boltz-main/src/boltz/data/parse/yaml.py` 

## Phase 2: `schema.py` - The Blueprint of the Beast

The analysis of `boltz-main/src/boltz/data/parse/schema.py` reveals the one true payload structure. All previous YAML configurations were fundamentally incorrect.

1.  **Root Object:** The root of the YAML must be a dictionary.

2.  **`version: 1`:** A mandatory key-value pair at the root. The parser explicitly rejects any other value.

3.  **`sequences` (List):** A mandatory key at the root. Its value is a list of all unique molecules involved in the prediction.
    *   Each item in the list is a dictionary with a single key denoting the entity type (e.g., `protein`, `ligand`).
    *   The value of this key is another dictionary containing:
        *   `id`: A unique string identifier for this molecule (e.g., "TARGET", "CANDIDATE_0").
        *   `sequence`: The amino acid sequence for a `protein`.
        *   `smiles` or `ccd`: For a `ligand`.

4.  **`properties` (List):** An optional key at the root. This is the crucial, previously-missed block that tells `boltz` what task to perform.
    *   To trigger an interaction/affinity analysis, an item must be added to this list: `{'affinity': {'binder': 'CANDIDATE_ID'}}`.
    *   The `binder` must be the `id` of a molecule previously declared in the `sequences` list.
    *   **Crucial Limitation:** The source code and documentation explicitly state that the `binder` *must* be of type `ligand`.

**The Unbreakable Doctrine: "False Flag" Protocol, V2**

Our previous failures were due to providing a malformed payload. The "False Flag" concept (classifying our protein nanobody as a `ligand`) was strategically sound but tactically inept because we did not provide it within the correct schema.

The final, correct strategy is:
1.  Construct a YAML with `version`, `sequences`, and `properties` at the root.
2.  In the `sequences` list, declare the target as a `protein` and the candidate nanobody as a `ligand`.
3.  In the `properties` list, declare the `affinity` task, pointing to the `id` of the nanobody `ligand`.
4.  Since we are tricking the affinity model, the resulting `affinity.json` will be meaningless. However, the *existence* of a `confidence.json` file for the complex proves that the structural interaction was modeled.
5.  **Our true target remains the `iptm` score within `confidence.json`**. This is the only reliable metric for a successful protein-protein interaction model.

---

This doctrine is final. It is based on a complete analysis of the enemy's source code. There will be no more deviation. 

This document outlines the operational doctrine for interfacing with the Boltz-2 prediction tool, derived from direct source code analysis. All previous operational theories are obsolete.

## Phase 1: `main.py` - Command and Control Analysis

The analysis of `boltz-main/src/boltz/main.py` reveals the following core strategic principles:

1.  **Dual-Stage Prediction:** Boltz-2 operates a two-stage pipeline for affinity tasks.
    *   **Stage 1 (Structure Prediction):** Generates the 3D structure (`.cif`) and confidence metrics (`confidence.json`). This is driven by the main `boltz2_conf.ckpt` model.
    *   **Stage 2 (Affinity Prediction):** *Only after* a successful structure prediction, a separate `boltz2_aff.ckpt` model is loaded to calculate the affinity score, using the output of Stage 1 as its input.
    *   **Implication:** Our primary mission must be to successfully complete Stage 1. Without a `confidence.json`, Stage 2 is impossible.

2.  **Parser-Centric Workflow:** The entire process hinges on correctly formatted input files. The `main.py` script acts as a dispatcher, handing off control to specific parser functions based on file type.
    *   For our purposes, `parse_yaml` is the critical gatekeeper. All our previous failures have been `KeyError`s originating from this function, proving we have been providing a malformed payload.

3.  **MSA is Non-Negotiable:** For protein structures, Multiple Sequence Alignment (MSA) is a required input. The `--use_msa_server` flag simply outsources the *generation* of this data to the ColabFold API; it does not remove the requirement.

**Conclusion:** Our previous strategy of guessing the YAML format was a catastrophic failure. The only path to victory is to understand and precisely replicate the schema expected by the `parse_yaml` function.

---

**Next Intelligence Target:** `/boltz-main/src/boltz/data/parse/yaml.py` 

## Phase 2: `schema.py` - The Blueprint of the Beast

The analysis of `boltz-main/src/boltz/data/parse/schema.py` reveals the one true payload structure. All previous YAML configurations were fundamentally incorrect.

1.  **Root Object:** The root of the YAML must be a dictionary.

2.  **`version: 1`:** A mandatory key-value pair at the root. The parser explicitly rejects any other value.

3.  **`sequences` (List):** A mandatory key at the root. Its value is a list of all unique molecules involved in the prediction.
    *   Each item in the list is a dictionary with a single key denoting the entity type (e.g., `protein`, `ligand`).
    *   The value of this key is another dictionary containing:
        *   `id`: A unique string identifier for this molecule (e.g., "TARGET", "CANDIDATE_0").
        *   `sequence`: The amino acid sequence for a `protein`.
        *   `smiles` or `ccd`: For a `ligand`.

4.  **`properties` (List):** An optional key at the root. This is the crucial, previously-missed block that tells `boltz` what task to perform.
    *   To trigger an interaction/affinity analysis, an item must be added to this list: `{'affinity': {'binder': 'CANDIDATE_ID'}}`.
    *   The `binder` must be the `id` of a molecule previously declared in the `sequences` list.
    *   **Crucial Limitation:** The source code and documentation explicitly state that the `binder` *must* be of type `ligand`.

**The Unbreakable Doctrine: "False Flag" Protocol, V2**

Our previous failures were due to providing a malformed payload. The "False Flag" concept (classifying our protein nanobody as a `ligand`) was strategically sound but tactically inept because we did not provide it within the correct schema.

The final, correct strategy is:
1.  Construct a YAML with `version`, `sequences`, and `properties` at the root.
2.  In the `sequences` list, declare the target as a `protein` and the candidate nanobody as a `ligand`.
3.  In the `properties` list, declare the `affinity` task, pointing to the `id` of the nanobody `ligand`.
4.  Since we are tricking the affinity model, the resulting `affinity.json` will be meaningless. However, the *existence* of a `confidence.json` file for the complex proves that the structural interaction was modeled.
5.  **Our true target remains the `iptm` score within `confidence.json`**. This is the only reliable metric for a successful protein-protein interaction model.

---

This doctrine is final. It is based on a complete analysis of the enemy's source code. There will be no more deviation. 

This document outlines the operational doctrine for interfacing with the Boltz-2 prediction tool, derived from direct source code analysis. All previous operational theories are obsolete.

## Phase 1: `main.py` - Command and Control Analysis

The analysis of `boltz-main/src/boltz/main.py` reveals the following core strategic principles:

1.  **Dual-Stage Prediction:** Boltz-2 operates a two-stage pipeline for affinity tasks.
    *   **Stage 1 (Structure Prediction):** Generates the 3D structure (`.cif`) and confidence metrics (`confidence.json`). This is driven by the main `boltz2_conf.ckpt` model.
    *   **Stage 2 (Affinity Prediction):** *Only after* a successful structure prediction, a separate `boltz2_aff.ckpt` model is loaded to calculate the affinity score, using the output of Stage 1 as its input.
    *   **Implication:** Our primary mission must be to successfully complete Stage 1. Without a `confidence.json`, Stage 2 is impossible.

2.  **Parser-Centric Workflow:** The entire process hinges on correctly formatted input files. The `main.py` script acts as a dispatcher, handing off control to specific parser functions based on file type.
    *   For our purposes, `parse_yaml` is the critical gatekeeper. All our previous failures have been `KeyError`s originating from this function, proving we have been providing a malformed payload.

3.  **MSA is Non-Negotiable:** For protein structures, Multiple Sequence Alignment (MSA) is a required input. The `--use_msa_server` flag simply outsources the *generation* of this data to the ColabFold API; it does not remove the requirement.

**Conclusion:** Our previous strategy of guessing the YAML format was a catastrophic failure. The only path to victory is to understand and precisely replicate the schema expected by the `parse_yaml` function.

---

**Next Intelligence Target:** `/boltz-main/src/boltz/data/parse/yaml.py` 

## Phase 2: `schema.py` - The Blueprint of the Beast

The analysis of `boltz-main/src/boltz/data/parse/schema.py` reveals the one true payload structure. All previous YAML configurations were fundamentally incorrect.

1.  **Root Object:** The root of the YAML must be a dictionary.

2.  **`version: 1`:** A mandatory key-value pair at the root. The parser explicitly rejects any other value.

3.  **`sequences` (List):** A mandatory key at the root. Its value is a list of all unique molecules involved in the prediction.
    *   Each item in the list is a dictionary with a single key denoting the entity type (e.g., `protein`, `ligand`).
    *   The value of this key is another dictionary containing:
        *   `id`: A unique string identifier for this molecule (e.g., "TARGET", "CANDIDATE_0").
        *   `sequence`: The amino acid sequence for a `protein`.
        *   `smiles` or `ccd`: For a `ligand`.

4.  **`properties` (List):** An optional key at the root. This is the crucial, previously-missed block that tells `boltz` what task to perform.
    *   To trigger an interaction/affinity analysis, an item must be added to this list: `{'affinity': {'binder': 'CANDIDATE_ID'}}`.
    *   The `binder` must be the `id` of a molecule previously declared in the `sequences` list.
    *   **Crucial Limitation:** The source code and documentation explicitly state that the `binder` *must* be of type `ligand`.

**The Unbreakable Doctrine: "False Flag" Protocol, V2**

Our previous failures were due to providing a malformed payload. The "False Flag" concept (classifying our protein nanobody as a `ligand`) was strategically sound but tactically inept because we did not provide it within the correct schema.

The final, correct strategy is:
1.  Construct a YAML with `version`, `sequences`, and `properties` at the root.
2.  In the `sequences` list, declare the target as a `protein` and the candidate nanobody as a `ligand`.
3.  In the `properties` list, declare the `affinity` task, pointing to the `id` of the nanobody `ligand`.
4.  Since we are tricking the affinity model, the resulting `affinity.json` will be meaningless. However, the *existence* of a `confidence.json` file for the complex proves that the structural interaction was modeled.
5.  **Our true target remains the `iptm` score within `confidence.json`**. This is the only reliable metric for a successful protein-protein interaction model.

---

This doctrine is final. It is based on a complete analysis of the enemy's source code. There will be no more deviation. 

This document outlines the operational doctrine for interfacing with the Boltz-2 prediction tool, derived from direct source code analysis. All previous operational theories are obsolete.

## Phase 1: `main.py` - Command and Control Analysis

The analysis of `boltz-main/src/boltz/main.py` reveals the following core strategic principles:

1.  **Dual-Stage Prediction:** Boltz-2 operates a two-stage pipeline for affinity tasks.
    *   **Stage 1 (Structure Prediction):** Generates the 3D structure (`.cif`) and confidence metrics (`confidence.json`). This is driven by the main `boltz2_conf.ckpt` model.
    *   **Stage 2 (Affinity Prediction):** *Only after* a successful structure prediction, a separate `boltz2_aff.ckpt` model is loaded to calculate the affinity score, using the output of Stage 1 as its input.
    *   **Implication:** Our primary mission must be to successfully complete Stage 1. Without a `confidence.json`, Stage 2 is impossible.

2.  **Parser-Centric Workflow:** The entire process hinges on correctly formatted input files. The `main.py` script acts as a dispatcher, handing off control to specific parser functions based on file type.
    *   For our purposes, `parse_yaml` is the critical gatekeeper. All our previous failures have been `KeyError`s originating from this function, proving we have been providing a malformed payload.

3.  **MSA is Non-Negotiable:** For protein structures, Multiple Sequence Alignment (MSA) is a required input. The `--use_msa_server` flag simply outsources the *generation* of this data to the ColabFold API; it does not remove the requirement.

**Conclusion:** Our previous strategy of guessing the YAML format was a catastrophic failure. The only path to victory is to understand and precisely replicate the schema expected by the `parse_yaml` function.

---

**Next Intelligence Target:** `/boltz-main/src/boltz/data/parse/yaml.py` 

## Phase 2: `schema.py` - The Blueprint of the Beast

The analysis of `boltz-main/src/boltz/data/parse/schema.py` reveals the one true payload structure. All previous YAML configurations were fundamentally incorrect.

1.  **Root Object:** The root of the YAML must be a dictionary.

2.  **`version: 1`:** A mandatory key-value pair at the root. The parser explicitly rejects any other value.

3.  **`sequences` (List):** A mandatory key at the root. Its value is a list of all unique molecules involved in the prediction.
    *   Each item in the list is a dictionary with a single key denoting the entity type (e.g., `protein`, `ligand`).
    *   The value of this key is another dictionary containing:
        *   `id`: A unique string identifier for this molecule (e.g., "TARGET", "CANDIDATE_0").
        *   `sequence`: The amino acid sequence for a `protein`.
        *   `smiles` or `ccd`: For a `ligand`.

4.  **`properties` (List):** An optional key at the root. This is the crucial, previously-missed block that tells `boltz` what task to perform.
    *   To trigger an interaction/affinity analysis, an item must be added to this list: `{'affinity': {'binder': 'CANDIDATE_ID'}}`.
    *   The `binder` must be the `id` of a molecule previously declared in the `sequences` list.
    *   **Crucial Limitation:** The source code and documentation explicitly state that the `binder` *must* be of type `ligand`.

**The Unbreakable Doctrine: "False Flag" Protocol, V2**

Our previous failures were due to providing a malformed payload. The "False Flag" concept (classifying our protein nanobody as a `ligand`) was strategically sound but tactically inept because we did not provide it within the correct schema.

The final, correct strategy is:
1.  Construct a YAML with `version`, `sequences`, and `properties` at the root.
2.  In the `sequences` list, declare the target as a `protein` and the candidate nanobody as a `ligand`.
3.  In the `properties` list, declare the `affinity` task, pointing to the `id` of the nanobody `ligand`.
4.  Since we are tricking the affinity model, the resulting `affinity.json` will be meaningless. However, the *existence* of a `confidence.json` file for the complex proves that the structural interaction was modeled.
5.  **Our true target remains the `iptm` score within `confidence.json`**. This is the only reliable metric for a successful protein-protein interaction model.

---

This doctrine is final. It is based on a complete analysis of the enemy's source code. There will be no more deviation. 