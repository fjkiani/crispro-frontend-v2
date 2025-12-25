# 05 — Mechanism-Targeted Therapy (Pathway Bin Decisioning)

**Purpose:** If DDR_bin is not the only escape signal, choose therapy based on *which* pathway is rising.

---

## Inputs

- DDR_bin, MAPK_bin, PI3K_bin, Efflux_bin (as available)
- Current therapy class
- Eligibility constraints (comorbidities, prior therapies)

---

## Rule (Template)

- If **DDR_bin↓** and other bins stable → “HR restoration” playbook
- If **MAPK_bin↑** → MAPK bypass playbook (MEK/RAF strategies)
- If **PI3K_bin↑** → PI3K bypass playbook
- If **Efflux_bin↑** → switch away from exported drugs / consider non-PARP strategy

---

## Outputs

- Mechanism hypothesis
- Tiered actions + evidence tiers
- Next tests targeted to confirm (e.g., RAD51C/BRCA reversions, KRAS/NF1)

