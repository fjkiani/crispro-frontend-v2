# AYESHA CASE: VALIDATED vs NOT VALIDATED (Agent Ground Truth)

**Date:** January 27, 2026  
**Purpose:** Single source of truth for all agents on what is VALIDATED vs ASSUMED  
**Approach:** No-assumption approach - explicitly mark confidence levels  

---

## ✅ VALIDATED (High Confidence)

### Patient Data (Ground Truth)
| Item | Value | Source | Confidence |
|------|-------|--------|------------|
| **MBD4** | p.K431Nfs*54 (homozygous frameshift) | Germline test | ✅ 100% |
| **TP53** | Mutant type (IHC positive) | Pathology | ✅ 100% |
| **PDGFRA** | p.S755P (VUS) | Germline test | ✅ 100% |
| **PD-L1 CPS** | 10 (POSITIVE) | IHC 22C3 | ✅ 100% |
| **MSI Status** | MSS (stable) | IHC | ✅ 100% |
| **Stage** | IVB HGSOC with carcinomatosis | Imaging/Pathology | ✅ 100% |
| **ER/PR** | ER weak+, PR- | Pathology | ✅ 100% |

### PDGFRA VUS Coordinates (Resolved)
| Item | Value | Source | Confidence |
|------|-------|--------|------------|
| **Chromosome** | 4 | Ensembl VEP | ✅ 100% |
| **Position** | 54280422 | Ensembl VEP | ✅ 100% |
| **Ref/Alt** | T>C | Ensembl VEP | ✅ 100% |
| **Assembly** | GRCh38 | Ensembl VEP | ✅ 100% |
| **VEP Impact** | MODERATE (missense) | Ensembl VEP | ✅ 100% |
| **ClinVar** | Uncertain Significance | ClinVar | ✅ 100% |

### Code Implementation (Verified)
| Item | Location | Validation | Confidence |
|------|----------|------------|------------|
| **SL Service Exists** | `api/services/synthetic_lethality/` | Code inspection | ✅ 100% |
| **BER→HR/PARP Mapping** | `constants.py` line 100-102 | Code inspection | ✅ 100% |
| **SL Pathway Accuracy** | 62.2% on GDSC2 | Paper validation | ✅ 100% |
| **S/P/E Weights** | 0.35/0.35/0.30 | `config.py` | ✅ 100% |
| **VUS Router Exists** | `api/routers/vus.py` | Code inspection | ✅ 100% |
| **Evo2 Service Exists** | `api/routers/evo.py` | Code inspection | ✅ 100% |

### Biological Mechanisms (Literature-Supported)
| Item | Evidence | Source | Confidence |
|------|----------|--------|------------|
| **MBD4 is BER gene** | DNA glycosylase | Literature | ✅ 100% |
| **MBD4 loss → BER deficiency** | Established biology | Literature | ✅ 100% |
| **BER deficiency ≠ HRD** | Different pathways | Literature | ✅ 100% |
| **PDGFRA is RTK** | Receptor tyrosine kinase | Literature | ✅ 100% |

---

## ⚠️ MECHANISM VALIDATED, CLINICAL OUTCOMES NOT VALIDATED

### MBD4 → PARP Synthetic Lethality
| Aspect | Status | Evidence |
|--------|--------|----------|
| **Mechanism (BER→HR dependency)** | ✅ VALIDATED | SL service implements, 62.2% accuracy |
| **Code Implementation** | ✅ VALIDATED | `constants.py` maps BER→HR/PARP |
| **Clinical Outcomes** | ⚠️ NOT VALIDATED | No patient outcome data (MBD4 too rare) |
| **Confidence for PARP Drugs** | ⚠️ ASSUMED 70% | Not validated against clinical trials |

### Synthetic Lethality Framework
| Aspect | Status | Evidence |
|--------|--------|----------|
| **SL Service Exists** | ✅ VALIDATED | Code inspection |
| **Pathway Mapping Works** | ✅ VALIDATED | 62.2% accuracy on GDSC2 |
| **Drug Recommendations** | ✅ VALIDATED | Code implements PARP for BER deficiency |
| **Patient Outcomes** | ⚠️ NOT VALIDATED | No clinical trial data for MBD4→PARP |

---

## ❌ NOT VALIDATED (Assumed - Need Verification)

### PDGFRA VUS Analysis
| Item | Assumed Value | Actual Status |
|------|---------------|---------------|
| **Evo2 Delta Score** | -0.08 (assumed) | ❌ NOT RUN - Cache hit invalidated |
| **Functionality Score** | 0.60 (assumed) | ❌ NOT RUN |
| **Essentiality Score** | 0.35 (assumed) | ❌ NOT RUN |
| **PDGFRA→VEGF Mapping** | Weight 0.3 (assumed) | ❌ NOT IN CODE - PDGFRA not in `get_pathway_weights_for_gene()` |
| **AlphaMissense Score** | 0.5-0.7 (assumed) | ❌ NOT RUN |
| **VUS Classification** | "LIKELY PATHOGENIC" (claimed) | ❌ INVALIDATED - Based on cache hit |

### Drug Confidence Scores
| Drug | Assumed Confidence | Validation Status |
|------|-------------------|-------------------|
| **Olaparib** | 70% | ⚠️ Mechanism validated, confidence not validated |
| **Niraparib** | 65% | ⚠️ Mechanism validated, confidence not validated |
| **Pembrolizumab** | 65% | ❌ TMB not confirmed by NGS |
| **Bevacizumab** | 60% + 5% boost | ❌ Boost not validated |
| **Imatinib** | 35% | ❌ No ovarian cancer evidence |
| **Sunitinib** | 40% | ❌ No ovarian cancer evidence |
| **Regorafenib** | 40% | ❌ No ovarian cancer evidence |

### Missing Data
| Item | Status | Impact |
|------|--------|--------|
| **HRD Score** | NOT TESTED | Cannot determine HRD status |
| **TMB** | NOT CONFIRMED | IO eligibility assumed from MBD4 |
| **CA-125** | NOT IN PROFILE | Missing monitoring data |
| **Expression Data** | NOT AVAILABLE | Cannot enhance pathway analysis |
| **NGS Somatic Panel** | NOT AVAILABLE | TP53 variant inferred from IHC |

---

## 📊 CONFIDENCE SCALE FOR AGENTS

### When Making Claims:

| Confidence Level | When to Use | Example |
|-----------------|-------------|---------|
| **✅ VALIDATED** | Code implemented + accuracy measured | "SL service exists (62.2% accuracy)" |
| **✅ LITERATURE** | Published evidence supports | "MBD4 is a BER gene" |
| **⚠️ MECHANISM VALIDATED** | Code works, outcomes unknown | "MBD4→PARP via SL (mechanism validated)" |
| **⚠️ ASSUMED** | Reasonable inference, not validated | "TMB likely HIGH from MBD4" |
| **❌ NOT RUN** | Claimed but not actually executed | "Evo2 delta = -0.08 (NOT RUN)" |
| **❌ NOT IN CODE** | Claimed mapping doesn't exist | "PDGFRA→VEGF (NOT IN CODE)" |
| **❌ CACHE HIT** | Stale data, not fresh computation | "Evo2 result (CACHE HIT - INVALID)" |

---

## 🔴 CRITICAL CHECKS FOR ALL AGENTS

### Before Claiming Evo2/API Results:
```python
# ALWAYS check cache_hit in response
response = await call_evo2(...)
if response.get("provenance", {}).get("cache_hit"):
    # THIS IS STALE DATA - DO NOT CLAIM AS VALIDATED
    return "❌ CACHE HIT - NOT VALIDATED"
```

### Before Claiming Pathway Mapping:
```python
# VERIFY gene is in actual mapping function
from drug_mapping import get_pathway_weights_for_gene
weights = get_pathway_weights_for_gene("PDGFRA")
if not weights:
    # GENE NOT IN MAPPING - DO NOT CLAIM
    return "❌ NOT IN CODE"
```

### Before Claiming Clinical Outcomes:
```python
# CHECK if clinical trial data exists
# MBD4 is too rare for clinical trials
# Therefore: mechanism ✅ validated, outcomes ❌ not validated
```

---

## 📋 WHAT AGENTS CAN CONFIDENTLY SAY

### ✅ CAN SAY (Validated):
- "SL service implements BER→PARP pathway (62.2% accuracy)"
- "MBD4 is a BER gene (literature)"
- "PDGFRA VUS coordinates are chr4:54280422 T>C (GRCh38)"
- "ClinVar classifies PDGFRA S755P as VUS"
- "S/P/E weights are 0.35/0.35/0.30 (code verified)"

### ⚠️ CAN SAY WITH CAVEAT:
- "MBD4→PARP via SL (mechanism validated, clinical outcomes not validated)"
- "TMB likely HIGH from MBD4 hypermutation phenotype (assumed, not NGS confirmed)"
- "PARP drugs may benefit MBD4 patients (mechanism validated, confidence score assumed)"

### ❌ CANNOT SAY (Not Validated):
- "PDGFRA S755P is LIKELY PATHOGENIC" (cache hit - invalid)
- "Evo2 delta = -0.08 for PDGFRA" (not actually run)
- "PDGFRA maps to VEGF pathway with weight 0.3" (not in code)
- "TKI drugs have 35-40% confidence for PDGFRA VUS" (assumed)
- "Bevacizumab gets +5% boost from PDGFRA" (not validated)

---

## 🔧 ACTIONS NEEDED TO VALIDATE ASSUMPTIONS

### P0 (Critical):
| # | Item | How to Validate |
|---|------|-----------------|
| 1 | PDGFRA Evo2 Score | Run fresh Evo2 with `EVO_CACHE_DISABLE=1` |
| 2 | PDGFRA Pathway | Add to `get_pathway_weights_for_gene()` if appropriate |
| 3 | S/P/E Documentation | Fix docs (says 30/40/30, code is 35/35/30) |
| 4 | MBD4 in SAE | Verify MBD4 in DDR_bin calculation |

### P1 (High):
| # | Item | How to Validate |
|---|------|-----------------|
| 5 | TKI Ovarian Evidence | Literature search for PDGFRA TKIs in ovarian |
| 6 | PARP Confidence | Define confidence calibration for SL-based drugs |
| 7 | TMB Calculation | Wait for NGS or build TMB estimator |
| 8 | VUS Threshold | Run fresh Evo2, calculate percentile |

---

**For Agent X and all agents: This is the ground truth. Before claiming anything, check this document for validation status.**
