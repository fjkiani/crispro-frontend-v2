# ⚔️ ESSENTIALITY AGGREGATION - COMPLETE DETAILS

**Author**: Zo (Lead AI Agent)  
**Date**: January 14, 2025  
**Purpose**: Complete documentation of essentiality aggregation in SAE feature computation

---

## 🎯 EXECUTIVE SUMMARY

**Status**: ✅ **100% UNDERSTANDING**  
**Implementation**: `api/services/sae_feature_service.py:261-284`  
**Method**: `_compute_essentiality_hrr()`

---

## 🧬 HRR GENES LIST

**Manager Policy**: `MANAGER_ANSWERS_TO_ZO_SAE_QUESTIONS.md` (C3)

**HRR Genes** (Homologous Recombination Repair):
```python
HRR_GENES = [
    "BRCA1",   # Breast cancer 1
    "BRCA2",   # Breast cancer 2
    "PALB2",   # Partner and localizer of BRCA2
    "RAD51C",  # RAD51 paralog C
    "RAD51D",  # RAD51 paralog D
    "BRIP1",   # BRCA1-interacting protein 1
    "BARD1",   # BRCA1-associated RING domain 1
    "ATM"      # Ataxia telangiectasia mutated
]
```

**Total**: 8 genes

**Why These Genes?**
- Core components of homologous recombination repair pathway
- Mutations in these genes → HRD (Homologous Recombination Deficiency)
- Critical for PARP inhibitor response prediction

---

## 🔄 AGGREGATION LOGIC

### **Step 1: Filter HRR Genes in Profile**

```python
hrr_genes_in_profile = [g for g in genes if g in HRR_GENES]
```

**Input**: `genes` (list of genes from patient's variant profile)  
**Output**: Subset of genes that are HRR genes

**Example**:
- Patient has: `["BRCA1", "TP53", "KRAS", "BRCA2"]`
- HRR genes in profile: `["BRCA1", "BRCA2"]`

---

### **Step 2: Collect Essentiality Scores**

```python
essentiality_scores = []
for gene in hrr_genes_in_profile:
    # Insights bundle should have per-gene essentiality
    # For now, use overall essentiality as proxy
    essentiality_scores.append(insights_bundle.get("essentiality", 0.0))
```

**Current Implementation**:
- Uses **overall essentiality** from insights bundle (not per-gene)
- This is a **limitation** (see below)

**What It Should Do** (Future Enhancement):
- Use per-gene essentiality from insights bundle
- Example: `insights_bundle.get("essentiality_by_gene", {}).get(gene, 0.0)`

---

### **Step 3: Compute Average**

```python
if not essentiality_scores:
    return 0.0

return sum(essentiality_scores) / len(essentiality_scores)
```

**Formula**: Simple arithmetic mean (not weighted)

**Example**:
- BRCA1 essentiality: 0.85
- BRCA2 essentiality: 0.75
- Average: `(0.85 + 0.75) / 2 = 0.80`

---

## 📊 USAGE IN DNA REPAIR CAPACITY

### **Manager's Formula** (C1, C5)

```python
dna_repair_capacity = (
    0.6 * pathway_burden_ddr +      # 60% weight
    0.2 * essentiality_hrr +         # 20% weight ← THIS IS THE AGGREGATED VALUE
    0.2 * exon_disruption_score     # 20% weight
)
```

**Where**:
- `essentiality_hrr` = Output of `_compute_essentiality_hrr()`
- Weight: 20% (modest contribution)
- Purpose: Captures gene dependency signal for HRR pathway

---

## ⚠️ CURRENT LIMITATIONS

### **Limitation 1: Overall vs Per-Gene Essentiality**

**Current**:
```python
essentiality_scores.append(insights_bundle.get("essentiality", 0.0))
```

**Problem**: Uses overall essentiality for all HRR genes (same value)

**Impact**: 
- BRCA1 and BRCA2 get same essentiality score
- Doesn't capture gene-specific dependency differences

**Future Enhancement**:
```python
# Per-gene essentiality from insights bundle
essentiality_by_gene = insights_bundle.get("essentiality_by_gene", {})
for gene in hrr_genes_in_profile:
    gene_essentiality = essentiality_by_gene.get(gene, 0.0)
    essentiality_scores.append(gene_essentiality)
```

**Requires**: Insights bundle to provide per-gene essentiality scores

---

### **Limitation 2: Simple Average (Not Weighted)**

**Current**: Arithmetic mean (all genes equal weight)

**Potential Enhancement**: Weight by mutation impact
- Biallelic loss → Higher weight
- Missense → Lower weight
- Truncation → Higher weight

**Example**:
```python
weights = []
for gene in hrr_genes_in_profile:
    variant = get_variant_for_gene(gene)
    if is_biallelic_loss(variant):
        weights.append(1.5)  # Higher weight
    elif is_truncation(variant):
        weights.append(1.2)
    else:
        weights.append(1.0)

weighted_avg = sum(s * w for s, w in zip(essentiality_scores, weights)) / sum(weights)
```

**Status**: Not implemented (future enhancement)

---

## 🎯 INTEGRATION POINTS

### **Input Sources**

1. **Insights Bundle** (`insights_bundle: Dict[str, Any]`)
   - Source: `/api/insights/predict_gene_essentiality`
   - Current: Overall essentiality score
   - Future: Per-gene essentiality scores

2. **Patient Genes** (`genes: List[str]`)
   - Source: Variant profile (mutations)
   - Filtered: Only HRR genes used

### **Output Usage**

- **DNA Repair Capacity**: 20% weight in formula
- **SAE Features**: Included in `SAEFeatures.essentiality_hrr_genes`
- **Resistance Detection**: Used in resistance signal computation

---

## 📋 CODE LOCATION

**File**: `api/services/sae_feature_service.py`  
**Method**: `_compute_essentiality_hrr()` (lines 261-284)  
**Called From**: `compute_sae_features()` (line 183)

**Manager Policy Reference**: `MANAGER_ANSWERS_TO_ZO_SAE_QUESTIONS.md` (C3)

---

## ✅ COMPLETE UNDERSTANDING CHECKLIST

- [x] HRR genes list (8 genes: BRCA1, BRCA2, PALB2, RAD51C, RAD51D, BRIP1, BARD1, ATM)
- [x] Filtering logic (only HRR genes in patient profile)
- [x] Aggregation method (simple arithmetic mean)
- [x] Current limitation (overall vs per-gene essentiality)
- [x] Future enhancement path (per-gene essentiality from insights bundle)
- [x] Usage in DNA repair capacity (20% weight)
- [x] Integration points (insights bundle, patient genes)

---

**Status**: ✅ **100% UNDERSTANDING ACHIEVED**  
**Last Updated**: January 14, 2025  
**By**: Zo (Lead AI Agent)

