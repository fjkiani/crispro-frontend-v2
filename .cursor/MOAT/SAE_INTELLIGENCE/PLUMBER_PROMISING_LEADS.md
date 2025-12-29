# PLUMBER - Promising External Ovarian Data Leads

**Date:** January 28, 2025  
**Status:** 🔍 **PROMISING LEADS IDENTIFIED - NEED VERIFICATION**

---

## Web Search Findings

### 1. PTRC-HGSOC Collection ⭐ **MOST PROMISING**
**Source:** The Cancer Imaging Archive (TCIA)  
**URL:** https://stage.cancerimagingarchive.net/collection/ptrc-hgsoc/

**Description:**
- Proteogenomic data from HGSOC tumor biopsies
- **Annotated for patient sensitivity to platinum chemotherapy** (refractory or sensitive)
- High-grade serous ovarian cancer (HGSOC)

**Critical Checks Needed:**
- ✅ Platinum response labels: **CONFIRMED** (refractory/sensitive)
- ❓ Variant data: Need to verify if WES/WGS variants included
- ❓ Patient count: Need to verify (need ≥50)
- ❓ Variant coordinates: Need to verify (chrom, pos, ref, alt)

**Action:** **HIGH PRIORITY** - Investigate this collection immediately

---

### 2. ProteomeXchange PXD028225 ⚠️ **NEEDS VERIFICATION**
**Source:** ProteomeXchange  
**URL:** https://proteomecentral.proteomexchange.org/dataset/PXD028225

**Description:**
- Proteomic profiles comparing platinum-sensitive vs platinum-resistant HGSOC patients
- **Has platinum response labels**

**Critical Checks Needed:**
- ✅ Platinum response labels: **CONFIRMED**
- ❓ Variant data: **CRITICAL** - Does it include WES/WGS variants or only proteomics?
- ❓ Patient count: Need to verify (need ≥50)
- ❓ Variant coordinates: Need to verify

**Action:** **HIGH PRIORITY** - Verify if variant data included

---

### 3. GEO GSE274006 ⚠️ **NEEDS VERIFICATION**
**Source:** Gene Expression Omnibus  
**URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274006

**Description:**
- Gene expression profiles from ovarian cancer cell lines
- Focus on cell-intrinsic platinum response
- Genetic and gene expression signatures

**Critical Checks Needed:**
- ❓ Platinum response: Cell lines (not patients) - may not be suitable
- ❓ Variant data: Need to verify
- ❓ Patient count: Cell lines, not patient cohort

**Action:** **LOW PRIORITY** - Cell lines may not meet patient cohort requirement

---

### 4. GEO GSE274659 ⚠️ **NEEDS VERIFICATION**
**Source:** Gene Expression Omnibus  
**URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274659

**Description:**
- Gene expression profiles of platinum-resistant HGSOC
- **Has platinum resistance labels**

**Critical Checks Needed:**
- ✅ Platinum response: **CONFIRMED** (resistant)
- ❓ Variant data: Expression only? Need to verify if variants included
- ❓ Patient count: Need to verify
- ❓ Variant coordinates: Need to verify

**Action:** **MEDIUM PRIORITY** - Verify if variant data included

---

## Priority Investigation Order

### Priority 1: PTRC-HGSOC Collection ⭐
**Why:** Proteogenomic = likely has both proteomics AND genomics  
**Action:**
1. Access TCIA portal
2. Check for WES/WGS variant data (VCF/MAF files)
3. Verify patient count
4. Verify platinum response labels format
5. Download if suitable

### Priority 2: ProteomeXchange PXD028225
**Why:** Confirmed platinum response labels  
**Action:**
1. Access PRIDE/ProteomeXchange portal
2. Check dataset contents
3. Verify if WES/WGS variants included
4. Check patient count
5. Download if suitable

### Priority 3: GEO GSE274659
**Why:** Has platinum resistance labels  
**Action:**
1. Check GEO dataset page
2. Verify if variant data available
3. Check patient count
4. Download if suitable

---

## Current Status Summary

### ✅ Extracted
- MSK Ovarian 2025: 27 patients, 116 mutations (coordinates), no platinum response

### ⏳ Needs Investigation
- PTRC-HGSOC: Platinum response ✅, variant data ❓
- PXD028225: Platinum response ✅, variant data ❓
- GSE274659: Platinum response ✅, variant data ❓

### ❌ Checked (No Platinum Response)
- All cBioPortal ovarian studies (patient-level and sample-level clinical data)
- Treatment endpoints (not available via API)

---

## Next Immediate Actions

1. **Investigate PTRC-HGSOC Collection** (Highest Priority)
   - Access TCIA portal
   - Check for variant data
   - Verify patient count
   - Extract if suitable

2. **Verify ProteomeXchange PXD028225**
   - Check if WES/WGS variants included
   - Verify patient count
   - Extract if suitable

3. **Check GEO GSE274659**
   - Verify variant data availability
   - Check patient count
   - Extract if suitable

---

**Created:** January 28, 2025  
**Status:** 🔍 **PROMISING LEADS IDENTIFIED**  
**Next:** Investigate PTRC-HGSOC and PXD028225

