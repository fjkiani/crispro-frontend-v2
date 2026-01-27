# 🔬 AYESHA ENGINE AUDIT - What's Built vs What's Displayed

**Last Updated:** January 26, 2026  
**Synced With:** `AYESHA_FRONTEND_PRODUCTION_AUDIT.md`  
**Purpose:** Identify built engines NOT connected to Ayesha's main dashboard

---

## 🚨 CRITICAL FINDING: MBD4 NOT IN DDR GENE LISTS

### The Problem
**MBD4 is NOT in the DDR engine's gene lists!**

From `api/services/resistance/config/ddr_config.py`:
```python
"ovary": {
    "core_hrr_genes": [
        "BRCA1", "BRCA2", "RAD51C", "RAD51D",
        "PALB2", "BARD1", "BRIP1"
    ],
    "extended_ddr_genes": [
        "ATM", "ATR", "CHEK1", "CHEK2",
        "FANCA", "FANCD2", "RAD50", "MRE11", "NBN", "POLQ"
    ],
}
```

**MBD4 is NOT in either list!**

### Impact on Ayesha
If we call the DDR endpoint with MBD4 mutation:
1. MBD4 won't match `core_hrr_genes` → `core_HRR_pathogenic = False`
2. MBD4 won't match `extended_ddr_genes` → `extended_DDR_pathogenic = False`  
3. No HRD assay data → `hrd_positive_inferred = False`
4. **Result: `DDR_bin_status = "unknown"`**

### Biological Reality
MBD4 is a **Base Excision Repair (BER)** gene, NOT an HRR gene. But:
- MBD4 loss → BER pathway loss → DDR deficiency
- MBD4 + TP53 → Synthetic lethality with PARP

### Solution Required
**Option A:** Add MBD4 to `extended_ddr_genes` in ddr_config.py
**Option B:** Create special handling for MBD4 in the scoring logic
**Option C:** Use the Synthetic Lethality engine (already detects MBD4+TP53!)

---

## 📋 EXECUTIVE SUMMARY - REVISED

### The Real Gaps (Priority Order)

| # | Gap | Impact | Solution |
|---|-----|--------|----------|
| **G1** | MBD4 not in DDR gene lists | DDR engine returns "unknown" for Ayesha | Add MBD4 to ddr_config.py OR use SL engine |
| **G2** | DDR components on separate page | User doesn't see DDR status | Import to AyeshaTrialExplorer |
| **G3** | WIWFM awaiting NGS | No drug predictions for Ayesha | Use germline fallback OR wait for NGS |
| **G4** | CA-125 null | CA-125 Intelligence blocked | Add entry form |

---

## ✓ VERIFIED: WHAT WORKS

### 1. Synthetic Lethality Engine ✅
- **DOES detect MBD4+TP53 combination**
- Shows on Tab 5 (Synthetic Lethality)
- Shows PARP as recommended therapy
- **This is the workaround for G1!**

### 2. DDR Components ✅
- `DDRStatusCard.jsx` - Shows DDR_bin_status
- `DDRTreatmentEligibility.jsx` - Shows "PARP ELIGIBLE"
- Exports correctly from `components/ddr/index.js`
- **But will show "unknown" without G1 fix**

### 3. API Endpoints ✅
- `POST /api/resistance/ddr-status` - Working
- Request/Response schemas verified
- **But MBD4 won't trigger DDR_defective**

---

## 🎯 REVISED DELIVERABLES

### D0: Add MBD4 to DDR Gene Config (30 minutes) - **NEW**

**Why:** MBD4 not in gene lists → DDR engine returns "unknown"

**File:** `api/services/resistance/config/ddr_config.py`

**Change:**
```python
"ovary": {
    ...
    "extended_ddr_genes": [
        "ATM", "ATR", "CHEK1", "CHEK2",
        "FANCA", "FANCD2", "RAD50", "MRE11", "NBN", "POLQ",
        "MBD4"  # Add MBD4 - BER pathway gene
    ],
}
```

**Alternative:** Add to default config as well for pan-cancer support.

---

### D1: Add DDR Components to Ayesha's Page (4-6 hours)

**Why:** DDR status not visible on main dashboard

**Prerequisites:** 
- D0 must be done first (otherwise shows "unknown")
- OR use SL engine result instead of DDR engine

**Files to Modify:**
- `AyeshaTrialExplorer.jsx` - Add imports and useEffect

**Verified Request Format:**
```jsx
const requestData = {
  patient_id: 'AK',
  disease_site: 'ovary',
  tumor_subtype: 'HGSOC',
  mutations: [
    { gene_symbol: 'MBD4', variant_classification: 'pathogenic' }
  ]
};
```

**Verified Response Fields:**
```javascript
// ddrStatus object shape:
{
  DDR_bin_status: "DDR_defective" | "DDR_proficient" | "unknown",
  DDR_score: 1.0,  // Only 1.0 if extended_DDR_pathogenic
  extended_DDR_pathogenic: true,  // Will be true if MBD4 added to config
  HRD_status_inferred: "unknown",  // No HRD assay
  provenance: {...}
}
```

---

### D2: Treatment Options Summary (6-8 hours)

**Why:** No unified view of options

**Data Sources (Verified):**
- `socRecommendation` - From complete_care_v2
- `trials` - Array from complete_care_v2
- `wiwfm` - Drug efficacy (blocked by "awaiting_ngs")
- `slResult` - Synthetic lethality detection

---

### D3: CA-125 Entry Form (4 hours)

**Why:** CA-125 Intelligence blocked

---

## 📊 VERIFIED API SCHEMAS

### DDR Status Endpoint

**Request (`POST /api/resistance/ddr-status`):**
```python
class DDRStatusRequest(BaseModel):
    patient_id: str                           # Required
    disease_site: str                         # Required - "ovary"
    tumor_subtype: Optional[str] = None       # Optional - "HGSOC"
    mutations: List[MutationInputDDR] = []    # Required
    cna: Optional[List[CNAInputDDR]] = None   # Optional
    hrd_assay: Optional[HRDAssayInputDDR] = None  # Optional

class MutationInputDDR(BaseModel):
    gene_symbol: str           # "MBD4"
    variant_classification: str # "pathogenic"
    variant_type: Optional[str] = None
```

**Response:**
```python
class DDRStatusResponse(BaseModel):
    patient_id: str
    DDR_bin_status: str          # "DDR_defective" | "unknown"
    DDR_score: float             # 0.0 - ~8.5
    BRCA_pathogenic: bool
    core_HRR_pathogenic: bool
    extended_DDR_pathogenic: bool  # True if MBD4 added to config
    HRD_status_inferred: str
    DDR_features_used: Dict
    provenance: Dict
```

---

## 📊 DDR SCORING LOGIC (Verified)

From `ddr_bin_scoring.py` lines 278-324:

```python
# Priority-ordered rules:
1. BRCA_pathogenic → DDR_defective
2. HRD_positive_inferred → DDR_defective
3. core_HRR_pathogenic → DDR_defective
4. extended_DDR_pathogenic → DDR_defective  # ← MBD4 would trigger this (if added to config)
5. no DDR/HRD data → unknown
6. else → DDR_proficient
```

**Score Weights:**
```python
BRCA_pathogenic: 3.0
HRD_positive: 2.5
core_hrr_pathogenic: 2.0
extended_ddr_pathogenic: 1.0  # MBD4 would get this score
```

---

## ✅ ACTION PLAN (Honest Assessment)

### Before ANY Frontend Work:

**Step 1:** Add MBD4 to ddr_config.py (30 min)
- Without this, DDR components will show "unknown"
- Alternative: Use SL result which already detects MBD4+TP53

**Step 2:** Test DDR endpoint with MBD4 (15 min)
```bash
curl -X POST http://localhost:8000/api/resistance/ddr-status \
  -H 'Content-Type: application/json' \
  -d '{
    "patient_id": "AK",
    "disease_site": "ovary", 
    "tumor_subtype": "HGSOC",
    "mutations": [{"gene_symbol": "MBD4", "variant_classification": "pathogenic"}]
  }'
```

**Step 3:** Verify response shows `DDR_bin_status: "DDR_defective"`

### Then Proceed with Frontend:

**Step 4:** Import DDR components (2h)
**Step 5:** Wire up useEffect with verified request (2h)
**Step 6:** Test on Ayesha's page (1h)

---

## 📝 NOTES TO SELF

### ❌ DO NOT ASSUME:
1. ~~MBD4 triggers DDR_defective~~ - **WRONG!** Not in gene lists
2. Components will just work - Always verify imports
3. API responses match expectations - Test first

### ✅ VERIFIED FACTS:
1. DDR endpoint exists at `/api/resistance/ddr-status`
2. DDR components exist and export correctly
3. MBD4 is NOT in core_hrr_genes or extended_ddr_genes
4. SL engine DOES detect MBD4+TP53 → Can use as fallback
5. DDR score weights: extended_ddr_pathogenic = 1.0

### ⚠️ STILL NEED TO VERIFY:
1. Will adding MBD4 to config work without restart?
2. Does SL result have a field we can use instead of DDR?
3. What's the exact shape of `slResult` for PARP eligibility?

---

**SYNC STATUS:** Aligned with `AYESHA_FRONTEND_PRODUCTION_AUDIT.md`  
**Last Updated:** January 26, 2026
