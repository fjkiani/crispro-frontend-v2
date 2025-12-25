# VUS Coordinate Resolution - Complete Audit & Solution

## Executive Summary

**Problem**: PDGFRA p.S755P VUS resolution workflow is blocked because GRCh38 coordinates are required but not automatically resolved from HGVS protein notation.

**Root Cause**: Coordinate resolution capability exists but is **NOT integrated** into the VUS workflow.

**Solution**: Integrate existing coordinate resolution into VUS workflow via new backend endpoint and frontend integration.

---

## 🔍 Audit Results

### ✅ **EXISTING CAPABILITIES FOUND**

#### 1. **`src/tools/threat_assessor.py`** - VEP Annotation Method
**Location**: `src/tools/threat_assessor.py:214-282`

**Capability**:
```python
def _get_vep_annotation(self):
    """
    Gets variant annotation from Ensembl VEP using a two-step process:
    1. Find canonical transcript for the gene.
    2. Use transcript to construct a full HGVS notation for the VEP query.
    """
    # Step 1: Get canonical transcript
    lookup_server = "https://rest.ensembl.org"
    lookup_ext = f"/lookup/symbol/homo_sapiens/{self.gene}"
    # ... gets canonical_transcript
    
    # Step 2: Query VEP with transcript:HGVS
    self.hgvs_notation = f"{transcript_id}:{self.protein_change}"
    vep_ext = "/vep/human/hgvs"
    # ... returns genomic coordinates
```

**What it does**:
- ✅ Gets canonical transcript for gene
- ✅ Constructs transcript-specific HGVS notation
- ✅ Queries Ensembl VEP API
- ✅ Returns genomic coordinates (chrom, pos, ref, alt)

**Status**: ✅ **WORKING** but **NOT EXPOSED** as API endpoint

---

#### 2. **`src/tools/coordinate_handler.py`** - Coordinate Utilities
**Location**: `src/tools/coordinate_handler.py`

**Capabilities**:
- `GeneCoordinateConverter` - Converts gene symbols to coordinates
- `_fetch_gene_coordinates()` - Fetches from Ensembl API
- Coordinate format conversions (UCSC, Ensembl, IGV)

**Status**: ✅ **AVAILABLE** but **NOT INTEGRATED** into VUS workflow

---

#### 3. **Backend Endpoints - Partial Support**

**`/api/safety/ensembl_context`**:
- ✅ Takes coordinates (chrom, pos, ref, alt)
- ❌ Does NOT resolve from HGVS
- ❌ Requires coordinates as input

**`/api/evidence/deep_analysis`**:
- ✅ Can accept gene + hgvs_p
- ❌ Still requires coordinates for full functionality
- ⚠️ Partial support only

**Status**: ⚠️ **PARTIAL** - endpoints exist but don't resolve coordinates

---

### ❌ **MISSING INTEGRATIONS**

#### 1. **No Backend Endpoint for HGVS→Coordinates**
**Gap**: No `/api/vus/resolve_coordinates` or similar endpoint

**Impact**: Frontend cannot resolve coordinates automatically

**Required**:
```python
POST /api/vus/resolve_coordinates
{
  "gene": "PDGFRA",
  "hgvs_p": "S755P",
  "c_dna": "c.2263T>C"  # optional
}

Response:
{
  "chrom": "4",
  "pos": 55152000,  # example
  "ref": "T",
  "alt": "C",
  "transcript": "ENST00000257290",
  "provenance": {...}
}
```

---

#### 2. **VUS Components Don't Auto-Resolve**
**Gap**: `AnalysisResults.jsx` expects coordinates but doesn't resolve them

**Current Code** (`AnalysisResults.jsx:66-69`):
```javascript
const chrom = activeMutation?.chrom;
const pos = activeMutation?.pos;
const ref = activeMutation?.ref;
const alt = activeMutation?.alt;
```

**Problem**: If `activeMutation` only has `gene` and `hgvs_p`, coordinates are `undefined`

**Impact**: 
- ❌ Insights bundle fails (requires coordinates)
- ❌ ClinVar lookup fails (requires coordinates)
- ❌ Fusion coverage check fails (requires coordinates)
- ❌ Regulatory/chromatin analysis fails (requires coordinates)

---

#### 3. **VUS Master Doctrine Doesn't Specify Resolution**
**Gap**: `.cursor/rules/specialized_systems/vus_master.mdc` Step 1 mentions "Input Formats: HGVS, VCF, coordinates" but doesn't specify HOW to convert HGVS to coordinates.

**Missing**: 
- Step 1.5: "Coordinate Resolution (if HGVS provided)"
- Integration point for `threat_assessor.py` logic
- Frontend hook for coordinate resolution

---

## 🎯 SOLUTION ARCHITECTURE

### **Phase 1: Backend Endpoint** (Priority: HIGH)

**Create**: `oncology-coPilot/oncology-backend-minimal/api/routers/vus.py`

**Endpoint**: `POST /api/vus/resolve_coordinates`

**Implementation**:
```python
@router.post("/resolve_coordinates")
async def resolve_coordinates(request: Dict[str, Any]):
    """
    Resolve GRCh38 coordinates from HGVS protein notation.
    
    Input: { gene, hgvs_p, c_dna? }
    Output: { chrom, pos, ref, alt, transcript, provenance }
    """
    # Use threat_assessor.py logic:
    # 1. Get canonical transcript
    # 2. Construct transcript:HGVS notation
    # 3. Query VEP API
    # 4. Extract coordinates
    # 5. Return normalized GRCh38 coordinates
```

**Dependencies**:
- Adapt `threat_assessor.py` logic
- Use Ensembl VEP API
- Handle errors gracefully

---

### **Phase 2: Frontend Hook** (Priority: HIGH)

**Create**: `oncology-coPilot/oncology-frontend/src/hooks/useCoordinateResolution.js`

**Implementation**:
```javascript
export function useCoordinateResolution({ gene, hgvs_p, c_dna }) {
  // Calls /api/vus/resolve_coordinates
  // Returns { chrom, pos, ref, alt, loading, error }
  // Caches results by {gene, hgvs_p} key
}
```

**Integration**: Update `AnalysisResults.jsx` to use hook when coordinates missing

---

### **Phase 3: VUS Workflow Integration** (Priority: MEDIUM)

**Update**: `.cursor/rules/specialized_systems/vus_master.mdc`

**Add Step 1.5**:
```markdown
#### **Step 1.5: Coordinate Resolution (if HGVS provided)**
- **Input**: Gene + HGVS protein notation
- **Process**: Call `/api/vus/resolve_coordinates`
- **Output**: GRCh38 coordinates (chrom, pos, ref, alt)
- **Fallback**: If resolution fails, continue with gene-level analysis only
```

**Update**: `resolve_pdgra_vus.py` to use new endpoint

---

## 📋 IMPLEMENTATION CHECKLIST

### Backend
- [ ] Create `api/routers/vus.py`
- [ ] Implement `resolve_coordinates` endpoint
- [ ] Adapt `threat_assessor.py` logic
- [ ] Add error handling
- [ ] Add caching (optional)
- [ ] Register router in `main.py`
- [ ] Test with PDGFRA p.S755P

### Frontend
- [ ] Create `hooks/useCoordinateResolution.js`
- [ ] Update `AnalysisResults.jsx` to use hook
- [ ] Add loading states
- [ ] Add error handling
- [ ] Add retry logic
- [ ] Test with PDGFRA p.S755P

### Documentation
- [ ] Update `vus_master.mdc` with Step 1.5
- [ ] Update `PDGFRA_VUS_RESOLUTION_PLAN.md`
- [ ] Add API documentation
- [ ] Add usage examples

### Testing
- [ ] Test coordinate resolution for PDGFRA p.S755P
- [ ] Test error cases (invalid gene, invalid HGVS)
- [ ] Test caching behavior
- [ ] Test full VUS workflow with resolved coordinates

---

## 🚀 QUICK START SOLUTION

### Immediate Fix for PDGFRA p.S755P

**Option 1: Use Existing Script** (Manual)
```bash
# Adapt threat_assessor.py logic
python3 -c "
from src.tools.threat_assessor import ThreatAssessor
ta = ThreatAssessor('PDGFRA', 'p.S755P')
ta._get_vep_annotation()
print(ta.vep_data)
"
```

**Option 2: Create Standalone Script** (Automated)
```python
# resolve_pdgra_coordinates.py
import requests

gene = "PDGFRA"
hgvs_p = "p.S755P"

# Step 1: Get canonical transcript
lookup_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}"
r = requests.get(lookup_url)
transcript = r.json()['canonical_transcript']

# Step 2: Query VEP
vep_url = "https://rest.ensembl.org/vep/human/hgvs"
vep_data = {"hgvs_notations": [f"{transcript}:{hgvs_p}"]}
r = requests.post(vep_url, json=vep_data)
coords = r.json()[0]

print(f"Chrom: {coords['seq_region_name']}")
print(f"Pos: {coords['start']}")
print(f"Ref: {coords.get('allele_string', '').split('/')[0]}")
print(f"Alt: {coords.get('allele_string', '').split('/')[1]}")
```

---

## 📊 GAP ANALYSIS SUMMARY

| Component | Status | Integration | Priority |
|-----------|--------|-------------|----------|
| `threat_assessor.py` | ✅ Working | ❌ Not exposed | HIGH |
| `coordinate_handler.py` | ✅ Available | ❌ Not used | MEDIUM |
| Backend endpoint | ❌ Missing | N/A | HIGH |
| Frontend hook | ❌ Missing | N/A | HIGH |
| VUS workflow | ⚠️ Partial | ❌ No auto-resolve | HIGH |
| Documentation | ⚠️ Incomplete | N/A | MEDIUM |

---

## 🎯 RECOMMENDED ACTION PLAN

### **Immediate (Today)**
1. ✅ Create backend endpoint `/api/vus/resolve_coordinates`
2. ✅ Create frontend hook `useCoordinateResolution`
3. ✅ Update `AnalysisResults.jsx` to auto-resolve coordinates
4. ✅ Test with PDGFRA p.S755P

### **Short-term (This Week)**
1. Update `vus_master.mdc` with Step 1.5
2. Update `resolve_pdgra_vus.py` to use new endpoint
3. Add error handling and retry logic
4. Add caching for resolved coordinates

### **Long-term (Next Sprint)**
1. Add coordinate resolution to `GenomicQueryPanel.jsx`
2. Add batch coordinate resolution
3. Add coordinate validation
4. Add coordinate normalization

---

## ✅ CONCLUSION

**We HAVE the capability** (`threat_assessor.py`) but it's **NOT INTEGRATED** into the VUS workflow.

**Solution**: Create backend endpoint + frontend hook + integrate into VUS components.

**Timeline**: Can be implemented today (2-3 hours) to unblock PDGFRA VUS resolution.

**Next Step**: Implement `/api/vus/resolve_coordinates` endpoint.




## Executive Summary

**Problem**: PDGFRA p.S755P VUS resolution workflow is blocked because GRCh38 coordinates are required but not automatically resolved from HGVS protein notation.

**Root Cause**: Coordinate resolution capability exists but is **NOT integrated** into the VUS workflow.

**Solution**: Integrate existing coordinate resolution into VUS workflow via new backend endpoint and frontend integration.

---

## 🔍 Audit Results

### ✅ **EXISTING CAPABILITIES FOUND**

#### 1. **`src/tools/threat_assessor.py`** - VEP Annotation Method
**Location**: `src/tools/threat_assessor.py:214-282`

**Capability**:
```python
def _get_vep_annotation(self):
    """
    Gets variant annotation from Ensembl VEP using a two-step process:
    1. Find canonical transcript for the gene.
    2. Use transcript to construct a full HGVS notation for the VEP query.
    """
    # Step 1: Get canonical transcript
    lookup_server = "https://rest.ensembl.org"
    lookup_ext = f"/lookup/symbol/homo_sapiens/{self.gene}"
    # ... gets canonical_transcript
    
    # Step 2: Query VEP with transcript:HGVS
    self.hgvs_notation = f"{transcript_id}:{self.protein_change}"
    vep_ext = "/vep/human/hgvs"
    # ... returns genomic coordinates
```

**What it does**:
- ✅ Gets canonical transcript for gene
- ✅ Constructs transcript-specific HGVS notation
- ✅ Queries Ensembl VEP API
- ✅ Returns genomic coordinates (chrom, pos, ref, alt)

**Status**: ✅ **WORKING** but **NOT EXPOSED** as API endpoint

---

#### 2. **`src/tools/coordinate_handler.py`** - Coordinate Utilities
**Location**: `src/tools/coordinate_handler.py`

**Capabilities**:
- `GeneCoordinateConverter` - Converts gene symbols to coordinates
- `_fetch_gene_coordinates()` - Fetches from Ensembl API
- Coordinate format conversions (UCSC, Ensembl, IGV)

**Status**: ✅ **AVAILABLE** but **NOT INTEGRATED** into VUS workflow

---

#### 3. **Backend Endpoints - Partial Support**

**`/api/safety/ensembl_context`**:
- ✅ Takes coordinates (chrom, pos, ref, alt)
- ❌ Does NOT resolve from HGVS
- ❌ Requires coordinates as input

**`/api/evidence/deep_analysis`**:
- ✅ Can accept gene + hgvs_p
- ❌ Still requires coordinates for full functionality
- ⚠️ Partial support only

**Status**: ⚠️ **PARTIAL** - endpoints exist but don't resolve coordinates

---

### ❌ **MISSING INTEGRATIONS**

#### 1. **No Backend Endpoint for HGVS→Coordinates**
**Gap**: No `/api/vus/resolve_coordinates` or similar endpoint

**Impact**: Frontend cannot resolve coordinates automatically

**Required**:
```python
POST /api/vus/resolve_coordinates
{
  "gene": "PDGFRA",
  "hgvs_p": "S755P",
  "c_dna": "c.2263T>C"  # optional
}

Response:
{
  "chrom": "4",
  "pos": 55152000,  # example
  "ref": "T",
  "alt": "C",
  "transcript": "ENST00000257290",
  "provenance": {...}
}
```

---

#### 2. **VUS Components Don't Auto-Resolve**
**Gap**: `AnalysisResults.jsx` expects coordinates but doesn't resolve them

**Current Code** (`AnalysisResults.jsx:66-69`):
```javascript
const chrom = activeMutation?.chrom;
const pos = activeMutation?.pos;
const ref = activeMutation?.ref;
const alt = activeMutation?.alt;
```

**Problem**: If `activeMutation` only has `gene` and `hgvs_p`, coordinates are `undefined`

**Impact**: 
- ❌ Insights bundle fails (requires coordinates)
- ❌ ClinVar lookup fails (requires coordinates)
- ❌ Fusion coverage check fails (requires coordinates)
- ❌ Regulatory/chromatin analysis fails (requires coordinates)

---

#### 3. **VUS Master Doctrine Doesn't Specify Resolution**
**Gap**: `.cursor/rules/specialized_systems/vus_master.mdc` Step 1 mentions "Input Formats: HGVS, VCF, coordinates" but doesn't specify HOW to convert HGVS to coordinates.

**Missing**: 
- Step 1.5: "Coordinate Resolution (if HGVS provided)"
- Integration point for `threat_assessor.py` logic
- Frontend hook for coordinate resolution

---

## 🎯 SOLUTION ARCHITECTURE

### **Phase 1: Backend Endpoint** (Priority: HIGH)

**Create**: `oncology-coPilot/oncology-backend-minimal/api/routers/vus.py`

**Endpoint**: `POST /api/vus/resolve_coordinates`

**Implementation**:
```python
@router.post("/resolve_coordinates")
async def resolve_coordinates(request: Dict[str, Any]):
    """
    Resolve GRCh38 coordinates from HGVS protein notation.
    
    Input: { gene, hgvs_p, c_dna? }
    Output: { chrom, pos, ref, alt, transcript, provenance }
    """
    # Use threat_assessor.py logic:
    # 1. Get canonical transcript
    # 2. Construct transcript:HGVS notation
    # 3. Query VEP API
    # 4. Extract coordinates
    # 5. Return normalized GRCh38 coordinates
```

**Dependencies**:
- Adapt `threat_assessor.py` logic
- Use Ensembl VEP API
- Handle errors gracefully

---

### **Phase 2: Frontend Hook** (Priority: HIGH)

**Create**: `oncology-coPilot/oncology-frontend/src/hooks/useCoordinateResolution.js`

**Implementation**:
```javascript
export function useCoordinateResolution({ gene, hgvs_p, c_dna }) {
  // Calls /api/vus/resolve_coordinates
  // Returns { chrom, pos, ref, alt, loading, error }
  // Caches results by {gene, hgvs_p} key
}
```

**Integration**: Update `AnalysisResults.jsx` to use hook when coordinates missing

---

### **Phase 3: VUS Workflow Integration** (Priority: MEDIUM)

**Update**: `.cursor/rules/specialized_systems/vus_master.mdc`

**Add Step 1.5**:
```markdown
#### **Step 1.5: Coordinate Resolution (if HGVS provided)**
- **Input**: Gene + HGVS protein notation
- **Process**: Call `/api/vus/resolve_coordinates`
- **Output**: GRCh38 coordinates (chrom, pos, ref, alt)
- **Fallback**: If resolution fails, continue with gene-level analysis only
```

**Update**: `resolve_pdgra_vus.py` to use new endpoint

---

## 📋 IMPLEMENTATION CHECKLIST

### Backend
- [ ] Create `api/routers/vus.py`
- [ ] Implement `resolve_coordinates` endpoint
- [ ] Adapt `threat_assessor.py` logic
- [ ] Add error handling
- [ ] Add caching (optional)
- [ ] Register router in `main.py`
- [ ] Test with PDGFRA p.S755P

### Frontend
- [ ] Create `hooks/useCoordinateResolution.js`
- [ ] Update `AnalysisResults.jsx` to use hook
- [ ] Add loading states
- [ ] Add error handling
- [ ] Add retry logic
- [ ] Test with PDGFRA p.S755P

### Documentation
- [ ] Update `vus_master.mdc` with Step 1.5
- [ ] Update `PDGFRA_VUS_RESOLUTION_PLAN.md`
- [ ] Add API documentation
- [ ] Add usage examples

### Testing
- [ ] Test coordinate resolution for PDGFRA p.S755P
- [ ] Test error cases (invalid gene, invalid HGVS)
- [ ] Test caching behavior
- [ ] Test full VUS workflow with resolved coordinates

---

## 🚀 QUICK START SOLUTION

### Immediate Fix for PDGFRA p.S755P

**Option 1: Use Existing Script** (Manual)
```bash
# Adapt threat_assessor.py logic
python3 -c "
from src.tools.threat_assessor import ThreatAssessor
ta = ThreatAssessor('PDGFRA', 'p.S755P')
ta._get_vep_annotation()
print(ta.vep_data)
"
```

**Option 2: Create Standalone Script** (Automated)
```python
# resolve_pdgra_coordinates.py
import requests

gene = "PDGFRA"
hgvs_p = "p.S755P"

# Step 1: Get canonical transcript
lookup_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}"
r = requests.get(lookup_url)
transcript = r.json()['canonical_transcript']

# Step 2: Query VEP
vep_url = "https://rest.ensembl.org/vep/human/hgvs"
vep_data = {"hgvs_notations": [f"{transcript}:{hgvs_p}"]}
r = requests.post(vep_url, json=vep_data)
coords = r.json()[0]

print(f"Chrom: {coords['seq_region_name']}")
print(f"Pos: {coords['start']}")
print(f"Ref: {coords.get('allele_string', '').split('/')[0]}")
print(f"Alt: {coords.get('allele_string', '').split('/')[1]}")
```

---

## 📊 GAP ANALYSIS SUMMARY

| Component | Status | Integration | Priority |
|-----------|--------|-------------|----------|
| `threat_assessor.py` | ✅ Working | ❌ Not exposed | HIGH |
| `coordinate_handler.py` | ✅ Available | ❌ Not used | MEDIUM |
| Backend endpoint | ❌ Missing | N/A | HIGH |
| Frontend hook | ❌ Missing | N/A | HIGH |
| VUS workflow | ⚠️ Partial | ❌ No auto-resolve | HIGH |
| Documentation | ⚠️ Incomplete | N/A | MEDIUM |

---

## 🎯 RECOMMENDED ACTION PLAN

### **Immediate (Today)**
1. ✅ Create backend endpoint `/api/vus/resolve_coordinates`
2. ✅ Create frontend hook `useCoordinateResolution`
3. ✅ Update `AnalysisResults.jsx` to auto-resolve coordinates
4. ✅ Test with PDGFRA p.S755P

### **Short-term (This Week)**
1. Update `vus_master.mdc` with Step 1.5
2. Update `resolve_pdgra_vus.py` to use new endpoint
3. Add error handling and retry logic
4. Add caching for resolved coordinates

### **Long-term (Next Sprint)**
1. Add coordinate resolution to `GenomicQueryPanel.jsx`
2. Add batch coordinate resolution
3. Add coordinate validation
4. Add coordinate normalization

---

## ✅ CONCLUSION

**We HAVE the capability** (`threat_assessor.py`) but it's **NOT INTEGRATED** into the VUS workflow.

**Solution**: Create backend endpoint + frontend hook + integrate into VUS components.

**Timeline**: Can be implemented today (2-3 hours) to unblock PDGFRA VUS resolution.

**Next Step**: Implement `/api/vus/resolve_coordinates` endpoint.









