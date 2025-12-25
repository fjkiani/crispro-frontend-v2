# PDGFRA p.S755P VUS - TRUE STATUS

## 🎯 **THE VUS**

| Field | Value |
|-------|-------|
| **Gene** | PDGFRA |
| **Protein** | p.S755P |
| **cDNA** | c.2263T>C |
| **Coordinates** | chr4:54280422 (GRCh38) |
| **ClinVar** | Uncertain Significance (rs780701113) |
| **Status** | **VUS - NOT RESOLVED** |

---

## 📖 **HOW TO RESOLVE VUS (from Evo2 Paper)**

### **Method 1: Zero-Shot Delta Likelihood (Primary)**

From the Evo2 paper, VUS resolution works by:

1. **Fetch 8,192 nt window** centered on variant position
2. **Create sequences**:
   - `ref_sequence`: Original sequence with T at position
   - `alt_sequence`: Mutated sequence with C at position
3. **Score with Evo2**: 
   ```python
   delta_score = model.score_sequences([alt_sequence])[0] - model.score_sequences([ref_sequence])[0]
   ```
4. **Apply threshold** (from BRCA1 calibration):
   - `delta < -0.0009178519` → **Likely Pathogenic**
   - `delta >= -0.0009178519` → **Likely Benign**
5. **Calculate confidence**:
   - LOF std: 0.0015140239
   - FUNC std: 0.0009016589

### **Method 2: Supervised Classifier (For BRCA1/BRCA2)**

For high-stakes genes like BRCA1/BRCA2 (AUROC 0.95):

1. Extract embeddings from Evo2-40B block 20 (pre-normalization layer)
2. Average 128-nt window around variant for ref and alt
3. Concatenate embeddings (16,384 dimensions)
4. MLP classifier: Input → 512 → 128 → 32 → 1 (sigmoid)
5. Output: `prob_pathogenic`

### **Key Insight from Paper**

> "Evo2 achieves 0.94 AUROC on BRCA1 variants without task-specific training"

Negative delta = more deleterious  
The model learns that pathogenic variants make the genome "less likely"

---

## 🏗️ **THE REAL INFRASTRUCTURE**

### **Evo Service (Modal H100:2)**
- **URL**: `https://crispro--evo-service-evoservice-api.modal.run`
- **Endpoints**:
  - `/score_delta` - Raw sequence scoring
  - `/score_variant` - Coordinate-based (fetches from Ensembl)
  - `/score_variant_multi` - Multi-window scoring

### **Genesis Engine (Modal H100)**
- **Code**: `src/services/genesis_engine/main.py`
- **Endpoint**: `/analyze_single_variant`
- **Uses**: Evo2-40B model with ClinVar enrichment

### **Proper API Call**
```bash
curl -X POST https://crispro--evo-service-evoservice-api.modal.run/score_variant \
  -H "Content-Type: application/json" \
  -d '{
    "assembly": "GRCh38",
    "chrom": "4",
    "pos": 54280422,
    "ref": "T",
    "alt": "C",
    "window": 8192
  }'
```

**Expected Response**:
```json
{
  "ref_likelihood": -1234.56,
  "alt_likelihood": -1250.90,
  "delta_score": -16.34,
  "window_start": 54276326,
  "window_end": 54284518,
  "variant_index": 4096
}
```

---

## ❌ **CURRENT BLOCKER**

**Evo Service is returning CUDA errors**:
```
CUDA error: an illegal memory access was encountered
```

This indicates the Modal deployment needs attention:
1. GPU memory issue
2. Model not properly loaded
3. Container restart needed

---

## 🔧 **TO COMPLETE VUS RESOLUTION**

### **Option 1: Fix Evo Service**
```bash
cd src/services/evo_service
modal deploy main.py
```

### **Option 2: Use Genesis Engine**
```bash
cd src/services/genesis_engine
modal deploy main.py
# Then call /analyze_single_variant
```

### **Option 3: Local Evo2 (if GPU available)**
```python
from evo2 import Evo2
model = Evo2('evo2_7b')  # or evo2_40b

# Fetch sequence from Ensembl
import requests
url = "https://rest.ensembl.org/sequence/region/human/4:54276326-54284518:1?content-type=text/plain;coord_system_version=GRCh38"
ref_seq = requests.get(url).text.strip().upper()

# Create alt sequence
idx = 54280422 - 54276326  # Position in window
alt_seq = ref_seq[:idx] + 'C' + ref_seq[idx+1:]

# Score
scores = model.score_sequences([ref_seq, alt_seq])
delta = scores[1] - scores[0]

# Classify
threshold = -0.0009178519
if delta < threshold:
    print(f"Likely Pathogenic (delta={delta:.6f})")
else:
    print(f"Likely Benign (delta={delta:.6f})")
```

---

## 📊 **WHAT WE KNOW SO FAR**

| Step | Status | Data |
|------|--------|------|
| Coordinate Resolution | ✅ | chr4:54280422 T>C |
| ClinVar Lookup | ✅ | Uncertain Significance (rs780701113) |
| Pathway Mapping | ✅ | RTK, MAPK, PI3K |
| Triumvirate Gate | ✅ | PASS (not truncating) |
| **Evo2 Delta Score** | ❌ | **SERVICE DOWN** |
| Classification | ❌ | **CANNOT COMPLETE** |

---

## 🎯 **BOTTOM LINE**

**VUS Status**: PDGFRA p.S755P remains **UNRESOLVED**

**Blocker**: Evo Service CUDA error

**To resolve**: 
1. Redeploy Evo Service: `modal deploy src/services/evo_service/main.py`
2. Call `/score_variant` with chr4:54280422 T>C
3. Get delta_score → classify

**The framework exists. The infrastructure exists. The service just needs to be running.**



## 🎯 **THE VUS**

| Field | Value |
|-------|-------|
| **Gene** | PDGFRA |
| **Protein** | p.S755P |
| **cDNA** | c.2263T>C |
| **Coordinates** | chr4:54280422 (GRCh38) |
| **ClinVar** | Uncertain Significance (rs780701113) |
| **Status** | **VUS - NOT RESOLVED** |

---

## 📖 **HOW TO RESOLVE VUS (from Evo2 Paper)**

### **Method 1: Zero-Shot Delta Likelihood (Primary)**

From the Evo2 paper, VUS resolution works by:

1. **Fetch 8,192 nt window** centered on variant position
2. **Create sequences**:
   - `ref_sequence`: Original sequence with T at position
   - `alt_sequence`: Mutated sequence with C at position
3. **Score with Evo2**: 
   ```python
   delta_score = model.score_sequences([alt_sequence])[0] - model.score_sequences([ref_sequence])[0]
   ```
4. **Apply threshold** (from BRCA1 calibration):
   - `delta < -0.0009178519` → **Likely Pathogenic**
   - `delta >= -0.0009178519` → **Likely Benign**
5. **Calculate confidence**:
   - LOF std: 0.0015140239
   - FUNC std: 0.0009016589

### **Method 2: Supervised Classifier (For BRCA1/BRCA2)**

For high-stakes genes like BRCA1/BRCA2 (AUROC 0.95):

1. Extract embeddings from Evo2-40B block 20 (pre-normalization layer)
2. Average 128-nt window around variant for ref and alt
3. Concatenate embeddings (16,384 dimensions)
4. MLP classifier: Input → 512 → 128 → 32 → 1 (sigmoid)
5. Output: `prob_pathogenic`

### **Key Insight from Paper**

> "Evo2 achieves 0.94 AUROC on BRCA1 variants without task-specific training"

Negative delta = more deleterious  
The model learns that pathogenic variants make the genome "less likely"

---

## 🏗️ **THE REAL INFRASTRUCTURE**

### **Evo Service (Modal H100:2)**
- **URL**: `https://crispro--evo-service-evoservice-api.modal.run`
- **Endpoints**:
  - `/score_delta` - Raw sequence scoring
  - `/score_variant` - Coordinate-based (fetches from Ensembl)
  - `/score_variant_multi` - Multi-window scoring

### **Genesis Engine (Modal H100)**
- **Code**: `src/services/genesis_engine/main.py`
- **Endpoint**: `/analyze_single_variant`
- **Uses**: Evo2-40B model with ClinVar enrichment

### **Proper API Call**
```bash
curl -X POST https://crispro--evo-service-evoservice-api.modal.run/score_variant \
  -H "Content-Type: application/json" \
  -d '{
    "assembly": "GRCh38",
    "chrom": "4",
    "pos": 54280422,
    "ref": "T",
    "alt": "C",
    "window": 8192
  }'
```

**Expected Response**:
```json
{
  "ref_likelihood": -1234.56,
  "alt_likelihood": -1250.90,
  "delta_score": -16.34,
  "window_start": 54276326,
  "window_end": 54284518,
  "variant_index": 4096
}
```

---

## ❌ **CURRENT BLOCKER**

**Evo Service is returning CUDA errors**:
```
CUDA error: an illegal memory access was encountered
```

This indicates the Modal deployment needs attention:
1. GPU memory issue
2. Model not properly loaded
3. Container restart needed

---

## 🔧 **TO COMPLETE VUS RESOLUTION**

### **Option 1: Fix Evo Service**
```bash
cd src/services/evo_service
modal deploy main.py
```

### **Option 2: Use Genesis Engine**
```bash
cd src/services/genesis_engine
modal deploy main.py
# Then call /analyze_single_variant
```

### **Option 3: Local Evo2 (if GPU available)**
```python
from evo2 import Evo2
model = Evo2('evo2_7b')  # or evo2_40b

# Fetch sequence from Ensembl
import requests
url = "https://rest.ensembl.org/sequence/region/human/4:54276326-54284518:1?content-type=text/plain;coord_system_version=GRCh38"
ref_seq = requests.get(url).text.strip().upper()

# Create alt sequence
idx = 54280422 - 54276326  # Position in window
alt_seq = ref_seq[:idx] + 'C' + ref_seq[idx+1:]

# Score
scores = model.score_sequences([ref_seq, alt_seq])
delta = scores[1] - scores[0]

# Classify
threshold = -0.0009178519
if delta < threshold:
    print(f"Likely Pathogenic (delta={delta:.6f})")
else:
    print(f"Likely Benign (delta={delta:.6f})")
```

---

## 📊 **WHAT WE KNOW SO FAR**

| Step | Status | Data |
|------|--------|------|
| Coordinate Resolution | ✅ | chr4:54280422 T>C |
| ClinVar Lookup | ✅ | Uncertain Significance (rs780701113) |
| Pathway Mapping | ✅ | RTK, MAPK, PI3K |
| Triumvirate Gate | ✅ | PASS (not truncating) |
| **Evo2 Delta Score** | ❌ | **SERVICE DOWN** |
| Classification | ❌ | **CANNOT COMPLETE** |

---

## 🎯 **BOTTOM LINE**

**VUS Status**: PDGFRA p.S755P remains **UNRESOLVED**

**Blocker**: Evo Service CUDA error

**To resolve**: 
1. Redeploy Evo Service: `modal deploy src/services/evo_service/main.py`
2. Call `/score_variant` with chr4:54280422 T>C
3. Get delta_score → classify

**The framework exists. The infrastructure exists. The service just needs to be running.**








