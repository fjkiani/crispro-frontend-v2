# 📊 ITERATION 4: AI SERVICES & MODEL INTEGRATION

**Status**: ✅ **COMPLETE**  
**Duration**: 3-4 hours  
**Created**: January 14, 2025

---

### **4.1 MODAL DEPLOYMENT ARCHITECTURE**

#### **4.1.1 Core Pattern: HTTP-Based Communication**
**Critical Doctrine**: Backend communicates with Modal services via **HTTP**, NOT Modal SDK
- **Why**: Decouples backend from Modal SDK dependencies
- **Pattern**: Direct HTTP calls to Modal service URLs
- **Example** (`config.py:139-142`):
  ```python
  EVO_SERVICE_URL = "https://crispro--evo-service-evoservice-api.modal.run"
  EVO_URL_1B = "https://crispro--evo-service-evoservice1b-api-1b.modal.run"
  EVO_URL_7B = "https://crispro--evo-service-evoservice7b-api-7b.modal.run"
  EVO_URL_40B = "https://crispro--evo-service-evoservice-api.modal.run"
  ```

#### **4.1.2 Modal Service URL Naming Convention**
**Pattern**: `https://{workspace}--{app-name}-{class-name}-api.modal.run`
- **Workspace**: `crispro`
- **App Name**: Service name (e.g., `evo-service`, `zeta-oracle`)
- **Class Name**: Modal class name (e.g., `evoservice`, `evoservice1b`)
- **Suffix**: `-api` (Modal ASGI app endpoint)

#### **4.1.3 Dynamic Model URL Fallback** (`config.py:146-161`):
```python
def get_model_url(model_id: str) -> str:
    """Get URL for model with fallback logic, evaluated at runtime"""
    url_1b = os.getenv("EVO_URL_1B", "")
    url_7b = os.getenv("EVO_URL_7B", "...")
    url_40b = os.getenv("EVO_URL_40B", "...")
    
    if model_id == "evo2_1b":
        return url_1b or url_7b or url_40b  # Fallback chain
    elif model_id == "evo2_7b":
        return url_7b or url_40b
    elif model_id == "evo2_40b":
        return url_40b
    else:
        return url_1b or url_7b or url_40b
```

**Fallback Strategy**:
- **1B request** → Try 1B, fallback to 7B, then 40B
- **7B request** → Try 7B, fallback to 40B
- **40B request** → Only 40B (no fallback)

**Why**: Cost optimization (smaller models first) + resilience (fallback if unavailable)

---

### **4.2 EVO2 MODAL SERVICE ARCHITECTURE**

#### **4.2.1 Service Structure** (`src/services/evo_service/main.py`):

**Main Service Class** (`EvoService`):
- **GPU**: `H100:2` (2x H100 GPUs)
- **Volumes**: Model cache volume (`evo-model-cache`)
- **Scaledown**: 300 seconds (5 min idle before shutdown)
- **Timeout**: 1800 seconds (30 min max request time)
- **Default Model**: `evo2_1b_base` (can be overridden via `EVO_MODEL_ID`)

**Additional Services**:
- **EvoService7B**: `H100:1`, model `evo2_7b` (if `ENABLE_EVO_7B=1`)
- **EvoService1B**: `H100:1`, model `evo2_1b_base` (dedicated 1B service)

#### **4.2.2 Image Definition** (`evo2_image`):
```python
evo2_image = (
    modal.Image.from_registry("nvidia/cuda:12.4.0-devel-ubuntu22.04", add_python="3.11")
    .apt_install(["build-essential", "cmake", "ninja-build", "libcudnn8", "libcudnn8-dev", "git", "gcc", "g++"])
    .env({"CC": "/usr/bin/gcc", "CXX": "/usr/bin/g++", "PYTORCH_CUDA_ALLOC_CONF": "expandable_segments:True"})
    .run_commands("mkdir -p /tools/llvm/bin", "ln -s /usr/bin/ar /tools/llvm/bin/llvm-ar")
    .run_commands("git clone --recurse-submodules https://github.com/ArcInstitute/evo2.git && cd evo2 && pip install .")
    .run_commands("pip uninstall -y transformer-engine transformer_engine")
    .run_commands("pip install 'transformer_engine[pytorch]==1.13' --no-build-isolation")
    .pip_install("fastapi", "uvicorn[standard]", "loguru", "pydantic", "numpy==1.22.0", "httpx")
)
```

**Key Dependencies**:
- **CUDA 12.4.0**: For GPU acceleration
- **Evo2**: Cloned from ArcInstitute GitHub
- **Transformer Engine 1.13**: Required for Evo2
- **Numpy 1.22.0**: Pinned version (build compatibility)

#### **4.2.3 Endpoints** (`EvoService`):

**1. `/score_delta`** (POST):
- **Input**: `{"ref_sequence": str, "alt_sequence": str}`
- **Output**: `{"ref_likelihood": float, "alt_likelihood": float, "delta_score": float}`
- **Use**: Direct sequence scoring (no Ensembl fetch)

**2. `/score_batch`** (POST):
- **Input**: `{"pairs": [{"ref_sequence": str, "alt_sequence": str}, ...]}`
- **Output**: `{"results": [{"ref_likelihood": float, "alt_likelihood": float, "delta_score": float}, ...]}`
- **Use**: Batch sequence scoring

**3. `/score_variant`** (POST):
- **Input**: `{"assembly": "GRCh38", "chrom": "7", "pos": 140453136, "ref": "A", "alt": "T", "window": 8192}`
- **Output**: Same as `/score_delta` + `{"window_start": int, "window_end": int, "variant_index": int}`
- **Use**: Single-window variant scoring (fetches from Ensembl)

**4. `/score_variant_multi`** (POST):
- **Input**: `{"assembly": "GRCh38", "chrom": "7", "pos": 140453136, "ref": "A", "alt": "T", "windows": [1024, 2048, 4096, 8192]}`
- **Output**: `{"deltas": [{"window": int, "delta": float}, ...], "min_delta": float, "window_used": int}`
- **Use**: Multi-window variant scoring (picks most negative delta)

**5. `/score_variant_exon`** (POST):
- **Input**: `{"assembly": "GRCh38", "chrom": "7", "pos": 140453136, "ref": "A", "alt": "T", "flank": 600}`
- **Output**: `{"exon_delta": float, "window_used": int}`
- **Use**: Tight-window (exon-context) scoring

**6. `/score_variant_profile`** (POST):
- **Input**: `{"assembly": "GRCh38", "chrom": "7", "pos": 140453136, "ref": "A", "alt": "T", "flank": 600, "radius": 100}`
- **Output**: `{"profile": [{"offset": int, "delta": float}, ...], "peak_delta": float, "peak_offset": int}`
- **Use**: Local delta profile (applies ALT at offsets ±100 bp)

**7. `/score_variant_probe`** (POST):
- **Input**: `{"assembly": "GRCh38", "chrom": "7", "pos": 140453136, "ref": "A"}`
- **Output**: `{"probes": [{"alt": str, "delta": float}, ...], "top_alt": str, "top_delta": float}`
- **Use**: 3-alt sensitivity probe (scores ref→A/C/G/T)

**8. `/generate`** (POST):
- **Input**: `{"prompt": str, "n_tokens": int}`
- **Output**: `{"job_id": str}` (async job submission)
- **Use**: Sequence generation (background job)

**9. `/status/{job_id}`** (GET):
- **Output**: `{"status": "pending"|"running"|"complete"|"failed", "result": {...}, "error": str}`
- **Use**: Poll generation job status

#### **4.2.4 Ensembl Integration**:
- **API**: `https://rest.ensembl.org/sequence/region/human/{region}?content-type=text/plain;coord_system_version={asm}`
- **Region Format**: `{chrom}:{start}-{end}:1`
- **Assembly Mapping**: `GRCh38`/`hg38` → `GRCh38`, else → `GRCh37`
- **Validation**: Checks fetched base matches provided REF (allows `N`)

#### **4.2.5 Shared State**:
- **Job Status**: `modal.Dict.from_name("evo-job-status")` (shared across containers)
- **Model Cache**: `modal.Volume.from_name("evo-model-cache")` (persistent model weights)

---

### **4.3 OTHER MODAL SERVICES**

#### **4.3.1 Zeta Oracle** (`src/services/oracle/main.py`):
- **Purpose**: Variant scoring using AlphaMissense + ESM + Evo2 fusion
- **GPU**: `H100:2`
- **Image**: Same `evo2_image` (Evo2 dependencies)
- **Endpoints**:
  - `/score_variants`: Fused scoring (AlphaMissense + ESM + Evo2)
  - `/invoke`: General-purpose variant analysis
- **URL**: `https://crispro--zeta-oracle-v2-zetaoracle-api.modal.run`

#### **4.3.2 Zeta Forge** (`src/services/forge/main.py`):
- **Purpose**: Generate novel biologic inhibitors using Evo2
- **GPU**: `H100:2`
- **Image**: Same `evo2_image`
- **Endpoints**:
  - `/generate_inhibitor`: Submit generation job (async)
  - `/status/{job_id}`: Poll job status
- **Doctrine**: "Ambush" - reverse complement of target motif
- **URL**: `https://crispro--zeta-forge-zetaforge-api.modal.run`

#### **4.3.3 Fusion Engine** (`src/services/fusion_engine/main.py`):
- **Purpose**: Lightweight variant scoring (AlphaMissense + ESM, NO Evo2)
- **GPU**: `any` (generic GPU sufficient)
- **Image**: Lightweight (no Evo2 dependencies)
- **Volume**: `alphamissense-data` (AlphaMissense parquet file)
- **Endpoints**:
  - `/score_variants`: Fused scoring (AlphaMissense + ESM only)
- **Why**: Fast scoring without heavy Evo2 dependencies
- **URL**: `https://crispro--fusion-engine-fusionengine-api.modal.run`

#### **4.3.4 Boltz Service** (`src/services/boltz_service/main.py`):
- **Purpose**: Protein structure prediction using Boltz-2
- **GPU**: `H100`
- **Model**: Boltz-2 from Hugging Face (`boltz-community/boltz-2`)
- **Volume**: `boltz-models` (persistent model weights)
- **Endpoints**:
  - `/v1/predict_structure`: Structure prediction
  - `/v1/predict_interaction`: Binding affinity prediction
- **Features**:
  - **Fast mode**: 2-5 min predictions
  - **iPTM scores**: Binding affinity (>0.7 = high confidence)
  - **SMILES handling**: Direct compound analysis
- **URL**: `https://crispro--boltz-service-boltzservice-api.modal.run`

#### **4.3.5 Genesis Engine** (`src/services/genesis_engine/main.py`):
- **Purpose**: Variant analysis using Evo2 40B
- **GPU**: `H100`
- **Model**: `evo2_40b`
- **Endpoints**:
  - `/analyze_single_variant`: Single variant analysis with enrichment
- **URL**: `https://crispro--genesis-engine-genesisengine-api.modal.run`

---

### **4.4 ALPHAFOLD 3 INTEGRATION STATUS**

#### **4.4.1 Current State**:
- **Status**: ❌ **NOT AUTOMATED** - Manual workflow only
- **Approach**: JSON input → AlphaFold Server → JSON output → Review
- **Why Manual**: Faster than automated deployment, already validated (15/15 guides, 100% pass rate)

#### **4.4.2 Manual Workflow**:
1. **Generate JSON**: Structure description (protein/DNA/RNA chains)
2. **Submit to AF Server**: Web interface or API
3. **Download Results**: pLDDT/PAE metrics, structure files
4. **Review**: Validate quality metrics

#### **4.4.3 AF Server Capabilities**:
- **Molecule Types**: `proteinChain`, `dnaSequence`, `rnaSequence`, `ligand`, `ion`
- **Complexes**: Protein-DNA complexes (e.g., Cas9:gRNA:target DNA)
- **Multiple Chains**: Full CRISPR complexes
- **PTMs**: Post-translational modifications
- **DNA Modifications**: 5mC, 8-oxoG, etc.
- **Template Control**: `useStructureTemplate: false` to disable PDB templates

#### **4.4.4 Boltz as Alternative**:
- **Status**: ✅ **FULLY AUTOMATED** on Modal
- **Speed**: 2-5 min per structure (vs 60+ min for ColabFold)
- **Use Cases**: Fast structure validation, binding affinity prediction
- **Limitation**: Does not model complexes/interfaces (single-chain only)

#### **4.4.5 Future Roadmap**:
- **Phase 1**: ESMFold for fast single-chain predictions (1-2 min/structure)
- **Phase 2**: AlphaFold3 integration (pending Google DeepMind weights approval)
- **Phase 3**: Boltz for binding affinity (NO DiffDock needed)

---

### **4.5 BACKEND → MODAL COMMUNICATION PATTERNS**

#### **4.5.1 Request Flow** (`routers/evo.py`):
```python
# 1. Get model URL with fallback
target_url = get_model_url(model_id)

# 2. Make HTTP request
async with httpx.AsyncClient(timeout=EVO_TIMEOUT) as client:
    response = await client.post(
        f"{target_url}/score_variant_multi",
        json=payload,
        headers={"Content-Type": "application/json"}
    )
    response.raise_for_status()
    return response.json()
```

#### **4.5.2 Fallback Logic** (`routers/evo.py:235-236`):
```python
# If 1B or 7B fails, try falling back to 40B (if allowed)
if (not DISABLE_EVO_FALLBACK) and model_id in ["evo2_1b", "evo2_7b"] and EVO_URL_40B:
    # Retry with 40B
    fallback_request = {**request, "model_id": "evo2_40b"}
    response = await client.post(f"{EVO_URL_40B}/generate", json=fallback_request)
    result["warning"] = f"Requested {model_id} failed, fell back to evo2_40b"
```

#### **4.5.3 Timeout Configuration** (`config.py:143`):
```python
EVO_TIMEOUT = Timeout(60.0, connect=10.0)  # 60s total, 10s connect
```

#### **4.5.4 Health Checks** (`routers/evo.py:383-436`):
- **Endpoint**: `/health` (GET)
- **Purpose**: Check service availability
- **Response**: `{"status": "healthy"|"unhealthy", "url": str}`

---

### **4.6 MODAL DEPLOYMENT PATTERNS**

#### **4.6.1 Image Definition Pattern**:
1. **Base Image**: CUDA image for GPU services, Debian slim for CPU services
2. **System Dependencies**: `apt_install` for build tools, libraries
3. **Python Dependencies**: `pip_install` for Python packages
4. **Local Code**: `add_local_dir` for project code
5. **Environment Variables**: `.env()` for configuration

#### **4.6.2 Service Class Pattern**:
```python
@app.cls(
    gpu="H100:2",           # GPU allocation
    volumes={...},          # Persistent storage
    scaledown_window=300,   # Idle timeout
    timeout=1800            # Max request time
)
class ServiceName:
    @modal.enter()
    def load_model(self):
        # Load model on container startup
        self.model = Model()
    
    @modal.asgi_app()
    def api(self):
        # FastAPI app
        return fastapi_app
```

#### **4.6.3 Shared State Pattern**:
- **Job Status**: `modal.Dict.from_name("service-jobs")` (shared across containers)
- **Model Cache**: `modal.Volume.from_name("service-models")` (persistent storage)

#### **4.6.4 Async Job Pattern**:
```python
# 1. Submit job (returns job_id)
job_id = str(uuid.uuid4())
job_status_dict[job_id] = {"status": "pending"}
background_task.spawn(job_id, request)
return {"job_id": job_id}

# 2. Poll status
@fastapi_app.get("/status/{job_id}")
def get_status(job_id: str):
    return job_status_dict[job_id]
```

---

### **4.7 KEY INSIGHTS**

#### **Modal Service Architecture**:
1. **HTTP-Based**: Backend uses HTTP, not Modal SDK (decoupling)
2. **Fallback Chains**: Model URL fallback (1B → 7B → 40B)
3. **Shared State**: `modal.Dict` for job status, `modal.Volume` for model cache
4. **Scaledown**: 300s idle timeout (cost optimization)

#### **Evo2 Service**:
1. **Multi-Model**: 1B, 7B, 40B variants (separate services)
2. **Ensembl Integration**: Fetches reference sequences from Ensembl API
3. **Multi-Window**: Tests multiple window sizes, picks best delta
4. **Async Jobs**: Generation uses background jobs with status polling

#### **Other Services**:
1. **Oracle**: Fused scoring (AlphaMissense + ESM + Evo2)
2. **Forge**: Sequence generation (inhibitor design)
3. **Fusion**: Lightweight scoring (AlphaMissense + ESM only)
4. **Boltz**: Structure prediction (2-5 min, automated)

#### **AlphaFold 3**:
1. **Manual Workflow**: JSON → AF Server → Review (faster than automated)
2. **Boltz Alternative**: Automated structure prediction (single-chain only)
3. **Future**: ESMFold for fast validation, AF3 for complexes (pending weights)

---

**Status**: 🔄 **ITERATION 5 IN PROGRESS** - Frontend Architecture & User Experience  
**Next**: Complete I5 documentation, then move to I6 (Clinical Systems & Workflows)

---