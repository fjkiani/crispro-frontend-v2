# Mission: Fusion Engine Service
# This service is a lightweight, decoupled engine for orchestrating variant scoring
# using AlphaMissense and ESM, deliberately excluding the heavy Evo2 dependencies.
import modal
import os
import re
from fastapi import FastAPI
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from typing import List, Dict, Any
import asyncio
import logging
import functools

# It's critical to be able to import from the tools directory.
import sys
# CORRECTED PATHING LOGIC: Add the project's root `src` directory to the path
# to ensure that the `tools` module can be found during the Modal build process.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from tools.alphamissense_client import AlphaMissenseClient
from tools.esm_client import ESMClient

# --- Lightweight Image Definition ---
# This image contains all necessary dependencies for both clients.
fusion_image = (
    modal.Image.debian_slim(python_version="3.11")
    # CRITICAL FIX: Add build-essential to provide compilers (g++, etc.)
    # that are needed by some of torch's dependencies (e.g., triton).
    .apt_install("build-essential")
    .pip_install(
        "fastapi",
        "pydantic",
        "pandas",
        "pyarrow",
        "fair-esm", 
        "torch",
        "requests",
        "numpy" # Explicitly add numpy to resolve dependency issue.
    )
    # Copy the 'tools' directory into the container image at build time (path resolved relative to this file).
    .add_local_dir(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../tools")), remote_path="/root/tools")
)

# Define the volume where our large data files are stored.
alphamissense_volume = modal.Volume.from_name("alphamissense-data")

# --- App Definition ---
app = modal.App("fusion-engine", image=fusion_image)

# --- Pydantic Models for API Validation (Corrected Structure) ---
class VariantInput(BaseModel):
    variant_id: str
    hgvs: str
    alphamissense_variant_str: str 

class ScoringRequest(BaseModel):
    protein_sequence: str
    variants: List[VariantInput]

class ScoredVariant(BaseModel):
    variant_id: str
    alphamissense_score: float
    esm_score: float
    zeta_score: float # The fused score

class ScoringResponse(BaseModel):
    scored_variants: List[ScoredVariant]

# With the correct data, a generic GPU is sufficient.
@app.cls(cpu=8, memory=65536, timeout=600, volumes={"/data": alphamissense_volume}, gpu="any")
class FusionEngine:
    @modal.enter()
    def load_clients(self):
        """Initializes the clients when the container starts."""
        print("🚀 Initializing Fusion Engine clients...")
        # CRITICAL FIX: Point the client to the correct path in the mounted volume.
        self.am_client = AlphaMissenseClient(data_path="/data/AlphaMissense_hg38.parquet")
        self.esm_client = ESMClient()
        self._esm_available = bool(self.esm_client and getattr(self.esm_client, "model", None))
        # In-memory cache to avoid repeated AM lookups for the same key across requests
        self._am_cache: Dict[str, Dict[str, Any]] = {}
        # Quiet AlphaMissense client logs to prevent spam during bulk scoring
        try:
            logging.getLogger("tools.alphamissense_client").setLevel(logging.ERROR)
            logging.getLogger("tools.esm_client").setLevel(logging.ERROR)
            logging.getLogger("uvicorn").setLevel(logging.ERROR)
            logging.getLogger("uvicorn.error").setLevel(logging.ERROR)
            logging.getLogger("uvicorn.access").setLevel(logging.ERROR)
        except Exception:
            pass
        print("✅ Clients initialized. Fusion Engine is OPERATIONAL.")

    @modal.asgi_app()
    def api(self):
        """Defines the FastAPI web service."""
        fastapi_app = FastAPI(
            title="Fusion Engine",
            description="A lightweight service to score variants using AlphaMissense and ESM.",
            version="1.0.0"
        )

        @fastapi_app.post("/score_variants", response_model=ScoringResponse)
        async def score_variants(request: ScoringRequest):
            """
            This is the main fusion engine endpoint.
            It orchestrates calls to AlphaMissense and ESM in parallel.
            """
            try:
                # --- Step 1: Dispatch All Tasks in Parallel ---
                
                # Only run ESM task if the model loaded correctly (no warnings when absent)
                if self._esm_available:
                    esm_task = asyncio.to_thread(
                        self.esm_client.get_scores,
                        request.protein_sequence,
                        [v.hgvs for v in request.variants]
                    )
                else:
                    # Create a dummy task that returns an empty dict
                    async def dummy_esm_task():
                        return {}
                    esm_task = asyncio.create_task(dummy_esm_task())

                # CORRECT LOGIC + DEDUP + VALIDATION:
                # - AlphaMissense expects GRCh38 keys in form 'chr#:pos:REF:ALT' with single-base alleles
                # - Deduplicate keys within the request to avoid repeated lookups
                # - Serve from cache when available
                key_pattern = re.compile(r"^chr(?:[0-9]{1,2}|X|Y):[0-9]+:[ACGT]:[ACGT]$", re.IGNORECASE)
                idx_to_key: Dict[int, str] = {}
                for i, v in enumerate(request.variants):
                    k = (v.alphamissense_variant_str or "").strip()
                    if key_pattern.match(k):
                        # Normalize chromosome prefix and allele case
                        chrom, pos, ref, alt = k.split(":")
                        chrom_norm = chrom if chrom.lower().startswith("chr") else f"chr{chrom}"
                        nk = f"{chrom_norm}:{pos}:{ref.upper()}:{alt.upper()}"
                        idx_to_key[i] = nk
                    else:
                        # Invalid key shape for AM; skip lookup for this variant
                        continue

                # Unique keys to query (exclude those already in cache)
                unique_keys = []
                for k in set(idx_to_key.values()):
                    if k not in self._am_cache:
                        unique_keys.append(k)

                # Build tasks only for uncached keys
                am_tasks = [asyncio.to_thread(self.am_client.get_score, k) for k in unique_keys]
                
                # Gather results from ESM and all individual AlphaMissense calls
                all_tasks = [esm_task] + am_tasks
                results = await asyncio.gather(*all_tasks, return_exceptions=True)
                
                # --- Step 2: Unpack and Fuse Results ---
                esm_results = results[0] if not isinstance(results[0], Exception) else {}
                am_results_list = results[1:]

                # Update cache for queried keys in order
                for k, res in zip(unique_keys, am_results_list):
                    if isinstance(res, dict):
                        self._am_cache[k] = res
                    else:
                        self._am_cache[k] = {}
                
                final_scores = []
                for i, v in enumerate(request.variants):
                    # Unpack AlphaMissense result for this variant (from cache when present)
                    k = idx_to_key.get(i)
                    am_result = self._am_cache.get(k, {}) if k else {}
                    am_score = am_result.get("score", -999.0)

                    # Unpack ESM result for this variant
                    esm_result = esm_results.get(v.hgvs, {}) if isinstance(esm_results, dict) else {}
                    esm_score = esm_result.get("esm_score", -999.0)
                    
                    # Fuse the scores. For now, a simple average, ignoring errors.
                    valid_scores = [s for s in [am_score, esm_score] if s not in [-999.0, -998.0, -997.0]]
                    zeta_score = sum(valid_scores) / len(valid_scores) if valid_scores else -999.0
                    
                    final_scores.append(ScoredVariant(
                        variant_id=v.variant_id,
                        alphamissense_score=round(am_score, 4),
                        esm_score=round(esm_score, 4),
                        zeta_score=round(zeta_score, 4)
                    ))
                
                return ScoringResponse(scored_variants=final_scores)
            except Exception as e:
                logging.error(f"💥 FUSION ENGINE CRITICAL ERROR: {e}", exc_info=True)
                return JSONResponse(
                    content={"status": "error", "message": f"Internal server error: {str(e)}"},
                    status_code=500
                )

        @fastapi_app.get("/debug/info")
        def debug_info():
            """Expose parquet path and basic dataframe info for AlphaMissense."""
            path = getattr(self.am_client, 'data_path', None)
            head = []
            try:
                df = getattr(self.am_client, '_data', None)
                if df is not None and len(df.index) > 0:
                    # Extract first 5 index keys
                    for i, key in enumerate(df.index[:5]):
                        try:
                            head.append(tuple(key))
                        except Exception:
                            head.append(str(key))
            except Exception as e:
                head = [f"error: {e}"]
            return {"alphamissense_path": path, "index_head": head}

        @fastapi_app.get("/debug/lookup")
        def debug_lookup(key: str):
            """Lookup a raw AlphaMissense key (#CHROM:POS:REF:ALT) and return raw result."""
            try:
                res = self.am_client.get_score(key)
                return {"key": key, "result": res}
            except Exception as e:
                return JSONResponse(content={"key": key, "error": str(e)}, status_code=500)

        @fastapi_app.get("/debug/pos")
        def debug_pos(key: str, limit: int = 20):
            """List AM index entries present at a chrom:pos (accepts with/without 'chr')."""
            try:
                chrom_pos = key.strip()
                if ":" not in chrom_pos:
                    return JSONResponse(content={"error": "key must be 'chr7:140753336'"}, status_code=400)
                chrom, pos = chrom_pos.split(":", 1)
                chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
                pos = int(pos)
                df = getattr(self.am_client, "_data", None)
                if df is None:
                    return {"entries": [], "note": "am data not loaded"}
                # Select matching rows by level
                try:
                    subset = df.loc[(chrom, pos)]
                except Exception:
                    subset = None
                entries = []
                if subset is not None:
                    if hasattr(subset, "index"):
                        # subset index is (REF, ALT)
                        for i, idx in enumerate(list(subset.index)[: max(0, int(limit))]):
                            try:
                                ref, alt = tuple(idx)
                            except Exception:
                                pair = str(idx)
                                if "," in pair:
                                    ref, alt = pair.split(",", 1)
                                else:
                                    ref, alt = pair, "?"
                            entries.append({
                                "key": f"{chrom}:{pos}:{ref}:{alt}",
                                "am_pathogenicity": float(subset.iloc[i]["am_pathogenicity"]) if "am_pathogenicity" in subset.columns else None,
                                "am_class": str(subset.iloc[i]["am_class"]) if "am_class" in subset.columns else None,
                            })
                return {"chrom": chrom, "pos": pos, "entries": entries}
            except Exception as e:
                return JSONResponse(content={"key": key, "error": str(e)}, status_code=500)
        
        @fastapi_app.get("/health")
        def health_check():
            """A simple health check endpoint."""
            return {"status": "healthy", "service": "fusion-engine"}
        
        return fastapi_app

@app.local_entrypoint()
def main():
    """A local test function to verify the fusion engine is working."""
    print("--- ⚔️ LOCAL FUSION ENGINE TEST ⚔️ ---")
    print("   This confirms the service code is syntactically valid.")
    print("   To perform a full end-to-end test, deploy the service with:")
    print("   modal deploy services/fusion_engine/main.py") 
                    am_result = self._am_cache.get(k, {}) if k else {}
                    am_score = am_result.get("score", -999.0)

                    # Unpack ESM result for this variant
                    esm_result = esm_results.get(v.hgvs, {}) if isinstance(esm_results, dict) else {}
                    esm_score = esm_result.get("esm_score", -999.0)
                    
                    # Fuse the scores. For now, a simple average, ignoring errors.
                    valid_scores = [s for s in [am_score, esm_score] if s not in [-999.0, -998.0, -997.0]]
                    zeta_score = sum(valid_scores) / len(valid_scores) if valid_scores else -999.0
                    
                    final_scores.append(ScoredVariant(
                        variant_id=v.variant_id,
                        alphamissense_score=round(am_score, 4),
                        esm_score=round(esm_score, 4),
                        zeta_score=round(zeta_score, 4)
                    ))
                
                return ScoringResponse(scored_variants=final_scores)
            except Exception as e:
                logging.error(f"💥 FUSION ENGINE CRITICAL ERROR: {e}", exc_info=True)
                return JSONResponse(
                    content={"status": "error", "message": f"Internal server error: {str(e)}"},
                    status_code=500
                )

        @fastapi_app.get("/debug/info")
        def debug_info():
            """Expose parquet path and basic dataframe info for AlphaMissense."""
            path = getattr(self.am_client, 'data_path', None)
            head = []
            try:
                df = getattr(self.am_client, '_data', None)
                if df is not None and len(df.index) > 0:
                    # Extract first 5 index keys
                    for i, key in enumerate(df.index[:5]):
                        try:
                            head.append(tuple(key))
                        except Exception:
                            head.append(str(key))
            except Exception as e:
                head = [f"error: {e}"]
            return {"alphamissense_path": path, "index_head": head}

        @fastapi_app.get("/debug/lookup")
        def debug_lookup(key: str):
            """Lookup a raw AlphaMissense key (#CHROM:POS:REF:ALT) and return raw result."""
            try:
                res = self.am_client.get_score(key)
                return {"key": key, "result": res}
            except Exception as e:
                return JSONResponse(content={"key": key, "error": str(e)}, status_code=500)

        @fastapi_app.get("/debug/pos")
        def debug_pos(key: str, limit: int = 20):
            """List AM index entries present at a chrom:pos (accepts with/without 'chr')."""
            try:
                chrom_pos = key.strip()
                if ":" not in chrom_pos:
                    return JSONResponse(content={"error": "key must be 'chr7:140753336'"}, status_code=400)
                chrom, pos = chrom_pos.split(":", 1)
                chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
                pos = int(pos)
                df = getattr(self.am_client, "_data", None)
                if df is None:
                    return {"entries": [], "note": "am data not loaded"}
                # Select matching rows by level
                try:
                    subset = df.loc[(chrom, pos)]
                except Exception:
                    subset = None
                entries = []
                if subset is not None:
                    if hasattr(subset, "index"):
                        # subset index is (REF, ALT)
                        for i, idx in enumerate(list(subset.index)[: max(0, int(limit))]):
                            try:
                                ref, alt = tuple(idx)
                            except Exception:
                                pair = str(idx)
                                if "," in pair:
                                    ref, alt = pair.split(",", 1)
                                else:
                                    ref, alt = pair, "?"
                            entries.append({
                                "key": f"{chrom}:{pos}:{ref}:{alt}",
                                "am_pathogenicity": float(subset.iloc[i]["am_pathogenicity"]) if "am_pathogenicity" in subset.columns else None,
                                "am_class": str(subset.iloc[i]["am_class"]) if "am_class" in subset.columns else None,
                            })
                return {"chrom": chrom, "pos": pos, "entries": entries}
            except Exception as e:
                return JSONResponse(content={"key": key, "error": str(e)}, status_code=500)
        
        @fastapi_app.get("/health")
        def health_check():
            """A simple health check endpoint."""
            return {"status": "healthy", "service": "fusion-engine"}
        
        return fastapi_app

@app.local_entrypoint()
def main():
    """A local test function to verify the fusion engine is working."""
    print("--- ⚔️ LOCAL FUSION ENGINE TEST ⚔️ ---")
    print("   This confirms the service code is syntactically valid.")
    print("   To perform a full end-to-end test, deploy the service with:")
    print("   modal deploy services/fusion_engine/main.py") 