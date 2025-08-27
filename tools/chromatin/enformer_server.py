import hashlib
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel


class PredictRequest(BaseModel):
    chrom: str
    start: int
    end: int


app = FastAPI(title="Enformer Proxy (local)", version="0.1.0")


@app.get("/health")
def health():
    return {"ok": True}


@app.post("/predict")
def predict(req: PredictRequest):
    try:
        if req.end <= req.start:
            raise HTTPException(status_code=400, detail="end must be > start")
        # Deterministic pseudo-score in [0,1] based on coordinates as a placeholder
        key = f"{req.chrom}:{req.start}-{req.end}".encode()
        h = hashlib.sha256(key).hexdigest()
        # use first 8 hex digits to derive a float in [0,1]
        v = int(h[:8], 16) / 0xFFFFFFFF
        return {
            "accessibility_score": round(v, 4),
            "provenance": {"model": "enformer-proxy-local", "note": "deterministic placeholder"},
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))



from fastapi import FastAPI, HTTPException
from pydantic import BaseModel


class PredictRequest(BaseModel):
    chrom: str
    start: int
    end: int


app = FastAPI(title="Enformer Proxy (local)", version="0.1.0")


@app.get("/health")
def health():
    return {"ok": True}


@app.post("/predict")
def predict(req: PredictRequest):
    try:
        if req.end <= req.start:
            raise HTTPException(status_code=400, detail="end must be > start")
        # Deterministic pseudo-score in [0,1] based on coordinates as a placeholder
        key = f"{req.chrom}:{req.start}-{req.end}".encode()
        h = hashlib.sha256(key).hexdigest()
        # use first 8 hex digits to derive a float in [0,1]
        v = int(h[:8], 16) / 0xFFFFFFFF
        return {
            "accessibility_score": round(v, 4),
            "provenance": {"model": "enformer-proxy-local", "note": "deterministic placeholder"},
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


