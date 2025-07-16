from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import uuid
import time

app = FastAPI(title="Mock Evo Service")

# --- In-memory "database" for job statuses ---
mock_jobs = {}

class GuideRNARequest(BaseModel):
    target_sequence: str
    num_guides: int = 5

class JobSubmitResponse(BaseModel):
    job_id: str

@app.post("/generate", status_code=202, response_model=JobSubmitResponse)
async def submit_generation(request: GuideRNARequest):
    """
    Accepts a sequence and simulates submitting a background job.
    """
    job_id = str(uuid.uuid4())
    print(f"[Mock Evo Service] Received generation request. Assigning job ID: {job_id}")
    
    # Simulate the job being processed
    mock_jobs[job_id] = {
        "status": "running",
        "request_time": time.time(),
        "num_guides": request.num_guides
    }
    return {"job_id": job_id}

@app.get("/status/{job_id}")
async def get_status(job_id: str):
    """
    Polls for the status and result of a guide RNA generation job.
    """
    print(f"[Mock Evo Service] Received status query for job ID: {job_id}")
    job = mock_jobs.get(job_id)

    if not job:
        raise HTTPException(status_code=404, detail="Job ID not found.")

    # Simulate completion after a short delay
    if time.time() - job["request_time"] > 2:
        if job["status"] != "complete":
             # Generate plausible fake guides
            fake_guides = [
                f"FAKEGUIDE{i+1}ACGTACGTACGT" for i in range(job["num_guides"])
            ]
            job["status"] = "complete"
            job["result"] = {
                "guides": [
                    {
                        "guide_sequence": g,
                        "assassin_score": 0.85 + (i * 0.01), # some variance
                        "on_target_score": 0.9,
                        "off_target_score": 0.95,
                        "immunogenicity_score": 0.1
                    } for i, g in enumerate(fake_guides)
                ]
            }
            print(f"[Mock Evo Service] Job {job_id} is now complete.")

    return job

@app.get("/")
def read_root():
    return {"message": "Mock Evo Service is running"}

# To run this mock service:
# uvicorn tests.mocks.mock_evo_service:app --port 8001 