"""Pydantic models for Evo2 service API."""
from pydantic import BaseModel, Field


class JobSubmitResponse(BaseModel):
    """Response model for job submission."""
    job_id: str


class GenerationRequest(BaseModel):
    """Request model for sequence generation."""
    prompt: str = Field(..., description="A detailed natural language prompt for the generative model.")
    n_tokens: int = Field(100, description="The number of tokens to generate.")
