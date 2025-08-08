from pydantic import BaseModel, Field
from typing import List, Optional

class RankedCandidate(BaseModel):
    """
    Represents a single candidate nanobody scored by Boltz-2.
    """
    sequence: str
    binding_affinity: float
    commentary: str

class InteractionRequest(BaseModel):
    """
    The request model for the Boltz-2 interaction service.
    It takes a single target sequence and a list of candidate nanobodies to screen against it.
    """
    target_sequence: str
    candidate_sequences: List[str]
    job_id: str

class InteractionResponse(BaseModel):
    """
    The response model for the Bol-2 interaction service.
    It returns a ranked list of the candidates by their predicted binding affinity.
    """
    job_id: str
    ranked_candidates: List[RankedCandidate]
    error: Optional[str] = None

# --- New Models for Structural Integrity (pLDDT) ---

class StructuralIntegrityRequest(BaseModel):
    """
    The request model for the new single-protein structural integrity check.
    """
    protein_sequence: str = Field(..., description="The single protein sequence to be validated.")
    job_id: str = Field(..., description="A unique identifier for this validation job.")

class StructuralIntegrityResponse(BaseModel):
    """
    The response model for the structural integrity check.
    Returns the pLDDT score.
    """
    job_id: str
    plddt_score: Optional[float] = None
    ptm_score: Optional[float] = None
    fraction_disordered: Optional[float] = None
    status: str
    error: Optional[str] = None

# --- New Models for Structural Integrity (pLDDT) ---