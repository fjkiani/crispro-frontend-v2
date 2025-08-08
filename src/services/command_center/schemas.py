
from pydantic import BaseModel, Field
from typing import Optional, Dict, Any, List

class PatientDataPacket(BaseModel):
    patient_identifier: str
    gene: str
    mutation_hgvs_c: str
    mutation_hgvs_p: str
    transcript_id: str
    bam_file_path: str
    protein_sequence: str
    sequence_for_perplexity: str

class BattlePlan(BaseModel):
    id: int
    patient_identifier: str
    gene: str
    mutation_hgvs_p: str
    status: str
    ensemble_scores: Optional[Dict[str, Any]] = None
    perplexity_score: Optional[float] = None
    structural_variants: Optional[Dict[str, Any]] = None
    evidence_features: Optional[Dict[str, Any]] = None
    proposed_interventions: Optional[Dict[str, Any]] = None

    class Config:
        from_attributes = True

class BattlePlanResponse(BaseModel):
    plan_id: str
    message: str

class GuideDesignResponse(BaseModel):
    battle_plan_id: int
    message: str
    generated_guides: list

# --- DOCTRINE: SUBDOMAIN HUNTER ---
# Pydantic models for the new Reconnaissance feature.

class SubdomainReconRequest(BaseModel):
    gene_symbol: str

class Domain(BaseModel):
    name: str
    start: int
    end: int
    # --- DOCTRINE: ACTIONABLE INTELLIGENCE ---
    # The data contract must include the priority and the justification for it.
    priority: str
    justification: str
    zeta_threat_assessment: Optional[float] = None # Placeholder for ZetaOracle score

class SubdomainReconResponse(BaseModel):
    gene_symbol: str
    domains: List[Domain]

# --- DOCTRINE: ACTIONABLE INTELLIGENCE ---
# Pydantic models for exposing the Threat Matrix.
class ThreatMatrixResponse(BaseModel):
    keywords: List[str]
    cdd_identifiers: dict[str, str]

# --- DOCTRINE: HUNTER-KILLER ---
# Pydantic models for the new Targeted Strike Plan workflow.

class StrikePlanRequest(BaseModel):
    gene_symbol: str

class TargetedStrikePlan(BaseModel):
    gene_symbol: str
    mission_briefing: str
    prioritized_targets: List[Domain]

# --- DOCTRINE: UI UPLINK ---
# Pydantic models for the new integrated Hunter-Analyst endpoint.

class HunterAnalystRequest(BaseModel):
    gene_symbol: str

class HunterAnalystResponse(BaseModel):
    briefing: str
    primary_target: Optional[Domain] = None
    full_dna_sequence: Optional[str] = None

# --- NEW DOCTRINE: INHIBITOR VALIDATION (THE CONFIRMATION KILL) ---
# THESE SCHEMAS ARE BEING DELETED. THEY BELONG IN THE ZETA ORACLE SERVICE.
# class InhibitorValidationRequest(BaseModel):
#     target_dna_sequence: str = Field(..., description="The short DNA motif of the target.")
#     candidate_inhibitors: list[str] = Field(..., description="List of raw DNA inhibitor candidates from the ZetaForge.")
#
# class ValidatedCandidate(BaseModel):
#     inhibitor_sequence: str
#     stability_score: float
#     commentary: str
#
# class InhibitorValidationResponse(BaseModel):
#     ranked_candidates: list[ValidatedCandidate]

# --- DOCTRINE: TRIUMVIRATE THREAT ASSESSMENT ---

class ThreatAssessmentRequest(BaseModel):
    gene_symbol: str
    protein_change: str # e.g., p.Ala2fs, p.V600E

class Assessment(BaseModel):
    assessment_source: str # e.g., "Truncation Sieve", "Zeta Oracle"
    is_pathogenic: bool
    zeta_score: float
    details: str
    error: Optional[str] = None

class ThreatAssessmentResponse(BaseModel):
    assessment: Assessment 

class JobResponse(BaseModel):
    job_id: str
    status: str
    message: str
    result: List[Dict[str, Any]] | None = None

# --- DOCTRINE: WORKFLOW ORCHESTRATION ---
# Pydantic models for managing long-running workflows and their status.

class WorkflowRequest(BaseModel):
    # The primary identifier for the workflow target.
    # This can be a gene symbol, a specific variant, or a patient ID.
    target_gene_symbol: Optional[str] = None
    # An optional bait sequence for generative tasks.
    bait_sequence: Optional[str] = None

class WorkflowStatus(BaseModel):
    workflow_id: str
    status: str = Field(..., description="Current status of the workflow (e.g., 'pending', 'in_progress', 'completed', 'failed').")
    message: str = Field(..., description="A human-readable message describing the current state or outcome.")
    
class JobStatusResponse(BaseModel):
    job_id: str
    status: str
    message: str
    result: Optional[Dict[str, Any]] = None