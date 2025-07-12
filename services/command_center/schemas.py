from pydantic import BaseModel, Field
from typing import Optional, Dict, Any

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