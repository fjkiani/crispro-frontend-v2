"""
Treatment History Schema for Treatment Line Integration

Defines the data model for prior therapy tracking and treatment line context.
"""

from typing import List, Optional, Dict, Any
from pydantic import BaseModel, Field


class TreatmentHistory(BaseModel):
    """
    Treatment history for a patient, used to compute line appropriateness and cross-resistance.
    
    Attributes:
        current_line: Current line of therapy (1 = first-line, 2 = second-line, etc.)
        prior_therapies: List of prior drug/regimen names (e.g., ["carboplatin", "paclitaxel"])
        outcomes: Optional list of outcome dicts with keys like pfs_months, response_type (added in P1)
    """
    current_line: int = Field(
        ...,
        ge=1,
        le=10,
        description="Current treatment line (1-10)"
    )
    prior_therapies: List[str] = Field(
        default_factory=list,
        description="List of prior drug/regimen names"
    )
    outcomes: Optional[List[Dict[str, Any]]] = Field(
        default=None,
        description="Optional outcome data for each prior therapy (P1)"
    )


class TreatmentLineProvenance(BaseModel):
    """
    Provenance tracking for treatment line calculations.
    
    Stored in efficacy response provenance.treatment_line
    """
    current_line: int
    prior_therapies: List[str]
    line_appropriateness: float = Field(
        ge=0.0,
        le=1.0,
        description="How appropriate is this drug for the current treatment line"
    )
    cross_resistance_risk: float = Field(
        ge=0.0,
        le=1.0,
        description="Risk of cross-resistance with prior therapies"
    )
    sequencing_fitness: float = Field(
        ge=0.0,
        le=1.0,
        description="Overall sequencing fitness score"
    )
    nccn_category: Optional[str] = Field(
        default=None,
        description="NCCN category (1, 2A, 2B, 3) if available"
    )
    rationale: Optional[str] = Field(
        default=None,
        description="Human-readable explanation of treatment line calculations"
    )

