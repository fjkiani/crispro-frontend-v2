import datetime
from sqlalchemy import (
    Column,
    Integer,
    String,
    DateTime,
    ForeignKey,
    JSON,
    Float
)
from sqlalchemy.orm import relationship
from .database import Base

class Patient(Base):
    __tablename__ = "patients"
    id = Column(Integer, primary_key=True, index=True)
    patient_identifier = Column(String, unique=True, index=True)
    mutations = relationship("Mutation", back_populates="patient")

class Mutation(Base):
    __tablename__ = "mutations"
    id = Column(Integer, primary_key=True, index=True)
    patient_id = Column(Integer, ForeignKey("patients.id"))
    gene = Column(String)
    hgvs_c = Column(String)
    hgvs_p = Column(String)
    
    patient = relationship("Patient", back_populates="mutations")
    battle_plans = relationship("BattlePlan", back_populates="mutation")

class BattlePlan(Base):
    __tablename__ = "battle_plans"
    id = Column(Integer, primary_key=True, index=True)
    mutation_id = Column(Integer, ForeignKey("mutations.id"), nullable=True)  # Make optional for now
    created_at = Column(DateTime, default=datetime.datetime.utcnow)

    # Flattened fields for direct API usage
    patient_identifier = Column(String, index=True)
    gene = Column(String)
    mutation_hgvs_p = Column(String)
    status = Column(String, default="pending_review")

    # Legacy field for general data
    plan_data = Column(JSON)

    # New fields for Guardian Protocol v2
    structural_variants = Column(JSON)
    ensemble_scores = Column(JSON)
    perplexity_score = Column(Float)
    evidence_features = Column(JSON)
    proposed_interventions = Column(JSON) # To store generated guide RNAs

    mutation = relationship("Mutation", back_populates="battle_plans") 