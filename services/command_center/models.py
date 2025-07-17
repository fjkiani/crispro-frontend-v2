from sqlalchemy import (
    Column,
    Integer,
    String,
    Float,
    ForeignKey,
    DateTime,
    JSON,
    Text
)
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
from datetime import datetime

Base = declarative_base()

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
    hgvs_p = Column(String, unique=True) # Assuming protein change is unique for simplicity
    
    patient = relationship("Patient", back_populates="mutations")
    battle_plans = relationship("BattlePlan", back_populates="mutation")

class BattlePlan(Base):
    __tablename__ = "battle_plans"
    id = Column(Integer, primary_key=True, index=True)
    mutation_id = Column(Integer, ForeignKey("mutations.id"), nullable=False)
    created_at = Column(DateTime, default=datetime.utcnow)
    status = Column(String, default="pending_design")
    target_sequence = Column(Text, nullable=True)
    results = Column(JSON, default={})
    
    mutation = relationship("Mutation", back_populates="battle_plans")

    def __repr__(self):
        return f"<BattlePlan(id={self.id}, mutation_id='{self.mutation_id}', status='{self.status}')>" 