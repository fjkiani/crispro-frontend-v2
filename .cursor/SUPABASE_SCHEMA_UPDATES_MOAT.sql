-- ============================================================================
-- SUPABASE SCHEMA UPDATES - MOAT ORCHESTRATOR & CAPABILITIES
-- ============================================================================
-- Purpose: Add tables for MOAT orchestrator pipeline, SAE features, resistance playbook
-- Based on: Frontend MOAT Capabilities Audit (.cursor/MOAT/FRONTEND_MOAT_CAPABILITIES_AUDIT.md)
-- Run this AFTER running SUPABASE_SCHEMA_UPDATES.sql
-- ============================================================================

-- ============================================================================
-- 1. ORCHESTRATOR PIPELINE STATE TABLE
-- ============================================================================
-- Tracks the 7-phase orchestrator pipeline execution state
CREATE TABLE IF NOT EXISTS orchestrator_pipeline_state (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES user_profiles(id) ON DELETE CASCADE,
    patient_id TEXT NOT NULL REFERENCES patient_profiles(patient_id) ON DELETE CASCADE,
    
    -- Pipeline phase tracking
    current_phase TEXT NOT NULL DEFAULT 'initialized' CHECK (current_phase IN (
        'initialized', 'extracting', 'analyzing', 'ranking', 
        'matching', 'planning', 'monitoring', 'complete', 'error'
    )),
    progress_percentage FLOAT DEFAULT 0.0 CHECK (progress_percentage >= 0.0 AND progress_percentage <= 100.0),
    
    -- Disease context
    disease TEXT, -- "ovarian", "myeloma", etc.
    
    -- Agent outputs (stored as JSONB for flexibility)
    patient_profile JSONB, -- From 01_DATA_EXTRACTION
    biomarker_profile JSONB, -- From 02_BIOMARKER
    resistance_prediction JSONB, -- From 03_RESISTANCE
    drug_ranking JSONB, -- From 04_DRUG_EFFICACY (S/P/E)
    synthetic_lethality_result JSONB, -- From 14_SYNTHETIC_LETHALITY
    trial_matches JSONB, -- From 05_TRIAL_MATCHING
    nutrition_plan JSONB, -- From 06_NUTRITION
    care_plan JSONB, -- From 07_CARE_PLAN
    monitoring_config JSONB, -- From 08_MONITORING
    
    -- Derived data
    mechanism_vector JSONB, -- 7D vector: [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
    mutations JSONB, -- Array of mutation objects
    
    -- Quality & alerts
    data_quality_flags TEXT[],
    alerts JSONB, -- Array of alert objects
    
    -- Status tracking
    is_complete BOOLEAN DEFAULT FALSE,
    has_errors BOOLEAN DEFAULT FALSE,
    error_message TEXT,
    
    -- Classification
    data_classification TEXT DEFAULT 'PHI' CHECK (data_classification IN ('PHI', 'NON_PHI', 'SENSITIVE', 'PUBLIC')),
    
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    completed_at TIMESTAMPTZ,
    
    -- Unique constraint: one pipeline state per patient
    UNIQUE(user_id, patient_id)
);

CREATE INDEX IF NOT EXISTS idx_orchestrator_state_user_id ON orchestrator_pipeline_state(user_id);
CREATE INDEX IF NOT EXISTS idx_orchestrator_state_patient_id ON orchestrator_pipeline_state(patient_id);
CREATE INDEX IF NOT EXISTS idx_orchestrator_state_phase ON orchestrator_pipeline_state(current_phase);
CREATE INDEX IF NOT EXISTS idx_orchestrator_state_complete ON orchestrator_pipeline_state(is_complete) WHERE is_complete = TRUE;
CREATE INDEX IF NOT EXISTS idx_orchestrator_state_data_classification ON orchestrator_pipeline_state(data_classification);

-- ============================================================================
-- 2. AGENT EXECUTION HISTORY TABLE
-- ============================================================================
-- Tracks individual agent executions within the pipeline
CREATE TABLE IF NOT EXISTS agent_executions (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES user_profiles(id) ON DELETE CASCADE,
    patient_id TEXT NOT NULL REFERENCES patient_profiles(patient_id) ON DELETE CASCADE,
    pipeline_state_id UUID REFERENCES orchestrator_pipeline_state(id) ON DELETE CASCADE,
    
    -- Agent identification
    agent_id TEXT NOT NULL, -- "01_DATA_EXTRACTION", "02_BIOMARKER", etc.
    agent_name TEXT NOT NULL, -- Human-readable name
    
    -- Execution status
    status TEXT NOT NULL DEFAULT 'pending' CHECK (status IN ('pending', 'running', 'completed', 'failed', 'skipped')),
    phase TEXT NOT NULL, -- Which pipeline phase this agent belongs to
    
    -- Timing
    started_at TIMESTAMPTZ,
    completed_at TIMESTAMPTZ,
    duration_seconds FLOAT,
    
    -- Results
    output_data JSONB, -- Agent output (stored in pipeline_state, but also here for history)
    error_message TEXT,
    error_traceback TEXT,
    
    -- Metadata
    input_data_snapshot JSONB, -- Snapshot of inputs used
    dependencies TEXT[], -- Which agents this depends on
    
    -- Classification
    data_classification TEXT DEFAULT 'PHI' CHECK (data_classification IN ('PHI', 'NON_PHI', 'SENSITIVE', 'PUBLIC')),
    
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_agent_executions_user_id ON agent_executions(user_id);
CREATE INDEX IF NOT EXISTS idx_agent_executions_patient_id ON agent_executions(patient_id);
CREATE INDEX IF NOT EXISTS idx_agent_executions_pipeline_state_id ON agent_executions(pipeline_state_id);
CREATE INDEX IF NOT EXISTS idx_agent_executions_agent_id ON agent_executions(agent_id);
CREATE INDEX IF NOT EXISTS idx_agent_executions_status ON agent_executions(status);
CREATE INDEX IF NOT EXISTS idx_agent_executions_phase ON agent_executions(phase);
CREATE INDEX IF NOT EXISTS idx_agent_executions_data_classification ON agent_executions(data_classification);

-- ============================================================================
-- 3. SAE FEATURES TABLE
-- ============================================================================
-- Stores SAE (Sparse Autoencoder) features and diagnostics
-- Supports both TRUE SAE and PROXY SAE
CREATE TABLE IF NOT EXISTS sae_features (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES user_profiles(id) ON DELETE CASCADE,
    patient_id TEXT NOT NULL REFERENCES patient_profiles(patient_id) ON DELETE CASCADE,
    pipeline_state_id UUID REFERENCES orchestrator_pipeline_state(id) ON DELETE SET NULL,
    
    -- SAE Provenance
    sae_provenance TEXT NOT NULL CHECK (sae_provenance IN ('true_sae', 'proxy', 'proxy+true')),
    sae_version TEXT, -- Version of SAE model used
    mapping_version TEXT, -- Version of feature→pathway mapping used
    
    -- DNA Repair Capacity
    dna_repair_capacity FLOAT CHECK (dna_repair_capacity >= 0.0 AND dna_repair_capacity <= 1.0),
    
    -- Pathway Burden Scores
    pathway_burden JSONB, -- {"ddr": 0.85, "mapk": 0.12, "pi3k": 0.08, ...}
    
    -- TRUE SAE Diagnostics (if TRUE SAE used)
    ddr_bin_score FLOAT, -- DDR pathway bin score (0-1)
    ddr_sae_score FLOAT, -- DDR SAE score
    mapk_sae_score FLOAT, -- MAPK SAE score
    pi3k_sae_score FLOAT, -- PI3K SAE score
    other_pathway_scores JSONB, -- Other pathway scores
    
    -- Mechanism Vector (7D)
    mechanism_vector JSONB, -- [DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux]
    
    -- Resistance Signals
    resistance_signals JSONB, -- Detected resistance patterns
    
    -- Validation Metrics (for TRUE SAE)
    validation_metrics JSONB, -- {"auroc": 0.783, "p_value": 0.0020, "diamond_features": [27607, ...]}
    
    -- Classification
    data_classification TEXT DEFAULT 'PHI' CHECK (data_classification IN ('PHI', 'NON_PHI', 'SENSITIVE', 'PUBLIC')),
    
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    
    -- Unique constraint: one SAE feature set per patient (latest)
    UNIQUE(user_id, patient_id)
);

CREATE INDEX IF NOT EXISTS idx_sae_features_user_id ON sae_features(user_id);
CREATE INDEX IF NOT EXISTS idx_sae_features_patient_id ON sae_features(patient_id);
CREATE INDEX IF NOT EXISTS idx_sae_features_pipeline_state_id ON sae_features(pipeline_state_id);
CREATE INDEX IF NOT EXISTS idx_sae_features_provenance ON sae_features(sae_provenance);
CREATE INDEX IF NOT EXISTS idx_sae_features_ddr_bin_score ON sae_features(ddr_bin_score) WHERE ddr_bin_score IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_sae_features_data_classification ON sae_features(data_classification);

-- ============================================================================
-- 4. RESISTANCE PLAYBOOK TABLE
-- ============================================================================
-- Stores resistance playbook recommendations and alternative treatments
CREATE TABLE IF NOT EXISTS resistance_playbook (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES user_profiles(id) ON DELETE CASCADE,
    patient_id TEXT NOT NULL REFERENCES patient_profiles(patient_id) ON DELETE CASCADE,
    pipeline_state_id UUID REFERENCES orchestrator_pipeline_state(id) ON DELETE SET NULL,
    
    -- Resistance Detection
    resistance_detected BOOLEAN NOT NULL DEFAULT FALSE,
    resistance_risk_level TEXT CHECK (resistance_risk_level IN ('HIGH', 'MODERATE', 'LOW', 'NONE')),
    resistance_confidence FLOAT CHECK (resistance_confidence >= 0.0 AND resistance_confidence <= 1.0),
    
    -- Resistance Pathways (JSONB array of pathway objects with mechanisms)
    -- Example: [{"pathway": "DDR", "resistance_mechanism": "TP53 mutation", "detected_genes": ["TP53", "BRCA1"], "confidence": 0.85}, ...]
    resistance_pathways JSONB,
    
    -- Alternative Treatment Options (JSONB array of alternative drug recommendations)
    -- Example: [{"drug_name": "Olaparib", "drug_class": "PARP_inhibitor", "rationale": "Alternative to platinum-based therapy", "evidence_tier": "A", "priority": 1, "source_gene": "TP53"}, ...]
    alternative_drugs JSONB,
    
    -- Recommendations (JSONB array of action recommendations)
    -- Example: [{"action": "Switch to PARP inhibitor", "priority": "high", "rationale": "TP53 mutation detected, platinum resistance likely", "evidence_level": "A"}, ...]
    recommendations JSONB,
    
    -- Cross-Resistance Analysis (JSONB array of cross-resistance patterns)
    -- Example: [{"current_drug": "carboplatin", "potential_drug": "cisplatin", "resistance_risk": 0.8, "mechanism": "TP53 mutation", "evidence_level": "A"}, ...]
    cross_resistance_patterns JSONB,
    
    -- Monitoring Updates (JSONB array of recommended monitoring changes)
    -- Example: [{"biomarker": "CA-125", "frequency": "monthly", "rationale": "Increased monitoring due to resistance risk"}, ...]
    monitoring_updates JSONB,
    
    -- Classification
    data_classification TEXT DEFAULT 'PHI' CHECK (data_classification IN ('PHI', 'NON_PHI', 'SENSITIVE', 'PUBLIC')),
    
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    
    -- Unique constraint: one resistance playbook per patient (latest)
    UNIQUE(user_id, patient_id)
);

CREATE INDEX IF NOT EXISTS idx_resistance_playbook_user_id ON resistance_playbook(user_id);
CREATE INDEX IF NOT EXISTS idx_resistance_playbook_patient_id ON resistance_playbook(patient_id);
CREATE INDEX IF NOT EXISTS idx_resistance_playbook_pipeline_state_id ON resistance_playbook(pipeline_state_id);
CREATE INDEX IF NOT EXISTS idx_resistance_playbook_resistance_detected ON resistance_playbook(resistance_detected) WHERE resistance_detected = TRUE;
CREATE INDEX IF NOT EXISTS idx_resistance_playbook_risk_level ON resistance_playbook(resistance_risk_level);
CREATE INDEX IF NOT EXISTS idx_resistance_playbook_data_classification ON resistance_playbook(data_classification);

-- ============================================================================
-- 5. FILE UPLOADS TABLE
-- ============================================================================
-- Tracks VCF/PDF/MAF file uploads for data extraction
CREATE TABLE IF NOT EXISTS file_uploads (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES user_profiles(id) ON DELETE CASCADE,
    patient_id TEXT REFERENCES patient_profiles(patient_id) ON DELETE SET NULL,
    pipeline_state_id UUID REFERENCES orchestrator_pipeline_state(id) ON DELETE SET NULL,
    
    -- File metadata
    file_name TEXT NOT NULL,
    file_type TEXT NOT NULL CHECK (file_type IN ('VCF', 'PDF', 'MAF', 'OTHER')),
    file_size_bytes BIGINT NOT NULL,
    mime_type TEXT,
    
    -- Storage
    storage_path TEXT, -- Path in Supabase Storage or S3
    storage_url TEXT, -- Public URL if available
    
    -- Processing status
    upload_status TEXT NOT NULL DEFAULT 'uploading' CHECK (upload_status IN (
        'uploading', 'uploaded', 'processing', 'processed', 'failed', 'deleted'
    )),
    processing_status TEXT CHECK (processing_status IN (
        'pending', 'extracting', 'parsing', 'validating', 'completed', 'failed'
    )),
    
    -- Extraction results
    extracted_data JSONB, -- Extracted mutations, patient data, etc.
    extraction_errors TEXT[],
    extraction_warnings TEXT[],
    
    -- Validation
    is_valid BOOLEAN,
    validation_errors TEXT[],
    validation_warnings TEXT[],
    
    -- Metadata
    uploaded_by TEXT DEFAULT 'user',
    processing_started_at TIMESTAMPTZ,
    processing_completed_at TIMESTAMPTZ,
    processing_duration_seconds FLOAT,
    
    -- Classification
    data_classification TEXT DEFAULT 'PHI' CHECK (data_classification IN ('PHI', 'NON_PHI', 'SENSITIVE', 'PUBLIC')),
    
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    deleted_at TIMESTAMPTZ
);

CREATE INDEX IF NOT EXISTS idx_file_uploads_user_id ON file_uploads(user_id);
CREATE INDEX IF NOT EXISTS idx_file_uploads_patient_id ON file_uploads(patient_id);
CREATE INDEX IF NOT EXISTS idx_file_uploads_pipeline_state_id ON file_uploads(pipeline_state_id);
CREATE INDEX IF NOT EXISTS idx_file_uploads_file_type ON file_uploads(file_type);
CREATE INDEX IF NOT EXISTS idx_file_uploads_upload_status ON file_uploads(upload_status);
CREATE INDEX IF NOT EXISTS idx_file_uploads_processing_status ON file_uploads(processing_status);
CREATE INDEX IF NOT EXISTS idx_file_uploads_data_classification ON file_uploads(data_classification);

-- ============================================================================
-- 6. TOXICITY RISK ASSESSMENTS TABLE
-- ============================================================================
-- Stores toxicity risk assessments from MOAT toxicity risk service
CREATE TABLE IF NOT EXISTS toxicity_risk_assessments (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES user_profiles(id) ON DELETE CASCADE,
    patient_id TEXT NOT NULL REFERENCES patient_profiles(patient_id) ON DELETE CASCADE,
    pipeline_state_id UUID REFERENCES orchestrator_pipeline_state(id) ON DELETE SET NULL,
    
    -- Risk Assessment
    risk_score FLOAT NOT NULL CHECK (risk_score >= 0.0 AND risk_score <= 1.0),
    risk_level TEXT NOT NULL CHECK (risk_level IN ('HIGH', 'MODERATE', 'LOW', 'NONE')),
    confidence FLOAT CHECK (confidence >= 0.0 AND confidence <= 1.0),
    reason TEXT, -- Human-readable reason for risk level
    
    -- Contributing Factors (JSONB array of factor objects)
    -- Example: [{"type": "pharmacogene", "detail": "DPYD *2A/*2A", "weight": 0.6, "confidence": 0.95}, ...]
    contributing_factors JSONB,
    
    -- Therapeutic Candidate (JSONB object with drug/candidate info)
    -- Example: {"type": "drug", "name": "5-fluorouracil", "moa": "antimetabolite"}
    therapeutic_candidate JSONB,
    
    -- Mitigating Foods/Supplements (JSONB array of food/supplement recommendations)
    -- Example: [{"food": "NAC", "rationale": "Supports DNA repair pathways", "evidence_level": "B"}, ...]
    mitigating_foods JSONB,
    
    -- Pathway Overlap Analysis (JSONB object with pathway overlap scores)
    -- Example: {"DDR": 0.9, "MAPK": 0.1}
    pathway_overlap JSONB,
    
    -- Classification
    data_classification TEXT DEFAULT 'PHI' CHECK (data_classification IN ('PHI', 'NON_PHI', 'SENSITIVE', 'PUBLIC')),
    
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    
    -- Unique constraint: one assessment per patient/candidate combination (latest)
    UNIQUE(user_id, patient_id, therapeutic_candidate)
);

CREATE INDEX IF NOT EXISTS idx_toxicity_risk_user_id ON toxicity_risk_assessments(user_id);
CREATE INDEX IF NOT EXISTS idx_toxicity_risk_patient_id ON toxicity_risk_assessments(patient_id);
CREATE INDEX IF NOT EXISTS idx_toxicity_risk_pipeline_state_id ON toxicity_risk_assessments(pipeline_state_id);
CREATE INDEX IF NOT EXISTS idx_toxicity_risk_risk_level ON toxicity_risk_assessments(risk_level);
CREATE INDEX IF NOT EXISTS idx_toxicity_risk_risk_score ON toxicity_risk_assessments(risk_score);
CREATE INDEX IF NOT EXISTS idx_toxicity_risk_data_classification ON toxicity_risk_assessments(data_classification);

-- ============================================================================
-- 7. DOSING GUIDANCE TABLE
-- ============================================================================
-- Stores dosing guidance recommendations from MOAT dosing guidance service
CREATE TABLE IF NOT EXISTS dosing_guidance (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES user_profiles(id) ON DELETE CASCADE,
    patient_id TEXT NOT NULL REFERENCES patient_profiles(patient_id) ON DELETE CASCADE,
    pipeline_state_id UUID REFERENCES orchestrator_pipeline_state(id) ON DELETE SET NULL,
    
    -- Pharmacogene Information
    gene TEXT NOT NULL, -- "DPYD", "UGT1A1", "TPMT", etc.
    variant TEXT, -- HGVS notation (e.g., "c.1905+1G>A")
    phenotype TEXT, -- "Poor Metabolizer", "Intermediate Metabolizer", etc.
    
    -- Drug Information
    drug TEXT NOT NULL, -- Drug name
    standard_dose TEXT, -- Standard dose (e.g., "1000 mg/m²")
    
    -- Dosing Recommendation
    adjustment_type TEXT NOT NULL CHECK (adjustment_type IN (
        'CONTRAINDICATED', 'REDUCE', 'INCREASE', 'MONITOR', 'NO_CHANGE'
    )),
    adjustment_factor FLOAT, -- Percentage adjustment (e.g., 0.5 for 50% reduction)
    recommended_dose TEXT, -- Recommended dose after adjustment
    
    -- CPIC Evidence Level
    cpic_level TEXT CHECK (cpic_level IN ('A', 'B', 'C', 'D')),
    
    -- Recommendation Details
    recommendation TEXT NOT NULL, -- Plain English recommendation
    rationale TEXT, -- Scientific rationale
    monitoring TEXT[], -- Monitoring recommendations
    
    -- Alternatives
    alternatives TEXT[], -- Alternative drugs if contraindicated
    
    -- Cumulative Toxicity Alert
    cumulative_toxicity_alert TEXT, -- Alert message if cumulative toxicity risk
    
    -- Treatment Context
    treatment_line INTEGER,
    prior_therapies TEXT[],
    disease TEXT,
    
    -- Classification
    data_classification TEXT DEFAULT 'PHI' CHECK (data_classification IN ('PHI', 'NON_PHI', 'SENSITIVE', 'PUBLIC')),
    
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    
    -- Unique constraint: one guidance per patient/gene/drug combination (latest)
    UNIQUE(user_id, patient_id, gene, drug)
);

CREATE INDEX IF NOT EXISTS idx_dosing_guidance_user_id ON dosing_guidance(user_id);
CREATE INDEX IF NOT EXISTS idx_dosing_guidance_patient_id ON dosing_guidance(patient_id);
CREATE INDEX IF NOT EXISTS idx_dosing_guidance_pipeline_state_id ON dosing_guidance(pipeline_state_id);
CREATE INDEX IF NOT EXISTS idx_dosing_guidance_gene ON dosing_guidance(gene);
CREATE INDEX IF NOT EXISTS idx_dosing_guidance_drug ON dosing_guidance(drug);
CREATE INDEX IF NOT EXISTS idx_dosing_guidance_adjustment_type ON dosing_guidance(adjustment_type);
CREATE INDEX IF NOT EXISTS idx_dosing_guidance_cpic_level ON dosing_guidance(cpic_level);
CREATE INDEX IF NOT EXISTS idx_dosing_guidance_data_classification ON dosing_guidance(data_classification);

-- ============================================================================
-- 8. ROW-LEVEL SECURITY (RLS) POLICIES
-- ============================================================================

-- Enable RLS on new tables
ALTER TABLE orchestrator_pipeline_state ENABLE ROW LEVEL SECURITY;
ALTER TABLE agent_executions ENABLE ROW LEVEL SECURITY;
ALTER TABLE sae_features ENABLE ROW LEVEL SECURITY;
ALTER TABLE resistance_playbook ENABLE ROW LEVEL SECURITY;
ALTER TABLE file_uploads ENABLE ROW LEVEL SECURITY;
ALTER TABLE toxicity_risk_assessments ENABLE ROW LEVEL SECURITY;
ALTER TABLE dosing_guidance ENABLE ROW LEVEL SECURITY;

-- Policy: Users can only see their own orchestrator pipeline states
CREATE POLICY "Users can view own pipeline states"
    ON orchestrator_pipeline_state FOR SELECT
    USING (auth.uid() = user_id);

CREATE POLICY "Users can insert own pipeline states"
    ON orchestrator_pipeline_state FOR INSERT
    WITH CHECK (auth.uid() = user_id);

CREATE POLICY "Users can update own pipeline states"
    ON orchestrator_pipeline_state FOR UPDATE
    USING (auth.uid() = user_id);

-- Policy: Users can only see their own agent executions
CREATE POLICY "Users can view own agent executions"
    ON agent_executions FOR SELECT
    USING (auth.uid() = user_id);

CREATE POLICY "Users can insert own agent executions"
    ON agent_executions FOR INSERT
    WITH CHECK (auth.uid() = user_id);

CREATE POLICY "Users can update own agent executions"
    ON agent_executions FOR UPDATE
    USING (auth.uid() = user_id);

-- Policy: Users can only see their own SAE features
CREATE POLICY "Users can view own SAE features"
    ON sae_features FOR SELECT
    USING (auth.uid() = user_id);

CREATE POLICY "Users can insert own SAE features"
    ON sae_features FOR INSERT
    WITH CHECK (auth.uid() = user_id);

CREATE POLICY "Users can update own SAE features"
    ON sae_features FOR UPDATE
    USING (auth.uid() = user_id);

-- Policy: Users can only see their own resistance playbooks
CREATE POLICY "Users can view own resistance playbooks"
    ON resistance_playbook FOR SELECT
    USING (auth.uid() = user_id);

CREATE POLICY "Users can insert own resistance playbooks"
    ON resistance_playbook FOR INSERT
    WITH CHECK (auth.uid() = user_id);

CREATE POLICY "Users can update own resistance playbooks"
    ON resistance_playbook FOR UPDATE
    USING (auth.uid() = user_id);

-- Policy: Users can only see their own file uploads
CREATE POLICY "Users can view own file uploads"
    ON file_uploads FOR SELECT
    USING (auth.uid() = user_id);

CREATE POLICY "Users can insert own file uploads"
    ON file_uploads FOR INSERT
    WITH CHECK (auth.uid() = user_id);

CREATE POLICY "Users can update own file uploads"
    ON file_uploads FOR UPDATE
    USING (auth.uid() = user_id);

CREATE POLICY "Users can delete own file uploads"
    ON file_uploads FOR DELETE
    USING (auth.uid() = user_id);

-- Policy: Users can only see their own toxicity risk assessments
CREATE POLICY "Users can view own toxicity risk assessments"
    ON toxicity_risk_assessments FOR SELECT
    USING (auth.uid() = user_id);

CREATE POLICY "Users can insert own toxicity risk assessments"
    ON toxicity_risk_assessments FOR INSERT
    WITH CHECK (auth.uid() = user_id);

CREATE POLICY "Users can update own toxicity risk assessments"
    ON toxicity_risk_assessments FOR UPDATE
    USING (auth.uid() = user_id);

-- Policy: Users can only see their own dosing guidance
CREATE POLICY "Users can view own dosing guidance"
    ON dosing_guidance FOR SELECT
    USING (auth.uid() = user_id);

CREATE POLICY "Users can insert own dosing guidance"
    ON dosing_guidance FOR INSERT
    WITH CHECK (auth.uid() = user_id);

CREATE POLICY "Users can update own dosing guidance"
    ON dosing_guidance FOR UPDATE
    USING (auth.uid() = user_id);

-- ============================================================================
-- 9. HELPER FUNCTIONS
-- ============================================================================

-- Function to get latest pipeline state for a patient
CREATE OR REPLACE FUNCTION get_latest_pipeline_state(
    p_user_id UUID,
    p_patient_id TEXT
) RETURNS orchestrator_pipeline_state AS $$
BEGIN
    RETURN (
        SELECT * FROM orchestrator_pipeline_state
        WHERE user_id = p_user_id AND patient_id = p_patient_id
        ORDER BY updated_at DESC
        LIMIT 1
    );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Function to update pipeline progress
CREATE OR REPLACE FUNCTION update_pipeline_progress(
    p_pipeline_state_id UUID,
    p_phase TEXT,
    p_progress_percentage FLOAT
) RETURNS VOID AS $$
BEGIN
    UPDATE orchestrator_pipeline_state
    SET 
        current_phase = p_phase,
        progress_percentage = p_progress_percentage,
        updated_at = NOW()
    WHERE id = p_pipeline_state_id;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Function to mark pipeline as complete
CREATE OR REPLACE FUNCTION mark_pipeline_complete(
    p_pipeline_state_id UUID
) RETURNS VOID AS $$
BEGIN
    UPDATE orchestrator_pipeline_state
    SET 
        current_phase = 'complete',
        progress_percentage = 100.0,
        is_complete = TRUE,
        completed_at = NOW(),
        updated_at = NOW()
    WHERE id = p_pipeline_state_id;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- ============================================================================
-- DONE! ✅
-- ============================================================================
-- New tables added for MOAT Orchestrator & Capabilities:
-- 1. orchestrator_pipeline_state - Track 7-phase pipeline execution
-- 2. agent_executions - Track individual agent execution history
-- 3. sae_features - Store SAE features, DDR_bin scores, TRUE SAE provenance
-- 4. resistance_playbook - Store resistance playbook recommendations
-- 5. file_uploads - Track VCF/PDF/MAF file uploads
-- 6. toxicity_risk_assessments - Store toxicity risk assessments
-- 7. dosing_guidance - Store dosing guidance recommendations
--
-- These tables support:
-- ✅ Orchestrator pipeline integration (Priority 1)
-- ✅ SAE Features component (Priority 2)
-- ✅ Resistance Playbook component (Priority 2)
-- ✅ File upload capability (Priority 1)
-- ✅ Pipeline status tracking (Priority 1)
-- ✅ All MOAT capabilities storage
--
-- Next steps:
-- 1. Update backend code to use these tables
-- 2. Migrate existing orchestrator state (if any) to database
-- 3. Test RLS policies
-- 4. Update frontend to query these tables
-- ============================================================================

