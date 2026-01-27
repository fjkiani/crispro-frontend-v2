-- ============================================================================
-- SUPABASE SCHEMA UPDATES - Based on Recent Work
-- ============================================================================
-- Purpose: Add tables for trial matching, research intelligence, biomarker tracking
-- Run this AFTER running SUPABASE_COMPLETE_SETUP.sql
-- ============================================================================

-- ============================================================================
-- 1. TRIAL MATCHES TABLE
-- ============================================================================
-- Stores trial matching results with scores, reasoning, and MoA vectors
CREATE TABLE IF NOT EXISTS trial_matches (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES user_profiles(id) ON DELETE CASCADE,
    patient_id TEXT REFERENCES patient_profiles(patient_id),
    
    -- Trial identification
    nct_id TEXT NOT NULL,
    title TEXT NOT NULL,
    phase TEXT,
    status TEXT, -- Recruiting, Active, etc.
    
    -- Matching scores
    mechanism_fit_score FLOAT DEFAULT 0.0 CHECK (mechanism_fit_score >= 0.0 AND mechanism_fit_score <= 1.0),
    eligibility_score FLOAT DEFAULT 0.0 CHECK (eligibility_score >= 0.0 AND eligibility_score <= 1.0),
    combined_score FLOAT DEFAULT 0.0 CHECK (combined_score >= 0.0 AND combined_score <= 1.0),
    
    -- Mechanism of Action (7D vector: DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux)
    trial_moa_vector JSONB, -- {"ddr": 0.8, "mapk": 0.2, ...}
    patient_mechanism_vector JSONB, -- Patient's mechanism vector used for matching
    
    -- Reasoning and eligibility
    reasoning JSONB, -- Structured reasoning object (why_eligible, why_good_fit, etc.)
    eligibility_criteria JSONB, -- Meets, excludes, uncertain criteria
    
    -- Trial details
    interventions TEXT[],
    locations JSONB, -- Array of location objects
    contact_info JSONB, -- Contact name, phone, email
    source_url TEXT, -- ClinicalTrials.gov URL
    
    -- Freshness metadata
    last_refreshed_at TIMESTAMPTZ,
    stale BOOLEAN DEFAULT FALSE,
    staleness_reason TEXT,
    refresh_needed BOOLEAN DEFAULT FALSE,
    refresh_sla_hours INTEGER, -- Hours until refresh needed
    
    -- Provenance
    query_matched TEXT, -- Which search query found this trial
    tag_confidence FLOAT, -- Confidence of MoA tag (0-1)
    tag_provider TEXT, -- "cohere", "gemini", etc.
    tag_model TEXT, -- Model used for tagging
    
    -- User interaction
    pinned BOOLEAN DEFAULT FALSE,
    viewed_at TIMESTAMPTZ,
    viewed_count INTEGER DEFAULT 0,
    
    -- Classification
    data_classification TEXT DEFAULT 'NON_PHI' CHECK (data_classification IN ('PHI', 'NON_PHI', 'SENSITIVE', 'PUBLIC')),
    
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    
    -- Unique constraint: one match per user/patient/trial combination
    UNIQUE(user_id, patient_id, nct_id)
);

CREATE INDEX IF NOT EXISTS idx_trial_matches_user_id ON trial_matches(user_id);
CREATE INDEX IF NOT EXISTS idx_trial_matches_patient_id ON trial_matches(patient_id);
CREATE INDEX IF NOT EXISTS idx_trial_matches_nct_id ON trial_matches(nct_id);
CREATE INDEX IF NOT EXISTS idx_trial_matches_combined_score ON trial_matches(combined_score DESC);
CREATE INDEX IF NOT EXISTS idx_trial_matches_pinned ON trial_matches(pinned) WHERE pinned = TRUE;
CREATE INDEX IF NOT EXISTS idx_trial_matches_stale ON trial_matches(stale) WHERE stale = TRUE;
CREATE INDEX IF NOT EXISTS idx_trial_matches_data_classification ON trial_matches(data_classification);

-- ============================================================================
-- 2. RESEARCH INTELLIGENCE QUERIES TABLE
-- ============================================================================
-- Stores research intelligence queries and their results
CREATE TABLE IF NOT EXISTS research_intelligence_queries (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES user_profiles(id) ON DELETE CASCADE,
    patient_id TEXT REFERENCES patient_profiles(patient_id),
    
    -- Query details
    question TEXT NOT NULL,
    research_plan JSONB, -- Structured research plan (entities, sub-questions, etc.)
    context JSONB, -- Patient context (disease, treatment_line, etc.)
    
    -- Results
    portal_results JSONB, -- PubMed, GDC, Project Data Sphere results
    parsed_content JSONB, -- Deep-parsed article content
    synthesized_findings JSONB, -- LLM-synthesized findings (mechanisms, dosage, safety, outcomes)
    moat_analysis JSONB, -- MOAT pathway analysis
    article_summaries JSONB, -- Per-article summaries
    sub_question_answers JSONB, -- Individual sub-question answers
    
    -- Provenance
    methods_used TEXT[], -- List of methods used (llm_deep_research, diffbot_extraction, etc.)
    run_id TEXT, -- Unique run identifier
    llm_provider TEXT, -- "cohere", "gemini", etc.
    llm_model TEXT, -- Model used
    
    -- Metadata
    articles_parsed INTEGER DEFAULT 0,
    mechanisms_found INTEGER DEFAULT 0,
    sub_questions_answered INTEGER DEFAULT 0,
    overall_confidence FLOAT CHECK (overall_confidence >= 0.0 AND overall_confidence <= 1.0),
    
    -- Classification
    data_classification TEXT DEFAULT 'NON_PHI' CHECK (data_classification IN ('PHI', 'NON_PHI', 'SENSITIVE', 'PUBLIC')),
    
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_research_queries_user_id ON research_intelligence_queries(user_id);
CREATE INDEX IF NOT EXISTS idx_research_queries_patient_id ON research_intelligence_queries(patient_id);
CREATE INDEX IF NOT EXISTS idx_research_queries_created_at ON research_intelligence_queries(created_at DESC);
CREATE INDEX IF NOT EXISTS idx_research_queries_data_classification ON research_intelligence_queries(data_classification);

-- ============================================================================
-- 3. BIOMARKER HISTORY TABLE
-- ============================================================================
-- Tracks biomarker values over time (CA-125, PSA, CEA, etc.)
CREATE TABLE IF NOT EXISTS biomarker_history (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES user_profiles(id) ON DELETE CASCADE,
    patient_id TEXT REFERENCES patient_profiles(patient_id),
    
    -- Biomarker details
    biomarker_name TEXT NOT NULL, -- "CA-125", "PSA", "CEA", etc.
    biomarker_value FLOAT NOT NULL,
    biomarker_unit TEXT, -- "U/mL", "ng/mL", etc.
    
    -- Clinical context
    measurement_date DATE NOT NULL,
    treatment_line INTEGER, -- Treatment line at time of measurement
    treatment_context TEXT, -- What treatment was ongoing
    
    -- Analysis results (from biomarker intelligence services)
    burden_assessment JSONB, -- High/Medium/Low burden assessment
    response_forecast JSONB, -- Forecasted response trajectory
    resistance_signals JSONB, -- Resistance detection signals
    
    -- Classification
    data_classification TEXT DEFAULT 'PHI' CHECK (data_classification IN ('PHI', 'NON_PHI', 'SENSITIVE', 'PUBLIC')),
    
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_biomarker_history_user_id ON biomarker_history(user_id);
CREATE INDEX IF NOT EXISTS idx_biomarker_history_patient_id ON biomarker_history(patient_id);
CREATE INDEX IF NOT EXISTS idx_biomarker_history_biomarker_name ON biomarker_history(biomarker_name);
CREATE INDEX IF NOT EXISTS idx_biomarker_history_measurement_date ON biomarker_history(measurement_date DESC);
CREATE INDEX IF NOT EXISTS idx_biomarker_history_data_classification ON biomarker_history(data_classification);

-- ============================================================================
-- 4. CARE PLAN VERSIONS TABLE
-- ============================================================================
-- Tracks care plan history and versions
CREATE TABLE IF NOT EXISTS care_plan_versions (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES user_profiles(id) ON DELETE CASCADE,
    patient_id TEXT REFERENCES patient_profiles(patient_id),
    
    -- Version details
    version_number INTEGER NOT NULL,
    care_plan_data JSONB NOT NULL, -- Full care plan JSON
    
    -- What changed
    change_summary TEXT, -- Human-readable summary of changes
    changed_sections TEXT[], -- List of sections that changed
    
    -- Context
    created_by TEXT DEFAULT 'system', -- "system", "user", "agent"
    trigger_event TEXT, -- What triggered this version (e.g., "new_ca125_value", "trial_match_update")
    
    -- Classification
    data_classification TEXT DEFAULT 'PHI' CHECK (data_classification IN ('PHI', 'NON_PHI', 'SENSITIVE', 'PUBLIC')),
    
    created_at TIMESTAMPTZ DEFAULT NOW(),
    
    -- Unique constraint: one version number per user/patient
    UNIQUE(user_id, patient_id, version_number)
);

CREATE INDEX IF NOT EXISTS idx_care_plan_versions_user_id ON care_plan_versions(user_id);
CREATE INDEX IF NOT EXISTS idx_care_plan_versions_patient_id ON care_plan_versions(patient_id);
CREATE INDEX IF NOT EXISTS idx_care_plan_versions_version_number ON care_plan_versions(version_number DESC);
CREATE INDEX IF NOT EXISTS idx_care_plan_versions_data_classification ON care_plan_versions(data_classification);

-- ============================================================================
-- 5. UPDATE PATIENT PROFILES TABLE
-- ============================================================================
-- Add mechanism_vector column to patient_profiles if it doesn't exist
DO $$ 
BEGIN
    IF NOT EXISTS (
        SELECT 1 FROM information_schema.columns 
        WHERE table_name = 'patient_profiles' 
        AND column_name = 'mechanism_vector'
    ) THEN
        ALTER TABLE patient_profiles 
        ADD COLUMN mechanism_vector JSONB;
        
        COMMENT ON COLUMN patient_profiles.mechanism_vector IS 
        '7D mechanism vector: DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux';
    END IF;
END $$;

-- ============================================================================
-- 6. ROW-LEVEL SECURITY (RLS) POLICIES
-- ============================================================================

-- Enable RLS on new tables
ALTER TABLE trial_matches ENABLE ROW LEVEL SECURITY;
ALTER TABLE research_intelligence_queries ENABLE ROW LEVEL SECURITY;
ALTER TABLE biomarker_history ENABLE ROW LEVEL SECURITY;
ALTER TABLE care_plan_versions ENABLE ROW LEVEL SECURITY;

-- Policy: Users can only see their own trial matches
CREATE POLICY "Users can view own trial matches"
    ON trial_matches FOR SELECT
    USING (auth.uid() = user_id);

CREATE POLICY "Users can insert own trial matches"
    ON trial_matches FOR INSERT
    WITH CHECK (auth.uid() = user_id);

CREATE POLICY "Users can update own trial matches"
    ON trial_matches FOR UPDATE
    USING (auth.uid() = user_id);

-- Policy: Users can only see their own research queries
CREATE POLICY "Users can view own research queries"
    ON research_intelligence_queries FOR SELECT
    USING (auth.uid() = user_id);

CREATE POLICY "Users can insert own research queries"
    ON research_intelligence_queries FOR INSERT
    WITH CHECK (auth.uid() = user_id);

-- Policy: Users can only see their own biomarker history
CREATE POLICY "Users can view own biomarker history"
    ON biomarker_history FOR SELECT
    USING (auth.uid() = user_id);

CREATE POLICY "Users can insert own biomarker history"
    ON biomarker_history FOR INSERT
    WITH CHECK (auth.uid() = user_id);

CREATE POLICY "Users can update own biomarker history"
    ON biomarker_history FOR UPDATE
    USING (auth.uid() = user_id);

-- Policy: Users can only see their own care plan versions
CREATE POLICY "Users can view own care plan versions"
    ON care_plan_versions FOR SELECT
    USING (auth.uid() = user_id);

CREATE POLICY "Users can insert own care plan versions"
    ON care_plan_versions FOR INSERT
    WITH CHECK (auth.uid() = user_id);

-- ============================================================================
-- 7. HELPER FUNCTIONS
-- ============================================================================

-- Function to get latest care plan version number for a patient
CREATE OR REPLACE FUNCTION get_latest_care_plan_version(
    p_user_id UUID,
    p_patient_id TEXT
) RETURNS INTEGER AS $$
BEGIN
    RETURN COALESCE(
        (SELECT MAX(version_number) 
         FROM care_plan_versions 
         WHERE user_id = p_user_id AND patient_id = p_patient_id),
        0
    );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Function to auto-increment care plan version
CREATE OR REPLACE FUNCTION auto_increment_care_plan_version()
RETURNS TRIGGER AS $$
BEGIN
    NEW.version_number := get_latest_care_plan_version(NEW.user_id, NEW.patient_id) + 1;
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

-- Trigger to auto-increment version number
DROP TRIGGER IF EXISTS auto_increment_care_plan_version_trigger ON care_plan_versions;
CREATE TRIGGER auto_increment_care_plan_version_trigger
    BEFORE INSERT ON care_plan_versions
    FOR EACH ROW
    WHEN (NEW.version_number IS NULL OR NEW.version_number = 0)
    EXECUTE FUNCTION auto_increment_care_plan_version();

-- ============================================================================
-- DONE! ✅
-- ============================================================================
-- New tables added:
-- 1. trial_matches - Store trial matching results with scores and reasoning
-- 2. research_intelligence_queries - Store research queries and results
-- 3. biomarker_history - Track biomarker values over time
-- 4. care_plan_versions - Track care plan history
--
-- Updates:
-- - patient_profiles.mechanism_vector - Added 7D mechanism vector column
--
-- Next steps:
-- 1. Update backend code to use these tables
-- 2. Migrate existing data if needed
-- 3. Test RLS policies
-- ============================================================================

