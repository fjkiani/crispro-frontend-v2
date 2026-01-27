-- ============================================================================
-- COMPLETE SUPABASE SETUP - Run ALL of this in Supabase SQL Editor
-- ============================================================================
-- Purpose: Create all required tables for CrisPRO Oncology application
-- Run this in: Supabase Dashboard → SQL Editor → New Query
-- ============================================================================

-- ============================================================================
-- 0. EXTENSIONS (Required for gen_random_uuid)
-- ============================================================================
CREATE EXTENSION IF NOT EXISTS pgcrypto;


-- ============================================================================
-- 1. USER PROFILES TABLE
-- ============================================================================
CREATE TABLE IF NOT EXISTS user_profiles (
    id UUID PRIMARY KEY REFERENCES auth.users(id) ON DELETE CASCADE,
    email TEXT UNIQUE NOT NULL,
    role TEXT DEFAULT 'patient' CHECK (role IN ('patient', 'researcher', 'admin')),
    tier TEXT DEFAULT 'free' CHECK (tier IN ('free', 'pro', 'enterprise')),
    full_name TEXT,
    institution TEXT,
    country TEXT,
    timezone TEXT,
    data_classification TEXT DEFAULT 'NON_PHI' CHECK (data_classification IN ('PHI', 'NON_PHI', 'SENSITIVE', 'PUBLIC')),
    mfa_enabled BOOLEAN DEFAULT FALSE,
    mfa_secret TEXT,
    mfa_verified_at TIMESTAMPTZ,
    mfa_backup_codes TEXT[],
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- Indexes for user_profiles
CREATE INDEX IF NOT EXISTS idx_user_profiles_email ON user_profiles(email);
CREATE INDEX IF NOT EXISTS idx_user_profiles_role ON user_profiles(role);
CREATE INDEX IF NOT EXISTS idx_user_profiles_tier ON user_profiles(tier);
CREATE INDEX IF NOT EXISTS idx_user_profiles_mfa_enabled ON user_profiles(mfa_enabled) WHERE mfa_enabled = TRUE;

-- ============================================================================
-- 2. USER QUOTAS TABLE
-- ============================================================================
CREATE TABLE IF NOT EXISTS user_quotas (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES user_profiles(id) ON DELETE CASCADE,
    tier TEXT DEFAULT 'free' CHECK (tier IN ('free', 'pro', 'enterprise')),
    variant_analyses_limit INTEGER DEFAULT 10,
    variant_analyses_used INTEGER DEFAULT 0,
    drug_queries_limit INTEGER DEFAULT 5,
    drug_queries_used INTEGER DEFAULT 0,
    food_queries_limit INTEGER DEFAULT 3,
    food_queries_used INTEGER DEFAULT 0,
    clinical_trials_limit INTEGER DEFAULT 0,
    clinical_trials_used INTEGER DEFAULT 0,
    period_start TIMESTAMPTZ DEFAULT NOW(),
    period_end TIMESTAMPTZ DEFAULT (NOW() + INTERVAL '1 month'),
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    UNIQUE(user_id)
);

CREATE INDEX IF NOT EXISTS idx_user_quotas_user_id ON user_quotas(user_id);

-- ============================================================================
-- 3. PATIENT PROFILES TABLE
-- ============================================================================
CREATE TABLE IF NOT EXISTS patient_profiles (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES user_profiles(id) ON DELETE CASCADE,
    patient_id TEXT UNIQUE,
    disease TEXT,
    treatment_line INTEGER DEFAULT 1,
    ca125_value FLOAT,
    ca125_baseline FLOAT,
    ca125_last_updated TIMESTAMPTZ,
    stage TEXT,
    germline_status TEXT,
    has_ascites BOOLEAN DEFAULT FALSE,
    has_peritoneal_disease BOOLEAN DEFAULT FALSE,
    location_state TEXT,
    ecog_status INTEGER,
    tumor_context JSONB,
    treatment_history JSONB,
    data_classification TEXT DEFAULT 'PHI' CHECK (data_classification IN ('PHI', 'NON_PHI', 'SENSITIVE', 'PUBLIC')),
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_patient_profiles_user_id ON patient_profiles(user_id);
CREATE INDEX IF NOT EXISTS idx_patient_profiles_data_classification ON patient_profiles(data_classification);

-- ============================================================================
-- 4. PATIENT SESSIONS TABLE
-- ============================================================================
CREATE TABLE IF NOT EXISTS patient_sessions (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES user_profiles(id) ON DELETE CASCADE,
    patient_id TEXT REFERENCES patient_profiles(patient_id),
    session_type TEXT DEFAULT 'care_plan',
    session_data JSONB,
    ca125_value FLOAT,
    data_classification TEXT DEFAULT 'PHI' CHECK (data_classification IN ('PHI', 'NON_PHI', 'SENSITIVE', 'PUBLIC')),
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    expires_at TIMESTAMPTZ DEFAULT (NOW() + INTERVAL '24 hours')
);

CREATE INDEX IF NOT EXISTS idx_patient_sessions_user_id ON patient_sessions(user_id);
CREATE INDEX IF NOT EXISTS idx_patient_sessions_expires_at ON patient_sessions(expires_at);
CREATE INDEX IF NOT EXISTS idx_patient_sessions_data_classification ON patient_sessions(data_classification);

-- ============================================================================
-- 5. PATIENT CARE PLANS TABLE
-- ============================================================================
CREATE TABLE IF NOT EXISTS patient_care_plans (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES user_profiles(id) ON DELETE CASCADE,
    patient_id TEXT REFERENCES patient_profiles(patient_id),
    care_plan_data JSONB NOT NULL,
    data_classification TEXT DEFAULT 'PHI' CHECK (data_classification IN ('PHI', 'NON_PHI', 'SENSITIVE', 'PUBLIC')),
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_patient_care_plans_user_id ON patient_care_plans(user_id);
CREATE INDEX IF NOT EXISTS idx_patient_care_plans_data_classification ON patient_care_plans(data_classification);

-- ============================================================================
-- 6. ROW-LEVEL SECURITY (RLS) POLICIES
-- ============================================================================
Error: Failed to run sql query: ERROR: 42710: policy "Users can view own profile" for table "user_profiles" already exists



-- Enable RLS on all tables
ALTER TABLE user_profiles ENABLE ROW LEVEL SECURITY;
ALTER TABLE user_quotas ENABLE ROW LEVEL SECURITY;
ALTER TABLE patient_profiles ENABLE ROW LEVEL SECURITY;
ALTER TABLE patient_sessions ENABLE ROW LEVEL SECURITY;
ALTER TABLE patient_care_plans ENABLE ROW LEVEL SECURITY;

-- Policy: Users can only see their own profile
CREATE POLICY "Users can view own profile"
    ON user_profiles FOR SELECT
    USING (auth.uid() = id);

CREATE POLICY "Users can update own profile"
    ON user_profiles FOR UPDATE
    USING (auth.uid() = id);

-- Policy: Users can only see their own quotas
CREATE POLICY "Users can view own quotas"
    ON user_quotas FOR SELECT
    USING (auth.uid() = user_id);

-- Policy: Users can only see their own patient profiles
CREATE POLICY "Users can view own patient profiles"
    ON patient_profiles FOR SELECT
    USING (auth.uid() = user_id);

CREATE POLICY "Users can insert own patient profiles"
    ON patient_profiles FOR INSERT
    WITH CHECK (auth.uid() = user_id);

CREATE POLICY "Users can update own patient profiles"
    ON patient_profiles FOR UPDATE
    USING (auth.uid() = user_id);

-- Policy: Users can only see their own sessions
CREATE POLICY "Users can view own sessions"
    ON patient_sessions FOR SELECT
    USING (auth.uid() = user_id);

CREATE POLICY "Users can insert own sessions"
    ON patient_sessions FOR INSERT
    WITH CHECK (auth.uid() = user_id);

-- Function to create user profile when user signs up
CREATE OR REPLACE FUNCTION public.handle_new_user()
RETURNS TRIGGER AS $$
BEGIN
    INSERT INTO public.user_profiles (id, email, role, tier)
    VALUES (
        NEW.id,
        NEW.email,
        COALESCE(NEW.raw_user_meta_data->>'role', 'patient'),
        'free'
    );
    
    -- Create default quotas
    INSERT INTO public.user_quotas (user_id, tier)
    VALUES (NEW.id, 'free');
    
    RETURN NEW;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Trigger to call function on new user
DROP TRIGGER IF EXISTS on_auth_user_created ON auth.users;
CREATE TRIGGER on_auth_user_created
    AFTER INSERT ON auth.users
    FOR EACH ROW
    EXECUTE FUNCTION public.handle_new_user();


-- Policy: Users can only see their own care plans
-- Before (would fail if policy exists):
CREATE POLICY "Users can view own care plans"
    ON patient_care_plans FOR SELECT
    USING (auth.uid() = user_id);

-- After (safe to run multiple times):
DROP POLICY IF EXISTS "Users can view own care plans" ON patient_care_plans;

CREATE POLICY "Users can view own care plans"
    ON patient_care_plans FOR SELECT
    USING (auth.uid() = user_id);

-- ============================================================================
-- 7. TRIGGERS (Auto-create user profile on signup)
-- ============================================================================



-- ============================================================================
-- 8. PGX + ADVANCED CARE PLAN (PRODUCTION STORAGE)
-- ============================================================================
-- Purpose: Persist the artifacts we actually generate in production:
-- - Germline/PGx variants (input + extracted)
-- - Complete care plan run receipts (request/response/provenance)
-- - Optional: per-drug PGx dosing guidance run receipts

-- 8.1 User quotas: add PGx + care-plan run counters (keeps gating honest)
ALTER TABLE user_quotas ADD COLUMN IF NOT EXISTS pgx_queries_limit INTEGER DEFAULT 10;
ALTER TABLE user_quotas ADD COLUMN IF NOT EXISTS pgx_queries_used INTEGER DEFAULT 0;
ALTER TABLE user_quotas ADD COLUMN IF NOT EXISTS complete_care_limit INTEGER DEFAULT 5;
ALTER TABLE user_quotas ADD COLUMN IF NOT EXISTS complete_care_used INTEGER DEFAULT 0;

-- 8.2 Patient profile: store germline/PGx variants and current meds (do NOT store raw VCF)
ALTER TABLE patient_profiles ADD COLUMN IF NOT EXISTS germline_variants JSONB;
ALTER TABLE patient_profiles ADD COLUMN IF NOT EXISTS pgx_variants JSONB;
ALTER TABLE patient_profiles ADD COLUMN IF NOT EXISTS current_medications JSONB;
ALTER TABLE patient_profiles ADD COLUMN IF NOT EXISTS ca125_history JSONB;

CREATE INDEX IF NOT EXISTS idx_patient_profiles_germline_variants_gin ON patient_profiles USING GIN (germline_variants);
CREATE INDEX IF NOT EXISTS idx_patient_profiles_pgx_variants_gin ON patient_profiles USING GIN (pgx_variants);
CREATE INDEX IF NOT EXISTS idx_patient_profiles_current_meds_gin ON patient_profiles USING GIN (current_medications);

-- 8.3 Care plan runs: store request/response/provenance for /api/ayesha/complete_care_v2 and /api/complete_care/v2
CREATE TABLE IF NOT EXISTS patient_care_plan_runs (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES user_profiles(id) ON DELETE CASCADE,
    patient_id TEXT REFERENCES patient_profiles(patient_id),

    -- Which orchestrator produced the run
    orchestrator TEXT NOT NULL, -- e.g. "ayesha_complete_care_v2", "universal_complete_care_v2"

    -- Full request/response receipts (JSON) for reproducibility & auditing
    request_payload JSONB NOT NULL,
    response_payload JSONB NOT NULL,
    provenance JSONB,

    -- Run status & diagnostics
    status TEXT DEFAULT 'ok' CHECK (status IN ('ok', 'error')),
    error_message TEXT,

    -- Classification (default PHI because patient payloads can include tumor context)
    data_classification TEXT DEFAULT 'PHI' CHECK (data_classification IN ('PHI', 'NON_PHI', 'SENSITIVE', 'PUBLIC')),

    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_patient_care_plan_runs_user_id ON patient_care_plan_runs(user_id);
CREATE INDEX IF NOT EXISTS idx_patient_care_plan_runs_patient_id ON patient_care_plan_runs(patient_id);
CREATE INDEX IF NOT EXISTS idx_patient_care_plan_runs_created_at ON patient_care_plan_runs(created_at);
CREATE INDEX IF NOT EXISTS idx_patient_care_plan_runs_request_gin ON patient_care_plan_runs USING GIN (request_payload);
CREATE INDEX IF NOT EXISTS idx_patient_care_plan_runs_response_gin ON patient_care_plan_runs USING GIN (response_payload);

-- 8.4 PGx dosing guidance runs (optional audit trail per drug/gene)
CREATE TABLE IF NOT EXISTS pgx_dosing_guidance_runs (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES user_profiles(id) ON DELETE CASCADE,
    patient_id TEXT REFERENCES patient_profiles(patient_id),

    gene TEXT NOT NULL,
    variant TEXT,
    drug TEXT NOT NULL,

    toxicity_tier TEXT,
    adjustment_factor FLOAT,
    contraindicated BOOLEAN DEFAULT FALSE,
    cpic_level TEXT,

    recommendations JSONB,
    alerts JSONB,
    provenance JSONB,

    data_classification TEXT DEFAULT 'PHI' CHECK (data_classification IN ('PHI', 'NON_PHI', 'SENSITIVE', 'PUBLIC')),

    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_pgx_dosing_runs_user_id ON pgx_dosing_guidance_runs(user_id);
CREATE INDEX IF NOT EXISTS idx_pgx_dosing_runs_patient_id ON pgx_dosing_guidance_runs(patient_id);
CREATE INDEX IF NOT EXISTS idx_pgx_dosing_runs_gene ON pgx_dosing_guidance_runs(gene);
CREATE INDEX IF NOT EXISTS idx_pgx_dosing_runs_drug ON pgx_dosing_guidance_runs(drug);

-- ============================================================================
-- 9. UPDATED_AT TRIGGERS (Quality-of-life for UI + audit)
-- ============================================================================
CREATE OR REPLACE FUNCTION public.set_updated_at()
RETURNS TRIGGER AS $$
BEGIN
  NEW.updated_at = NOW();
  RETURN NEW;
END;
$$ LANGUAGE plpgsql;

DROP TRIGGER IF EXISTS set_user_profiles_updated_at ON user_profiles;
CREATE TRIGGER set_user_profiles_updated_at BEFORE UPDATE ON user_profiles
FOR EACH ROW EXECUTE FUNCTION public.set_updated_at();

DROP TRIGGER IF EXISTS set_user_quotas_updated_at ON user_quotas;
CREATE TRIGGER set_user_quotas_updated_at BEFORE UPDATE ON user_quotas
FOR EACH ROW EXECUTE FUNCTION public.set_updated_at();

DROP TRIGGER IF EXISTS set_patient_profiles_updated_at ON patient_profiles;
CREATE TRIGGER set_patient_profiles_updated_at BEFORE UPDATE ON patient_profiles
FOR EACH ROW EXECUTE FUNCTION public.set_updated_at();

DROP TRIGGER IF EXISTS set_patient_sessions_updated_at ON patient_sessions;
CREATE TRIGGER set_patient_sessions_updated_at BEFORE UPDATE ON patient_sessions
FOR EACH ROW EXECUTE FUNCTION public.set_updated_at();

DROP TRIGGER IF EXISTS set_patient_care_plans_updated_at ON patient_care_plans;
CREATE TRIGGER set_patient_care_plans_updated_at BEFORE UPDATE ON patient_care_plans
FOR EACH ROW EXECUTE FUNCTION public.set_updated_at();

DROP TRIGGER IF EXISTS set_patient_care_plan_runs_updated_at ON patient_care_plan_runs;
CREATE TRIGGER set_patient_care_plan_runs_updated_at BEFORE UPDATE ON patient_care_plan_runs
FOR EACH ROW EXECUTE FUNCTION public.set_updated_at();

DROP TRIGGER IF EXISTS set_pgx_dosing_guidance_runs_updated_at ON pgx_dosing_guidance_runs;
CREATE TRIGGER set_pgx_dosing_guidance_runs_updated_at BEFORE UPDATE ON pgx_dosing_guidance_runs
FOR EACH ROW EXECUTE FUNCTION public.set_updated_at();

-- ============================================================================
-- 10. RLS POLICIES (New tables)
-- ============================================================================
ALTER TABLE patient_care_plan_runs ENABLE ROW LEVEL SECURITY;
ALTER TABLE pgx_dosing_guidance_runs ENABLE ROW LEVEL SECURITY;

CREATE POLICY "Users can view own care plan runs"
    ON patient_care_plan_runs FOR SELECT
    USING (auth.uid() = user_id);

CREATE POLICY "Users can insert own care plan runs"
    ON patient_care_plan_runs FOR INSERT
    WITH CHECK (auth.uid() = user_id);

CREATE POLICY "Users can view own PGx dosing runs"
    ON pgx_dosing_guidance_runs FOR SELECT
    USING (auth.uid() = user_id);

CREATE POLICY "Users can insert own PGx dosing runs"
    ON pgx_dosing_guidance_runs FOR INSERT
    WITH CHECK (auth.uid() = user_id);

-- Optional: allow updates for user's own rows (needed if you upsert / patch runs)
CREATE POLICY "Users can update own care plan runs"
    ON patient_care_plan_runs FOR UPDATE
    USING (auth.uid() = user_id);

CREATE POLICY "Users can update own PGx dosing runs"
    ON pgx_dosing_guidance_runs FOR UPDATE
    USING (auth.uid() = user_id);

-- ============================================================================
-- DONE! ✅
-- ============================================================================
-- After running this:
-- 1. Go to Authentication → Users → Create user (ak@ak.com / 786)
-- 2. Configure .env files with Supabase credentials
-- 3. Restart servers
-- 4. Test login!
-- ============================================================================
