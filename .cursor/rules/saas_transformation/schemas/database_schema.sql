-- ⚔️ SAAS TRANSFORMATION DATABASE SCHEMA
-- PostgreSQL (Supabase) Schema for User Management, Feature Flags, Quotas, Sessions

-- ============================================================================
-- AUTHENTICATION (Managed by Supabase)
-- ============================================================================
-- Note: auth.users table is managed by Supabase Auth
-- We only create custom tables in public schema

-- ============================================================================
-- USER PROFILES (Custom Metadata)
-- ============================================================================

CREATE TABLE IF NOT EXISTS public.user_profiles (
    id UUID PRIMARY KEY REFERENCES auth.users(id) ON DELETE CASCADE,
    email VARCHAR(255) NOT NULL,
    full_name VARCHAR(255),
    institution VARCHAR(255),
    role VARCHAR(50) CHECK (role IN ('researcher', 'clinician', 'admin', 'enterprise')) DEFAULT 'researcher',
    tier VARCHAR(50) CHECK (tier IN ('free', 'pro', 'enterprise')) DEFAULT 'free',
    avatar_url TEXT,
    bio TEXT,
    country VARCHAR(100),
    timezone VARCHAR(50),
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    UNIQUE(email)
);

CREATE INDEX idx_user_profiles_email ON public.user_profiles(email);
CREATE INDEX idx_user_profiles_tier ON public.user_profiles(tier);

-- ============================================================================
-- USER SUBSCRIPTIONS
-- ============================================================================

CREATE TABLE IF NOT EXISTS public.user_subscriptions (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE,
    tier VARCHAR(50) NOT NULL CHECK (tier IN ('free', 'pro', 'enterprise')),
    status VARCHAR(50) CHECK (status IN ('active', 'canceled', 'past_due', 'trialing')) DEFAULT 'active',
    stripe_customer_id VARCHAR(255),
    stripe_subscription_id VARCHAR(255),
    current_period_start TIMESTAMPTZ,
    current_period_end TIMESTAMPTZ,
    trial_end TIMESTAMPTZ,
    cancel_at TIMESTAMPTZ,
    canceled_at TIMESTAMPTZ,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_user_subscriptions_user_id ON public.user_subscriptions(user_id);
CREATE INDEX idx_user_subscriptions_status ON public.user_subscriptions(status);

-- ============================================================================
-- USAGE QUOTAS
-- ============================================================================

CREATE TABLE IF NOT EXISTS public.user_quotas (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE,
    tier VARCHAR(50) NOT NULL,
    -- Quotas
    variant_analyses_limit INT DEFAULT 10,
    variant_analyses_used INT DEFAULT 0,
    drug_queries_limit INT DEFAULT 5,
    drug_queries_used INT DEFAULT 0,
    food_queries_limit INT DEFAULT 3,
    food_queries_used INT DEFAULT 0,
    clinical_trials_limit INT DEFAULT 0,
    clinical_trials_used INT DEFAULT 0,
    -- Reset tracking
    period_start TIMESTAMPTZ DEFAULT NOW(),
    period_end TIMESTAMPTZ DEFAULT NOW() + INTERVAL '1 month',
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    UNIQUE(user_id)
);

CREATE INDEX idx_user_quotas_user_id ON public.user_quotas(user_id);
CREATE INDEX idx_user_quotas_period_end ON public.user_quotas(period_end);

-- ============================================================================
-- FEATURE FLAGS
-- ============================================================================

CREATE TABLE IF NOT EXISTS public.user_feature_flags (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE,
    feature_name VARCHAR(100) NOT NULL,
    enabled BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    UNIQUE(user_id, feature_name)
);

CREATE INDEX idx_user_feature_flags_user_id ON public.user_feature_flags(user_id);
CREATE INDEX idx_user_feature_flags_feature ON public.user_feature_flags(feature_name);

-- Available Features Registry
CREATE TABLE IF NOT EXISTS public.features (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name VARCHAR(100) UNIQUE NOT NULL,
    display_name VARCHAR(255),
    description TEXT,
    tier_required VARCHAR(50) CHECK (tier_required IN ('free', 'pro', 'enterprise')),
    enabled_by_default BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_features_tier ON public.features(tier_required);

-- Insert default features
INSERT INTO public.features (name, display_name, tier_required, enabled_by_default) VALUES
('variant_analysis', 'Variant Analysis', 'free', TRUE),
('drug_efficacy', 'Drug Efficacy Prediction', 'free', TRUE),
('food_validator', 'Food/Supplement Validator', 'free', TRUE),
('sae_features', 'SAE Features (Treatment Line Intelligence)', 'pro', FALSE),
('clinical_trials', 'Clinical Trials Matching', 'pro', FALSE),
('fusion_engine', 'Fusion Engine (AlphaMissense)', 'pro', FALSE),
('cohort_lab', 'Cohort Lab & Benchmarking', 'enterprise', FALSE),
('crispr_design', 'CRISPR Guide Design', 'enterprise', FALSE),
('pdf_export', 'PDF Export', 'pro', FALSE),
('api_access', 'Programmatic API Access', 'enterprise', FALSE)
ON CONFLICT (name) DO NOTHING;

-- ============================================================================
-- SESSION & ANALYSIS PERSISTENCE
-- ============================================================================

CREATE TABLE IF NOT EXISTS public.user_sessions (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE,
    session_name VARCHAR(255),
    patient_context JSONB,  -- Disease, biomarkers, treatment history
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    last_accessed_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_user_sessions_user_id ON public.user_sessions(user_id);
CREATE INDEX idx_user_sessions_last_accessed ON public.user_sessions(last_accessed_at);

-- Saved Analyses
CREATE TABLE IF NOT EXISTS public.saved_analyses (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE,
    session_id UUID REFERENCES public.user_sessions(id) ON DELETE CASCADE,
    analysis_type VARCHAR(50) CHECK (analysis_type IN (
        'variant_analysis', 
        'drug_efficacy', 
        'food_validator', 
        'clinical_trials',
        'crispr_design'
    )),
    input_data JSONB,  -- Original request
    output_data JSONB,  -- API response
    provenance JSONB,  -- run_id, model, flags, timestamp
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_saved_analyses_user_id ON public.saved_analyses(user_id);
CREATE INDEX idx_saved_analyses_session_id ON public.saved_analyses(session_id);
CREATE INDEX idx_saved_analyses_type ON public.saved_analyses(analysis_type);
CREATE INDEX idx_saved_analyses_created_at ON public.saved_analyses(created_at);

-- ============================================================================
-- USAGE LOGS (for billing & analytics)
-- ============================================================================

CREATE TABLE IF NOT EXISTS public.usage_logs (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE,
    endpoint VARCHAR(255),
    analysis_type VARCHAR(50),
    request_data JSONB,
    response_time_ms INT,
    status_code INT,
    error_message TEXT,
    ip_address INET,
    user_agent TEXT,
    created_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_usage_logs_user_id ON public.usage_logs(user_id);
CREATE INDEX idx_usage_logs_created_at ON public.usage_logs(created_at);
CREATE INDEX idx_usage_logs_endpoint ON public.usage_logs(endpoint);

-- ============================================================================
-- ROW LEVEL SECURITY (RLS) POLICIES
-- ============================================================================

-- Enable RLS on all tables
ALTER TABLE public.user_profiles ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.user_subscriptions ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.user_quotas ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.user_feature_flags ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.user_sessions ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.saved_analyses ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.usage_logs ENABLE ROW LEVEL SECURITY;

-- Users can only read/write their own data
CREATE POLICY "Users can view own profile" ON public.user_profiles
    FOR SELECT USING (auth.uid() = id);

CREATE POLICY "Users can update own profile" ON public.user_profiles
    FOR UPDATE USING (auth.uid() = id);

-- Similar policies for other tables...

-- ============================================================================
-- TRIGGERS
-- ============================================================================

-- Auto-create user profile on signup
CREATE OR REPLACE FUNCTION public.handle_new_user()
RETURNS TRIGGER AS $$
BEGIN
    INSERT INTO public.user_profiles (id, email, tier)
    VALUES (NEW.id, NEW.email, 'free');
    
    INSERT INTO public.user_quotas (user_id, tier)
    VALUES (NEW.id, 'free');
    
    RETURN NEW;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

CREATE TRIGGER on_auth_user_created
    AFTER INSERT ON auth.users
    FOR EACH ROW EXECUTE FUNCTION public.handle_new_user();

-- Auto-update updated_at timestamp
CREATE OR REPLACE FUNCTION public.update_updated_at_column()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = NOW();
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER update_user_profiles_updated_at
    BEFORE UPDATE ON public.user_profiles
    FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();

CREATE TRIGGER update_user_subscriptions_updated_at
    BEFORE UPDATE ON public.user_subscriptions
    FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();

CREATE TRIGGER update_user_quotas_updated_at
    BEFORE UPDATE ON public.user_quotas
    FOR EACH ROW EXECUTE FUNCTION public.update_updated_at_column();








