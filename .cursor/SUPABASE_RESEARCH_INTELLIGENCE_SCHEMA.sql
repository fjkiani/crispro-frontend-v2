-- ============================================================================
-- RESEARCH INTELLIGENCE SCHEMA
-- ============================================================================
-- Purpose: Persist Research Intelligence queries and dossiers
-- Run this in: Supabase Dashboard → SQL Editor → New Query
-- ============================================================================

-- ============================================================================
-- 1. RESEARCH INTELLIGENCE QUERIES TABLE
-- ============================================================================
CREATE TABLE IF NOT EXISTS public.research_intelligence_queries (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE,
    session_id UUID,  -- Optional: Link to patient_sessions if needed (no FK to avoid dependency)
    
    -- Query Input
    question TEXT NOT NULL,
    context JSONB NOT NULL,  -- {disease, treatment_line, biomarkers}
    options JSONB,  -- {portals, synthesize, run_moat_analysis, persona}
    
    -- Results (Full API response)
    result JSONB NOT NULL,
    provenance JSONB,  -- {run_id, methods_used, timestamp}
    
    -- Metadata
    dossier_id UUID,  -- Link to generated dossier
    persona VARCHAR(20) DEFAULT 'patient' CHECK (persona IN ('patient', 'doctor', 'r&d')),
    
    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    last_accessed_at TIMESTAMPTZ DEFAULT NOW()
);

-- Indexes for performance
CREATE INDEX IF NOT EXISTS idx_ri_queries_user_id ON public.research_intelligence_queries(user_id);
CREATE INDEX IF NOT EXISTS idx_ri_queries_session_id ON public.research_intelligence_queries(session_id) WHERE session_id IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_ri_queries_created_at ON public.research_intelligence_queries(created_at DESC);
CREATE INDEX IF NOT EXISTS idx_ri_queries_question_search ON public.research_intelligence_queries USING gin(to_tsvector('english', question));
CREATE INDEX IF NOT EXISTS idx_ri_queries_persona ON public.research_intelligence_queries(persona);
CREATE INDEX IF NOT EXISTS idx_ri_queries_dossier_id ON public.research_intelligence_queries(dossier_id);

-- ============================================================================
-- 2. RESEARCH INTELLIGENCE DOSSIERS TABLE
-- ============================================================================
CREATE TABLE IF NOT EXISTS public.research_intelligence_dossiers (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    query_id UUID REFERENCES public.research_intelligence_queries(id) ON DELETE CASCADE NOT NULL,
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE,
    
    -- Dossier Content
    persona VARCHAR(20) NOT NULL CHECK (persona IN ('patient', 'doctor', 'r&d')),
    markdown TEXT NOT NULL,  -- Full markdown content
    pdf_path TEXT,  -- Path to generated PDF (if generated)
    
    -- Sharing
    shareable_link TEXT UNIQUE,  -- Unique shareable URL (e.g., uuid)
    shareable_expires_at TIMESTAMPTZ,  -- Optional expiration
    is_public BOOLEAN DEFAULT false,
    
    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- Indexes for performance
CREATE INDEX IF NOT EXISTS idx_ri_dossiers_query_id ON public.research_intelligence_dossiers(query_id);
CREATE INDEX IF NOT EXISTS idx_ri_dossiers_user_id ON public.research_intelligence_dossiers(user_id);
CREATE INDEX IF NOT EXISTS idx_ri_dossiers_shareable_link ON public.research_intelligence_dossiers(shareable_link);
CREATE INDEX IF NOT EXISTS idx_ri_dossiers_persona ON public.research_intelligence_dossiers(persona);

-- ============================================================================
-- 3. ROW-LEVEL SECURITY (RLS) POLICIES
-- ============================================================================

-- Enable RLS
ALTER TABLE public.research_intelligence_queries ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.research_intelligence_dossiers ENABLE ROW LEVEL SECURITY;

-- Drop existing policies if they exist (for idempotency)
DROP POLICY IF EXISTS "Users can view own queries" ON public.research_intelligence_queries;
DROP POLICY IF EXISTS "Users can insert own queries" ON public.research_intelligence_queries;
DROP POLICY IF EXISTS "Users can update own queries" ON public.research_intelligence_queries;
DROP POLICY IF EXISTS "Users can view own dossiers" ON public.research_intelligence_dossiers;
DROP POLICY IF EXISTS "Users can insert own dossiers" ON public.research_intelligence_dossiers;
DROP POLICY IF EXISTS "Users can update own dossiers" ON public.research_intelligence_dossiers;

-- Policy: Users can only see their own queries
CREATE POLICY "Users can view own queries" 
    ON public.research_intelligence_queries
    FOR SELECT 
    USING (auth.uid() = user_id);

CREATE POLICY "Users can insert own queries" 
    ON public.research_intelligence_queries
    FOR INSERT 
    WITH CHECK (auth.uid() = user_id);

CREATE POLICY "Users can update own queries" 
    ON public.research_intelligence_queries
    FOR UPDATE 
    USING (auth.uid() = user_id);

-- Policy: Users can only see their own dossiers (unless public)
CREATE POLICY "Users can view own dossiers" 
    ON public.research_intelligence_dossiers
    FOR SELECT 
    USING (auth.uid() = user_id OR is_public = true);

CREATE POLICY "Users can insert own dossiers" 
    ON public.research_intelligence_dossiers
    FOR INSERT 
    WITH CHECK (auth.uid() = user_id);

CREATE POLICY "Users can update own dossiers" 
    ON public.research_intelligence_dossiers
    FOR UPDATE 
    USING (auth.uid() = user_id);

-- ============================================================================
-- 4. HELPER FUNCTIONS
-- ============================================================================

-- Function to update updated_at timestamp
CREATE OR REPLACE FUNCTION update_updated_at_column()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = NOW();
    RETURN NEW;
END;
$$ language 'plpgsql';

-- Drop existing triggers if they exist (for idempotency)
-- Using DO block to safely drop triggers only if tables exist
DO $$
BEGIN
    IF EXISTS (SELECT 1 FROM information_schema.tables WHERE table_schema = 'public' AND table_name = 'research_intelligence_queries') THEN
        DROP TRIGGER IF EXISTS update_ri_queries_updated_at ON public.research_intelligence_queries;
    END IF;
    
    IF EXISTS (SELECT 1 FROM information_schema.tables WHERE table_schema = 'public' AND table_name = 'research_intelligence_dossiers') THEN
        DROP TRIGGER IF EXISTS update_ri_dossiers_updated_at ON public.research_intelligence_dossiers;
    END IF;
END $$;

-- Trigger for research_intelligence_queries
CREATE TRIGGER update_ri_queries_updated_at 
    BEFORE UPDATE ON public.research_intelligence_queries
    FOR EACH ROW 
    EXECUTE FUNCTION update_updated_at_column();

-- Trigger for research_intelligence_dossiers
CREATE TRIGGER update_ri_dossiers_updated_at 
    BEFORE UPDATE ON public.research_intelligence_dossiers
    FOR EACH ROW 
    EXECUTE FUNCTION update_updated_at_column();

-- ============================================================================
-- 5. VERIFICATION QUERIES
-- ============================================================================

-- Test insert (replace 'test-user-id' with actual user_id)
-- INSERT INTO public.research_intelligence_queries (user_id, question, context, result)
-- VALUES (
--     'test-user-id',
--     'How does curcumin help with cancer?',
--     '{"disease": "breast_cancer"}',
--     '{"test": "result"}'
-- );

-- Test select
-- SELECT * FROM public.research_intelligence_queries WHERE user_id = 'test-user-id';

