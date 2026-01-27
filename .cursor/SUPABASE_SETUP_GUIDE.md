# 🔐 SUPABASE SETUP GUIDE - Production Ready

**Purpose:** Configure Supabase for authentication and database access

---

## 📋 STEP 1: Create Supabase Project

1. Go to https://supabase.com
2. Sign up / Log in
3. Click "New Project"
4. Fill in:
   - **Project Name:** `crispro-oncology` (or your choice)
   - **Database Password:** (save this!)
   - **Region:** Choose closest to you
5. Wait for project to be created (~2 minutes)

---

## 📋 STEP 2: Get Supabase Credentials

Once project is created:

1. Go to **Settings** → **API**
2. Copy these values:
   - **Project URL** (e.g., `https://xxxxx.supabase.co`)
   - **anon/public key** (starts with `eyJ...`)
   - **service_role key** (starts with `eyJ...`) - **KEEP SECRET!**

---

## 📋 STEP 3: Configure Frontend (.env)

**File:** `oncology-coPilot/oncology-frontend/.env`

```bash
# Supabase Configuration
VITE_SUPABASE_URL=https://youid.supabase.co
VITE_SUPABASE_ANON_KEY=your-anon-key-here
```

**Example:**
```bash
VITE_SUPABASE_URL=https://abcdefghijklmnop.supabase.co
VITE_SUPABASE_ANON_KEY=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6ImFiY2RlZmdoaWprbG1ub3AiLCJyb2xlIjoiYW5vbiIsImlhdCI6MTYzODU2Nzg5MCwiZXhwIjoxOTU0MTQzODkwfQ.example
```

---

## 📋 STEP 4: Configure Backend (.env)

**File:** `oncology-coPilot/oncology-backend-minimal/.env`

```bash
# Supabase Configuration
SUPABASE_URL=https://your-project-id.supabase.co
SUPABASE_ANON_KEY=your-anon-key-here
SUPABASE_SERVICE_KEY=your-service-role-key-here

# Database (if using direct connection)
DATABASE_URL=postgresql://postgres:[YOUR-PASSWORD]@db.[PROJECT-REF].supabase.co:5432/postgres
```

**Example:**
```bash
SUPABASE_URL=https://abcdefghijklmnop.supabase.co
SUPABASE_ANON_KEY=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6ImFiY2RlZmdoaWprbG1ub3AiLCJyb2xlIjoiYW5vbiIsImlhdCI6MTYzODU2Nzg5MCwiZXhwIjoxOTU0MTQzODkwfQ.example
SUPABASE_SERVICE_Y=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6ImFiY2RlZmdoaWprbG1ub3AiLCJyb2xlIjoic2VydmljZV9yb2xlIiwiaWF0IjoxNjM4NTY3ODkwLCJleHAiOjE5NTQxNDM4OTB9.example
```

---

## 📋 STEP 5: Run Database Migrations

**In Supabase Dashboard:**

1. Go to **SQL Editor**
2. Run each migration file in order:

**Migration 1: Complete Base Setup**
```sql
-- Copy and run ALL contents of: .cursor/SUPABASE_COMPLETE_SETUP.sql
-- This creates: user_profiles, user_quotas, patient_profiles, patient_sessions, patient_care_plans
-- Includes: RLS policies, triggers, indexes
```

**Migration 2: Schema Updates (NEW - Jan 2025)**
```sql
-- Copy and run ALL contents of: .cursor/SUPABASE_SCHEMA_UPDATES.sql
-- This adds: trial_matches, research_intelligence_queries, biomarker_history, care_plan_versions
-- Updates: patient_profiles (adds mechanism_vector column)
```

**Or use psql:**
```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main
psql $DATABASE_URL < .cursor/SUPABASE_COMPLETE_SETUP.sql
psql $DATABASE_URL < .cursor/SUPABASE_SCHEMA_UPDATES.sql
```

**Note:** Run `SUPABASE_COMPLETE_SETUP.sql` FIRST, then `SUPABASE_SCHEMA_UPDATES.sql`

---

## 📋 STEP 5B: What Was Added (Schema Updates - Jan 2025)

### 🆕 New Tables Added:

1. **`trial_matches`** - Stores clinical trial matching results
   - **Purpose:** Persist trial matches with scores, reasoning, and MoA vectors
   - **Key Fields:**
     - Matching scores (mechanism_fit, eligibility, combined)
     - 7D MoA vectors (DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux)
     - Reasoning objects (why_eligible, why_good_fit, etc.)
     - Freshness metadata (last_refreshed_at, stale flags)
     - User interaction (pinned, viewed_count)
   - **Use Case:** Store trial matches from `/api/ayesha/trials/search` endpoint

2. **`research_intelligence_queries`** - Stores research queries and results
   - **Purpose:** Persist Research Intelligence framework queries and synthesized findings
   - **Key Fields:**
     - Research questions and structured plans
     - Portal results (PubMed, GDC, Project Data Sphere)
     - LLM-synthesized findings (mechanisms, dosage, safety, outcomes)
     - MOAT pathway analysis
     - Article summaries and sub-question answers
     - LLM provider metadata (Cohere/Gemini)
   - **Use Case:** Store queries from Research Intelligence orchestrator

3. **`biomarker_history`** - Tracks biomarker values over time
   - **Purpose:** Historical tracking of CA-125, PSA, CEA, and other biomarkers
   - **Key Fields:**
     - Biomarker name, value, unit
     - Measurement date and treatment context
     - Burden assessment, response forecasts, resistance signals
   - **Use Case:** Track biomarker trends for monitoring and resistance detection

4. **`care_plan_versions`** - Tracks care plan history
   - **Purpose:** Version control and audit trail for care plans
   - **Key Fields:**
     - Version numbers (auto-incremented)
     - Full care plan JSON
     - Change summaries and trigger events
   - **Use Case:** Track care plan changes over time, enable rollback

### 🔄 Updates to Existing Tables:

- **`patient_profiles.mechanism_vector`** - Added JSONB column
  - Stores 7D mechanism vector: `{"ddr": 0.8, "mapk": 0.2, "pi3k": 0.0, ...}`
  - Used for mechanism-based trial matching

### 🔒 Security Features:

- **Row-Level Security (RLS)** enabled on all new tables
- **Data Classification** support (PHI/NON_PHI/SENSITIVE/PUBLIC)
- **User isolation** - users can only see their own data

### 📊 Indexes Added:

- Performance indexes on foreign keys, scores, dates, and classification fields
- Composite indexes for common query patterns

---


## 📋 STEP 5C: PGx + Advanced Care Plan Storage (Jan 2025)

### 🆕 New Extension Required:

1. **`pgcrypto` Extension** - Required for UUID generation
   - **Purpose:** Enables `gen_random_uuid()` function for primary keys
   - **Location:** Section 0 in `SUPABASE_COMPLETE_SETUP.sql`
   - **Note:** Must be created BEFORE any tables that use `gen_random_uuid()`

### 🆕 New Columns on `patient_profiles`:

1. **`germline_variants JSONB`** - Stores extraine variants
   - **Purpose:** Persist germline variants from VCF files or structured profiles
   - **Format:** `[{"gene": "DPYD", "variant": "*2A/*2A", "hgvs": "c.1905+1G>A", ...}, ...]`
   - **Use Case:** Input for PGx screening service

2. **`pgx_variants JSONB`** - Stores PGx-relevant variants only
   - **Purpose:** Filtered subset of germline variants in known pharmacogenes (DPYD, TPMT, UGT1A1, CYP2D6, CYP2C19)
   - **Format:** Same as `germline_variants`, but filtered to PGx genes
   - **Use Case:** Direct input for `DosingGuidanceService` and `PGxScreeningService`

3. **`current_medications JSONB`** - Stores patient's current medications
   - **Purpose:** Track medications for PGx drug-drug interaction checks
   - **Format:** `[{"drug": "5-fluorouracil", "dose": "500mg/m²", "frequency": "daily"}, ...]`
   - **Use Case:** Input for cumulative toxicity assessment

4. **`ca125_history JSONB`** - Stores CA-125 biomarker measurements over time
   - **Purpose:** Historical tracking for resistance detection and monitoring
   - **Format:** `[{"value": 2842.0, "unit": "U/mL", "date": "2025-01-01", "cycle": 0}, ...]`
   - **Use Case:** Input for CA-125 intelligence service and resistance prophet

### 🆕 New Tables:

1. **`patient_care_plan_runs`** - Stores complete care plan orchestrator runs
   - **Purpose:** Persist request/response/provenance from:
     - `POST /api/ayesha/complete_care_v2` (Ayesha-specific orchestrator)
     - `POST /api/complete_care/v2` (Universal orchestrator)
   - **Key Fields:**
     - `request_payload JSONB` - Full request (mutations, disease, tumor_context, etc.)
     - `response_payload JSONB` - Full response (wiwfm, trials, soc, ca125, resistance, pgx_screening_summary)
     - `provenance JSONB` - Run ID, profile, flags, timestamps
     - `orchestrator VARCHAR(50)` - "ayesha" or "universal"
     - `status VARCHAR(20)` - "success", "error", "partial"
     - `error_message TEXT` - Error details if status = "error"
   - **Use Case:** 
     - Reproduce care plaruns for debugging
     - Audit trail for clinical decision support
     - Analytics on care plan usage patterns

2. **`pgx_dosing_guidance_runs`** - Stores per-drug PGx dosing guidance runs
   - **Purpose:** Persist individual drug PGx screening results
   - **Key Fields:**
     - `drug_name VARCHAR(255)` - Drug being screened
     - `patient_pgx_variants JSONB` - PGx variants used for screening
     - `metabolizer_status VARCHAR(50)` - "Poor", "Intermediate", "Normal", "Rapid", "Ultrarapid"
     - `toxicity_tier VARCHAR(20)` - "AVOID", "REDUCE_DOSE", "MONITOR", "SAFE"
     - `dose_adjustment_factor FLOAT` - Multiplier for dose adjustment (e.g., 0.5 = 50% dose)
     - `recommendation TEXT` - Human-readable dosing guidance
     - `evidence_tier VARCHAR(20)` - "CPIC_STANDARD", "CLINVAR_BRIDGE", "INSUFFICIENT"
   - **Use Case:**
     - Track PGx screening history per patient
     - Audit PGx recommendations
     - Analytics on PGx variant frequencies

### 🔄 Updates to Existing Tables:

- **`user_quotas.pgxueries_lifetime`** - Added INT column (default 0)
  - Tracks total PGx queries across all patients
  - Used for quota gating and analytics

- **`user_quotas.care_plan_runs_lifetime`** - Added INT column (default 0)
  - Tracks total care plan runs across all patients
  - Used for quota gating and analytics

### 🔒 Security Features:

- **Row-Level Security (RLS)** enabled on all new tables
- **User isolation** - users can only see their own care plan runs
- **Patient isolation** - care plan runs linked to specific patient profiles

### 📊 Indexes Added:

- Performance indexes on foreign keys (`user_id`, `patient_id`)
- Indexes on `orchestrator` and `status` for filtering
- Indexes on `created_at` for time-based queries

### 🔗 Integration Points:

**Backend Services That Write to These Tables:**
- `api/services/pgx_extraction_service.py` → writes to `patient_profiles.pgx_variants
- `api/routers/ayesha_orchestrator_v2.py` → writes to `patient_care_plan_runs` (when integrated)
- `api/routers/completversal.py` → writes to `patient_care_plan_runs` (when integrated)
- `api/services/pgx_screening_service.py` → writes to `pgx_dosing_guidance_runs` (when integrated)

**Frontend Components That Read from These Tables:**
- `UniversalCompleteCare.jsx` → reads `patient_care_plan_runs` for care plan history
- `AyeshaTrialExplorer.jsx` → reads `patient_care_plan_runs` for Ayesha-specific care plans
- `SafetyGateCard.jsx` → displays PGx safety from `pgx_dosing_guidance_runs`


## 📋 STEP 6: Create Test User (ak@ak.com)

**Option A: Via Supabase Dashboard**
1. Go to **Authentication** → **Users**
2. Click "Add User" → "Create new user"
3. Email: `ak@ak.com`
4. Password: `786`
5. Auto  ✅ (check this)

**Option B: Via API (Backend)**
```bash
curl -X POST http://localhost:8000/api/auth/signup \
  -H "Content-Type: application/json" \
  -d '{
    "email": "ak@ak.com",
    "password": "786",
    "role": "patient",
    "full_name": "Test User"
  }'
```

---

## 📋 STEP 7: Verify Setup

### Frontend Check:
1. Restart frontend: `npm run dev`
2. Go to `http://localhost:5173/login`
3. Should NOT see "Supabase not configured"
4. Try logging in with `ak@ak.com` / `786`

### Backend Check:
```bash
# Test health endpoint
curl http://localhost:8000/health

# Test signup (if user doesn't exist)
curl -X POST http://localhost:8000/api/auth/signup \
  -H "Content-Type: application/json" \
  -d '{"email": "test@test.com", "password": "test123", "role": "patient"}'

# Test login
curl -X POST http://localhost:8000/api/auth/login \
  -H "Content-Type: application/json" \
  -d '{"email": "ak@ak.com", "password": "786"}'
```

---

## 📋 STEP 8: Test New Schema Tables

### Test 1: Verify Tables Exist

**In Supabase Dashboard:**
1. Go to **Table Editor**
2. Verify these tables exist:
   - ✅ `trial_matches`
   - ✅ `research_intelligence_queries`
   - ✅ `biomarker_history`
   - ✅ `care_plan_versions`
   - ✅ `patient_care_plan_runs` (NEW - Jan 2025)
   - ✅ `pgx_dosing_guidance_runs` (NEW - Jan 2025)
   - ✅ `patient_profiles.pgx_variants column (NEW - Jan 2025)
   - ✅ `patient_profiles.germline_variants` column (NEW - Jan 2025)
   - ✅ `patient_profiles.current_medications` column (NEW - Jan 2025)
   - ✅ `patient_profiles.ca125_history` column (NEW - Jan 2025)
   - ✅ `patient_profiles` (check for `mechanism_vector` column)

**Or via SQL:**
```sql
-- Check tables exist
SELECT table_name 
FROM information_schema.tables 
WHERE table_schema = 'public' 
AND table_name IN (
    'trial_matches',
    'research_intelligence_queries',
    'biomarker_history',
    'care_plan_versions'
);

-- Check mechanism_vector column exists
SELECT column_name, data_type 
FROM information_schema.columns 
WHERE table_name = 'patient_profiles' 
AND column_name = 'mechanism_vector';
```

### Test 2: Test RLS Policies

**In Supabase SQL Editor (as authenticated user):**
```sql
-- Test: Should only see your own trial matches
SELECT COUNT(*) FROM trial_matches WHERE user_id = auth.uid();

-- Test: Should NOT see other users' data
SELECT COUNT(*) FROM trial_matches WHERE user_id != auth.uid();
-- Should return 0 (RLS blocks cross-user access)
```

### Test 3: Test Trial Matches Insert

**Via Backend API (after backend integration):**
```bash
# Get auth token first
TOKEN=$(curl -X POST http://localhost:8000/api/auth/login \
  -H "Content-Type: application/json" \
  -d '{"email": "ak@ak.com", "password": "786"}' | jq -r '.access_token')

# Test trial match insert (when endpoint is ready)
curl -X POST http://localhost:8000/api/ayesha/trials/search \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "stage": "IVB",
    "treatment_line": "either",
    "germline_status": "negative",
    "location_state": "NY"
  }'
```

### Test 4: Test Biomarker History Insert

**Via SQL (manual test):**
```sql
-- Insert test biomarker record
INSERT INTO biomarker_history (
    user_id,
    patient_id,
    biomarker_name,
    biomarker_value,
    biomarker_unit,
    measurement_date
) VALUES (
    (SELECT id FROM user_profiles WHERE email = 'ak@ak.com'),
    'test-patient-001',
    'CA-125',
    2842.0,
    'U/mL',
    CURRENT_DATE
);

-- Verify insert
SELECT * FROM biomarker_history 
WHERE biomarker_name = 'CA-125' 
ORDER BY measurement_date DESC;
```

### Test 5: Test Care Plan Versions

**Via SQL (manual test):**
```sql
-- Insert test care plan version
INSERT INTO care_plan_versions (
    user_id,
    patient_id,
    care_plan_data,
    change_summary,
    trigger_event
) VALUES (
    (SELECT id FROM user_profiles WHERE email = 'ak@ak.com'),
    'test-patient-001',
    '{"trials": [], "soc": {}}'::jsonb,
    'Initial care plan creation',
    'system_init'
);

-- Verify version auto-increment
SELECT version_number, change_summary, created_at 
FROM care_plan_versions 
ORDER BY version_number DESC;
```

### Test 6: Test Mechanism Vector Column

**Via SQL:**
```sql
-- Update patient profile with mechanism vector
UPDATE patient_profiles 
SET mechanism_vector = '{"ddr": 0.8, "mapk": 0.2, "pi3k": 0.0, "vegf": 0.0, "her2": 0.0, "io": 0.1, "efflux": 0.0}'::jsonb
WHERE patient_id = 'test-patient-001';

-- Verify update
SELECT patient_id, mechanism_vector 
FROM patient_profiles 
WHERE mechanism_vector IS NOT NULL;
```

### Test 7: Test PGx Variants Storage

**Via SQL (manual test):**
```sql
-- Insert test PGx variant
UPDATE patient_profiles 
SET pgx_variants = '[{"gene": "DPYD", "variant": "*2A/*2A", "hgvs": "c.1905+1G>A", "chrom": "1", "pos": 97915614, "ref": "G", "alt": "A"}]'::jsonb
WHERE patient_id = 'test-patient-001';

-- Verify update
SELECT patient_id, pgx_variants 
FROM patient_profiles 
WHERE pgx_variants IS NOT NULL;
```

### Test 8: Test Care Plan Runs Storage

**Via SQL (manual test):**
```sql
-- Insert test care plan run
INSERT INTO patient_care_plan_runs (
    user_id,
    patient_id,
    orchestrator,
    request_payload,
    response_payload,
    provenance,
    status
) VALUES (
    (SELECT id FROM user_profiles WHERE email = 'ak@ak.com'),
    'test-patient-001',
    'ayesha',
    '{"mutations": [], "disease": "ovarian_cancer_hgs"}'::jsonb,
    '{"wiwfm": {"drugs": []}, "trials": []}'::jsonb,
    '{"run_id": "test-run-001", "profile": "baseline"}'::jsonb,
    'success'
);

-- Verify insert
SELECT orchestrator, status, created_at 
FROM patient_care_plan_runs 
ORDER BY created_at DESC;
```

### Test 9: Test PGx Dosing Guidance Runs Storage

**Via SQL (manual test):**
```sql
-- Insert test PGx dosing guidance run
INSERT INTO pgx_dosing_guidance_runs (
    user_id,
    patient_id,
    drug_name,
    patient_pgx_variants,
    metabolizer_status,
    toxicity_tier,
    dose_adjustment_factor,
    recommendation,
    evidence_tier
) VALUES (
    (SELECT id FROM user_profiles WHERE email = 'ak@ak.com'),
    'test-patient-001',
    '5-fluorouracil',
    '[{"gene": "DPYD", "variant": "*2A/*2A"}]'::jsonb,
    'Poor Metabolizer',
    'AVOID',
    0.0,
    'AVOID: DPYD *2A/*2A confers poor metabolizer status. Severe toxicity risk (grade 3-4).',
    'CPIC_STANDARD'
);

-- Verify insert
SELECT drug_name, toxicity_tier, metabolizer_status 
FROM pgx_dosing_guidance_runs 
ORDER BY created_at DESC;
```


---

## 📋 STEP 9: Backend Integration Plan

### Phase 1: Database Connection Setup
- [ ] Create Supabase client utility in backend
- [ ] Add connection pooling configuration
- [ ] Test connection with health check

### Phase 2: Trial Matches Integration
- [ ] Update `ayesha_trials.py` router to save matches to `trial_matches` table
- [ ] Add endpoint to retrieve saved matches: `GET /api/ayesha/trials/matches`
- [ ] Add endpoint to pin/unpin trials: `POST /api/ayesha/trials/pin`
- [ ] Add endpoint to refresh stale trials: `POST /api/ayesha/trials/refresh`

### Phase 3: Research Intelligence Integration
- [ ] Update Research Intelligence orchestrator to save queries to `research_intelligence_queries` table
- [ ] Add endpoint to retrieve query history: `GET /api/research/intelligence/history`
- [ ] Add endpoint to retrieve specific query result: `GET /api/research/intelligence/{query_id}`

### Phase 4: Biomarker History Integration
- [ ] Update CA-125 intelligence service to save to `biomarker_history` table
- [ ] Add endpoint to retrieve biomarker history: `GET /api/biomarkers/history`
- [ ] Add endpoint to add biomarker measurement: `POST /api/biomarkers/measurement`

### Phase 5: Care Plan Versions Integration
- [ ] Update care plan orchestrator to save versions to `care_plan_versions` table
- [ ] Add endpoint to retrieve care plan history: `GET /api/care-plans/versions`
- [ ] Add endpoint to rollback to previous version: `POST /api/care-plans/rollback`


### Phase 6: PGx + Advanced Care Plan Integration (NEW - Jan 2025)
- [ ] Update `pgx_extraction_service.py` to save extracted PGx variants to patient_profiles.pgx_variants
- [ ] Update `ayesha_orchestrator_v2.py` to save care plan runs to `patient_care_plan_runs` table
- [ ] Update `complete_care_universal.py` to save care plan runs to `patient_care_plan_runs` table
- [ ] Update `pgx_screening_service.py` to save dosing guidance runs to `pgx_dosing_guidance_runs` table
- [ ] Add endpoint to retrieve care plan history: `GET /api/care-plans/history?patient_id={id}`
- [ ] Add endpoint to retrieve PGx screening history: `GET /api/pgx/history?patient_id={id}`
- [ ] Update quota tracking to increment `pgx_queries_lifetime` and `care_plan_runs_lifetime`

### Testing Strategy:
1. **Unit Tests:** Test database operations in isolation
2. **Integration Tests:** Test API endpoints with real database
3. **E2E Tests:** Test full user flows (trial search → save → retrieve)
4. **RLS Tests:** Verify users can only access their own data

---

## 🔧 TROUBLESHOOTING

### "Supabase not configured" Error

**Frontend:**
- `.env` file exists in `oncology-frontend/`
- Check `VITE_SUPABASE_URL` and `VITE_SUPABASE_ANON_KEY` are set
- Restart dev server after changing `.env`
- Check browser console for errors

**Backend:**
- Check `.env` file exists in `oncology-backend-minimal/`
- Check `SUPABASE_URL` and `SUPABASE_ANON_KEY` are set
- Restart backend after changing `.env`

### "Invalid API key" Error
- Verify you copied the correct key (anon vs service_role)
- Check for extra spaces/newlines in `.env` file
- Ensure URL doesn't have trailing slash

### "User not found" Error
- Create user via Supabase Dashboard or signup endpoint
- Check user exists in Supabase → Authentication → Users

### "Database connection failed"
- Verify `DATABASE_URL` is correct
- Check password is correct (from project creation)
- Verify project is active in Supabase dashboard

---

## ✅ QUICK SETUP CHECKLIST

### Initial Setup:
- [ ] Supabase project created
- [ ] Credentials copied (URL, anon key, service key)
- [ ] Frontend `.env` configured (`VITE_SUPABASE_URL`, `VITE_SUPABASE_ANON_KEY`)
- [ ] Backend `.env` configured (`SUPABASE_URL`, `SUPABASE_ANON_KEY`, `SUPABASE_SERVICE_KEY`)

### Database Migrations:
- [ ] Base setup migration run (`SUPABASE_COMPLETE_SETUP.sql`)
  - [ ] Verify `pgcrypto` extension created (Section 0)
  - [ ] Verify Section 8 (PGX + Advanced Care Plan) tables created
- [ ] Schema updates migration run (`SUPABASE_SCHEMA_UPDATES.sql`)
- [ ] All tables verified in Supabase Table Editor

### User Setup:
- [ ] Test user created (`ak@ak.com` / `786`)
- [ ] User profile auto-created (via trigger)

### Verification:
- [ ] Frontend restarted
- [ ] Backend restarted
- [ ] Login works at `http://localhost:5173/login`
- [ ] RLS policies tested (can't see other users' data)
- [ ] New tables accessible (trial_matches, research_intelligence_queries, biomarker_history, care_plan_versions, patient_care_plan_runs, pgx_dosing_guidance_runs)

### Backend Integration (Next Steps):
- [ ] Supabase client configured in backend
- [ ] Trial matches saving implemented
- [ ] Research Intelligence queries saving implemented
- [ ] Biomarker history tracking implemented
- [ ] Care plan versions tracking implemented

---

## 🧪 TESTING PLAN SUMMARY

### What to Test:

1. **Schema Validation:**
   - ✅ All 4 new tables exist
   - ✅ `patient_profiles.mechanism_vector` column exists
   - ✅ All indexes created
   - ✅ RLS policies enabled

2. **Security (RLS):**
   - ✅ Users can only see their own data
   - ✅ Cross-user access blocked
   - ✅ Data classification respected

3. **Functionality:**
   - ✅ Trial matches can be inserted
   - ✅ Biomarker history can be tracked
   - ✅ Care plan versions auto-increment
   - ✅ Mechanism vectors can be stored

4. **Backend Integration (Future):**
   - ✅ API endpoints save to new tables
   - ✅ API endpoints retrieve from new tables
   - ✅ Data persists across sessions

### Test Scripts Location:
- **SQL Tests:** Run in Supabase SQL Editor
- **API Tests:** Use `curl` commands (see Step 8)
- **E2E Tests:** Create in `oncology-coPilot/oncology-backend-minimal/tests/`

---

**Ready to test!** 🚀

**Next:** Implement backend integration (see Step 9 for detailed plan)
