# 🔐 SUPABASE COMPLETE SETUP GUIDE - Single Source of Truth

**Date**: January 2025  
**Status**: ✅ **COMPLETE SETUP GUIDE**  
**Purpose**: Single source of truth for all Supabase setup, configuration, and schema management

---

## 📊 EXECUTIVE SUMMARY

This guide consolidates all Supabase setup documentation into one comprehensive resource covering:
- **Quick Setup** (5 minutes) - Get login working immediately
- **Complete Setup** (30 minutes) - Full database schema and configuration
- **Schema Updates** - All migrations and updates
- **Environment Configuration** - Frontend and backend `.env` setup
- **Troubleshooting** - Common issues and solutions

---

## ⚡ PART 1: QUICK SETUP (5 MINUTES) - GET LOGIN WORKING NOW

### Step 1: Get Supabase Credentials

1. **Go to**: https://supabase.com/dashboard
2. **Create/Select Project**
3. **Go to**: Settings → API
4. **Copy these 3 values:**
   - **Project URL** (looks like: `https://xxxxx.supabase.co`)
   - **anon public key** (long string starting with `eyJ...`)
   - **service_role key** (long string starting with `eyJ...`) - **KEEP SECRET!**

### Step 2: Add to Frontend .env

**File**: `oncology-coPilot/oncology-frontend/.env`

```bash
VITE_SUPABASE_URL=https://your-project-id.supabase.co
VITE_SUPABASE_ANON_KEY=your-anon-key-here
```

### Step 3: Add to Backend .env

**File**: `oncology-coPilot/oncology-backend-minimal/.env`

```bash
SUPABASE_URL=https://your-project-id.supabase.co
SUPABASE_ANON_KEY=your-anon-key-here
SUPABASE_SERVICE_KEY=your-service-role-key-here
```

### Step 4: Restart Servers

```bash
# Terminal 1 - Backend
cd oncology-coPilot/oncology-backend-minimal
uvicorn api.main:app --reload --host 0.0.0.0 --port 8000

# Terminal 2 - Frontend  
cd oncology-coPilot/oncology-frontend
npm run dev
```

### Step 5: Create Test User

**Option A: Via Supabase Dashboard (EASIEST)**
1. Go to: **Authentication** → **Users** → **Add User**
2. Email: `ak@ak.com`
3. Password: `786`
4. ✅ Check "Auto Confirm User"

**Option B: Via Backend API**
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

### Step 6: Test Login

1. Go to: `http://localhost:5173/login`
2. Should NOT see "Supabase not configured" anymore
3. Enter: `ak@ak.com` / `786`
4. Click "Sign In"

---

## 🚀 PART 2: COMPLETE DATABASE SETUP (30 MINUTES)

### Step 1: Create Supabase Project

1. Go to https://supabase.com
2. Sign up / Log in
3. Click "New Project"
4. Fill in:
   - **Project Name**: `crispro-oncology` (or your choice)
   - **Database Password**: ⚠️ **SAVE THIS!** You'll need it for `DATABASE_URL`
   - **Region**: Choose closest to you
5. Wait 2-3 minutes for project creation

### Step 2: Run Base Schema Migration

**In Supabase Dashboard:**
1. Go to **SQL Editor**
2. Click **"New query"**
3. **Copy and run ALL contents of**: `.cursor/SUPABASE_COMPLETE_SETUP.sql`

**This creates:**
- `user_profiles` table (with MFA columns, data classification)
- `user_quotas` table (with PGx and care plan quotas)
- `patient_profiles` table (with germline/PGx variants, current medications, CA-125 history)
- `patient_sessions` table
- `patient_care_plans` table
- `patient_care_plan_runs` table (stores complete care plan orchestrator runs)
- `pgx_dosing_guidance_runs` table (stores per-drug PGx screening results)
- RLS policies (Row-Level Security)
- Triggers (auto-create user profile on signup, updated_at timestamps)

### Step 3: Run Schema Updates Migration

**Still in SQL Editor:**
1. Click **"New query"**
2. **Copy and run ALL contents of**: `.cursor/SUPABASE_SCHEMA_UPDATES.sql`

**This adds:**
- `trial_matches` table (trial matching results with scores and reasoning)
- `research_intelligence_queries` table (research queries and results)
- `biomarker_history` table (biomarker values over time)
- `care_plan_versions` table (care plan history and versioning)
- `patient_profiles.mechanism_vector` column (7D mechanism vector)

### Step 4: Run MOAT Schema Updates

**Still in SQL Editor:**
1. Click **"New query"**
2. **Copy and run ALL contents of**: `.cursor/SUPABASE_SCHEMA_UPDATES_MOAT.sql`

**This adds:**
- `orchestrator_pipeline_state` table (7-phase pipeline execution state)
- `agent_executions` table (individual agent execution history)
- `sae_features` table (SAE features and diagnostics)
- `resistance_playbook` table (resistance playbook recommendations)
- `file_uploads` table (VCF/PDF/MAF file uploads)
- `toxicity_risk_assessments` table (toxicity risk assessments)
- `dosing_guidance` table (dosing guidance recommendations)

### Step 5: Run Research Intelligence Schema

**Still in SQL Editor:**
1. Click **"New query"**
2. **Copy and run ALL contents of**: `.cursor/SUPABASE_RESEARCH_INTELLIGENCE_SCHEMA.sql`

**This adds:**
- `research_intelligence_queries` table (if not already exists)
- `research_intelligence_dossiers` table (dossier content and sharing)

---

## 📋 PART 3: ENVIRONMENT CONFIGURATION

### Frontend .env

**File**: `oncology-coPilot/oncology-frontend/.env`

```bash
# Existing vars
VITE_API_ROOT=http://127.0.0.1:8000
VITE_WS_ROOT=ws://localhost:8000
GEMINI_API_KEY=AIzaSyDmPm3J2yqzJD1nXvd_5-8i6TX6rygwZ0Y

# Supabase Configuration (REQUIRED for authentication)
VITE_SUPABASE_URL=https://your-project-id.supabase.co
VITE_SUPABASE_ANON_KEY=your-anon-key-here
```

**Example with real values:**
```bash
VITE_SUPABASE_URL=https://abcdefghijklmnop.supabase.co
VITE_SUPABASE_ANON_KEY=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6ImFiY2RlZmdoaWprbG1ub3AiLCJyb2xlIjoiYW5vbiIsImlhdCI6MTYzODU2Nzg5MCwiZXhwIjoxOTU0MTQzODkwfQ.example
```

### Backend .env

**File**: `oncology-coPilot/oncology-backend-minimal/.env`

```bash
# Supabase Configuration (REQUIRED for authentication)
SUPABASE_URL=https://your-project-id.supabase.co
SUPABASE_ANON_KEY=your-anon-key-here
SUPABASE_SERVICE_KEY=your-service-role-key-here

# Database Connection (Optional, for direct DB access)
DATABASE_URL=postgresql://postgres:YourPassword@db.xxxxx.supabase.co:5432/postgres
```

**Example with real values:**
```bash
SUPABASE_URL=https://abcdefghijklmnop.supabase.co
SUPABASE_ANON_KEY=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6ImFiY2RlZmdoaWprbG1ub3AiLCJyb2xlIjoiYW5vbiIsImlhdCI6MTYzODU2Nzg5MCwiZXhwIjoxOTU0MTQzODkwfQ.example
SUPABASE_SERVICE_KEY=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6ImFiY2RlZmdoaWprbG1ub3AiLCJyb2xlIjoic2VydmljZV9yb2xlIiwiaWF0IjoxNjM4NTY3ODkwLCJleHAiOjE5NTQxNDM4OTB9.example
DATABASE_URL=postgresql://postgres:MySecurePass123!@db.abcdefghijklmnop.supabase.co:5432/postgres
```

### After Adding Credentials

1. **Restart frontend**: `npm run dev` (in `oncology-coPilot/oncology-frontend/`)
2. **Restart backend**: `uvicorn api.main:app --reload` (in `oncology-coPilot/oncology-backend-minimal/`)
3. **Test**: Go to http://localhost:5173/login - should NOT see "Supabase not configured"

---

## 🗄️ PART 4: DATABASE SCHEMA OVERVIEW

### Core Tables

1. **`user_profiles`**
   - User authentication and profile data
   - MFA columns: `mfa_enabled`, `mfa_secret`, `mfa_verified_at`, `mfa_backup_codes`
   - Data classification: `data_classification` (PHI/NON_PHI/SENSITIVE/PUBLIC)

2. **`user_quotas`**
   - User tier and quota limits
   - PGx quotas: `pgx_queries_limit`, `pgx_queries_used`
   - Care plan quotas: `complete_care_limit`, `complete_care_used`

3. **`patient_profiles`**
   - Patient medical data
   - PGx data: `germline_variants`, `pgx_variants`, `current_medications`, `ca125_history`
   - Mechanism vector: `mechanism_vector` (7D: DDR, MAPK, PI3K, VEGF, HER2, IO, Efflux)

4. **`patient_care_plan_runs`**
   - Complete care plan orchestrator runs
   - Stores: `request_payload`, `response_payload`, `provenance`
   - Tracks: `orchestrator` (ayesha/universal), `status`, `error_message`

5. **`pgx_dosing_guidance_runs`**
   - Per-drug PGx screening results
   - Stores: `gene`, `variant`, `drug`, `toxicity_tier`, `adjustment_factor`, `recommendation`

### Feature-Specific Tables

6. **`trial_matches`** - Clinical trial matching results
7. **`research_intelligence_queries`** - Research queries and results
8. **`research_intelligence_dossiers`** - Dossier content and sharing
9. **`biomarker_history`** - Biomarker values over time
10. **`care_plan_versions`** - Care plan history and versioning

### MOAT Orchestrator Tables

11. **`orchestrator_pipeline_state`** - 7-phase pipeline execution state
12. **`agent_executions`** - Individual agent execution history
13. **`sae_features`** - SAE features and diagnostics
14. **`resistance_playbook`** - Resistance playbook recommendations
15. **`file_uploads`** - VCF/PDF/MAF file uploads
16. **`toxicity_risk_assessments`** - Toxicity risk assessments
17. **`dosing_guidance`** - Dosing guidance recommendations

---

## 🧪 PART 5: VERIFICATION & TESTING

### Test 1: Verify Tables Exist

**In Supabase Dashboard:**
1. Go to **Table Editor**
2. Verify these tables exist:
   - ✅ `user_profiles`
   - ✅ `user_quotas`
   - ✅ `patient_profiles`
   - ✅ `patient_care_plan_runs`
   - ✅ `pgx_dosing_guidance_runs`
   - ✅ `trial_matches`
   - ✅ `research_intelligence_queries`
   - ✅ `biomarker_history`
   - ✅ `care_plan_versions`
   - ✅ `orchestrator_pipeline_state`
   - ✅ `agent_executions`
   - ✅ `sae_features`
   - ✅ `resistance_playbook`
   - ✅ `file_uploads`
   - ✅ `toxicity_risk_assessments`
   - ✅ `dosing_guidance`

### Test 2: Test Login

1. Go to: `http://localhost:5173/login`
2. Enter: `ak@ak.com` / `786`
3. Click "Sign In"
4. ✅ Should redirect to home/profile page

### Test 3: Backend Health Check

```bash
curl http://localhost:8000/health
# Should return: {"status": "healthy", "services": "operational"}
```

### Test 4: Verify User Profile Created

**In Supabase Dashboard:**
1. Go to **Table Editor** → **`user_profiles`**
2. ✅ Should see a row with `email = 'ak@ak.com'`
3. ✅ Should have columns: `id`, `email`, `role`, `data_classification`, `mfa_enabled`, etc.

---

## 🚨 TROUBLESHOOTING

### "Supabase not configured" Error

**Frontend:**
- ✅ Check `.env` file exists in `oncology-frontend/`
- ✅ Check variables start with `VITE_` (required for Vite)
- ✅ Restart frontend dev server (env vars only load on startup)
- ✅ Check browser console for errors

**Backend:**
- ✅ Check `.env` file exists in `oncology-backend-minimal/`
- ✅ Check `SUPABASE_URL` and `SUPABASE_ANON_KEY` are set
- ✅ Restart backend after changing `.env`

### "Invalid API key" Error

- ✅ Verify you copied the **anon key** (not service_role) for frontend
- ✅ Check for extra spaces/newlines in `.env` file
- ✅ Ensure URL doesn't have trailing slash

### "User not found" or Login Fails

- ✅ Create user via Supabase Dashboard (Authentication → Users → Add User)
- ✅ Check "Auto Confirm User" was checked
- ✅ Verify user exists: Supabase Dashboard → Authentication → Users
- ✅ Or use signup endpoint: `POST /api/auth/signup`

### Migration Errors

- ✅ If table doesn't exist, run `SUPABASE_COMPLETE_SETUP.sql` first
- ✅ If column already exists, migration will skip it (IF NOT EXISTS)
- ✅ Check SQL Editor for error messages

### Backend Can't Connect to Supabase

- ✅ Verify `SUPABASE_URL` and `SUPABASE_ANON_KEY` are set in backend `.env`
- ✅ Restart backend server
- ✅ Check backend logs for connection errors

---

## 📋 QUICK REFERENCE

**Supabase Dashboard**: https://supabase.com/dashboard  
**SQL Editor**: Dashboard → SQL Editor  
**Authentication**: Dashboard → Authentication → Users  
**Table Editor**: Dashboard → Table Editor  
**API Settings**: Dashboard → Settings → API

**Files to Edit:**
- `oncology-coPilot/oncology-frontend/.env`
- `oncology-coPilot/oncology-backend-minimal/.env`

**SQL Migration Files (Run in Order):**
1. `.cursor/SUPABASE_COMPLETE_SETUP.sql` (Base schema)
2. `.cursor/SUPABASE_SCHEMA_UPDATES.sql` (Trial matches, research intelligence, biomarkers)
3. `.cursor/SUPABASE_SCHEMA_UPDATES_MOAT.sql` (MOAT orchestrator tables)
4. `.cursor/SUPABASE_RESEARCH_INTELLIGENCE_SCHEMA.sql` (Research intelligence dossiers)

**Test Endpoints:**
- Frontend: `http://localhost:5173/login`
- Backend health: `http://localhost:8000/health`
- Backend login: `POST http://localhost:8000/api/auth/login`

---

## ✅ COMPLETE CHECKLIST

### Initial Setup:
- [ ] Supabase project created
- [ ] Credentials copied (URL, anon key, service key, database URL)
- [ ] Frontend `.env` configured (`VITE_SUPABASE_URL`, `VITE_SUPABASE_ANON_KEY`)
- [ ] Backend `.env` configured (`SUPABASE_URL`, `SUPABASE_ANON_KEY`, `SUPABASE_SERVICE_KEY`)

### Database Migrations:
- [ ] Base setup migration run (`SUPABASE_COMPLETE_SETUP.sql`)
- [ ] Schema updates migration run (`SUPABASE_SCHEMA_UPDATES.sql`)
- [ ] MOAT schema updates run (`SUPABASE_SCHEMA_UPDATES_MOAT.sql`)
- [ ] Research intelligence schema run (`SUPABASE_RESEARCH_INTELLIGENCE_SCHEMA.sql`)
- [ ] All tables verified in Supabase Table Editor

### User Setup:
- [ ] Test user created (`ak@ak.com` / `786`)
- [ ] User profile auto-created (via trigger)

### Verification:
- [ ] Frontend restarted
- [ ] Backend restarted
- [ ] Login works at `http://localhost:5173/login`
- [ ] RLS policies tested (can't see other users' data)
- [ ] All new tables accessible

---

**Last Updated**: January 2025  
**Status**: ✅ **COMPLETE SETUP GUIDE - SINGLE SOURCE OF TRUTH**
