# 📝 SUPABASE ENV TEMPLATE

## Quick Start: Get Credentials

### Step 1: Create Project (if needed)
1. Go to: **https://supabase.com/dashboard**
2. Click **"New Project"** (green button)
3. Fill in:
   - **Project Name:** `crispro-oncology` (or your choice)
   - **Database Password:** ⚠️ **SAVE THIS!** You'll need it for `DATABASE_URL`
   - **Region:** Choose closest to you (e.g., `US East (N. Virginia)`)
4. Wait 2-3 minutes for project creation

### Step 2: Get API Credentials
1. In project dashboard, click **Settings** (gear icon, bottom left)
2. Click **API** (under "Project Settings")
3. Copy these values:

#### a) Project URL
```
https://xfhiwodulrbbtfcqneqt.supabase.co
```
📋 **Copy this** → `VITE_SUPABASE_URL` / `SUPABASE_URL`

#### b) API Keys
- **anon public** (long JWT token starting with `eyJhbGci...`)
  📋 **Copy this** → `VITE_SUPABASE_ANON_KEY` / `SUPABASE_ANON_KEY`
  
- **service_role** (long JWT token starting with `eyJhbGci...`)
  📋 **Copy this** → `SUPABASE_SERVICE_KEY` (backend only, keep secret!)

#### c) Database Connection String
1. Still in Settings, click **Database** (left sidebar)
2. Scroll to **"Connection string"** section
3. Select **"URI"** tab
4. Copy connection string (looks like: `postgresql://postgres:[YOUR-PASSWORD]@db.xxxxx.supabase.co:5432/postgres`)
5. **Replace `[YOUR-PASSWORD]`** with password from Step 1
6. 📋 **Copy this** → `DATABASE_URL`

---

## Frontend .env
**Location:** `oncology-coPilot/oncology-frontend/.env`

**Edit the file and replace the placeholder values:**

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

**✅ Check when done:** File saved with actual values (not placeholders)

---

## Backend .env
**Location:** `oncology-coPilot/oncology-backend-minimal/.env`

**Add these lines (or edit existing):**

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

**✅ Check when done:** File saved with actual values (not placeholders)

---

## After Adding

1. **Restart frontend:** `npm run dev` (in `oncology-coPilot/oncology-frontend/`)
2. **Restart backend:** `uvicorn api.main:app --reload` (in `oncology-coPilot/oncology-backend-minimal/`)
3. **Test:** Go to http://localhost:5173/login - should NOT see "Supabase not configured"

---

## Quick Reference

**Dashboard:** https://supabase.com/dashboard  
**API Settings:** Dashboard → Settings → API  
**Database Settings:** Dashboard → Settings → Database  
**All 4 values needed:** Project URL, anon key, service_role key, database URL
