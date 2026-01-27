# ⚡ QUICK SUPABASE SETUP - Get Login Working NOW

## 🎯 IMMEDIATE FIX (5 minutes)

### Step 1: Get Supabase Credentials

1. **Go to:** https://supabase.com/dashboard
2. **Create/Select Project**
3. **Go to:** Settings → API
4. **Copy these 3 values:**
   - **Project URL** (looks like: `https://xxxxx.supabase.co`)
   - **anon public key** (long string starting with `eyJ...`)
   - **service_role key** (long string starting with `eyJ...`) - **KEEP SECRET!**

### Step 2: Add to Frontend .env

**File:** `oncology-coPilot/oncology-frontend/.env`

Add these lines (replace with YOUR values):
```bash
VITE_SUPABASE_URL=https://your-project-id.supabase.co
VITE_SUPABASE_ANON_KEY=your-anon-key-here
```

### Step 3: Add to Backend .env

**File:** `oncology-coPilot/oncology-backend-minimal/.env`

Add these lines (replace with YOUR values):
```bash
SUPABASE_URL=https://your-project-id.supabase.co
SUPANON_KEY=your-anon-key-here
SUPABASE_SERVICE_KEY=your-service-role-key-here
```

### Step 4: Restart Servers

```bash
# Terminal 1 - Backend
cd oncology-coPilot/oncology-backend-minimal
# Stop current server (Ctrl+C)
uvicorn api.main:app --reload --host 0.0.0.0 --port 8000

# Terminal 2 - Frontend  
cd oncology-coPilot/oncology-frontend
# Stop current server (Ctrl+C)
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

##IFICATION

### Check Frontend .env
```bash
cd oncology-coPilot/oncology-frontend
cat .env | grep VITE_SUPABASE
# Should show:
# VITE_SUPABASE_URL=https://...
# VITE_SUPABASE_ANON_KEY=eyJ...
```

### Check Backend .env
```bash
cd oncology-coPilot/oncology-backend-minimal
cat .env | grep SUPABASE
# Should show:
# SUPABASE_URL=https://...
# SUPABASE_ANON_KEY=eyJ...
# SUPABASE_SERVICE_KEY=eyJ...
```

### Test Backend Connection
```bash
curl http://localhost:8000/health
# Should return: {"status": "healthy", "services": "operational"}
```

---

## 🚨 TROUBLESHOOTING

### Still seeing "Supabase not configured"?
1. ✅ Check `.env` file exists in `oncology-frontend/`
2. ✅ Check variables start with `VITE_` (required for Vite)
3. ✅ Restart frontend dev server (env vars only load on startup)
4. ✅ Check browser console for errors

### "Invalid API key" error?
- Verify you copied the **anon key** (not service_role) for frontend
- Check for extra spaces/newlines in `.env` file
- Ensure URL doesn't have trailing# "User not found" or login fails?
- Create user via Supabase Dashboard (Option A above)
- Or use signup endpoint (Option B above)
- Check user exists: Supabase Dashboard → Authentication → Users

---

## 📋 WHERE TO GET CREDENTIALS

1. **Supabase Dashboard:** https://supabase.com/dashboard
2. **Select your project**
3. **Settings** (gear icon) → **API**
4. **Copy:**
   - **Project URL** (under "Project URL")
   - **anon public** key (under "Project API keys" → "anon public")
   - **service_role** key (under "Project API keys" → "service_role") - **SECRET!**

---

**After adding credentials and restarting, login should work!** 🚀
