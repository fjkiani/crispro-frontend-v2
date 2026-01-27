# 🚀 STARTUP GUIDE - Testing Ayesha Patient Journey

## Step 1: Start Backend Server

```bash
cd oncology-coPilot/oncology-backend-minimal

# Option 1: Using Makefile (if venv exists)
make backend

# Option 2: Direct uvicorn (if you have venv activated)
uvicorn api.main:app --reload --host 0.0.0.0 --port 8000

# Option 3: Using Python directly
python -m uvicorn api.main:app --reload --host 0.0.0.0 --port 8000
```

**Expected Output:**
```
INFO:     Uvicorn running on http://0.0.0.0:8000 (Press CTRL+C to quit)
INFO:     Started reloader process
INFO:     Started server process
INFO:     Waiting for application startup.
INFO:     Application startup complete.
```

**Verify Backend:**
```bash
curl http://localhost:8000/health
# Should return: {"status": "ok"}
```

---

## Step 2: Start Frontend Server

**Open a NEW terminal window/tab:**

```bash
cd oncology-coPilot/oncology-frontend

# Install dependencies (if not already done)
npm install

# Start dev server
m run dev
```

**Expected Output:**
```
  VITE v5.x.x  ready in xxx ms

  ➜  Local:   http://localhost:5173/
  ➜  Network: use --host to expose
```

**Frontend will be at:** `http://localhost:5173`

---

## Step 3: Create User Account (if needed)

### Option A: Use Signup Endpoint

```bash
curl -X POST http://localhost:8000/api/auth/signup \
  -H "Content-Type: application/json" \
  -d '{
    "email": "ak@ak.com",
    "password": "786",
    "role": "patient"
  }'
```

### Option B: Check if User Exists

```bash
# Check if user already exists (requires admin or direct DB access)
# If user doesn't exist, use signup endpoint above
```

---

## Step 4: Login Flow

1. **Open Browser:** `http://localhost:5173`

2. **Navigate to Login:**
   - Click "Login" in navigation
   - OR go directly to: `http://localhost:5173/login`

3. **Enter Credentials:**
   - Email: `ak@ak.com`
   - Password: `786`

4. **After Login:**
   - You'll be redirected to home or profile
   - Check browser console for any errors

---

## S5: Navigate Patient Pages

### Available Patient Routes:

1. **Patient Profile**
   - URL: `http://localhost:5173/patient/profile`
   - What it shows: Patient profile, medical history, genomic data

2. **Ayesha Complete Care Plan**
   - URL: `http://localhost:5173/ayesha-complete-care`
   - What it shows: Complete care plan with all insights

3. **Ayesha Trial Explorer**
   - URL: `http://localhost:5173/ayesha-trials`
   - What it shows: Clinical trials, resistance playbook, SAE features

4. **Ayesha Dossiers**
   - URL: `http://localhost:5173/ayesha-dossiers`
   - What it shows: Clinical trial dossiers browser

5. **Patient Onboarding** (if new patient)
   - URL: `http://localhost:5173/patient/onboarding`
   - What it shows: Onboarding flow for new patients

---

## Step 6: Testing Checklist

### Backend Health Checks:
- [ ] `curl http://localhost:8000/health` returns `{"status": "ok"}`
- [ ] `curl http://localhost:8000/api/health` works
- [ ] Backend logs show no errors

### Frontend Health Checks:
- [ ] Frontend loads at `http://localhost:5173`
- [ ] No console errors in browser
- [ ] Login page loads correctly

### Authentication:
- [ ] Signup endpoint works (if creating new user)
- [ ] Login endpoint works
- [ ] JWT token is received after login
- [ ] Token is stored in localStorage/sessionStorage

### Patient Pages:
- [ ] Can access `/patient/profile` (requires login)
- [ ] Can access `/ayesha-complete-care` (requires login)
- [ ] Can access `/ayesha-trials` (requires login)
- [ ] API calls to backend succeed (check Network tab)

---

## Troubleshooting

### Backend Won't Start:
```bash
# Check if port 8000 is in use
lsof -i :8000

# Kill process if needed
kill -9 <PID>

# Check Python environment
which python
python --version

# Check if dependencies installed
pip list | grep fastapi
pip list | grep uvicorn
```

### Frontend Won't Start:
```bash
# Check if port 5173 is in use
lsof -i :5173

# Clear node_modules and reinstall
rm -rf node_modules package-lock.json
npm install

# Check Node version
node --version  # Should be v16+
```

### Login Issues:
- Check browser console for errors
- Check Network tab for failed API calls
- Verify backend is running
- Check CORS settings in backend
- Verify user exists in database

### API Calls Failing:
- Check backend logs for errors
- Verify API endpoint URLs in frontend
- Check CORS configuration
- Verify authentication token is sent

---

## Quick Test Commands

```bash
# Test backend health
curl http://localhost:8000/health

# Test signup (if needed)
curl -X POST http://localhost:8000/api/auth/signup \
  -H "Content-Type: application/json" \
  -d '{"email": "ak@ak.com", "password": "786", "role": "patient"}'

# Test login
curl -X POST http://localhost:8000/api/auth/login \
  -H "Content-Type: application/json" \
  -d '{"email": "ak@ak.com", "password": "786"}'
```

---

## Expected Journey Flow

1. **Start Backend** → `http://localhost:8000`
2. **Start Frontend** → `http://localhost:5173`
3. **Go to Login** → `http://localhost:5173/login`
4. **Enter Creden* → `ak@ak.com` / `786`
5. **Login Success** → Redirected to home/profile
6. **Navigate to Patient Pages:**
   - `/patient/profile` - View profile
   - `/ayesha-complete-care` - Complete care plan
   - `/ayesha-trials` - Trial explorer
   - `/ayesha-dossiers` - Dossier browser

---

**Ready to test!** 🚀
