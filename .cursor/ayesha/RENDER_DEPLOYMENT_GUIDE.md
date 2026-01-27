# 🚀 RENDER DEPLOYMENT GUIDE - Backend

**Date**: January 29, 2025  
**Service**: Oncology CoPilot Backend (FastAPI)  
**Platform**: Render.com

---

## 📋 QUICK START

### **Option 1: Auto-Deploy with render.yaml (Recommended)**

1. **Push to GitHub**: Ensure your code is in a GitHub repository
2. **Connect to Render**: 
   - Go to https://dashboard.render.com
   - Click "New +" → "Blueprints"
   - Select your repository
   - Render will auto-detect `render.yaml` and create the service
3. **Set Environment Variables**: Add required env vars in Render dashboard
4. **Deploy**: Render will automatically build and deploy

---

### **Option 2: Manual Setup**

1. **Create Web Service**:
   - Go to https://dashboard.render.com
   - Click "New +" → "Web Service"
   - Connect your GitHub repository
   - Select the `oncology-coPilot/oncology-backend-minimal` directory

2. **Configure Service**:
   ```
   Name: oncology-backend
   Environment: Python 3
   Build Command: pip install -r requirements.txt
   Start Command: uvicorn main:app --host 0.0.0.0 --port $PORT
   ```

3. **Set Environment Variables** (see below)

4. **Deploy**: Click "Create Web Service"

---

## 🔧 REQUIRED ENVIRONMENT VARIABLES

Set these in Render dashboard → Your Service → Environment tab:

### **Essential (Required)**

```bash
# Python Version
PYTHON_VERSION=3.11.0

# Environment
ENVIRONMENT=production
```

### **Database (Supabase) - If using patient features**

```bash
# Supabase Connection
SUPABASE_URL=https://your-project.supabase.co
SUPABASE_KEY=your-anon-key
SUPABASE_SERVICE_ROLE_KEY=your-service-role-key
```

### **Optional API Keys**

```bash
# Cohere (for embeddings/vector search)
COHERE_API_KEY=your-cohere-api-key

# OpenAI (if used by services)
OPENAI_API_KEY=your-openai-api-key

# AstraDB (if using vector database)
ASTRA_DB_APPLICATION_TOKEN=your-astra-token
ASTRA_DB_API_ENDPOINT=https://your-db-id-us-east1.apps.astra.datastax.com
```

---

## 📁 FILES CREATED

✅ **`render.yaml`** - Render deployment configuration  
✅ **`Procfile`** - Already exists (used by Render)

---

## 🔍 VERIFICATION

After deployment, check:

1. **Health Endpoint**: `https://your-service.onrender.com/health`
2. **API Docs**: `https://your-service.onrender.com/docs`
3. **Root Endpoint**: `https://your-service.onrender.com/`

---

## ⚙️ RENDER SETTINGS

### **Recommended Settings**:

- **Plan**: Starter ($7/month) or Free (limited)
- **Region**: Oregon (or closest to your users)
- **Auto-Deploy**: Enable (deploys on git push)
- **Health Check**: `/api/health` (already configured in render.yaml)

### **Scaling**:

- **Free Plan**: Sleeps after 15 min inactivity, slow cold starts
- **Starter Plan**: Always on, faster response times

---

## 🐛 TROUBLESHOOTING

### **Build Fails**:
- Check `requirements.txt` exists
- Verify Python version matches (3.11.0)
- Check build logs in Render dashboard

### **Service Won't Start**:
- Verify `startCommand` is correct: `uvicorn main:app --host 0.0.0.0 --port $PORT`
- Check logs in Render dashboard → Logs tab
- Ensure `main.py` exists and `app` is defined

### **Health Check Failing**:
- Verify `/health` endpoint exists
- Check health check path in Render settings (should be `/health`)
- Review service logs for errors

### **Timeouts**:
- Free plan has timeout limits
- Consider upgrading to Starter plan for better performance
- Optimize slow endpoints (add caching, async processing)

---

## 🔗 POST-DEPLOYMENT

After successful deployment:

1. **Update Frontend API URL**:
   ```javascript
   // In frontend .env or config
   VITE_API_ROOT=https://your-service.onrender.com
   ```

2. **Test Endpoints**:
   ```bash
   curl https://your-service.onrender.com/health
   curl https://your-service.onrender.com/api/ayesha/complete_care_v2 \
     -X POST -H "Content-Type: application/json" -d '{}'
   ```

3. **Monitor Logs**: Render dashboard → Logs tab

---

## 📝 NOTES

- **Free Plan**: Services sleep after inactivity (15 min). First request after sleep will be slow (~30 seconds)
- **Starter Plan ($7/month)**: Always-on, no sleep, faster responses
- **Auto-Deploy**: Enabled by default. Pushes to main branch trigger automatic deployments
- **Environment Variables**: Sensitive values should be set in Render dashboard, not committed to git

---

**Last Updated**: January 29, 2025  
**Status**: ✅ **READY FOR DEPLOYMENT**
