# 🚀 CrisPro CRISPR Assistant Deployment Guide

This guide will help you deploy the CrisPro CRISPR Assistant backend and frontend as separate, independent services.

## **🏗️ Architecture Overview**

```
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│   Frontend      │    │    Backend      │    │  External       │
│   (Vercel)      │◄──►│   (Vercel)      │◄──►│  Services       │
│                 │    │                 │    │                 │
└─────────────────┘    └─────────────────┘    └─────────────────┘
```

**🎯 Recommended Setup:** Both services on Vercel for unified deployment experience

## **📁 Project Structure**

```
crispr-assistant-main/
├── src/                          # Main backend services
│   ├── services/
│   │   ├── command_center/       # Main orchestration service ⭐
│   │   ├── oracle/               # AI prediction service
│   │   ├── forge/                # Therapeutic generation
│   │   ├── boltz/                # Structural prediction
│   │   └── ...                   # Other specialized services
│   └── app/                      # Streamlit applications
├── oncology-coPilot/             # Oncology-specific frontend
│   ├── oncology-frontend/        # React/Vite frontend
│   └── oncology-backend/         # Legacy backend (deprecated)
└── tools/                        # Utility scripts and tools
```

## **🔧 Prerequisites**

- **Git** - for version control
- **Python 3.11+** - for backend development
- **Node.js 18+** - for frontend development
- **Vercel CLI** - for both frontend and backend deployment

## **📋 Step-by-Step Deployment**

### **Option 1: Unified Vercel Deployment (Recommended)**

#### **Step 1: Main Backend Deployment (Command Center)**

1. **Navigate to the main backend directory:**
   ```bash
   cd src/services/command_center
   ```

2. **Install Vercel CLI:**
   ```bash
   npm install -g vercel
   ```

3. **Deploy to Vercel:**
   ```bash
   chmod +x deploy-vercel.sh
   ./deploy-vercel.sh
   ```

4. **Note the Vercel URL** - you'll need this for the frontend configuration

#### **Step 2: Frontend Deployment (Oncology CoPilot)**

1. **Navigate to frontend directory:**
   ```bash
   cd oncology-coPilot/oncology-frontend
   ```

2. **Configure environment:**
   ```bash
   cp env.example .env
   # Edit .env with your backend Vercel URL
   ```

3. **Deploy to Vercel:**
   ```bash
   chmod +x deploy.sh
   ./deploy.sh
   ```

### **Option 2: Hybrid Deployment (Modal Backend + Vercel Frontend)**

#### **Step 1: Backend Deployment (Modal)**

1. **Navigate to the main backend directory:**
   ```bash
   cd src/services/command_center
   ```

2. **Install Modal CLI:**
   ```bash
   pip install modal
   ```

3. **Deploy to Modal:**
   ```bash
   modal deploy main.py
   ```

#### **Step 2: Frontend Deployment (Vercel)**

1. **Navigate to frontend directory:**
   ```bash
   cd oncology-coPilot/oncology-frontend
   ```

2. **Configure environment:**
   ```bash
   cp env.example .env
   # Edit .env with your backend Modal URL
   ```

3. **Deploy to Vercel:**
   ```bash
   chmod +x deploy.sh
   ./deploy.sh
   ```

## **🔗 Service Configuration**

### **Main Backend Environment Variables**

The main backend service (`src/services/command_center`) uses environment variables from the project root. Key ones include:

```bash
# API Keys
GOOGLE_API_KEY=your_actual_key
OPENAI_API_KEY=your_actual_key

# External Service URLs
UNIFIED_ORACLE_URL=your_oracle_service_url
BOLTZ_SERVICE_URL=your_boltz_service_url

# CORS (Update with your Vercel domain)
ALLOWED_ORIGINS=https://your-app.vercel.app
```

### **Frontend Environment Variables (.env)**

```bash
# Backend API (Update with your backend URL)
VITE_BACKEND_API_URL=https://your-backend-url.vercel.app
VITE_WEBSOCKET_URL=wss://your-backend-url.vercel.app/ws

# External Services
VITE_FUSION_ENGINE_URL=https://crispro--fusion-engine-v1-fusionengine-api.modal.run
VITE_EVO_SERVICE_URL=your_evo_service_url
VITE_ZETA_ORACLE_URL=your_zeta_oracle_url
```

## **🚀 Deployment Commands**

### **Quick Deploy Main Backend (Vercel)**
```bash
cd src/services/command_center
vercel --prod
```

### **Quick Deploy Main Backend (Modal)**
```bash
cd src/services/command_center
modal deploy main.py
```

### **Quick Deploy Frontend (Vercel)**
```bash
cd oncology-coPilot/oncology-frontend
vercel --prod
```

## **🔍 Post-Deployment Verification**

### **Main Backend Health Check (Vercel)**
```bash
curl https://your-backend-url.vercel.app/
```

### **Main Backend Health Check (Modal)**
```bash
curl https://your-backend-modal-url.modal.run/
```

### **Frontend Health Check (Vercel)**
```bash
curl https://your-frontend-url.vercel.app/
```

## **📊 Monitoring & Logs**

### **Backend Logs (Vercel)**
```bash
vercel logs
```

### **Backend Logs (Modal)**
```bash
modal logs command-center
```

### **Frontend Logs (Vercel)**
```bash
vercel logs
```

## **🔄 Updating Services**

### **Main Backend Updates (Vercel)**
1. Make your code changes
2. Run: `vercel --prod`
3. Vercel will automatically update the service

### **Main Backend Updates (Modal)**
1. Make your code changes
2. Run: `modal deploy main.py`
3. Modal will automatically update the service

### **Frontend Updates (Vercel)**
1. Make your code changes
2. Run: `vercel --prod`
3. Vercel will automatically update the service

## **🚨 Troubleshooting**

### **Common Backend Issues (Vercel)**
- **Import errors:** Check `vercel_main.py` path configuration
- **Environment variables:** Ensure all required env vars are set
- **Dependencies:** Verify all packages are in `requirements-vercel.txt`
- **Timeout errors:** Vercel has 10-second timeout limit for serverless functions

### **Common Backend Issues (Modal)**
- **Import errors:** Check `main.py` path configuration
- **Environment variables:** Ensure all required env vars are set
- **Dependencies:** Verify all packages are in `requirements.txt`

### **Common Frontend Issues**
- **Build failures:** Check for syntax errors and missing dependencies
- **Environment variables:** Ensure all `VITE_*` variables are set
- **API calls:** Verify backend URL is correct and accessible

### **CORS Issues**
If you encounter CORS errors:
1. Update `ALLOWED_ORIGINS` in backend configuration
2. Ensure frontend domain is included
3. Redeploy backend

## **🔐 Security Considerations**

- **Environment variables:** Never commit `.env` files to git
- **API keys:** Use environment variables for all sensitive data
- **CORS:** Restrict `ALLOWED_ORIGINS` to production domains only
- **HTTPS:** Both Vercel and Modal provide HTTPS by default

## **💰 Cost Optimization**

### **Vercel (Both Services)**
- **Backend:** Pay-per-use for serverless functions
- **Frontend:** Generous free tier + pay-per-use for bandwidth
- **Benefits:** Unified billing, better integration, global CDN

### **Modal (Backend Only)**
- Pay-per-use pricing
- Auto-scaling based on demand
- Better for long-running AI/ML workloads

## **⚡ Performance Considerations**

### **Vercel Backend Limitations**
- **10-second timeout** for serverless functions
- **Cold starts** for infrequently used endpoints
- **Memory limits** for large operations

### **Modal Backend Advantages**
- **No timeout limits** for long-running operations
- **Persistent memory** between requests
- **Better for AI/ML workloads**

## **🎯 Service-Specific Notes**

### **Command Center Service**
- **Location:** `src/services/command_center/`
- **Main File:** `main.py` (Modal) / `vercel_main.py` (Vercel)
- **Purpose:** Main orchestration service for all workflows
- **Key Endpoints:** `/workflow/execute`, `/workflow/assess_threat`

### **Oracle Service**
- **Location:** `src/services/oracle/`
- **Purpose:** AI prediction and scoring
- **Deployment:** Can be deployed separately or as part of main service

### **Forge Service**
- **Location:** `src/services/forge/`
- **Purpose:** Therapeutic generation and design
- **Deployment:** Can be deployed separately or as part of main service

## **📞 Support**

- **Vercel Documentation:** https://vercel.com/docs
- **Modal Documentation:** https://modal.com/docs
- **Project Issues:** Check the main repository issues

---

**🎯 Deployment Status:** Ready for deployment
**📅 Last Updated:** 2024
**👥 Maintainer:** Development Team
**💡 Recommendation:** Start with unified Vercel deployment for simplicity
**⭐ Main Service:** Command Center (`src/services/command_center/`) 