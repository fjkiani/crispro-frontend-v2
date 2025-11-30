# Modal Service URLs - Deployment Configuration

**Date**: January 20, 2025  
**Status**: ✅ Services Deployed

---

## 🔗 Service URLs

### SAE Service
- **URL**: `https://crispro--sae-service-saeservice-api.modal.run`
- **Purpose**: Extract SAE features from Evo2 layer 26 activations
- **Model**: evo2_7b (4096-dim activations)
- **SAE Weights**: Trained weights (4096×32768) from Goodfire/Evo-2-Layer-26-Mixed

### Evo2 7B Service
- **URL**: `https://crispro--evo-service-evoservice7b-api-7b.modal.run`
- **Purpose**: Evo2 7B model for variant scoring and activation extraction
- **Model**: evo2_7b_base

---

## 🔧 Environment Variables

**Backend Configuration**:
```bash
export SAE_SERVICE_URL="https://crispro--sae-service-saeservice-api.modal.run"
export EVO_URL_7B="https://crispro--evo-service-evoservice7b-api-7b.modal.run"
export ENABLE_EVO2_SAE=1
export ENABLE_TRUE_SAE=1
```

**Or in `.env` file**:
```
SAE_SERVICE_URL=https://crispro--sae-service-saeservice-api.modal.run
EVO_URL_7B=https://crispro--evo-service-evoservice7b-api-7b.modal.run
ENABLE_EVO2_SAE=1
ENABLE_TRUE_SAE=1
```

---

## 🧪 Health Check Endpoints

**SAE Service**:
```bash
curl https://crispro--sae-service-saeservice-api.modal.run/health
```

**Evo2 7B Service**:
```bash
curl https://crispro--evo-service-evoservice7b-api-7b.modal.run/health
```

**Backend Health Check** (via router):
```bash
curl http://localhost:8000/api/sae/health
```

---

## ⚠️ Cost Control

**User Requirement**: "keep it tight where we always test to make sure we dont burn through credits"

**Mitigation**:
1. ✅ Feature flags (`ENABLE_TRUE_SAE`) gate SAE extraction
2. ✅ Circuit breaker implemented (max 50 errors per 100 requests)
3. ✅ Test with small batches first (1-2 patients)
4. ✅ Monitor Modal dashboard for costs
5. ⏸️ Wait for user approval before full cohort extraction

**Testing Strategy**:
- Light health checks only (no expensive operations)
- Small batch tests (1-2 patients) before full extraction
- User approval required for full cohort re-extraction

---

## 📝 Notes

- **Health checks may fail on cold start** (services need to wake up)
- **Services are deployed** (user confirmed)
- **Ready for testing** when user approves

