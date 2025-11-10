# 🚀 Metastasis Assessment - Quick Start Guide

## Start Backend

```bash
cd oncology-coPilot/oncology-backend-minimal
PYTHONPATH=. venv/bin/uvicorn api.main:app --host 127.0.0.1 --port 8000
```

## Start Frontend

```bash
cd oncology-coPilot/oncology-frontend
npm run dev
# Navigate to: http://localhost:5173/vus-explorer
```

## Test Backend

```bash
# Health check
curl http://localhost:8000/api/metastasis/health

# Test assessment
curl -X POST http://localhost:8000/api/metastasis/assess \
  -H "Content-Type: application/json" \
  -d '{"mutations":[{"gene":"BRAF","hgvs_p":"p.Val600Glu","chrom":"7","pos":140453136,"ref":"T","alt":"A"}]}'
```

## Run Tests

```bash
cd oncology-backend-minimal
PYTHONPATH=. venv/bin/pytest tests/metastasis/ -v
# Expected: 15 passed in ~50s
```

## Files Modified/Created

**Backend:**
- ✅ `api/config/metastasis_rules.json`
- ✅ `api/schemas/metastasis.py`
- ✅ `api/services/metastasis_service.py`
- ✅ `api/routers/metastasis.py`
- ✅ `api/main.py` (updated)
- ✅ `tests/metastasis/` (3 test files)

**Frontend:**
- ✅ `src/hooks/useMetastasis.js`
- ✅ `src/components/metastasis/MetastasisReport.jsx`
- ✅ `src/components/vus/AnalysisResults.jsx` (updated)

## Key Settings

- **Default Model:** `evo2_1b`
- **Timeout:** 60s per insight
- **Retry:** 1-2 attempts with exponential backoff
- **Cohort Lifts:** Capped at ≤0.05 per step (stubbed in v1)
- **Cache TTL:** 10 minutes (frontend)

## Troubleshooting

**Backend won't start:**
```bash
# Check imports
cd oncology-backend-minimal
venv/bin/python -c "from api.main import app; print('OK')"
```

**Tests fail:**
```bash
# Ensure PYTHONPATH is set
cd oncology-backend-minimal
PYTHONPATH=. venv/bin/pytest tests/metastasis/ -v
```

**Frontend error:**
```bash
# Check API_ROOT in .env.local
echo "VITE_API_ROOT=http://localhost:8000" > .env.local
```

## Support

- Implementation: `.cursor/rules/use-cases/metastatic-intervention.md`
- Summary: `.cursor/rules/use-cases/METASTASIS_COMPLETE.md`
- This guide: `.cursor/rules/use-cases/METASTASIS_QUICKSTART.md`

