# Test LLM Food Enhancement Service

## ✅ Service is Ready

The LLM Food Enhancement Service has been fixed and is ready to test.

## Quick Test Command

```bash
cd /Users/fahadkiani/Desktop/development/crispr-assistant-main/oncology-coPilot/oncology-backend-minimal
python3 test_env_and_llm.py
```

Or use the simpler test:
```bash
python3 simple_llm_test.py
```

## What to Expect

The test will:

1. **Load .env file** from `/Users/fahadkiani/Desktop/development/crispr-assistant-main/.env`
2. **Check API key** - Verify `GEMINI_API_KEY` is loaded
3. **Initialize service** - Create the LLM enhancement service
4. **Test LLM call** - Generate a personalized rationale for Vitamin D in ovarian cancer

## Expected Output

```
✅ Loaded .env from: /Users/fahadkiani/Desktop/development/crispr-assistant-main/.env
🔑 GEMINI_API_KEY set: True
   Key length: XX chars

================================================================================
Testing LLM Service
================================================================================

✅ Service initialized
   Model: gemini-2.5-pro
   Provider: gemini
   LLM Available: True

🧪 Testing personalized rationale generation...
   (This may take 10-30 seconds)

================================================================================
📝 PERSONALIZED RATIONALE:
================================================================================
Vitamin D is recommended for High-Grade Serous Ovarian Cancer during 
first-line chemotherapy (L1) for patients with HRD+ (Homologous Recombination 
Deficiency). This compound targets DNA repair pathways that are compromised in 
HRD+ patients, supporting homologous recombination repair mechanisms critical 
during platinum-based chemotherapy. The evidence is moderate (15 studies, 2 RCTs), 
with treatment line appropriateness of 90%, making it well-suited for first-line 
therapy where DNA repair support is most critical.
================================================================================

✅ Test completed successfully!
```

## If You See Errors

### "LLM not available"
- Check that `GEMINI_API_KEY` is in your `.env` file
- Verify the .env file path is correct
- Make sure `python-dotenv` is installed: `pip install python-dotenv`

### "Model not found"
- The model is set to `gemini-2.5-pro`
- If that doesn't work, try changing to `gemini-1.5-pro` in the service file

### Import errors
- Make sure you're in the `oncology-backend-minimal` directory
- Check that `src/tools/llm_api.py` exists

## Service Features

The service provides 4 LLM-enhanced features:

1. **Personalized Rationale** - Treatment line + biomarker specific
2. **Mechanism Synthesis** - Beyond keyword matching
3. **Evidence Interpretation** - Treatment line context
4. **Patient Recommendations** - Timing, monitoring, safety

## Files Created

- `test_env_and_llm.py` - Full test with .env loading
- `simple_llm_test.py` - Simplified test
- `run_llm_test.py` - Alternative test script

All tests are in `oncology-coPilot/oncology-backend-minimal/`

## Next Steps

Once the test works:
1. Integrate into `validate_food_dynamic` endpoint
2. Add LLM-enhanced fields to response
3. Test with real patient scenarios
