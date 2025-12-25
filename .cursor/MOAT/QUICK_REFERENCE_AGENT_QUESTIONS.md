# ⚡ Quick Reference: Agent Implementation Questions

**Last Updated:** January 28, 2025

---

## 🎯 The 3 Questions Every Agent Asks

### Q1: What exactly should I build?

**Answer:** Check your assignment type

| If you have... | Then build... |
|---------------|--------------|
| `AGENT_XX_NAME` in registry | ✅ Full orchestrator module (Type A) |
| `XX_MODULE.mdc` file | ✅ Full orchestrator module (Type A) |
| Task list with file paths | ✅ Feature enhancement (Type B) |
| Both above | ✅ Both - module integrates enhancements |

**Example (JR Agent G - Nutrition):**
- ✅ Build `AGENT_06_NUTRITION` (orchestrator module)
- ✅ Integrate MOAT toxicity features (your module uses them)
- ❌ Don't build MOAT in isolation

---

### Q2: Should I call API endpoints or use services directly?

**Answer:** ✅ **Always use services directly**

```python
# ✅ CORRECT
from ..my_service import MyService
service = MyService()
result = await service.process(data)

# ❌ WRONG
import httpx
response = await httpx.post("http://localhost:8000/api/...")
```

**Why?**
- ⚡ Faster (no HTTP overhead)
- 🔒 Type-safe (IDE support)
- 🎯 Consistent (matches all other agents)
- 🧪 Testable (no server needed)

---

### Q3: Should I modify the care plan agent?

**Answer:** ❌ **No - it auto-consumes your data**

**Pattern:**
```python
# You: Populate state
state.nutrition_plan = result

# Care Plan: Reads from state
state.nutrition_plan.supplements  # ← Automatic
```

**Your job:**
1. ✅ Return correct dict structure (from MDC)
2. ✅ Wire to orchestrator (`state.your_field = result`)
3. ❌ Don't touch care plan agent

---

## 📝 The Pattern (Copy This)

```python
# In orchestrator.py
async def _run_your_agent(self, state: PatientState) -> Dict:
    """Run the YOUR_NAME agent."""
    execution = state.start_agent('your_name')
    
    try:
        # 1. Import service
        from ..your_module import YourAgent
        
        # 2. Build request from state
        agent = YourAgent()
        
        # 3. Call service
        result = await agent.process(
            patient_id=state.patient_id,
            mutations=state.mutations,
            # ... extract what you need from state
        )
        
        # 4. Complete execution
        execution.complete(result)
        return result
        
    except Exception as e:
        execution.fail(str(e))
        raise
```

**That's it! Follow this pattern for every agent.**

---

## 🔍 Where to Look

| Question | Check Here |
|----------|-----------|
| What's my assignment? | `.cursor/MOAT/AGENT_OWNERSHIP_REGISTRY.md` |
| What should I build? | `.cursor/MOAT/orchestration/XX_MODULE.mdc` |
| How do I integrate? | `api/services/orchestrator/orchestrator.py` |
| What's the output format? | MDC file → `## 📤 OUTPUTS` section |
| How do other agents do it? | Search `_run_biomarker_agent`, `_run_resistance_agent`, `_run_trial_matching_agent` |

---

## ✅ Done When...

- [ ] Returns dict matching MDC schema
- [ ] Uses direct service imports (no HTTP)
- [ ] Wired to orchestrator (replaces TODO)
- [ ] Populates `state.your_field`
- [ ] Has unit tests
- [ ] Has integration test with orchestrator
- [ ] Registry status updated

---

## 📖 Full Guide

For detailed explanations and code examples, see:
**`.cursor/MOAT/AGENT_IMPLEMENTATION_GUIDE.md`**

---

**Quick Reference Status:** ✅ ACTIVE  
**For:** All JR Agents building MOAT modules





