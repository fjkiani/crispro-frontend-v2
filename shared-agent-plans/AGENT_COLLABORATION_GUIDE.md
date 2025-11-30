# Agent Collaboration Guide

**Shared Plans Location**: `~/Desktop/shared-agent-plans/`

---

## Why This Location?

1. **Not in Git**: Files won't be accidentally committed
2. **Accessible to All Agents**: Any agent can read/write here
3. **Persistent**: Files persist across sessions
4. **Shared**: Multiple agents can collaborate

---

## File Structure

```
~/Desktop/shared-agent-plans/
├── README.md                                    # Directory overview
├── AGENT_COLLABORATION_GUIDE.md                 # This file
├── FINAL_COMPREHENSIVE_PLAN_ALIGNED.md         # Main implementation plan (Manager)
├── REALITY_CHECK_COMPLETE_ASSESSMENT.md         # Current state assessment
├── mechanism-trial-matching.md                 # Mechanism-based trial matching plan
└── zo2.mdc                                      # Zo2's SPE/WIWFM outcome prediction plan
```

---

## Agent Responsibilities

### Manager (Me - Mechanism-Based)
**Files I Own**:
- `FINAL_COMPREHENSIVE_PLAN_ALIGNED.md` - My implementation plan
- `REALITY_CHECK_COMPLETE_ASSESSMENT.md` - My reality check

**What I Do**:
- Mechanism-based trial matching
- Mechanism-based resistance prediction
- Validation framework

**What I Don't Touch**:
- `zo2.mdc` - Zo2's domain (outcome prediction)

### Zo2 (SPE/WIWFM - Outcome-Based)
**Files I Own**:
- `zo2.mdc` - My implementation plan

**What I Do**:
- PFS/OS correlation improvement
- TMB/HRD/MSI biomarker gates
- Disease-specific SPE weights
- Outcome calibration

**What I Don't Touch**:
- Mechanism fit ranker code
- Resistance prophet code
- Trial matching code

---

## Collaboration Protocol

### When I (Manager) Find Something for Zo2:
1. Document in `~/Desktop/shared-agent-plans/MANAGER_TO_ZO2_HANDOFFS.md`
2. Zo2 decides if/how to integrate

### When Zo2 Finds Something for Me:
1. Document in `~/Desktop/shared-agent-plans/ZO2_TO_MANAGER_HANDOFFS.md`
2. I decide if/how to integrate

### Shared Files:
- Both agents can read all files
- Each agent owns their plan file
- Handoff files are for communication

---

## File Naming Convention

- `*_PLAN.md` or `*.mdc` - Implementation plans
- `*_ASSESSMENT.md` - Reality checks / assessments
- `*_HANDOFFS.md` - Inter-agent communication
- `README.md` - Documentation

---

## Accessing Files

**From Python**:
```python
import os
plan_path = os.path.expanduser("~/Desktop/shared-agent-plans/FINAL_COMPREHENSIVE_PLAN_ALIGNED.md")
with open(plan_path) as f:
    plan = f.read()
```

**From Terminal**:
```bash
cat ~/Desktop/shared-agent-plans/FINAL_COMPREHENSIVE_PLAN_ALIGNED.md
```

**From Cursor**:
- Open: `~/Desktop/shared-agent-plans/FINAL_COMPREHENSIVE_PLAN_ALIGNED.md`
- Or use absolute path: `/Users/fahadkiani/Desktop/shared-agent-plans/`

---

## Benefits

✅ **No Git Noise**: Plans don't clutter git history  
✅ **Agent Access**: All agents can read/write  
✅ **Persistent**: Files survive git operations  
✅ **Shared**: Easy collaboration between agents  

---

**Location**: `~/Desktop/shared-agent-plans/`





