# Shared Plans Location

**IMPORTANT**: Plan documents are now saved in a shared location outside git repos.

## Location

**Shared Directory**: `~/Desktop/shared-agent-plans/`

**Full Path**: `/Users/fahadkiani/Desktop/shared-agent-plans/`

## Why This Location?

1. **Not in Git**: Files won't be tracked/committed accidentally
2. **Accessible to All Agents**: Any agent can read/write here
3. **Persistent**: Files persist across sessions
4. **Shared**: Multiple agents can collaborate on same files

## Current Files

- `FINAL_COMPREHENSIVE_PLAN_ALIGNED.md` - Main implementation plan
- `REALITY_CHECK_COMPLETE_ASSESSMENT.md` - Current state assessment  
- `mechanism-trial-matching.md` - Mechanism-based trial matching plan
- `zo2.mdc` - Zo2's SPE/WIWFM outcome prediction plan
- `README.md` - Directory documentation

## Usage

**To Read**:
```python
with open('~/Desktop/shared-agent-plans/FINAL_COMPREHENSIVE_PLAN_ALIGNED.md') as f:
    plan = f.read()
```

**To Write**:
```python
with open('~/Desktop/shared-agent-plans/my_plan.md', 'w') as f:
    f.write(plan_content)
```

## Migration

Old files in `.cursor/plans/` can remain for reference, but new plans should go to:
`~/Desktop/shared-agent-plans/`





