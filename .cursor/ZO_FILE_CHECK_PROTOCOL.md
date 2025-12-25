# ZO FILE CHECK PROTOCOL

## CRITICAL RULE
**NEVER claim a file doesn't exist without FIRST checking via terminal.**

## When Cursor File Tools Fail:
1. **IMMEDIATELY** run: `ls -la /path/to/file` or `test -f /path/to/file && echo "EXISTS" || echo "NOT FOUND"`
2. **ONLY THEN** make claims about file existence
3. **ALWAYS** use absolute paths when checking main repo vs worktree

## Worktree Sandbox Detection:
- If file tool returns: `ENOENT: ...worktrees/crispr-assistant-main/apd/...`
- This means: Cursor is redirecting to worktree
- **ACTION**: Use terminal commands to access main repo directly

## File Existence Claims:
- ❌ WRONG: "File does not exist" (without terminal check)
- ✅ RIGHT: "File tool failed, checking via terminal..." → then verify

## Sync Script Usage:
- `./sync_worktree_to_main.sh to-main` - Sync worktree → main
- `./sync_worktree_to_main.sh to-worktree` - Sync main → worktree
