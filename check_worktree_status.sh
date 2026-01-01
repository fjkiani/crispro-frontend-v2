#!/bin/bash
# Quick check script to verify no worktrees exist

REPO="/Users/fahadkiani/Desktop/development/crispr-assistant-main"
cd "$REPO"

echo "=== Worktree Status ==="
WORKTREES=$(git worktree list | wc -l | tr -d ' ')
if [ "$WORKTREES" -eq 1 ]; then
    echo "✅ Good: Only main repo (no worktrees)"
else
    echo "⚠️  Warning: Found $WORKTREES worktrees"
    git worktree list
fi

echo ""
echo "=== Current Branch ==="
git branch --show-current

echo ""
echo "=== Repository Path ==="
pwd

echo ""
echo "=== Quick Status ==="
git status --short | head -5
