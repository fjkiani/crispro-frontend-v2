#!/bin/bash
# Sync script: Worktree <-> Main Repo (FULL REPO SYNC)
# Usage: ./sync_worktree_to_main.sh [direction]
# direction: 'to-main' (default) or 'to-worktree'

WORKTREE_ROOT="/Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/apd"
MAIN_ROOT="/Users/fahadkiani/Desktop/development/crispr-assistant-main"
DIRECTION=${1:-to-main}

echo "🔄 Full Repo Sync: $DIRECTION"
echo "Worktree: $WORKTREE_ROOT"
echo "Main: $MAIN_ROOT"
echo "---"

# Exclusions: build artifacts, dependencies, git internals
EXCLUDE_ARGS="--exclude='node_modules' --exclude='venv' --exclude='__pycache__' --exclude='.git' --exclude='dist' --exclude='build' --exclude='.next' --exclude='.cache' --exclude='*.pyc' --exclude='.DS_Store'"

if [ "$DIRECTION" = "to-main" ]; then
    echo "Syncing FROM worktree TO main repo (FULL SYNC)..."
    rsync -av $EXCLUDE_ARGS \
        "$WORKTREE_ROOT/" \
        "$MAIN_ROOT/"
    echo "✅ Full sync complete: worktree -> main"
elif [ "$DIRECTION" = "to-worktree" ]; then
    echo "Syncing FROM main repo TO worktree (FULL SYNC)..."
    rsync -av $EXCLUDE_ARGS \
        "$MAIN_ROOT/" \
        "$WORKTREE_ROOT/"
    echo "✅ Full sync complete: main -> worktree"
else
    echo "❌ Invalid direction. Use 'to-main' or 'to-worktree'"
    exit 1
fi

echo ""
echo "📊 Sync Summary:"
echo "  - All source files synced"
echo "  - Excluded: node_modules, venv, __pycache__, .git, build artifacts"
