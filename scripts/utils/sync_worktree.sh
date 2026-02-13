#!/bin/bash
# Sync script to keep main repo and worktree in sync
# Usage: ./scripts/sync_worktree.sh [direction]
#   direction: "to-worktree" (default) or "to-main" or "both"

set -e

MAIN_REPO="/Users/fahadkiani/Desktop/development/crispr-assistant-main"
WORKTREE="/Users/fahadkiani/.cursor/worktrees/crispr-assistant-main/nau"

# Files to sync (relative to publications/sporadic_cancer/)
FILES=(
    "tables.md"
    "supplement.md"
    "submission_aacr/TABLES.md"
    "submission_aacr/SUPPLEMENT.md"
    "submission_aacr/FIGURES_TABLES_LIST.md"
    "submission_aacr/MANUSCRIPT_DRAFT.md"
    "submission_aacr/AUTHOR_CONTRIBUTIONS.md"
    "submission_aacr/COVER_LETTER.md"
)

DIRECTION="${1:-to-worktree}"

sync_file() {
    local file="$1"
    local src="$2"
    local dst="$3"
    
    if [ -f "$src" ]; then
        cp "$src" "$dst"
        echo "Synced: $file"
    else
        echo "Missing: $file"
    fi
}

case "$DIRECTION" in
    "to-worktree")
        echo "=== Syncing Main Repo to Worktree ==="
        for file in "${FILES[@]}"; do
            src="$MAIN_REPO/publications/sporadic_cancer/$file"
            dst="$WORKTREE/publications/sporadic_cancer/$file"
            mkdir -p "$(dirname "$dst")"
            sync_file "$file" "$src" "$dst"
        done
        ;;
    "to-main")
        echo "=== Syncing Worktree to Main Repo ==="
        for file in "${FILES[@]}"; do
            src="$WORKTREE/publications/sporadic_cancer/$file"
            dst="$MAIN_REPO/publications/sporadic_cancer/$file"
            mkdir -p "$(dirname "$dst")"
            sync_file "$file" "$src" "$dst"
        done
        ;;
    "both")
        echo "=== Syncing Both Directions ==="
        echo "Step 1: Main to Worktree"
        "$0" to-worktree
        echo ""
        echo "Step 2: Worktree to Main"
        "$0" to-main
        ;;
    *)
        echo "Usage: $0 [to-worktree|to-main|both]"
        exit 1
        ;;
esac

echo ""
echo "=== Sync Complete ==="
