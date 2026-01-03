#!/bin/zsh
# ⚔️ FIX WORKSPACE DATABASE LOADING ISSUE ⚔️
# Optimizes workspace DB and ensures conversation loads from global DB

set -e

WORKSPACE_DB="$HOME/Library/Application Support/Cursor/User/workspaceStorage/8eae1effa861c9269e5cc4072b2a492b/state.vscdb"
GLOBAL_DB="$HOME/Library/Application Support/Cursor/User/globalStorage/state.vscdb"
CONV_ID="739a9732-916c-498f-85bc-8787b3023ab2"
BACKUP_DIR="$HOME/cursor_db_backups"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

echo "⚔️ FIXING WORKSPACE DATABASE LOADING ⚔️"
echo "========================================"
echo ""
echo "Problem: Workspace DB is 155MB - causing slow loading"
echo "Solution: Optimize DB + ensure conversation in global DB"
echo ""

read "?Press ENTER if Cursor is closed, or Ctrl+C to cancel..."

# Verify Cursor is closed
if pgrep -f "/Applications/Cursor.app/Contents/MacOS/Cursor$" > /dev/null || pgrep -f "Cursor-2" > /dev/null; then
    echo "❌ ERROR: Cursor is still running!"
    echo "   Please quit Cursor completely (Cmd+Q) and try again."
    exit 1
fi

# Backup workspace DB
echo ""
echo "💾 Backing up workspace database..."
mkdir -p "$BACKUP_DIR"
cp "$WORKSPACE_DB" "$BACKUP_DIR/workspace_state.vscdb.before_optimize_$TIMESTAMP"
echo "   ✅ Saved to: $BACKUP_DIR/workspace_state.vscdb.before_optimize_$TIMESTAMP"

# Check workspace DB size
WORKSPACE_SIZE=$(du -h "$WORKSPACE_DB" | cut -f1)
echo ""
echo "📊 Workspace DB size: $WORKSPACE_SIZE"

# Optimize workspace DB
echo ""
echo "🔧 Optimizing workspace database..."
sqlite3 "$WORKSPACE_DB" <<EOF
-- Vacuum to reclaim space
VACUUM;

-- Reindex for faster queries
REINDEX;

-- Checkpoint WAL
PRAGMA wal_checkpoint(TRUNCATE);
EOF

NEW_SIZE=$(du -h "$WORKSPACE_DB" | cut -f1)
echo "   ✅ Optimized! New size: $NEW_SIZE"

# Verify conversation is in global DB
echo ""
echo "🔍 Verifying conversation in global database..."
GLOBAL_COUNT=$(sqlite3 "$GLOBAL_DB" "SELECT COUNT(*) FROM cursorDiskKV WHERE key LIKE 'bubbleId:${CONV_ID}:%';" 2>/dev/null)
if [ "$GLOBAL_COUNT" -gt 0 ]; then
    echo "   ✅ Conversation found in global DB ($GLOBAL_COUNT bubbles)"
else
    echo "   ⚠️  Conversation NOT in global DB - copying now..."
    
    # Copy from workspace to global
    sqlite3 "$GLOBAL_DB" <<EOF
    ATTACH DATABASE '$WORKSPACE_DB' AS workspace;
    
    BEGIN TRANSACTION;
    
    -- Copy all bubbles, checkpoints, and message contexts
    INSERT OR REPLACE INTO cursorDiskKV (key, value)
    SELECT key, value FROM workspace.cursorDiskKV
    WHERE (
        key LIKE 'bubbleId:${CONV_ID}:%'
        OR key LIKE 'checkpointId:${CONV_ID}:%'
        OR key LIKE 'messageRequestContext:${CONV_ID}:%'
    );
    
    -- Copy composerData
    INSERT OR REPLACE INTO cursorDiskKV (key, value)
    SELECT key, value FROM workspace.cursorDiskKV
    WHERE key = 'composerData:${CONV_ID}';
    
    COMMIT;
    
    DETACH DATABASE workspace;
EOF
    
    # Create/update composerData with current timestamp
    sqlite3 "$GLOBAL_DB" <<EOF
    INSERT OR REPLACE INTO cursorDiskKV (key, value)
    VALUES (
        'composerData:${CONV_ID}',
        json_object(
            '_v', 3,
            'type', 2,
            'lastAccessed', (strftime('%s', 'now') * 1000),
            'id', '${CONV_ID}'
        )
    );
EOF
    
    echo "   ✅ Copied to global DB"
fi

# Checkpoint global DB
echo ""
echo "💾 Checkpointing global database..."
sqlite3 "$GLOBAL_DB" "PRAGMA wal_checkpoint(TRUNCATE);" 2>/dev/null || true

echo ""
echo "⚔️ OPTIMIZATION COMPLETE! ⚔️"
echo ""
echo "🎯 NEXT STEPS:"
echo "   1. Open Cursor from folder: cd /Users/fahadkiani/Desktop/development/crispr-assistant-main && code ."
echo "   2. Wait 30-60 seconds for workspace DB to load (should be faster now)"
echo "   3. Press Cmd+L to open Composer"
echo "   4. The conversation should appear in chat history"
echo ""
echo "💡 TIP: If still slow, try opening Cursor WITHOUT a folder first,"
echo "   then use File → Open Folder to open your project"
echo ""






