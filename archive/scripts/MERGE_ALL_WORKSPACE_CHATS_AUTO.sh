#!/bin/zsh
# ⚔️ MERGE ALL WORKSPACE CHAT CONVERSATIONS INTO GLOBAL DATABASE ⚔️
# AUTO-RUN VERSION (no prompts)

set -e

CURSOR_DB="$HOME/Library/Application Support/Cursor/User/globalStorage/state.vscdb"
WORKSPACE_STORAGE="$HOME/Library/Application Support/Cursor/User/workspaceStorage"
BACKUP_DIR="$HOME/cursor_db_backups"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

echo "⚔️ MERGE ALL WORKSPACE CHAT CONVERSATIONS ⚔️"
echo "=============================================="
echo ""

# 1. Kill all Cursor processes
echo "🛑 Closing Cursor completely..."
pkill -9 -f "Cursor" || true
sleep 3

# Verify Cursor is closed
if pgrep -f "Cursor" > /dev/null; then
    echo "❌ ERROR: Cursor is still running! Please close it manually and try again."
    exit 1
fi

# 2. Backup current database
echo ""
echo "💾 Backing up current global database..."
mkdir -p "$BACKUP_DIR"
BACKUP_FILE="$BACKUP_DIR/state.vscdb.before_merge_$TIMESTAMP"
cp "$CURSOR_DB" "$BACKUP_FILE"
echo "   ✅ Saved to: $BACKUP_FILE"

# 3. Find all workspace databases with conversations
echo ""
echo "🔍 Finding all workspace databases with conversations..."

WORKSPACE_COUNT=0
TOTAL_CONVERSATIONS=0
WORKSPACE_DBS=()

for workspace_db in "$WORKSPACE_STORAGE"/*/state.vscdb; do
    if [ ! -f "$workspace_db" ]; then
        continue
    fi
    
    workspace_id=$(basename $(dirname "$workspace_db"))
    conv_count=$(sqlite3 "$workspace_db" "SELECT COUNT(DISTINCT substr(key, 10, 36)) FROM cursorDiskKV WHERE key LIKE 'bubbleId:%';" 2>/dev/null || echo "0")
    
    if [ "$conv_count" -gt 0 ]; then
        echo "   📁 Workspace $workspace_id: $conv_count conversations"
        WORKSPACE_DBS+=("$workspace_db")
        WORKSPACE_COUNT=$((WORKSPACE_COUNT + 1))
        TOTAL_CONVERSATIONS=$((TOTAL_CONVERSATIONS + conv_count))
    fi
done

echo ""
echo "📊 Found $WORKSPACE_COUNT workspace databases with $TOTAL_CONVERSATIONS total conversations"

if [ $WORKSPACE_COUNT -eq 0 ]; then
    echo "⚠️  No workspace databases with conversations found!"
    exit 1
fi

# 4. Merge into global database using ATTACH
echo ""
echo "⚔️ MERGING ALL CONVERSATIONS INTO GLOBAL DATABASE..."
echo "   This may take a minute..."

# Create a temporary SQL script for merging
TEMP_SQL="/tmp/merge_workspace_chats_$TIMESTAMP.sql"
> "$TEMP_SQL"

for workspace_db in "${WORKSPACE_DBS[@]}"; do
    workspace_id=$(basename $(dirname "$workspace_db"))
    echo "-- Merging workspace: $workspace_id" >> "$TEMP_SQL"
    echo "ATTACH DATABASE '$workspace_db' AS workspace;" >> "$TEMP_SQL"
    
    # Copy all chat-related entries from workspace to global
    cat >> "$TEMP_SQL" <<'EOFMERGE'
-- Copy bubbles
INSERT OR REPLACE INTO cursorDiskKV (key, value)
SELECT key, value FROM workspace.cursorDiskKV
WHERE key LIKE 'bubbleId:%' OR key LIKE 'checkpointId:%' OR key LIKE 'messageRequestContext:%' OR key LIKE 'composerData:%';

-- Copy composer UI state from ItemTable
INSERT OR REPLACE INTO ItemTable (key, value)
SELECT key, value FROM workspace.ItemTable
WHERE key LIKE '%composerChatViewPane%' OR key LIKE '%composerData%';

DETACH DATABASE workspace;
EOFMERGE
done

# Execute the merge
sqlite3 "$CURSOR_DB" < "$TEMP_SQL" 2>&1 | head -20

# 5. Create composerData entries for all conversations
echo ""
echo "📝 Creating composerData entries for UI visibility..."

sqlite3 "$CURSOR_DB" <<EOF
-- Get all unique conversation IDs
CREATE TEMP TABLE IF NOT EXISTS all_conv_ids AS
SELECT DISTINCT substr(key, 10, 36) as conv_id
FROM cursorDiskKV
WHERE key LIKE 'bubbleId:%';

-- Create composerData entries for each conversation
INSERT OR REPLACE INTO cursorDiskKV (key, value)
SELECT 
    'composerData:' || conv_id,
    json_object(
        'id', conv_id,
        'title', 'Conversation ' || substr(conv_id, 1, 8),
        'created', (SELECT strftime('%s', 'now')),
        'updated', (SELECT strftime('%s', 'now'))
    )
FROM all_conv_ids
WHERE NOT EXISTS (
    SELECT 1 FROM cursorDiskKV 
    WHERE key = 'composerData:' || conv_id
);

DROP TABLE all_conv_ids;
EOF

# 6. Checkpoint WAL
echo ""
echo "💾 Checkpointing WAL..."
sqlite3 "$CURSOR_DB" "PRAGMA wal_checkpoint(TRUNCATE);" 2>/dev/null || true

# 7. Verify
echo ""
echo "✅ VERIFICATION:"
GLOBAL_CONV_COUNT=$(sqlite3 "$CURSOR_DB" "SELECT COUNT(DISTINCT substr(key, 10, 36)) FROM cursorDiskKV WHERE key LIKE 'bubbleId:%';" 2>/dev/null)
COMPOSER_COUNT=$(sqlite3 "$CURSOR_DB" "SELECT COUNT(*) FROM cursorDiskKV WHERE key LIKE 'composerData:%';" 2>/dev/null)
echo "   📊 Conversations in global DB: $GLOBAL_CONV_COUNT"
echo "   📝 ComposerData entries: $COMPOSER_COUNT"

# 8. Cleanup
rm -f "$TEMP_SQL"

# 9. Reopen Cursor
echo ""
echo "🚀 Opening Cursor..."
open -a "Cursor-2" || open -a "Cursor"
sleep 3

echo ""
echo "⚔️ MERGE COMPLETE! ⚔️"
echo ""
echo "🎯 Next steps:"
echo "1. Wait for Cursor to fully load (30-60 seconds)"
echo "2. Click on 'Previous Chats' or check the chat history panel"
echo "3. Your conversations should now be visible!"
echo ""
echo "📊 Summary:"
echo "   - Merged from $WORKSPACE_COUNT workspace databases"
echo "   - Total conversations: $GLOBAL_CONV_COUNT"
echo "   - Backup saved to: $BACKUP_FILE"
echo ""
echo "⚔️ CONQUEST ACHIEVED! ⚔️"

