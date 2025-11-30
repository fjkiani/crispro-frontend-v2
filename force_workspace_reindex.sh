#!/bin/zsh

# Force Cursor to re-index workspace chat history
# This addresses the UI bug where chats exist but don't display

set -e

MAIN_WS_ID="8eae1effa861c9269e5cc4072b2a492b"
WS_DB="$HOME/Library/Application Support/Cursor/User/workspaceStorage/$MAIN_WS_ID/state.vscdb"
CACHED_DATA_DIR="$HOME/Library/Application Support/Cursor/CachedData"
BACKUP_DIR="$HOME/cursor_db_backups"

echo "🔄 FORCING WORKSPACE RE-INDEX"
echo ""

if pgrep -x "Cursor" > /dev/null; then
  echo "❌ Cursor is running! Please quit completely (Cmd+Q) first."
  exit 1
fi

# Backup
mkdir -p "$BACKUP_DIR"
BACKUP_FILE="$BACKUP_DIR/workspace_${MAIN_WS_ID}_before_reindex_$(date +%Y%m%d_%H%M%S)"
cp "$WS_DB" "$BACKUP_FILE"
echo "✅ Backed up to: $BACKUP_FILE"
echo ""

# Step 1: Update ALL composerData entries with fresh timestamps
echo "1️⃣ Updating composerData timestamps..."
sqlite3 "$WS_DB" << 'SQL_EOF'
BEGIN TRANSACTION;

-- Update all composerData entries with current timestamp
UPDATE cursorDiskKV
SET value = '{"_v":3,"type":2,"lastAccessed":' || (strftime('%s', 'now') * 1000) || '}'
WHERE key LIKE 'composerData:%';

-- Ensure composerData exists for all conversations
INSERT OR REPLACE INTO cursorDiskKV (key, value)
SELECT
  'composerData:' || substr(key, 10, 36) as new_key,
  '{"_v":3,"type":2,"lastAccessed":' || (strftime('%s', 'now') * 1000) || '}' as new_value
FROM cursorDiskKV
WHERE key LIKE 'bubbleId:%'
AND NOT EXISTS (
  SELECT 1 FROM cursorDiskKV AS T2 
  WHERE T2.key = 'composerData:' || substr(cursorDiskKV.key, 10, 36)
)
GROUP BY substr(key, 10, 36);

COMMIT;
SQL_EOF
echo "   ✅ ComposerData timestamps updated"
echo ""

# Step 2: Update conversation bubbleId entries (touch them)
echo "2️⃣ Touching conversation entries..."
sqlite3 "$WS_DB" << 'SQL_EOF'
BEGIN TRANSACTION;

-- Update the main bubbleId entry for each conversation (the one without :message: suffix)
UPDATE cursorDiskKV
SET value = value
WHERE key LIKE 'bubbleId:%'
AND key NOT LIKE 'bubbleId:%:%'
AND rowid IN (
  SELECT MIN(rowid)
  FROM cursorDiskKV
  WHERE key LIKE 'bubbleId:%'
  AND key NOT LIKE 'bubbleId:%:%'
  GROUP BY substr(key, 10, 36)
);

COMMIT;
SQL_EOF
echo "   ✅ Conversation entries touched"
echo ""

# Step 3: Checkpoint WAL
echo "3️⃣ Checkpointing WAL..."
sqlite3 "$WS_DB" "PRAGMA wal_checkpoint(TRUNCATE);"
echo "   ✅ WAL checkpointed"
echo ""

# Step 4: Touch database file
echo "4️⃣ Touching database file..."
touch "$WS_DB"
echo "   ✅ Database touched"
echo ""

# Step 5: Clear CachedData
echo "5️⃣ Clearing CachedData directory..."
if [ -d "$CACHED_DATA_DIR" ]; then
  rm -rf "$CACHED_DATA_DIR"
  echo "   ✅ CachedData cleared"
else
  echo "   ℹ️  CachedData directory not found"
fi
echo ""

# Step 6: Verify
FINAL_CONV=$(sqlite3 "$WS_DB" "SELECT COUNT(DISTINCT substr(key, 10, 36)) FROM cursorDiskKV WHERE key LIKE 'bubbleId:%';" 2>/dev/null || echo "0")
FINAL_COMPOSER=$(sqlite3 "$WS_DB" "SELECT COUNT(*) FROM cursorDiskKV WHERE key LIKE 'composerData:%';" 2>/dev/null || echo "0")

echo "🎯 RE-INDEX COMPLETE!"
echo ""
echo "📊 Final Stats:"
echo "   Conversations: $FINAL_CONV"
echo "   ComposerData entries: $FINAL_COMPOSER"
echo ""
echo "💾 Backup saved to: $BACKUP_FILE"
echo ""
echo "📋 NEXT STEPS:"
echo "   1. Open Cursor"
echo "   2. Open crispr-assistant-main project directory"
echo "   3. Wait 3-5 minutes for full re-indexing"
echo "   4. Try these actions to force UI refresh:"
echo "      - Click the chat history icon (clock) in sidebar"
echo "      - Press Cmd+Shift+P → 'Reload Window'"
echo "      - Close and reopen Cursor completely"
echo "   5. If still not visible, check if chats appear after searching (Cmd+F in chat panel)"
echo ""
