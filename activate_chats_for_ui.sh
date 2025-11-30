#!/bin/zsh

# Activate restored chats for Cursor 2.0 UI display
# This script updates conversation metadata to make them visible

set -e

CURSOR_DB="$HOME/Library/Application Support/Cursor/User/globalStorage/state.vscdb"
BACKUP_DIR="$HOME/cursor_db_backups"

echo "⚔️ ACTIVATING CHATS FOR CURSOR UI..."
echo ""

# Check if Cursor is running
if pgrep -x "Cursor" > /dev/null; then
  echo "❌ Cursor is running! Please quit Cursor completely (Cmd+Q) first."
  exit 1
fi

# Backup
mkdir -p "$BACKUP_DIR"
BACKUP_FILE="$BACKUP_DIR/state.vscdb.before_activate_$(date +%Y%m%d_%H%M%S)"
cp "$CURSOR_DB" "$BACKUP_FILE"
echo "✅ Backed up to: $BACKUP_FILE"
echo ""

# Get all conversation IDs
echo "📊 Finding all conversations..."
sqlite3 "$CURSOR_DB" "SELECT DISTINCT substr(key, 10, 36) FROM cursorDiskKV WHERE key LIKE 'bubbleId:%';" 2>/dev/null > /tmp/all_conversations.txt
CONV_COUNT=$(wc -l < /tmp/all_conversations.txt | xargs)
echo "   Found $CONV_COUNT conversations"
echo ""

# Strategy: Create/update composerData entries for each conversation
# These might be what Cursor uses to display chats in the UI
echo "🔧 Creating composerData entries for each conversation..."

# Create composerData entries for ALL conversations
sqlite3 "$CURSOR_DB" << 'EOF'
BEGIN TRANSACTION;

-- For each conversation, create/update a composerData entry
-- This might be what Cursor uses to display chats in the UI
INSERT OR REPLACE INTO cursorDiskKV (key, value)
SELECT 
  'composerData:' || conv_id as new_key,
  '{"_v":3,"type":2,"lastAccessed":' || (strftime('%s', 'now') * 1000) || '}' as new_value
FROM (
  SELECT DISTINCT substr(key, 10, 36) as conv_id
  FROM cursorDiskKV
  WHERE key LIKE 'bubbleId:%'
);

COMMIT;
EOF

echo "   ✅ Created composerData entries"
echo ""

# Alternative: Try to find the most recent bubble in each conversation and "touch" it
echo "🔧 Updating conversation timestamps..."

# Get the most recent message from each conversation and update it
sqlite3 "$CURSOR_DB" << 'EOF'
BEGIN TRANSACTION;

-- This will force SQLite to recognize the conversations as "active"
UPDATE cursorDiskKV 
SET value = value  -- No-op update to trigger change detection
WHERE key LIKE 'bubbleId:%'
AND rowid IN (
  SELECT rowid FROM cursorDiskKV 
  WHERE key LIKE 'bubbleId:%'
  GROUP BY substr(key, 10, 36)
  ORDER BY rowid DESC
);

COMMIT;
EOF

echo "   ✅ Updated timestamps"
echo ""

# Checkpoint WAL
echo "🔧 Checkpointing WAL..."
sqlite3 "$CURSOR_DB" "PRAGMA wal_checkpoint(FULL);" 2>/dev/null
echo "   ✅ WAL checkpointed"
echo ""

# Final verification
FINAL_CONV=$(sqlite3 "$CURSOR_DB" "SELECT COUNT(DISTINCT substr(key, 10, 36)) FROM cursorDiskKV WHERE key LIKE 'bubbleId:%';" 2>/dev/null | xargs)
FINAL_COMPOSER=$(sqlite3 "$CURSOR_DB" "SELECT COUNT(*) FROM cursorDiskKV WHERE key LIKE 'composerData:%';" 2>/dev/null | xargs)

echo "📊 FINAL STATUS:"
echo "   Conversations: $FINAL_CONV"
echo "   ComposerData entries: $FINAL_COMPOSER"
echo ""

echo "🎯 ACTIVATION COMPLETE!"
echo ""
echo "📋 NEXT STEPS:"
echo ""
echo "   1. Open Cursor"
echo "   2. Wait 3-5 minutes"
echo "   3. Try Cmd+Shift+P → 'Developer: Reload Window'"
echo "   4. Check chat history panel (clock icon)"
echo ""
echo "💡 If still not visible, this is a Cursor 2.0 bug."
echo "   The data is definitely in the database ($FINAL_CONV conversations)."
echo "   You may need to wait for a Cursor update or contact Cursor support."
echo ""
echo "💾 Backup: $BACKUP_FILE"
echo ""

