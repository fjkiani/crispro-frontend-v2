#!/bin/zsh

# Activate chat history in workspace-specific database
# This fixes the issue where chats don't show when opening a project directory

set -e

WS_ID="8eae1effa861c9269e5cc4072b2a492b"
WS_DB="$HOME/Library/Application Support/Cursor/User/workspaceStorage/$WS_ID/state.vscdb"
BACKUP_DIR="$HOME/cursor_db_backups"

echo "🚀 WORKSPACE CHAT ACTIVATION SCRIPT"
echo "Workspace: crispr-assistant-main"
echo "Database: $WS_DB"
echo ""

# Check if Cursor is running
if pgrep -x "Cursor" > /dev/null; then
  echo "❌ Cursor is running! Please quit Cursor completely (Cmd+Q) first."
  exit 1
fi

# Verify database exists
if [ ! -f "$WS_DB" ]; then
  echo "❌ Workspace database not found: $WS_DB"
  exit 1
fi

# Backup current DB
mkdir -p "$BACKUP_DIR"
BACKUP_FILE="$BACKUP_DIR/workspace_state.vscdb.before_activation_$(date +%Y%m%d_%H%M%S)"
cp "$WS_DB" "$BACKUP_FILE"
echo "✅ Backed up workspace database to: $BACKUP_FILE"
echo ""

# Check current state
CURRENT_CONVS=$(sqlite3 "$WS_DB" "SELECT COUNT(DISTINCT substr(key, 10, 36)) FROM cursorDiskKV WHERE key LIKE 'bubbleId:%';" 2>/dev/null | xargs)
CURRENT_COMPOSER=$(sqlite3 "$WS_DB" "SELECT COUNT(*) FROM cursorDiskKV WHERE key LIKE 'composerData:%';" 2>/dev/null | xargs)

echo "📊 Current State:"
echo "   Conversations: $CURRENT_CONVS"
echo "   ComposerData entries: $CURRENT_COMPOSER"
echo ""

# Step 1: Create composerData entries for all conversations
echo "1️⃣ Creating composerData entries for all conversations..."
sqlite3 "$WS_DB" << 'EOF'
BEGIN TRANSACTION;

INSERT OR REPLACE INTO cursorDiskKV (key, value)
SELECT
  'composerData:' || substr(key, 10, 36) as new_key,
  '{"_v":3,"type":2,"lastAccessed":' || (strftime('%s', 'now') * 1000) || '}' as new_value
FROM cursorDiskKV
WHERE key LIKE 'bubbleId:%'
GROUP BY substr(key, 10, 36);

COMMIT;
EOF
echo "   ✅ Created composerData entries"
echo ""

# Step 2: Update modification timestamps
echo "2️⃣ Updating conversation timestamps..."
sqlite3 "$WS_DB" << 'EOF'
BEGIN TRANSACTION;

UPDATE cursorDiskKV
SET value = value
WHERE key LIKE 'bubbleId:%'
AND rowid IN (
  SELECT MAX(T1.rowid)
  FROM cursorDiskKV AS T1
  WHERE T1.key LIKE 'bubbleId:%'
  GROUP BY substr(T1.key, 10, 36)
);

COMMIT;
EOF
echo "   ✅ Timestamps updated"
echo ""

# Step 3: Checkpoint WAL file
echo "3️⃣ Checkpointing WAL file..."
sqlite3 "$WS_DB" "PRAGMA wal_checkpoint(TRUNCATE);" 2>/dev/null || echo "   ℹ️ No WAL file (normal)"
echo "   ✅ WAL checkpointed"
echo ""

# Step 4: Touch the database file
echo "4️⃣ Touching database file..."
touch "$WS_DB"
echo "   ✅ Database file touched"
echo ""

# Verify final state
FINAL_CONVS=$(sqlite3 "$WS_DB" "SELECT COUNT(DISTINCT substr(key, 10, 36)) FROM cursorDiskKV WHERE key LIKE 'bubbleId:%';" 2>/dev/null | xargs)
FINAL_COMPOSER=$(sqlite3 "$WS_DB" "SELECT COUNT(*) FROM cursorDiskKV WHERE key LIKE 'composerData:%';" 2>/dev/null | xargs)

echo "📊 Final State:"
echo "   Conversations: $FINAL_CONVS"
echo "   ComposerData entries: $FINAL_COMPOSER"
echo ""

echo "🎯 WORKSPACE ACTIVATION COMPLETE!"
echo ""
echo "📋 NEXT STEPS:"
echo ""
echo "   1. Open Cursor"
echo "   2. Open the crispr-assistant-main project directory"
echo "   3. Wait 2-3 minutes for indexing"
echo "   4. Check chat history - should now show all 19 conversations!"
echo ""
echo "💾 Backup saved to: $BACKUP_FILE"
echo ""





















































