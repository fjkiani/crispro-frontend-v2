#!/bin/zsh

# Activate chats in ALL workspaces for crispr-assistant-main
# This ensures chats show regardless of which workspace Cursor uses

set -e

PROJECT_PATH="file:///Users/fahadkiani/Desktop/development/crispr-assistant-main"
GLOBAL_DB="$HOME/Library/Application Support/Cursor/User/globalStorage/state.vscdb"
BACKUP_DIR="$HOME/cursor_db_backups"

echo "🔄 ACTIVATING CHATS IN ALL WORKSPACES"
echo ""

# Check if Cursor is running
if pgrep -x "Cursor" > /dev/null; then
  echo "❌ Cursor is running! Please quit completely (Cmd+Q) first."
  exit 1
fi

mkdir -p "$BACKUP_DIR"

# Find all workspaces for crispr-assistant-main
echo "1️⃣ Finding all workspaces for crispr-assistant-main..."
find ~/Library/Application\ Support/Cursor/User/workspaceStorage -name "workspace.json" -exec grep -l "crispr-assistant-main" {} \; 2>/dev/null | while read ws_file; do
  ws_dir=$(dirname "$ws_file")
  ws_id=$(basename "$ws_dir")
  ws_db="$ws_dir/state.vscdb"
  
  if [ ! -f "$ws_db" ]; then
    echo "   ⚠️ Workspace $ws_id: No database file, skipping"
    continue
  fi
  
  echo ""
  echo "=== Processing Workspace: $ws_id ==="
  
  # Backup
  BACKUP_FILE="$BACKUP_DIR/workspace_${ws_id}_before_activation_$(date +%Y%m%d_%H%M%S)"
  cp "$ws_db" "$BACKUP_FILE"
  echo "   ✅ Backed up to: $BACKUP_FILE"
  
  # Get conversation count
  CONV_COUNT=$(sqlite3 "$ws_db" "SELECT COUNT(DISTINCT substr(key, 10, 36)) FROM cursorDiskKV WHERE key LIKE 'bubbleId:%';" 2>/dev/null || echo "0")
  COMPOSER_COUNT=$(sqlite3 "$ws_db" "SELECT COUNT(*) FROM cursorDiskKV WHERE key LIKE 'composerData:%';" 2>/dev/null || echo "0")
  
  echo "   Current: $CONV_COUNT conversations, $COMPOSER_COUNT composerData entries"
  
  # Sync from global if needed
  if [ "$CONV_COUNT" -eq 0 ]; then
    echo "   📥 No conversations found, syncing from global database..."
    
    # Get conversations from global
    sqlite3 "$GLOBAL_DB" "SELECT DISTINCT substr(key, 10, 36) FROM cursorDiskKV WHERE key LIKE 'bubbleId:%';" 2>/dev/null | while read conv_id; do
      # Copy conversation entries
      sqlite3 "$GLOBAL_DB" "SELECT 'INSERT OR REPLACE INTO cursorDiskKV (key, value) VALUES (''' || replace(key, '''', '''''') || ''', ''' || replace(value, '''', '''''') || ''');' FROM cursorDiskKV WHERE key LIKE 'bubbleId:$conv_id:%';" 2>/dev/null | sqlite3 "$ws_db"
      
      # Copy composerData
      sqlite3 "$GLOBAL_DB" "SELECT 'INSERT OR REPLACE INTO cursorDiskKV (key, value) VALUES (''' || replace(key, '''', '''''') || ''', ''' || replace(value, '''', '''''') || ''');' FROM cursorDiskKV WHERE key = 'composerData:$conv_id';" 2>/dev/null | sqlite3 "$ws_db"
    done
    
    echo "   ✅ Synced from global database"
  fi
  
  # Ensure composerData for all conversations
  echo "   🔧 Ensuring composerData entries..."
  sqlite3 "$ws_db" << 'SQL_EOF'
BEGIN TRANSACTION;

INSERT OR REPLACE INTO cursorDiskKV (key, value)
SELECT
  'composerData:' || substr(key, 10, 36) as new_key,
  '{"_v":3,"type":2,"lastAccessed":' || (strftime('%s', 'now') * 1000) || '}' as new_value
FROM cursorDiskKV
WHERE key LIKE 'bubbleId:%'
GROUP BY substr(key, 10, 36);

COMMIT;
SQL_EOF
  
  # Checkpoint WAL
  sqlite3 "$ws_db" "PRAGMA wal_checkpoint(TRUNCATE);"
  
  # Touch database
  touch "$ws_db"
  
  # Final stats
  FINAL_CONV=$(sqlite3 "$ws_db" "SELECT COUNT(DISTINCT substr(key, 10, 36)) FROM cursorDiskKV WHERE key LIKE 'bubbleId:%';" 2>/dev/null || echo "0")
  FINAL_COMPOSER=$(sqlite3 "$ws_db" "SELECT COUNT(*) FROM cursorDiskKV WHERE key LIKE 'composerData:%';" 2>/dev/null || echo "0")
  
  echo "   ✅ Final: $FINAL_CONV conversations, $FINAL_COMPOSER composerData entries"
done

echo ""
echo "🎯 ACTIVATION COMPLETE!"
echo ""
echo "📋 NEXT STEPS:"
echo "   1. Open Cursor"
echo "   2. Open crispr-assistant-main project directory"
echo "   3. Wait 2-3 minutes for indexing"
echo "   4. Check chat history - should now show all conversations!"
echo ""
