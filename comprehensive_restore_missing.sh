#!/bin/zsh
# ⚔️ COMPREHENSIVE CHAT RESTORATION - ALL ENTRIES ⚔️
# Restores missing conversations with ALL related entries (bubbleId, checkpointId, composerData, etc.)

set -e

MAIN_WS_ID="8eae1effa861c9269e5cc4072b2a492b"
WS_DB="$HOME/Library/Application Support/Cursor/User/workspaceStorage/$MAIN_WS_ID/state.vscdb"
GLOBAL_DB="$HOME/Library/Application Support/Cursor/User/globalStorage/state.vscdb"
CURSOR_SUPPORT="$HOME/Library/Application Support/Cursor"
BACKUP_DIR="$HOME/cursor_db_backups"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

echo "⚔️ COMPREHENSIVE CHAT RESTORATION ⚔️"
echo "===================================="
echo ""

# Check Cursor is not running
if pgrep -x "Cursor" > /dev/null; then
  echo "❌ Cursor is running! Please quit completely (Cmd+Q) first."
  echo ""
  echo "   Force quit? (y/n): "
  read -q "REPLY"
  echo ""
  if [[ "$REPLY" == "y" || "$REPLY" == "Y" ]]; then
    echo "   ⚔️ Force quitting Cursor..."
    killall -9 Cursor 2>/dev/null || true
    sleep 2
  else
    exit 1
  fi
fi

# Verify databases exist
if [ ! -f "$WS_DB" ]; then
  echo "❌ Workspace database not found: $WS_DB"
  exit 1
fi

if [ ! -f "$GLOBAL_DB" ]; then
  echo "❌ Global database not found: $GLOBAL_DB"
  exit 1
fi

# Backup global database
mkdir -p "$BACKUP_DIR"
GLOBAL_BACKUP="$BACKUP_DIR/global_before_comprehensive_restore_${TIMESTAMP}"
cp "$GLOBAL_DB" "$GLOBAL_BACKUP"
echo "✅ Backed up global database to: $(basename "$GLOBAL_BACKUP")"
echo ""

# Find missing conversations from ALL workspaceStorage directories
echo "1️⃣ Finding conversations in ALL workspaceStorage databases but not in global..."

# Check all workspaceStorage directories
WS_STORAGE_DIR="$CURSOR_SUPPORT/User/workspaceStorage"
ALL_WS_CONVS="/tmp/all_ws_convs.txt"
> "$ALL_WS_CONVS"

if [ -d "$WS_STORAGE_DIR" ]; then
  for ws_dir in "$WS_STORAGE_DIR"/*/; do
    if [ -d "$ws_dir" ]; then
      WS_DB_FILE="$ws_dir/state.vscdb"
      if [ -f "$WS_DB_FILE" ]; then
        WS_ID=$(basename "$ws_dir")
        echo "   Checking workspace: $WS_ID"
        sqlite3 "$WS_DB_FILE" "SELECT DISTINCT substr(key, 10, 36) FROM cursorDiskKV WHERE key LIKE 'bubbleId:%';" 2>/dev/null >> "$ALL_WS_CONVS"
      fi
    fi
  done
fi

# Also check main workspace database if it exists
if [ -f "$WS_DB" ]; then
  echo "   Checking main workspace database"
  sqlite3 "$WS_DB" "SELECT DISTINCT substr(key, 10, 36) FROM cursorDiskKV WHERE key LIKE 'bubbleId:%';" 2>/dev/null >> "$ALL_WS_CONVS"
fi

sort -u "$ALL_WS_CONVS" > /tmp/ws_convs_sorted.txt
sqlite3 "$GLOBAL_DB" "SELECT DISTINCT substr(key, 10, 36) FROM cursorDiskKV WHERE key LIKE 'bubbleId:%';" 2>/dev/null | sort > /tmp/global_convs.txt
comm -13 /tmp/global_convs.txt /tmp/ws_convs_sorted.txt > /tmp/missing_convs.txt
MISSING_COUNT=$(wc -l < /tmp/missing_convs.txt | xargs)

if [ "$MISSING_COUNT" -eq 0 ]; then
  echo "   ✅ No missing conversations found!"
  echo "   All workspace conversations are already in global database."
  exit 0
fi

echo "   📊 Found $MISSING_COUNT missing conversations"
echo ""

# Show missing conversation IDs
echo "2️⃣ Missing conversation IDs:"
head -10 /tmp/missing_convs.txt
if [ "$MISSING_COUNT" -gt 10 ]; then
  echo "   ... and $((MISSING_COUNT - 10)) more"
fi
echo ""

# Create comprehensive SQL script
echo "3️⃣ Creating restoration SQL script..."
echo "BEGIN TRANSACTION;" > /tmp/restore_all.sql

RESTORED_COUNT=0
while IFS= read -r conv_id; do
  echo "-- Conversation: $conv_id" >> /tmp/restore_all.sql
  
  # Find which workspace database has this conversation
  SOURCE_DB=""
  if [ -d "$WS_STORAGE_DIR" ]; then
    for ws_dir in "$WS_STORAGE_DIR"/*/; do
      if [ -d "$ws_dir" ]; then
        WS_DB_FILE="$ws_dir/state.vscdb"
        if [ -f "$WS_DB_FILE" ]; then
          EXISTS=$(sqlite3 "$WS_DB_FILE" "SELECT 1 FROM cursorDiskKV WHERE key LIKE 'bubbleId:${conv_id}:%' LIMIT 1;" 2>/dev/null || echo "")
          if [ -n "$EXISTS" ]; then
            SOURCE_DB="$WS_DB_FILE"
            break
          fi
        fi
      fi
    done
  fi
  
  # Fallback to main workspace database
  if [ -z "$SOURCE_DB" ] && [ -f "$WS_DB" ]; then
    EXISTS=$(sqlite3 "$WS_DB" "SELECT 1 FROM cursorDiskKV WHERE key LIKE 'bubbleId:${conv_id}:%' LIMIT 1;" 2>/dev/null || echo "")
    if [ -n "$EXISTS" ]; then
      SOURCE_DB="$WS_DB"
    fi
  fi
  
  if [ -z "$SOURCE_DB" ]; then
    echo "-- WARNING: Could not find source database for conversation $conv_id" >> /tmp/restore_all.sql
    continue
  fi
  
  # Copy ALL entry types for this conversation from the source database
  # 1. bubbleId entries (messages)
  sqlite3 "$SOURCE_DB" "SELECT 'INSERT OR REPLACE INTO cursorDiskKV (key, value) VALUES (' || quote(key) || ', ' || quote(value) || ');' FROM cursorDiskKV WHERE key LIKE 'bubbleId:${conv_id}:%';" 2>/dev/null >> /tmp/restore_all.sql
  
  # 2. checkpointId entries (conversation checkpoints)
  sqlite3 "$SOURCE_DB" "SELECT 'INSERT OR REPLACE INTO cursorDiskKV (key, value) VALUES (' || quote(key) || ', ' || quote(value) || ');' FROM cursorDiskKV WHERE key LIKE 'checkpointId:${conv_id}:%';" 2>/dev/null >> /tmp/restore_all.sql
  
  # 3. messageRequestContext entries (message context)
  sqlite3 "$SOURCE_DB" "SELECT 'INSERT OR REPLACE INTO cursorDiskKV (key, value) VALUES (' || quote(key) || ', ' || quote(value) || ');' FROM cursorDiskKV WHERE key LIKE 'messageRequestContext:${conv_id}:%';" 2>/dev/null >> /tmp/restore_all.sql
  
  # 4. composerData entry (composer metadata)
  sqlite3 "$SOURCE_DB" "SELECT 'INSERT OR REPLACE INTO cursorDiskKV (key, value) VALUES (' || quote(key) || ', ' || quote(value) || ');' FROM cursorDiskKV WHERE key = 'composerData:${conv_id}';" 2>/dev/null >> /tmp/restore_all.sql
  
  # 5. workbench.panel.composerChatViewPane entries (UI state)
  sqlite3 "$SOURCE_DB" "SELECT 'INSERT OR REPLACE INTO cursorDiskKV (key, value) VALUES (' || quote(key) || ', ' || quote(value) || ');' FROM cursorDiskKV WHERE key LIKE 'workbench.panel.composerChatViewPane%${conv_id}%';" 2>/dev/null >> /tmp/restore_all.sql
  
  # 6. Any other keys containing this conversation ID
  sqlite3 "$SOURCE_DB" "SELECT 'INSERT OR REPLACE INTO cursorDiskKV (key, value) VALUES (' || quote(key) || ', ' || quote(value) || ');' FROM cursorDiskKV WHERE key LIKE '%${conv_id}%' AND key NOT LIKE 'bubbleId:${conv_id}:%' AND key NOT LIKE 'checkpointId:${conv_id}:%' AND key NOT LIKE 'messageRequestContext:${conv_id}:%' AND key != 'composerData:${conv_id}' AND key NOT LIKE 'workbench.panel.composerChatViewPane%${conv_id}%';" 2>/dev/null >> /tmp/restore_all.sql
  
  RESTORED_COUNT=$((RESTORED_COUNT + 1))
done < /tmp/missing_convs.txt

echo "COMMIT;" >> /tmp/restore_all.sql

SQL_LINES=$(wc -l < /tmp/restore_all.sql | xargs)
echo "   ✅ Created SQL script with $SQL_LINES lines"
echo ""

# Execute restoration
echo "4️⃣ Executing restoration..."
sqlite3 "$GLOBAL_DB" < /tmp/restore_all.sql

if [ $? -eq 0 ]; then
  echo "   ✅ Restoration executed successfully!"
else
  echo "   ❌ Restoration failed! Restoring backup..."
  cp "$GLOBAL_BACKUP" "$GLOBAL_DB"
  exit 1
fi
echo ""

# Verify restoration
echo "5️⃣ Verifying restoration..."
FINAL_COUNT=$(sqlite3 "$GLOBAL_DB" "SELECT COUNT(DISTINCT substr(key, 10, 36)) FROM cursorDiskKV WHERE key LIKE 'bubbleId:%';" 2>/dev/null || echo "0")
echo "   📊 Final conversation count in global database: $FINAL_COUNT"
echo ""

# Ensure composerData entries exist
echo "6️⃣ Ensuring composerData entries..."
sqlite3 "$GLOBAL_DB" << 'SQL_EOF'
BEGIN TRANSACTION;

INSERT OR REPLACE INTO cursorDiskKV (key, value)
SELECT
  'composerData:' || substr(key, 10, 36) as new_key,
  '{"_v":3,"type":2,"lastAccessed":' || (strftime('%s', 'now') * 1000) || '}' as new_value
FROM cursorDiskKV
WHERE key LIKE 'bubbleId:%'
  AND NOT EXISTS (
    SELECT 1 FROM cursorDiskKV c2 
    WHERE c2.key = 'composerData:' || substr(cursorDiskKV.key, 10, 36)
  )
GROUP BY substr(key, 10, 36);

COMMIT;
SQL_EOF

COMPOSER_COUNT=$(sqlite3 "$GLOBAL_DB" "SELECT COUNT(*) FROM cursorDiskKV WHERE key LIKE 'composerData:%';" 2>/dev/null || echo "0")
echo "   ✅ ComposerData entries: $COMPOSER_COUNT"
echo ""

# Set permissions
chmod 664 "$GLOBAL_DB"
echo "✅ Database permissions set"
echo ""

echo "⚔️ RESTORATION COMPLETE! ⚔️"
echo "════════════════════════════════════════════════════════════"
echo ""
echo "📊 RESULTS:"
echo "   - Conversations restored: $RESTORED_COUNT"
echo "   - Final global database count: $FINAL_COUNT"
echo "   - ComposerData entries: $COMPOSER_COUNT"
echo ""
echo "💾 Backup saved:"
echo "   $GLOBAL_BACKUP"
echo ""
echo "🎯 NEXT STEPS:"
echo "   1. Open Cursor"
echo "   2. Wait 3-5 minutes for full re-indexing"
echo "   3. Try these actions to force UI refresh:"
echo "      - Click the chat history icon (clock) in sidebar"
echo "      - Press Cmd+Shift+P → 'Reload Window'"
echo "      - Close and reopen Cursor completely"
echo "   4. If still not visible, check if chats appear after searching (Cmd+F in chat panel)"
echo ""
echo "⚔️ CONQUEST ACHIEVED! ⚔️"

