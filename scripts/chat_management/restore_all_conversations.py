#!/usr/bin/env python3
"""
Restore ALL conversations from backup to Cursor's current database.
This will restore all 180 conversations with 153,516 messages.
"""

import json
import sqlite3
import sys
from pathlib import Path
from datetime import datetime
from collections import defaultdict

def get_all_conversation_ids(backup_db):
    """Extract all unique conversation IDs from backup."""
    conn = sqlite3.connect(backup_db)
    cursor = conn.cursor()
    
    print("🔍 Extracting all conversation IDs from backup...\n")
    cursor.execute("SELECT key FROM cursorDiskKV WHERE key LIKE 'bubbleId:%';")
    
    conv_ids = set()
    for (key,) in cursor.fetchall():
        parts = key.split(':')
        if len(parts) >= 3:
            conv_id = parts[1]
            if len(conv_id) == 36:  # UUID length
                conv_ids.add(conv_id)
    
    conn.close()
    return list(conv_ids)

def get_conversation_bubbles(backup_db, conv_id):
    """Extract all bubbles for a specific conversation."""
    conn = sqlite3.connect(backup_db)
    cursor = conn.cursor()
    
    cursor.execute(f"SELECT key, value FROM cursorDiskKV WHERE key LIKE 'bubbleId:{conv_id}:%';")
    bubbles = [(row[0], row[1]) for row in cursor.fetchall()]
    
    conn.close()
    return bubbles

def restore_all_conversations():
    """
    Restore ALL conversations from backup to current Cursor database.
    """
    # Paths
    backup_db = "/Users/fahadkiani/Library/Application Support/Cursor/Backups/state.vscdb.backup"
    current_db = "/Users/fahadkiani/Library/Application Support/Cursor/User/globalStorage/state.vscdb"
    
    print(f"🔄 Restoring ALL conversations to Cursor...\n")
    print(f"📁 Backup: {backup_db}")
    print(f"📁 Current DB: {current_db}\n")
    
    if not Path(backup_db).exists():
        print(f"❌ Backup DB not found: {backup_db}")
        return False
    if not Path(current_db).exists():
        print(f"❌ Current DB not found: {current_db}")
        return False

    # Safety check
    response = input("⚠️  Have you closed Cursor completely? (yes/no): ").strip().lower()
    if response != 'yes':
        print("❌ Please close Cursor first using: pkill -f 'Cursor'")
        return False

    # Get all conversation IDs
    conv_ids = get_all_conversation_ids(backup_db)
    print(f"📊 Found {len(conv_ids)} conversations to restore\n")

    conn_current = None
    try:
        conn_current = sqlite3.connect(current_db)
        cursor_current = conn_current.cursor()
        
        # Get existing composerData
        cursor_current.execute("SELECT value FROM ItemTable WHERE key = 'composer.composerData';")
        composer_data_row = cursor_current.fetchone()
        
        composer_data = {}
        if composer_data_row:
            composer_data = json.loads(composer_data_row[0])
        
        if 'allComposers' not in composer_data:
            composer_data['allComposers'] = []
        
        # Track existing conversations
        existing_composer_ids = {c.get('id') or c.get('composerId') for c in composer_data['allComposers']}
        
        # Statistics
        total_restored = 0
        total_skipped = 0
        total_bubbles = 0
        errors = []
        
        print("🚀 Starting batch restoration...\n")
        print("=" * 80)
        
        for i, conv_id in enumerate(conv_ids, 1):
            try:
                # Check if conversation already exists
                cursor_current.execute(f"SELECT COUNT(*) FROM cursorDiskKV WHERE key LIKE 'bubbleId:{conv_id}:%';")
                existing_bubbles = cursor_current.fetchone()[0]
                
                if existing_bubbles > 0:
                    print(f"[{i}/{len(conv_ids)}] ⏭️  Skipping {conv_id[:8]}... (already has {existing_bubbles} bubbles)")
                    total_skipped += 1
                    continue
                
                # Get bubbles from backup
                bubbles = get_conversation_bubbles(backup_db, conv_id)
                
                if not bubbles:
                    print(f"[{i}/{len(conv_ids)}] ⚠️  No bubbles found for {conv_id[:8]}...")
                    continue
                
                # Insert bubbles
                for key, value in bubbles:
                    cursor_current.execute(
                        "INSERT OR REPLACE INTO cursorDiskKV (key, value) VALUES (?, ?);",
                        (key, value)
                    )
                
                conn_current.commit()
                total_bubbles += len(bubbles)
                
                # Add to composerData if not already there
                if conv_id not in existing_composer_ids:
                    new_composer_entry = {
                        "id": conv_id,
                        "composerId": conv_id,
                        "title": f"Restored Conversation {i}",
                        "lastSendTime": int(datetime.now().timestamp() * 1000),
                        "messages": [],
                        "model": "gpt-4",
                        "context": []
                    }
                    composer_data['allComposers'].append(new_composer_entry)
                    existing_composer_ids.add(conv_id)
                
                total_restored += 1
                
                # Progress update every 10 conversations
                if i % 10 == 0:
                    print(f"[{i}/{len(conv_ids)}] ✅ Restored {total_restored} conversations ({total_bubbles} bubbles total)")
                else:
                    print(f"[{i}/{len(conv_ids)}] ✅ {conv_id[:8]}... ({len(bubbles)} bubbles)")
                
            except Exception as e:
                error_msg = f"Error restoring {conv_id[:8]}...: {e}"
                print(f"[{i}/{len(conv_ids)}] ❌ {error_msg}")
                errors.append(error_msg)
                continue
        
        # Update composerData
        print(f"\n📝 Updating composerData index...")
        cursor_current.execute(
            "INSERT OR REPLACE INTO ItemTable (key, value) VALUES (?, ?);",
            ('composer.composerData', json.dumps(composer_data))
        )
        conn_current.commit()
        
        print("=" * 80)
        print(f"\n🎉 Batch restoration complete!\n")
        print(f"📊 Statistics:")
        print(f"   ✅ Restored: {total_restored} conversations")
        print(f"   ⏭️  Skipped: {total_skipped} conversations (already exist)")
        print(f"   💬 Total bubbles: {total_bubbles:,} messages")
        print(f"   ❌ Errors: {len(errors)}")
        
        if errors:
            print(f"\n⚠️  Errors encountered:")
            for error in errors[:10]:  # Show first 10 errors
                print(f"   {error}")
        
        print(f"\n🔄 Restart Cursor to see all restored conversations in chat history.\n")
        return True

    except sqlite3.OperationalError as e:
        print(f"❌ SQLite error for {current_db}: {e}")
        print("   This might mean Cursor is still running. Please ensure Cursor is completely closed.")
        return False
    except Exception as e:
        print(f"❌ An unexpected error occurred during restoration: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if conn_current:
            conn_current.close()

if __name__ == "__main__":
    restore_all_conversations()

