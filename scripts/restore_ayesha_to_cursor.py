#!/usr/bin/env python3
"""
Restore the Ayesha conversation to Cursor's current database.
This will insert the conversation bubbles into the current state.vscdb.
"""

import json
import sqlite3
import sys
from pathlib import Path
from datetime import datetime

def restore_ayesha_conversation():
    """
    Restore Ayesha conversation from backup to current Cursor database.
    """
    # Paths
    backup_db = "/Users/fahadkiani/Library/Application Support/Cursor/Backups/state.vscdb.backup"
    current_db = "/Users/fahadkiani/Library/Application Support/Cursor/User/globalStorage/state.vscdb"
    ayesha_conv_id = "46d02246-1295-4920-80f6-7c62c2f2b158"
    
    print(f"🔄 Restoring Ayesha conversation to Cursor...\n")
    print(f"📁 Backup: {backup_db}")
    print(f"📁 Current DB: {current_db}")
    print(f"🆔 Conversation ID: {ayesha_conv_id}\n")
    
    # Check if conversation already exists in current DB
    try:
        conn_current = sqlite3.connect(current_db)
        cursor_current = conn_current.cursor()
        
        cursor_current.execute(f"SELECT COUNT(*) FROM cursorDiskKV WHERE key LIKE 'bubbleId:{ayesha_conv_id}:%';")
        existing_count = cursor_current.fetchone()[0]
        
        if existing_count > 0:
            print(f"⚠️  WARNING: Conversation already exists in current DB ({existing_count} bubbles)")
            response = input("Do you want to overwrite? (yes/no): ")
            if response.lower() != 'yes':
                print("❌ Cancelled")
                conn_current.close()
                return False
            else:
                # Delete existing bubbles
                print(f"🗑️  Deleting {existing_count} existing bubbles...")
                cursor_current.execute(f"DELETE FROM cursorDiskKV WHERE key LIKE 'bubbleId:{ayesha_conv_id}:%';")
                conn_current.commit()
        
        conn_current.close()
    except Exception as e:
        print(f"⚠️  Error checking current DB: {e}")
        print("   Continuing anyway...\n")
    
    # Extract bubbles from backup
    try:
        conn_backup = sqlite3.connect(backup_db)
        cursor_backup = conn_backup.cursor()
        
        print("📊 Extracting bubbles from backup...")
        cursor_backup.execute(f"SELECT key, value FROM cursorDiskKV WHERE key LIKE 'bubbleId:{ayesha_conv_id}:%';")
        
        bubbles_to_restore = []
        for key, value in cursor_backup:
            bubbles_to_restore.append((key, value))
        
        conn_backup.close()
        
        print(f"✅ Found {len(bubbles_to_restore)} bubbles to restore\n")
        
        if len(bubbles_to_restore) == 0:
            print("❌ No bubbles found in backup!")
            return False
        
        # Insert into current DB
        print("💾 Inserting bubbles into current database...")
        conn_current = sqlite3.connect(current_db)
        cursor_current = conn_current.cursor()
        
        # Use executemany for batch insert
        cursor_current.executemany(
            "INSERT OR REPLACE INTO cursorDiskKV (key, value) VALUES (?, ?);",
            bubbles_to_restore
        )
        
        conn_current.commit()
        conn_current.close()
        
        print(f"✅ Successfully restored {len(bubbles_to_restore)} bubbles!\n")
        
        # Verify
        print("🔍 Verifying restoration...")
        conn_current = sqlite3.connect(current_db)
        cursor_current = conn_current.cursor()
        
        cursor_current.execute(f"SELECT COUNT(*) FROM cursorDiskKV WHERE key LIKE 'bubbleId:{ayesha_conv_id}:%';")
        restored_count = cursor_current.fetchone()[0]
        
        conn_current.close()
        
        if restored_count == len(bubbles_to_restore):
            print(f"✅ Verification passed: {restored_count} bubbles in current DB\n")
            
            # Now add conversation to composerData so Cursor 2.0 can index it
            print("📝 Adding conversation to composerData index...")
            try:
                conn_current = sqlite3.connect(current_db)
                cursor_current = conn_current.cursor()
                
                # Get current composerData
                cursor_current.execute("SELECT value FROM ItemTable WHERE key = 'composer.composerData';")
                row = cursor_current.fetchone()
                
                if row:
                    composer_data = json.loads(row[0])
                    
                    # Check if conversation already in allComposers
                    all_composers = composer_data.get('allComposers', [])
                    existing = any(c.get('composerId') == ayesha_conv_id or c.get('id') == ayesha_conv_id for c in all_composers)
                    
                    if not existing:
                        # Add new composer entry
                        new_composer = {
                            'type': 'chat',
                            'composerId': ayesha_conv_id,
                            'id': ayesha_conv_id,
                            'createdAt': int(datetime.now().timestamp() * 1000),  # Current timestamp
                            'lastUpdatedAt': int(datetime.now().timestamp() * 1000),
                            'lastSendTime': int(datetime.now().timestamp() * 1000),
                            'unifiedMode': True,
                            'forceMode': None,
                            'hasUnreadMessages': False,
                            'contextUsagePercent': 0,
                            'totalLinesAdded': 0,
                            'totalLinesRemoved': 0,
                            'subtitle': 'Ayesha Conversation (Restored)',
                            'isArchived': False,
                            'isWorktree': False,
                            'isSpec': False,
                            'numSubComposers': 0,
                            'messages': []  # Messages are in bubbles, not here
                        }
                        
                        all_composers.append(new_composer)
                        composer_data['allComposers'] = all_composers
                        
                        # Update ItemTable
                        cursor_current.execute(
                            "UPDATE ItemTable SET value = ? WHERE key = 'composer.composerData';",
                            (json.dumps(composer_data),)
                        )
                        conn_current.commit()
                        
                        print(f"✅ Added conversation to composerData index\n")
                    else:
                        print(f"⚠️  Conversation already in composerData index\n")
                else:
                    print(f"⚠️  No composerData found - creating new entry...")
                    # Create new composerData
                    composer_data = {
                        'allComposers': [{
                            'type': 'chat',
                            'composerId': ayesha_conv_id,
                            'id': ayesha_conv_id,
                            'createdAt': int(datetime.now().timestamp() * 1000),
                            'lastUpdatedAt': int(datetime.now().timestamp() * 1000),
                            'lastSendTime': int(datetime.now().timestamp() * 1000),
                            'unifiedMode': True,
                            'subtitle': 'Ayesha Conversation (Restored)',
                            'isArchived': False,
                            'messages': []
                        }]
                    }
                    
                    cursor_current.execute(
                        "INSERT INTO ItemTable (key, value) VALUES ('composer.composerData', ?);",
                        (json.dumps(composer_data),)
                    )
                    conn_current.commit()
                    print(f"✅ Created new composerData entry\n")
                
                conn_current.close()
            except Exception as e:
                print(f"⚠️  Warning: Could not update composerData: {e}")
                print("   Bubbles are restored, but conversation may not appear in UI")
                print("   You may need to manually trigger a refresh\n")
            
            print("🎉 Ayesha conversation restored successfully!")
            print("   Restart Cursor to see the conversation in the chat history.\n")
            return True
        else:
            print(f"⚠️  Verification warning: Expected {len(bubbles_to_restore)}, found {restored_count}")
            return False
        
    except Exception as e:
        print(f"❌ Error during restoration: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("⚠️  IMPORTANT: Close Cursor completely before running this script!")
    print("   (Otherwise the database will be locked)\n")
    
    response = input("Have you closed Cursor? (yes/no): ")
    if response.lower() != 'yes':
        print("❌ Please close Cursor and try again")
        sys.exit(1)
    
    success = restore_ayesha_conversation()
    
    if success:
        print("\n✅ Next steps:")
        print("   1. Restart Cursor")
        print("   2. Check the chat history panel")
        print("   3. Look for conversation ID: 46d02246-1295-4920-80f6-7c62c2f2b158")
    else:
        print("\n❌ Restoration failed. Check errors above.")

