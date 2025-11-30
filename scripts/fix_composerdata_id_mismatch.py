#!/usr/bin/env python3
"""
Fix composerData ID mismatch - remove empty entries and add missing ones.
This ensures all conversations with bubbles have composerData entries.
"""

import json
import sqlite3
import sys
from pathlib import Path
from datetime import datetime

def fix_composerdata_mismatch():
    """
    Fix composerData by:
    1. Removing entries that have no bubbles (empty conversations)
    2. Adding entries for conversations that have bubbles but no composerData
    """
    global_db = "/Users/fahadkiani/Library/Application Support/Cursor/User/globalStorage/state.vscdb"
    
    print(f"🔧 Fixing composerData ID mismatch...\n")
    print(f"📁 Global DB: {global_db}\n")
    
    if not Path(global_db).exists():
        print(f"❌ Global DB not found: {global_db}")
        return False

    # Safety check
    response = input("⚠️  Have you closed Cursor completely? (yes/no): ")
    if response.lower() != 'yes':
        print("Please close Cursor completely before proceeding.")
        return False

    # Create a backup
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = Path(f"{global_db}.before_fix_mismatch_{timestamp}")
    try:
        import shutil
        shutil.copy(global_db, backup_path)
        print(f"✅ Created backup: {backup_path}\n")
    except Exception as e:
        print(f"❌ Error creating backup: {e}")
        return False

    conn = None
    try:
        conn = sqlite3.connect(str(global_db))
        cursor = conn.cursor()
        
        # 1. Get all conversation IDs from bubbles
        cursor.execute("SELECT DISTINCT key FROM cursorDiskKV WHERE key LIKE 'bubbleId:%';")
        all_bubble_keys = cursor.fetchall()
        
        conv_ids_from_bubbles = set()
        for (key,) in all_bubble_keys:
            parts = key.split(':')
            if len(parts) >= 2:
                conv_id = parts[1]  # Full UUID after bubbleId:
                conv_ids_from_bubbles.add(conv_id)
        
        print(f"📊 Found {len(conv_ids_from_bubbles)} conversations with bubbles\n")
        
        # 2. Get existing composerData
        cursor.execute("SELECT value FROM ItemTable WHERE key = 'composer.composerData';")
        composer_row = cursor.fetchone()
        
        composer_data = {"allComposers": []}
        if composer_row:
            composer_data = json.loads(composer_row[0])
        
        existing_composers = composer_data.get('allComposers', [])
        
        # 3. Build map of existing composerData entries
        existing_composer_map = {}
        for c in existing_composers:
            comp_id = c.get('composerId') or c.get('id', '')
            comp_id = comp_id.lstrip(':')  # Remove leading colon
            if comp_id:
                existing_composer_map[comp_id] = c
        
        print(f"📋 Found {len(existing_composer_map)} existing composerData entries\n")
        
        # 4. Keep only composerData entries that have bubbles
        valid_composers = []
        removed_count = 0
        
        for comp_id, composer_entry in existing_composer_map.items():
            if comp_id in conv_ids_from_bubbles:
                valid_composers.append(composer_entry)
            else:
                removed_count += 1
        
        print(f"✅ Kept {len(valid_composers)} valid entries (with bubbles)")
        print(f"🗑️  Removed {removed_count} empty entries (no bubbles)\n")
        
        # 5. Add composerData entries for conversations with bubbles but no entry
        added_count = 0
        existing_comp_ids = {c.get('composerId') or c.get('id', '').lstrip(':') for c in valid_composers}
        
        for conv_id in conv_ids_from_bubbles:
            if conv_id not in existing_comp_ids:
                # Create new composerData entry
                # Try to get first bubble to extract metadata
                cursor.execute("SELECT value FROM cursorDiskKV WHERE key LIKE ? LIMIT 1;", (f'bubbleId:{conv_id}:%',))
                bubble_row = cursor.fetchone()
                
                # Try to extract title from first message if available
                title = f"Restored Conversation {conv_id[:8]}..."
                if bubble_row:
                    try:
                        bubble_data = json.loads(bubble_row[0])
                        # Try to extract title from message content
                        if isinstance(bubble_data, dict):
                            content = bubble_data.get('content', '')
                            if content and len(content) > 0:
                                # Use first 50 chars as title
                                title = content[:50].replace('\n', ' ').strip()
                                if not title:
                                    title = f"Restored Conversation {conv_id[:8]}..."
                    except:
                        pass
                
                new_composer_entry = {
                    "id": conv_id,
                    "composerId": conv_id,
                    "name": title,
                    "type": "head",
                    "lastUpdatedAt": datetime.now().isoformat(),
                    "createdAt": datetime.now().isoformat(),
                    "isArchived": False,
                    "isDeleted": False,
                    "isDraft": False,
                    "isSpec": False,
                    "context": []
                }
                valid_composers.append(new_composer_entry)
                added_count += 1
        
        print(f"➕ Added {added_count} new composerData entries (had bubbles, missing from index)\n")
        
        # 6. Update composerData
        composer_data['allComposers'] = valid_composers
        
        cursor.execute(
            "INSERT OR REPLACE INTO ItemTable (key, value) VALUES (?, ?);",
            ('composer.composerData', json.dumps(composer_data))
        )
        conn.commit()
        
        print("🎉 ComposerData fixed successfully!")
        print(f"\n📊 Final Statistics:")
        print(f"   ✅ Total conversations: {len(valid_composers)}")
        print(f"   ✅ All have bubbles: {len(valid_composers)}")
        print(f"   🗑️  Removed empty: {removed_count}")
        print(f"   ➕ Added missing: {added_count}")
        print("\n🔄 Restart Cursor to see all conversations!\n")
        
        return True

    except sqlite3.Error as e:
        print(f"❌ SQLite error: {e}")
        return False
    except json.JSONDecodeError as e:
        print(f"❌ JSON decoding error: {e}")
        return False
    except Exception as e:
        print(f"❌ An unexpected error occurred: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if conn:
            conn.close()

if __name__ == "__main__":
    fix_composerdata_mismatch()







