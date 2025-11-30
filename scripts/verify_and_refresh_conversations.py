#!/usr/bin/env python3
"""
Verify restored conversations are in database and ensure they're properly indexed.
Also extracts real titles from first messages.
"""

import json
import sqlite3
from pathlib import Path
from datetime import datetime, timedelta

def verify_and_refresh():
    """
    Verify conversations are in DB and refresh their metadata.
    """
    backup_db = "/Users/fahadkiani/Library/Application Support/Cursor/Backups/state.vscdb.backup"
    current_db = "/Users/fahadkiani/Library/Application Support/Cursor/User/globalStorage/state.vscdb"
    
    print("🔍 Verifying restored conversations...\n")
    
    if not Path(current_db).exists():
        print(f"❌ Current DB not found!")
        return False
    
    conn = sqlite3.connect(current_db)
    cursor = conn.cursor()
    
    # Count conversations with bubbles
    cursor.execute("SELECT COUNT(DISTINCT substr(key, 9, 36)) FROM cursorDiskKV WHERE key LIKE 'bubbleId:%' AND length(key) > 45 AND value IS NOT NULL;")
    conv_count = cursor.fetchone()[0]
    print(f"✅ Conversations with bubbles: {conv_count}")
    
    # Get composerData
    cursor.execute("SELECT value FROM ItemTable WHERE key = 'composer.composerData';")
    composer_row = cursor.fetchone()
    
    if not composer_row:
        print("❌ No composerData found!")
        conn.close()
        return False
    
    composer_data = json.loads(composer_row[0])
    all_composers = composer_data.get('allComposers', [])
    print(f"✅ Conversations in composerData: {len(all_composers)}\n")
    
    # Find Ayesha conversation specifically
    ayesha_id = '46d02246-1295-4920-80f6-7c62c2f2b158'
    ayesha = next((c for c in all_composers if c.get('id') == ayesha_id or c.get('composerId') == ayesha_id), None)
    
    if ayesha:
        print(f"🎯 Ayesha conversation found:")
        print(f"   Title: {ayesha.get('title', 'No title')}")
        print(f"   ID: {ayesha.get('id', 'N/A')}")
        
        # Check bubbles
        cursor.execute(f"SELECT COUNT(*) FROM cursorDiskKV WHERE key LIKE 'bubbleId:{ayesha_id}:%' AND value IS NOT NULL;")
        bubble_count = cursor.fetchone()[0]
        print(f"   Bubbles: {bubble_count}")
        
        # Try to extract first user message for title
        cursor.execute(f"SELECT value FROM cursorDiskKV WHERE key LIKE 'bubbleId:{ayesha_id}:%' AND value IS NOT NULL ORDER BY key LIMIT 50;")
        for (value,) in cursor.fetchall():
            try:
                data = json.loads(value)
                if data.get('type') == 0 and data.get('text', '').strip():
                    first_msg = data.get('text', '').strip()[:80]
                    print(f"   First message: {first_msg}...")
                    
                    # Update title if it's still generic
                    if 'Restored' in ayesha.get('title', ''):
                        new_title = first_msg[:60].replace('\n', ' ').strip()
                        if len(new_title) > 60:
                            new_title = new_title[:60] + "..."
                        ayesha['title'] = new_title
                        print(f"   ✏️  Updated title to: {new_title}")
                    break
            except:
                continue
    else:
        print(f"⚠️  Ayesha conversation NOT in composerData!")
        print(f"   Adding it now...")
        
        # Check if bubbles exist
        cursor.execute(f"SELECT COUNT(*) FROM cursorDiskKV WHERE key LIKE 'bubbleId:{ayesha_id}:%' AND value IS NOT NULL;")
        bubble_count = cursor.fetchone()[0]
        
        if bubble_count > 0:
            # Add to composerData
            new_entry = {
                "id": ayesha_id,
                "composerId": ayesha_id,
                "title": "Ayesha Conversation (Restored)",
                "lastSendTime": int((datetime.now() - timedelta(days=120)).timestamp() * 1000),
                "messages": [],
                "model": "gpt-4",
                "context": [],
                "createdAt": int((datetime.now() - timedelta(days=120)).timestamp() * 1000)
            }
            all_composers.append(new_entry)
            print(f"   ✅ Added Ayesha conversation to composerData")
    
    # Update composerData
    composer_data['allComposers'] = all_composers
    cursor.execute(
        "INSERT OR REPLACE INTO ItemTable (key, value) VALUES (?, ?);",
        ('composer.composerData', json.dumps(composer_data))
    )
    conn.commit()
    conn.close()
    
    print(f"\n✅ Verification complete!")
    print(f"\n🔄 NEXT STEPS:")
    print(f"   1. Close Cursor completely (pkill -f 'Cursor')")
    print(f"   2. Restart Cursor")
    print(f"   3. Check chat history panel - conversations should appear")
    print(f"   4. If still not visible, try searching for 'Ayesha' in the search bar\n")
    
    return True

if __name__ == "__main__":
    verify_and_refresh()




