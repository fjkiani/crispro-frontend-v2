#!/usr/bin/env python3
"""
⚔️ ADD CONVERSATION TO CURSOR INDEX ⚔️

Adds imported conversation to Cursor's composer.composerData index so it appears in chat history.
"""

import sqlite3
import json
from pathlib import Path
from datetime import datetime

CURSOR_DB = Path.home() / "Library/Application Support/Cursor/User/globalStorage/state.vscdb"
CONVERSATION_ID = "5dac5408-f750-44f5-bbe6-713f7fb5e103"

def add_to_index():
    """Add conversation to composer.composerData index."""
    print("⚔️ ADDING CONVERSATION TO CURSOR INDEX ⚔️")
    print("=" * 60)
    
    if not CURSOR_DB.exists():
        print(f"❌ Database not found: {CURSOR_DB}")
        return False
    
    conn = sqlite3.connect(str(CURSOR_DB))
    cursor = conn.cursor()
    
    try:
        # Get existing composerData from ItemTable
        cursor.execute("SELECT value FROM ItemTable WHERE key = 'composer.composerData';")
        row = cursor.fetchone()
        
        if row:
            composer_data = json.loads(row[0])
        else:
            # Create new structure if it doesn't exist
            composer_data = {
                "_v": 10,
                "allComposers": []
            }
        
        # Check if conversation already in index
        existing = any(
            c.get('composerId') == CONVERSATION_ID or 
            c.get('id') == CONVERSATION_ID 
            for c in composer_data.get('allComposers', [])
        )
        
        if existing:
            print(f"⚠️  Conversation {CONVERSATION_ID} already in index")
            print("   🔄 Updating timestamp...")
            # Update lastUpdatedAt
            for c in composer_data.get('allComposers', []):
                if c.get('composerId') == CONVERSATION_ID or c.get('id') == CONVERSATION_ID:
                    c['lastUpdatedAt'] = int(datetime.now().timestamp() * 1000)
                    c['lastSendTime'] = int(datetime.now().timestamp() * 1000)
                    break
        else:
            # Get title from composerData in cursorDiskKV
            cursor.execute(
                "SELECT value FROM cursorDiskKV WHERE key = ?;",
                (f'composerData:{CONVERSATION_ID}',)
            )
            row = cursor.fetchone()
            title = "Designing a zero-tolerance content filter"
            if row:
                try:
                    data = json.loads(row[0])
                    title = data.get('title', title)
                except:
                    pass
            
            # Get first bubble timestamp for createdAt
            cursor.execute(
                "SELECT value FROM cursorDiskKV WHERE key LIKE ? ORDER BY key LIMIT 1;",
                (f'bubbleId:{CONVERSATION_ID}:%',)
            )
            row = cursor.fetchone()
            created_at = int(datetime.now().timestamp() * 1000)
            if row:
                try:
                    bubble_data = json.loads(row[0])
                    created_at = bubble_data.get('timestamp', created_at)
                except:
                    pass
            
            # Add new entry
            new_entry = {
                "type": "chat",
                "composerId": CONVERSATION_ID,
                "id": CONVERSATION_ID,
                "createdAt": created_at,
                "lastUpdatedAt": int(datetime.now().timestamp() * 1000),
                "lastSendTime": int(datetime.now().timestamp() * 1000),
                "unifiedMode": True,
                "forceMode": None,
                "hasUnreadMessages": False,
                "contextUsagePercent": 0,
                "totalLinesAdded": 0,
                "totalLinesRemoved": 0,
                "subtitle": "",
                "isArchived": False,
                "isWorktree": False,
                "isSpec": False,
                "numSubComposers": 0,
                "isDraft": False,
                "referencedPlans": [],
                "filesChangedCount": 0,
                "name": title,
                "hasBlockingPendingActions": False
            }
            
            if 'allComposers' not in composer_data:
                composer_data['allComposers'] = []
            
            composer_data['allComposers'].append(new_entry)
            print(f"✅ Added conversation to index: {title}")
        
        # Update ItemTable
        cursor.execute(
            "INSERT OR REPLACE INTO ItemTable (key, value) VALUES (?, ?);",
            ('composer.composerData', json.dumps(composer_data))
        )
        conn.commit()
        
        print(f"✅ Index updated successfully!")
        print(f"📊 Total conversations in index: {len(composer_data.get('allComposers', []))}")
        print("\n🎯 NEXT STEPS:")
        print("   1. Close Cursor completely: pkill -x Cursor")
        print("   2. Reopen Cursor: open -a Cursor")
        print("   3. Wait 30 seconds")
        print("   4. Check chat history panel")
        
        return True
        
    except Exception as e:
        conn.rollback()
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        conn.close()

if __name__ == "__main__":
    add_to_index()















