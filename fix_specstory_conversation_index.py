#!/usr/bin/env python3
"""
Fix SpecStory conversation indexing - add to composer.composerData properly.
This ensures the conversation shows up in Cursor's chat history.
"""

import json
import sqlite3
from pathlib import Path
from datetime import datetime, timedelta
import re

def extract_title_from_first_message(conn, conv_id):
    """Extract a meaningful title from the first user message."""
    cursor = conn.cursor()
    
    # Get first 20 bubbles to find a user message
    cursor.execute(f"SELECT value FROM cursorDiskKV WHERE key LIKE 'bubbleId:{conv_id}:%' AND value IS NOT NULL ORDER BY key LIMIT 20;")
    
    for (value,) in cursor.fetchall():
        try:
            data = json.loads(value)
            text = data.get('text', '').strip()
            msg_type = data.get('type', -1)
            
            # Type 0 = user message
            if msg_type == 0 and text and len(text) > 10:
                # Clean up the text for title
                title = text[:100].strip()
                # Remove markdown code blocks
                title = re.sub(r'```[\s\S]*?```', '', title)
                # Remove multiple spaces
                title = re.sub(r'\s+', ' ', title)
                # Take first sentence or first 60 chars
                if '.' in title:
                    title = title.split('.')[0]
                title = title[:60].strip()
                
                if title:
                    return title
        except:
            continue
    
    return f"Designing zero-tolerance content filter"

def fix_conversation_index(conv_id):
    """
    Add conversation to composer.composerData index in ItemTable.
    """
    db_path = Path.home() / "Library/Application Support/Cursor/User/globalStorage/state.vscdb"
    
    print(f"🔧 Fixing conversation index for: {conv_id}\n")
    
    if not db_path.exists():
        print(f"❌ Database not found: {db_path}")
        return False
    
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    
    # Check if conversation has bubbles
    cursor.execute(f"SELECT COUNT(*) FROM cursorDiskKV WHERE key LIKE 'bubbleId:{conv_id}:%';")
    bubble_count = cursor.fetchone()[0]
    
    if bubble_count == 0:
        print(f"❌ No bubbles found for conversation {conv_id}")
        conn.close()
        return False
    
    print(f"✅ Found {bubble_count} bubbles\n")
    
    # Get existing composerData from ItemTable
    cursor.execute("SELECT value FROM ItemTable WHERE key = 'composer.composerData';")
    row = cursor.fetchone()
    
    if row:
        composer_data = json.loads(row[0])
    else:
        print("⚠️  No composerData found - creating new structure...")
        composer_data = {'allComposers': []}
    
    if 'allComposers' not in composer_data:
        composer_data['allComposers'] = []
    
    # Check if conversation already in index
    existing_ids = {c.get('composerId'): c for c in composer_data['allComposers'] if c.get('composerId')}
    
    if conv_id in existing_ids:
        print(f"⚠️  Conversation already in index - updating timestamp...")
        existing = existing_ids[conv_id]
        existing['lastUpdatedAt'] = int(datetime.now().timestamp() * 1000)
        updated = True
    else:
        print(f"➕ Adding conversation to index...")
        # Extract title from first message
        title = extract_title_from_first_message(conn, conv_id)
        
        # Create entry matching Cursor's format
        now_ms = int(datetime.now().timestamp() * 1000)
        new_entry = {
            "type": "head",
            "composerId": conv_id,
            "lastUpdatedAt": now_ms,
            "createdAt": now_ms - 3600000,  # 1 hour before
            "unifiedMode": "agent",
            "forceMode": "edit",
            "hasUnreadMessages": False,
            "contextUsagePercent": 0.0,
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
            "hasBlockingPendingMessages": False
        }
        
        composer_data['allComposers'].append(new_entry)
        updated = False
    
    # Update ItemTable
    print("💾 Updating composer.composerData in ItemTable...")
    cursor.execute(
        "INSERT OR REPLACE INTO ItemTable (key, value) VALUES (?, ?);",
        ('composer.composerData', json.dumps(composer_data))
    )
    conn.commit()
    
    print(f"✅ Successfully {'updated' if updated else 'added'} conversation to index!")
    print(f"📊 Total conversations in index: {len(composer_data['allComposers'])}\n")
    
    conn.close()
    return True

if __name__ == "__main__":
    # The conversation ID from the SpecStory file
    conv_id = "5dac5408-f750-44f5-bbe6-713f7fb5e103"
    
    print("⚔️ FIX SPECSTORY CONVERSATION INDEX ⚔️")
    print("=" * 60)
    print()
    
    success = fix_conversation_index(conv_id)
    
    if success:
        print("=" * 60)
        print("✅ Fix complete!")
        print()
        print("🔄 Next steps:")
        print("   1. Close Cursor completely: pkill -f 'Cursor'")
        print("   2. Wait 5 seconds")
        print("   3. Restart Cursor")
        print("   4. Wait 1-2 minutes for indexing")
        print("   5. Check chat history panel")
    else:
        print("❌ Fix failed. Check errors above.")














