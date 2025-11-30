#!/usr/bin/env python3
"""
Fix restored conversations visibility by:
1. Extracting real titles from first user message
2. Setting proper timestamps (spread them out chronologically)
3. Ensuring all metadata is correct
"""

import json
import sqlite3
import re
from pathlib import Path
from datetime import datetime, timedelta

def extract_title_from_conversation(backup_db, conv_id):
    """Extract a meaningful title from the first user message."""
    conn = sqlite3.connect(backup_db)
    cursor = conn.cursor()
    
    # Get first few bubbles to find a user message
    cursor.execute(f"SELECT value FROM cursorDiskKV WHERE key LIKE 'bubbleId:{conv_id}:%' ORDER BY key LIMIT 20;")
    
    for (value,) in cursor.fetchall():
        try:
            data = json.loads(value)
            text = data.get('text', '').strip()
            msg_type = data.get('type', -1)
            
            # Type 0 = user message
            if msg_type == 0 and text and len(text) > 10:
                # Extract first sentence or first 50 chars
                first_line = text.split('\n')[0].strip()
                if len(first_line) > 60:
                    first_line = first_line[:60] + "..."
                
                # Clean up common prefixes
                first_line = re.sub(r'^(can you|could you|please|help me|i need|i want)\s+', '', first_line, flags=re.IGNORECASE)
                first_line = first_line.strip()
                
                if first_line:
                    conn.close()
                    return first_line
        except:
            continue
    
    conn.close()
    return f"Restored Conversation {conv_id[:8]}"

def fix_conversations_visibility():
    """
    Fix restored conversations to make them visible in Cursor UI.
    """
    backup_db = "/Users/fahadkiani/Library/Application Support/Cursor/Backups/state.vscdb.backup"
    current_db = "/Users/fahadkiani/Library/Application Support/Cursor/User/globalStorage/state.vscdb"
    
    print("🔧 Fixing restored conversations visibility...\n")
    
    if not Path(current_db).exists():
        print(f"❌ Current DB not found: {current_db}")
        return False
    
    conn = sqlite3.connect(current_db)
    cursor = conn.cursor()
    
    # Get composerData
    cursor.execute("SELECT value FROM ItemTable WHERE key = 'composer.composerData';")
    composer_row = cursor.fetchone()
    
    if not composer_row:
        print("❌ No composerData found!")
        conn.close()
        return False
    
    composer_data = json.loads(composer_row[0])
    all_composers = composer_data.get('allComposers', [])
    
    # Find all restored conversations
    restored = [c for c in all_composers if 'Restored' in c.get('title', '')]
    print(f"📊 Found {len(restored)} restored conversations to fix\n")
    
    # Extract real titles and update
    print("📝 Extracting real titles from conversations...\n")
    
    # Base date: 4 months ago (to match your "4mo ago" section)
    base_date = datetime.now() - timedelta(days=120)
    
    updated_count = 0
    for i, composer in enumerate(restored):
        conv_id = composer.get('id') or composer.get('composerId')
        if not conv_id:
            continue
        
        # Extract real title
        real_title = extract_title_from_conversation(backup_db, conv_id)
        
        # Set timestamp (spread them out over 4-5 months)
        days_ago = 120 - (i % 30)  # Spread over ~4 months
        timestamp = int((base_date + timedelta(days=days_ago)).timestamp() * 1000)
        
        # Update composer entry
        composer['title'] = real_title
        composer['lastSendTime'] = timestamp
        composer['createdAt'] = timestamp  # Add createdAt if missing
        
        updated_count += 1
        
        if (i + 1) % 20 == 0:
            print(f"   ✅ Updated {i + 1}/{len(restored)} conversations...")
    
    print(f"\n✅ Updated {updated_count} conversations with real titles and timestamps\n")
    
    # Save updated composerData
    composer_data['allComposers'] = all_composers
    cursor.execute(
        "INSERT OR REPLACE INTO ItemTable (key, value) VALUES (?, ?);",
        ('composer.composerData', json.dumps(composer_data))
    )
    conn.commit()
    conn.close()
    
    print("✅ composerData updated successfully!")
    print("\n🔄 Restart Cursor to see the conversations with proper titles and dates.\n")
    return True

if __name__ == "__main__":
    fix_conversations_visibility()




