#!/usr/bin/env python3
"""
Rebuild composerData from all conversations that have bubbles.
This ensures all conversations are properly indexed and visible in Cursor.
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
    
    return f"Conversation {conv_id[:8]}"

def rebuild_composerdata():
    """
    Rebuild composerData from all conversations that have bubbles.
    """
    current_db = "/Users/fahadkiani/Library/Application Support/Cursor/User/globalStorage/state.vscdb"
    backup_db = "/Users/fahadkiani/Library/Application Support/Cursor/Backups/state.vscdb.backup"
    
    print("🔧 Rebuilding composerData from all conversations with bubbles...\n")
    
    if not Path(current_db).exists():
        print(f"❌ Current DB not found!")
        return False
    
    conn = sqlite3.connect(current_db)
    cursor = conn.cursor()
    
    # Get all unique conversation IDs that have bubbles
    print("📊 Finding all conversations with bubbles...")
    cursor.execute("""
        SELECT DISTINCT substr(key, 9, 36) as conv_id 
        FROM cursorDiskKV 
        WHERE key LIKE 'bubbleId:%' 
        AND length(key) > 45 
        AND value IS NOT NULL
        ORDER BY conv_id;
    """)
    
    conv_ids_with_bubbles = [row[0] for row in cursor.fetchall()]
    print(f"✅ Found {len(conv_ids_with_bubbles)} conversations with bubbles\n")
    
    # Get existing composerData
    cursor.execute("SELECT value FROM ItemTable WHERE key = 'composer.composerData';")
    composer_row = cursor.fetchone()
    
    composer_data = {}
    if composer_row:
        composer_data = json.loads(composer_row[0])
    
    if 'allComposers' not in composer_data:
        composer_data['allComposers'] = []
    
    # Create lookup of existing entries (use composerId, not id)
    existing_ids = {c.get('composerId'): c for c in composer_data['allComposers'] if c.get('composerId')}
    
    print(f"📊 Existing entries in composerData: {len(existing_ids)}\n")
    
    # Connect to backup to extract titles
    conn_backup = None
    if Path(backup_db).exists():
        conn_backup = sqlite3.connect(backup_db)
        print("📁 Using backup DB to extract titles\n")
    else:
        print("⚠️  Backup DB not found, using generic titles\n")
    
    # Add missing conversations to composerData
    added_count = 0
    updated_count = 0
    
    # Spread timestamps over the last 4 months
    base_time = datetime.now() - timedelta(days=120)
    time_increment = timedelta(days=120 / max(len(conv_ids_with_bubbles), 1))
    
    for i, conv_id in enumerate(conv_ids_with_bubbles):
        if conv_id in existing_ids:
            # Update existing entry's lastUpdatedAt
            existing = existing_ids[conv_id]
            existing['lastUpdatedAt'] = int((base_time + time_increment * i).timestamp() * 1000)
            updated_count += 1
        else:
            # Extract title from backup or use generic
            if conn_backup:
                title = extract_title_from_first_message(conn_backup, conv_id)
            else:
                # Try to extract from current DB
                title = extract_title_from_first_message(conn, conv_id)
                if not title or title.startswith("Conversation"):
                    title = f"Conversation {conv_id[:8]}"
            
            # Create new entry matching existing structure
            now_ms = int((base_time + time_increment * i).timestamp() * 1000)
            new_entry = {
                "type": "head",
                "composerId": conv_id,
                "lastUpdatedAt": now_ms,
                "createdAt": now_ms - 3600000,  # 1 hour before lastUpdatedAt
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
                "hasBlockingPendingActions": False
            }
            
            composer_data['allComposers'].append(new_entry)
            added_count += 1
        
        if (i + 1) % 20 == 0:
            print(f"   Processed {i+1}/{len(conv_ids_with_bubbles)} conversations...")
    
    if conn_backup:
        conn_backup.close()
    
    print(f"\n✅ Added {added_count} new entries")
    print(f"✅ Updated {updated_count} existing entries")
    print(f"📊 Total entries in composerData: {len(composer_data['allComposers'])}\n")
    
    # Update database
    print("💾 Updating composerData in database...")
    cursor.execute(
        "INSERT OR REPLACE INTO ItemTable (key, value) VALUES (?, ?);",
        ('composer.composerData', json.dumps(composer_data))
    )
    conn.commit()
    
    print("✅ composerData rebuilt successfully!\n")
    print("🔄 Next steps:")
    print("   1. Close Cursor completely: pkill -f 'Cursor'")
    print("   2. Restart Cursor")
    print("   3. All conversations should now be visible!\n")
    
    conn.close()
    return True

if __name__ == "__main__":
    rebuild_composerdata()
