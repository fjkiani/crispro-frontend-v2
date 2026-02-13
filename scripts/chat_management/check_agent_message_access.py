#!/usr/bin/env python3
"""
Check if agent can access all messages in database, not just what UI shows.
This is critical for agent context recall.
"""
import sqlite3
import json

db = '/Users/fahadkiani/Library/Application Support/Cursor/User/globalStorage/state.vscdb'
conn = sqlite3.connect(db)
cursor = conn.cursor()

# Get the largest conversation
cursor.execute('''
    SELECT substr(key, 9, 36) as conv_id, COUNT(*) as msg_count
    FROM cursorDiskKV
    WHERE key LIKE 'bubbleId:%'
    AND length(key) > 45
    AND value IS NOT NULL
    GROUP BY conv_id
    ORDER BY msg_count DESC
    LIMIT 1;
''')

largest = cursor.fetchone()
if largest:
    conv_id = largest[0]
    msg_count = largest[1]
    print(f"📊 Testing conversation: {conv_id}")
    print(f"📊 Total messages in DB: {msg_count:,}\n")
    
    # Get first message (oldest)
    pattern = f'bubbleId:{conv_id}:%'
    query = "SELECT key, value FROM cursorDiskKV WHERE key LIKE ? AND value IS NOT NULL ORDER BY key ASC LIMIT 1;"
    cursor.execute(query, (pattern,))
    
    first = cursor.fetchone()
    if first:
        first_key = first[0]
        first_data = json.loads(first[1])
        print(f"✅ First message (oldest) found:")
        print(f"   Key: {first_key}")
        has_text = 'text' in first_data or 'content' in first_data
        print(f"   Has text: {has_text}")
        if 'text' in first_data:
            text_preview = str(first_data['text'])[:100]
            print(f"   Preview: {text_preview}...")
    
    # Get last message (newest)
    query = "SELECT key, value FROM cursorDiskKV WHERE key LIKE ? AND value IS NOT NULL ORDER BY key DESC LIMIT 1;"
    cursor.execute(query, (pattern,))
    
    last = cursor.fetchone()
    if last:
        last_key = last[0]
        last_data = json.loads(last[1])
        print(f"\n✅ Last message (newest) found:")
        print(f"   Key: {last_key}")
        has_text = 'text' in last_data or 'content' in last_data
        print(f"   Has text: {has_text}")
    
    # Check middle messages (to verify all are accessible)
    cursor.execute('''
        SELECT COUNT(*) 
        FROM cursorDiskKV
        WHERE key LIKE ?
        AND value IS NOT NULL;
    ''', (pattern,))
    
    total_in_db = cursor.fetchone()[0]
    print(f"\n📊 Messages accessible via database query: {total_in_db:,}")
    
    print(f"\n🔍 CRITICAL QUESTION:")
    print(f"   Database has: {total_in_db:,} messages")
    print(f"   UI shows: ~50-100 messages (estimated)")
    print(f"   Agent access: UNKNOWN - depends on how Cursor loads context")
    print(f"\n⚠️  If agent only sees UI-loaded messages, context is limited!")
    print(f"   If agent reads from database, it has full {total_in_db:,} message context")

conn.close()

