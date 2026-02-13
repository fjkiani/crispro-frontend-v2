#!/usr/bin/env python3
"""
Parse chat messages from SQL export and create readable text format.
Extracts actual conversation content from composerData and bubbleId entries.
"""

import json
import sqlite3
import sys
from pathlib import Path
from datetime import datetime

def parse_composer_data(value_str):
    """Extract readable content from composerData JSON."""
    try:
        data = json.loads(value_str)
        messages = []
        
        # Look for message arrays or conversation history
        if 'messages' in data:
            for msg in data['messages']:
                if isinstance(msg, dict):
                    role = msg.get('role', 'unknown')
                    content = msg.get('content', '')
                    if content:
                        messages.append(f"{role.upper()}: {content[:200]}")
        
        # Look for conversation threads
        if 'conversation' in data:
            for turn in data['conversation']:
                if isinstance(turn, dict):
                    role = turn.get('role', 'user')
                    text = turn.get('text', turn.get('content', ''))
                    if text:
                        messages.append(f"{role.upper()}: {text[:200]}")
        
        # Look for user/assistant fields
        if 'userMessage' in data:
            messages.append(f"USER: {data['userMessage'][:200]}")
        if 'assistantMessage' in data:
            messages.append(f"ASSISTANT: {data['assistantMessage'][:200]}")
        
        return '\n'.join(messages) if messages else None
    except Exception as e:
        return None

def extract_readable_chat(db_path, output_file=None):
    """Extract readable chat messages from database."""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Get composerData entries (these contain full conversations)
    query = """
    SELECT key, value 
    FROM cursorDiskKV 
    WHERE key LIKE 'composerData:%'
    ORDER BY key
    """
    
    cursor.execute(query)
    rows = cursor.fetchall()
    
    print(f"📊 Found {len(rows)} composer data entries")
    
    readable_chats = []
    for key, value in rows:
        if value:
            value_str = value.decode('utf-8') if isinstance(value, bytes) else str(value)
            parsed = parse_composer_data(value_str)
            if parsed:
                readable_chats.append(f"\n{'='*80}\n")
                readable_chats.append(f"Chat ID: {key}\n")
                readable_chats.append(parsed)
                readable_chats.append(f"\n")
    
    if output_file is None:
        output_file = Path(__file__).parent.parent / 'chat_history_readable.txt'
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(f"Chat History Export\n")
        f.write(f"Generated: {datetime.now().isoformat()}\n")
        f.write(f"Source: {db_path}\n")
        f.write(f"{'='*80}\n\n")
        f.write(''.join(readable_chats))
    
    print(f"✅ Readable chat history written to: {output_file}")
    print(f"📏 File size: {output_file.stat().st_size / 1024:.2f} KB")
    print(f"💬 Extracted {len(readable_chats)} conversation threads")
    
    conn.close()
    return output_file

if __name__ == '__main__':
    backup_db = "/Users/fahadkiani/Library/Application Support/Cursor/User/globalStorage/state.vscdb.backup"
    extract_readable_chat(backup_db)















