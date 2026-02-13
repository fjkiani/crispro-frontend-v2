#!/usr/bin/env python3
"""
Export ALL conversations (agents + regular) to text files.
This helps recover conversations that aren't visible in Cursor's UI.
"""

import sqlite3
import json
from pathlib import Path
from datetime import datetime
import re

def extract_text_from_bubble(bubble_data):
    """Extract readable text from a bubble, trying multiple fields."""
    if isinstance(bubble_data, str):
        try:
            bubble_data = json.loads(bubble_data)
        except:
            return bubble_data[:500] if len(bubble_data) > 500 else bubble_data
    
    if not isinstance(bubble_data, dict):
        return str(bubble_data)[:500]
    
    # Try direct text field
    if 'text' in bubble_data and isinstance(bubble_data['text'], str):
        return bubble_data['text']
    
    # Try content field (could be string or array)
    if 'content' in bubble_data:
        content = bubble_data['content']
        if isinstance(content, str):
            return content
        elif isinstance(content, list):
            text_parts = []
            for item in content:
                if isinstance(item, dict) and 'text' in item:
                    text_parts.append(item['text'])
                elif isinstance(item, str):
                    text_parts.append(item)
            if text_parts:
                return '\n'.join(text_parts)
    
    # Try message field
    if 'message' in bubble_data:
        msg = bubble_data['message']
        if isinstance(msg, str):
            return msg
        elif isinstance(msg, dict) and 'text' in msg:
            return msg['text']
    
    # Try to extract from codebaseContextChunks or other fields
    text_parts = []
    for field in ['codebaseContextChunks', 'attachedCodeChunks', 'toolResults']:
        if field in bubble_data and isinstance(bubble_data[field], list):
            for item in bubble_data[field][:5]:  # Limit to first 5
                if isinstance(item, dict):
                    for subfield in ['text', 'content', 'code', 'path']:
                        if subfield in item:
                            text_parts.append(f"[{field}.{subfield}]: {str(item[subfield])[:200]}")
    
    if text_parts:
        return '\n'.join(text_parts)
    
    return "[No readable text found in bubble]"

def export_all_conversations():
    """Export all conversations to text files."""
    current_db = '/Users/fahadkiani/Library/Application Support/Cursor/User/globalStorage/state.vscdb'
    output_dir = Path('/Users/fahadkiani/Desktop/development/crispr-assistant-main/exported_conversations')
    output_dir.mkdir(exist_ok=True)
    
    print('📤 EXPORTING ALL CONVERSATIONS TO TEXT FILES...\n')
    
    conn = sqlite3.connect(current_db)
    cursor = conn.cursor()
    
    # Get all composerData
    cursor.execute("SELECT value FROM ItemTable WHERE key = 'composer.composerData';")
    row = cursor.fetchone()
    composer_data = json.loads(row[0])
    all_composers = composer_data.get('allComposers', [])
    
    print(f'📊 Found {len(all_composers)} total conversations\n')
    
    exported_count = 0
    skipped_count = 0
    
    for i, conv in enumerate(all_composers, 1):
        conv_id = conv.get('composerId') or conv.get('id', 'unknown')
        clean_id = str(conv_id).lstrip(':')
        
        subtitle = conv.get('subtitle', conv.get('title', 'No title'))
        unified_mode = conv.get('unifiedMode', '')
        mode_str = 'agent' if unified_mode == 'agent' else 'regular'
        
        # Get all bubbles for this conversation
        cursor.execute("SELECT key, value FROM cursorDiskKV WHERE key LIKE ? ORDER BY key;", (f'bubbleId:{clean_id}%',))
        bubbles = cursor.fetchall()
        
        if not bubbles:
            skipped_count += 1
            continue
        
        # Sanitize filename
        clean_subtitle = re.sub(r'[^\w\s-]', '', subtitle)[:50] if subtitle != 'No title' else 'untitled'
        if not clean_subtitle:
            clean_subtitle = 'untitled'
        
        # Create output file
        filename = f"{i:03d}_{mode_str}_{clean_subtitle}_{clean_id[:8]}.txt"
        filepath = output_dir / filename
        
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(f"CONVERSATION EXPORT\n")
            f.write(f"{'='*80}\n\n")
            f.write(f"Conversation ID: {conv_id}\n")
            f.write(f"Type: {mode_str.upper()}\n")
            f.write(f"Title: {subtitle}\n")
            f.write(f"Total Messages: {len(bubbles)}\n")
            f.write(f"Export Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"\n{'='*80}\n\n")
            
            # Export each message
            for j, (key, value) in enumerate(bubbles, 1):
                f.write(f"\n{'─'*80}\n")
                f.write(f"MESSAGE {j}/{len(bubbles)}\n")
                f.write(f"Key: {key}\n")
                f.write(f"{'─'*80}\n\n")
                
                try:
                    bubble_data = json.loads(value)
                    text = extract_text_from_bubble(bubble_data)
                    f.write(text)
                    f.write('\n')
                except Exception as e:
                    f.write(f"[Error extracting message: {e}]\n")
                    f.write(f"Raw data preview: {str(value)[:500]}\n")
        
        exported_count += 1
        if i % 10 == 0:
            print(f'💾 Exported {i}/{len(all_composers)} conversations...')
    
    conn.close()
    
    print(f'\n✅ EXPORT COMPLETE!\n')
    print(f'📁 Files saved to: {output_dir}\n')
    print(f'📊 Statistics:')
    print(f'   ✅ Exported: {exported_count} conversations')
    print(f'   ⏭️  Skipped (no messages): {skipped_count} conversations\n')

if __name__ == "__main__":
    export_all_conversations()



