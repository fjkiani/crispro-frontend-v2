#!/usr/bin/env python3
"""
Export agent conversations that mention "Ayesha" to text files.
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

def find_ayesha_agent_conversations():
    """Find all agent conversations that mention Ayesha."""
    current_db = '/Users/fahadkiani/Library/Application Support/Cursor/User/globalStorage/state.vscdb'
    output_dir = Path('/Users/fahadkiani/Desktop/development/crispr-assistant-main/exported_conversations')
    output_dir.mkdir(exist_ok=True)
    
    print('🔍 SEARCHING FOR AGENT CONVERSATIONS MENTIONING "AYESHA"...\n')
    
    conn = sqlite3.connect(current_db)
    cursor = conn.cursor()
    
    # Get all composerData
    cursor.execute("SELECT value FROM ItemTable WHERE key = 'composer.composerData';")
    row = cursor.fetchone()
    composer_data = json.loads(row[0])
    all_composers = composer_data.get('allComposers', [])
    
    # Find agent conversations
    agent_conversations = []
    for conv in all_composers:
        if conv.get('unifiedMode') == 'agent':
            conv_id = conv.get('composerId') or conv.get('id', 'unknown')
            agent_conversations.append((conv_id, conv))
    
    print(f'📊 Found {len(agent_conversations)} agent conversations\n')
    
    # Search each agent conversation for "Ayesha"
    ayesha_conversations = []
    
    for conv_id, conv in agent_conversations:
        clean_id = str(conv_id).lstrip(':')
        
        # Get all bubbles for this conversation
        cursor.execute("SELECT key, value FROM cursorDiskKV WHERE key LIKE ? ORDER BY key;", (f'bubbleId:{clean_id}%',))
        bubbles = cursor.fetchall()
        
        # Search bubbles for "Ayesha" (case-insensitive)
        found_ayesha = False
        ayesha_count = 0
        
        for key, value in bubbles:
            try:
                bubble_data = json.loads(value)
                text = extract_text_from_bubble(bubble_data)
                
                # Search for Ayesha (case-insensitive)
                if re.search(r'ayesha', text, re.IGNORECASE):
                    found_ayesha = True
                    ayesha_count += text.lower().count('ayesha')
            except:
                pass
        
        if found_ayesha:
            subtitle = conv.get('subtitle', conv.get('title', 'No title'))
            ayesha_conversations.append({
                'id': clean_id,
                'subtitle': subtitle,
                'conv': conv,
                'bubble_count': len(bubbles),
                'ayesha_mentions': ayesha_count
            })
            print(f'✅ Found: {subtitle[:60]}... ({ayesha_count} mentions, {len(bubbles)} messages)')
    
    print(f'\n📊 Total agent conversations mentioning Ayesha: {len(ayesha_conversations)}\n')
    
    # Export each conversation to a text file
    for i, conv_info in enumerate(ayesha_conversations, 1):
        conv_id = conv_info['id']
        subtitle = conv_info['subtitle']
        clean_subtitle = re.sub(r'[^\w\s-]', '', subtitle)[:50]  # Sanitize filename
        
        # Get all bubbles for this conversation
        cursor.execute("SELECT key, value FROM cursorDiskKV WHERE key LIKE ? ORDER BY key;", (f'bubbleId:{conv_id}%',))
        bubbles = cursor.fetchall()
        
        # Create output file
        filename = f"{i:02d}_agent_ayesha_{clean_subtitle}_{conv_id[:8]}.txt"
        filepath = output_dir / filename
        
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(f"AGENT CONVERSATION EXPORT\n")
            f.write(f"{'='*80}\n\n")
            f.write(f"Conversation ID: {conv_id}\n")
            f.write(f"Title: {subtitle}\n")
            f.write(f"Total Messages: {len(bubbles)}\n")
            f.write(f"Ayesha Mentions: {conv_info['ayesha_mentions']}\n")
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
        
        print(f'💾 Exported: {filename}')
    
    conn.close()
    
    print(f'\n✅ EXPORT COMPLETE!\n')
    print(f'📁 Files saved to: {output_dir}\n')
    print(f'📊 Exported {len(ayesha_conversations)} agent conversations mentioning Ayesha\n')

if __name__ == "__main__":
    find_ayesha_agent_conversations()



