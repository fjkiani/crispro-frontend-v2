#!/usr/bin/env python3
"""
⚔️ IMPORT SPECSTORY MARKDOWN FILES TO CURSOR DATABASE ⚔️

Converts SpecStory markdown chat history files into Cursor's SQLite database format.
"""

import sqlite3
import json
import re
import sys
from pathlib import Path
from datetime import datetime
import uuid

# Database paths
GLOBAL_DB = Path.home() / "Library/Application Support/Cursor/User/globalStorage/state.vscdb"
BACKUP_DIR = Path.home() / "cursor_db_backups"

def extract_session_id(markdown_content):
    """Extract Cursor session ID from SpecStory markdown header."""
    match = re.search(r'<!-- cursor Session ([a-f0-9-]+)', markdown_content)
    if match:
        return match.group(1)
    return None

def parse_timestamp(timestamp_str):
    """Convert SpecStory timestamp (2025-11-19 22:05Z) to milliseconds since epoch."""
    try:
        # Parse format: "2025-11-19 22:05Z"
        dt = datetime.strptime(timestamp_str.replace('Z', ''), '%Y-%m-%d %H:%M')
        # Convert to UTC (Z means UTC)
        return int(dt.timestamp() * 1000)
    except:
        return None

def parse_specstory_markdown(file_path):
    """Parse SpecStory markdown file into conversation structure."""
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    session_id = extract_session_id(content)
    if not session_id:
        print(f"⚠️  Could not extract session ID from {file_path.name}")
        return None, []
    
    # Split content by message markers
    # Format: _**User (timestamp)**_ or _**Agent (model, mode)**_
    messages = []
    last_timestamp_ms = None
    
    # Split by message markers (including the marker line itself)
    parts = re.split(r'(?=_\*\*User \(|\*\*Agent \()', content)
    
    for part in parts:
        part = part.strip()
        if not part:
            continue
        
        # User message
        user_match = re.match(r'_\*\*User \(([^)]+)\)\*\*_\s*\n(.*?)(?=---|_\*\*|$)', part, re.DOTALL)
        if user_match:
            timestamp_str = user_match.group(1)
            text = user_match.group(2).strip()
            # Remove tool-use blocks and other markdown artifacts
            text = re.sub(r'<tool-use[^>]*>.*?</tool-use>', '', text, flags=re.DOTALL)
            text = re.sub(r'<details>.*?</details>', '', text, flags=re.DOTALL)
            text = re.sub(r'---+\s*', '', text)
            text = text.strip()
            if text:
                timestamp_ms = parse_timestamp(timestamp_str)
                if timestamp_ms:
                    last_timestamp_ms = timestamp_ms
                messages.append({
                    'type': 'user',
                    'timestamp_ms': timestamp_ms or (last_timestamp_ms + 1000 if last_timestamp_ms else int(datetime.now().timestamp() * 1000)),
                    'text': text
                })
            continue
        
        # Agent message
        agent_match = re.match(r'_\*\*Agent \(([^)]+)\)\*\*_\s*\n(.*?)(?=---|_\*\*|$)', part, re.DOTALL)
        if agent_match:
            agent_info = agent_match.group(1)
            text = agent_match.group(2).strip()
            # Remove tool-use blocks and other markdown artifacts
            text = re.sub(r'<tool-use[^>]*>.*?</tool-use>', '', text, flags=re.DOTALL)
            text = re.sub(r'<details>.*?</details>', '', text, flags=re.DOTALL)
            text = re.sub(r'---+\s*', '', text)
            text = text.strip()
            if text:
                parts_info = [p.strip() for p in agent_info.split(',')]
                # Agent messages don't have timestamps - infer from previous message + 1 second
                agent_timestamp_ms = (last_timestamp_ms + 1000) if last_timestamp_ms else int(datetime.now().timestamp() * 1000)
                last_timestamp_ms = agent_timestamp_ms
                messages.append({
                    'type': 'assistant',
                    'model': parts_info[0] if len(parts_info) > 0 else 'default',
                    'mode': parts_info[1] if len(parts_info) > 1 else 'default',
                    'timestamp_ms': agent_timestamp_ms,
                    'text': text
                })
            continue
    
    return session_id, messages

def create_bubble_entry(conversation_id, bubble_id, message, index):
    """Create a Cursor bubbleId entry from a message."""
    timestamp_ms = message.get('timestamp_ms', int(datetime.now().timestamp() * 1000) + index)
    
    if message['type'] == 'user':
        bubble_data = {
            "_v": 3,
            "type": 2,
            "text": message['text'],
            "bubbleId": bubble_id,
            "timestamp": timestamp_ms,
            "approximateLintErrors": [],
            "lints": [],
            "codebaseContextChunks": [],
            "commits": [],
            "pullRequests": [],
            "attachedCodeChunks": [],
            "assistantSuggestedDiffs": [],
            "gitDiffs": [],
            "interpreterResults": [],
            "images": [],
            "attachedFolders": [],
            "attachedFoldersNew": [],
            "userResponsesToSuggestedCodeBlocks": [],
            "suggestedCodeBlocks": [],
            "diffsForCompressingFiles": [],
            "relevantFiles": [],
            "toolResults": [],
            "notepads": [],
            "capabilities": [],
            "capabilityStatuses": {
                "mutate-request": [],
                "start-submit-chat": [],
                "before-submit-chat": [],
                "chat-stream-finished": [],
                "before-apply": [],
                "after-apply": []
            }
        }
    else:  # assistant
        bubble_data = {
            "_v": 3,
            "type": 2,
            "text": message['text'],
            "bubbleId": bubble_id,
            "timestamp": timestamp_ms,
            "approximateLintErrors": [],
            "lints": [],
            "codebaseContextChunks": [],
            "commits": [],
            "pullRequests": [],
            "attachedCodeChunks": [],
            "assistantSuggestedDiffs": [],
            "gitDiffs": [],
            "interpreterResults": [],
            "images": [],
            "attachedFolders": [],
            "attachedFoldersNew": [],
            "userResponsesToSuggestedCodeBlocks": [],
            "suggestedCodeBlocks": [],
            "diffsForCompressingFiles": [],
            "relevantFiles": [],
            "toolResults": [],
            "notepads": [],
            "capabilities": [],
            "capabilityStatuses": {
                "mutate-request": [],
                "start-submit-chat": [],
                "before-submit-chat": [],
                "chat-stream-finished": [],
                "before-apply": [],
                "after-apply": []
            }
        }
    
    return {
        'key': f'bubbleId:{conversation_id}:{bubble_id}',
        'value': json.dumps(bubble_data)
    }

def import_specstory_file(file_path, db_path):
    """Import a single SpecStory markdown file into Cursor database."""
    print(f"\n📄 Processing: {file_path.name}")
    
    session_id, messages = parse_specstory_markdown(file_path)
    if not session_id:
        print(f"   ❌ Failed to extract session ID")
        return False
    
    if not messages:
        print(f"   ⚠️  No messages found in file")
        return False
    
    print(f"   ✅ Session ID: {session_id}")
    print(f"   ✅ Found {len(messages)} messages")
    
    # Connect to database
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    
    try:
        # Check if conversation already exists
        existing = cursor.execute(
            "SELECT COUNT(*) FROM cursorDiskKV WHERE key LIKE ?",
            (f'bubbleId:{session_id}:%',)
        ).fetchone()[0]
        
        if existing > 0:
            print(f"   ⚠️  Conversation already exists with {existing} bubbles")
            print(f"   🔄 Overwriting with SpecStory version ({len(messages)} messages)...")
            # Delete existing conversation from THIS database
            cursor.execute(
                "DELETE FROM cursorDiskKV WHERE key LIKE ?",
                (f'bubbleId:{session_id}:%',)
            )
            cursor.execute(
                "DELETE FROM cursorDiskKV WHERE key = ?",
                (f'composerData:{session_id}',)
            )
        
        # Also delete from workspace databases
        workspace_storage = Path.home() / "Library/Application Support/Cursor/User/workspaceStorage"
        if workspace_storage.exists():
            for ws_dir in workspace_storage.iterdir():
                ws_db = ws_dir / "state.vscdb"
                if ws_db.exists():
                    try:
                        ws_conn = sqlite3.connect(str(ws_db))
                        ws_cursor = ws_conn.cursor()
                        ws_cursor.execute(
                            "DELETE FROM cursorDiskKV WHERE key LIKE ?",
                            (f'bubbleId:{session_id}:%',)
                        )
                        ws_cursor.execute(
                            "DELETE FROM cursorDiskKV WHERE key = ?",
                            (f'composerData:{session_id}',)
                        )
                        ws_conn.commit()
                        ws_conn.close()
                    except:
                        pass
        
        # Generate sequential bubble IDs (zero-padded for alphabetical ordering)
        bubble_entries = []
        max_digits = len(str(len(messages)))
        for i, message in enumerate(messages):
            # Use zero-padded sequential ID so they sort correctly
            bubble_id = f"{i:0{max_digits}d}"
            entry = create_bubble_entry(session_id, bubble_id, message, i)
            bubble_entries.append(entry)
        
        # Insert all bubbles
        print(f"   🔄 Inserting {len(bubble_entries)} bubbles...")
        for entry in bubble_entries:
            cursor.execute(
                "INSERT OR REPLACE INTO cursorDiskKV (key, value) VALUES (?, ?)",
                (entry['key'], entry['value'])
            )
        
        # Create composerData entry with title and preview
        # Extract title from first user message (first 50 chars)
        first_user_message = next((msg for msg in messages if msg.get('type') == 'user'), None)
        title = "Imported Conversation"
        preview = ""
        if first_user_message and first_user_message.get('text'):
            title_text = first_user_message['text'].strip()[:50]
            if title_text:
                title = title_text.replace('\n', ' ').replace('\r', ' ')
            preview = first_user_message['text'].strip()[:200].replace('\n', ' ').replace('\r', ' ')
        
        composer_data = {
            "_v": 3,
            "type": 2,
            "title": title,
            "preview": preview,
            "lastAccessed": int(datetime.now().timestamp() * 1000)
        }
        cursor.execute(
            "INSERT OR REPLACE INTO cursorDiskKV (key, value) VALUES (?, ?)",
            (f'composerData:{session_id}', json.dumps(composer_data))
        )
        
        conn.commit()
        print(f"   ✅ Successfully imported {len(bubble_entries)} messages to global database")
        
        # Also import to workspace databases
        import_to_workspace_databases(session_id, bubble_entries, composer_data)
        
        return True
        
    except Exception as e:
        conn.rollback()
        print(f"   ❌ Error: {e}")
        return False
    finally:
        conn.close()

def import_to_workspace_databases(session_id, bubble_entries, composer_data):
    """Import conversation to all workspace databases."""
    workspace_storage = Path.home() / "Library/Application Support/Cursor/User/workspaceStorage"
    if not workspace_storage.exists():
        return
    
    imported_count = 0
    for ws_dir in workspace_storage.iterdir():
        ws_db = ws_dir / "state.vscdb"
        if not ws_db.exists():
            continue
        
        try:
            ws_conn = sqlite3.connect(str(ws_db))
            ws_cursor = ws_conn.cursor()
            
            # Insert bubbles
            for entry in bubble_entries:
                ws_cursor.execute(
                    "INSERT OR REPLACE INTO cursorDiskKV (key, value) VALUES (?, ?)",
                    (entry['key'], entry['value'])
                )
            
            # Insert composerData
            ws_cursor.execute(
                "INSERT OR REPLACE INTO cursorDiskKV (key, value) VALUES (?, ?)",
                (f'composerData:{session_id}', json.dumps(composer_data))
            )
            
            ws_conn.commit()
            ws_conn.close()
            imported_count += 1
        except Exception as e:
            print(f"   ⚠️  Failed to import to workspace {ws_dir.name}: {e}")
    
    if imported_count > 0:
        print(f"   ✅ Also imported to {imported_count} workspace database(s)")

def main():
    """Main import function."""
    print("⚔️ IMPORT SPECSTORY FILES TO CURSOR DATABASE ⚔️")
    print("=" * 60)
    
    # Check if Cursor is running
    import subprocess
    try:
        result = subprocess.run(['pgrep', '-x', 'Cursor'], capture_output=True)
        if result.returncode == 0:
            print("❌ Cursor is running! Please quit completely (Cmd+Q) first.")
            sys.exit(1)
    except:
        pass
    
    # Verify database exists
    if not GLOBAL_DB.exists():
        print(f"❌ Database not found: {GLOBAL_DB}")
        sys.exit(1)
    
    # Backup database
    BACKUP_DIR.mkdir(exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = BACKUP_DIR / f"global_before_specstory_import_{timestamp}"
    import shutil
    shutil.copy2(GLOBAL_DB, backup_path)
    print(f"✅ Backed up database to: {backup_path.name}")
    
    # Files to import
    specstory_dir = Path(".specstory/history")
    files_to_import = [
        specstory_dir / "2025-11-19_17-41Z-designing-a-zero-tolerance-content-filter.md"
    ]
    
    print(f"\n📋 Files to import:")
    for f in files_to_import:
        if f.exists():
            print(f"   ✅ {f.name}")
        else:
            print(f"   ❌ {f.name} (not found)")
    
    print("\n🔄 Starting import...")
    
    success_count = 0
    for file_path in files_to_import:
        if file_path.exists():
            if import_specstory_file(file_path, GLOBAL_DB):
                success_count += 1
        else:
            print(f"\n❌ File not found: {file_path}")
    
    print("\n" + "=" * 60)
    print(f"⚔️ IMPORT COMPLETE! ⚔️")
    print(f"   ✅ Successfully imported: {success_count}/{len(files_to_import)} files")
    print(f"   💾 Backup saved: {backup_path.name}")
    print("\n🎯 NEXT STEPS:")
    print("   1. Open Cursor")
    print("   2. Wait 3-5 minutes for re-indexing")
    print("   3. Try: Cmd+Shift+P → 'Reload Window'")
    print("   4. Check chat history panel")

if __name__ == "__main__":
    main()

