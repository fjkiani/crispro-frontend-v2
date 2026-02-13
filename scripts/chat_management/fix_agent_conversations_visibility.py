#!/usr/bin/env python3
"""
Fix agent conversation visibility - ensure all agent conversations are properly indexed.
Agent conversations have unifiedMode: "agent" and may need special handling.
"""

import json
import sqlite3
import sys
from pathlib import Path
from datetime import datetime

def fix_agent_visibility():
    """
    Fix agent conversation visibility by ensuring all agent conversations
    with bubbles are properly indexed in composerData.
    """
    global_db = "/Users/fahadkiani/Library/Application Support/Cursor/User/globalStorage/state.vscdb"
    
    print(f"🔧 Fixing agent conversation visibility...\n")
    print(f"📁 Global DB: {global_db}\n")
    
    if not Path(global_db).exists():
        print(f"❌ Global DB not found: {global_db}")
        return False

    # Safety check
    response = input("⚠️  Have you closed Cursor completely? (yes/no): ")
    if response.lower() != 'yes':
        print("Please close Cursor completely before proceeding.")
        return False

    # Create a backup of the current global DB
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = Path(f"{global_db}.before_fix_agents_{timestamp}")
    try:
        Path(global_db).rename(backup_path)
        print(f"✅ Created backup of current global DB: {backup_path}")
        # Copy the backup back to the original location to work on it
        import shutil
        shutil.copy(backup_path, global_db)
    except Exception as e:
        print(f"❌ Error creating backup: {e}")
        return False

    conn = None
    try:
        conn = sqlite3.connect(str(global_db))
        cursor = conn.cursor()
        
        # 1. Get all unique conversation IDs from bubbles
        cursor.execute("SELECT DISTINCT key FROM cursorDiskKV WHERE key LIKE 'bubbleId:%';")
        all_bubble_keys = cursor.fetchall()
        
        conv_ids_from_bubbles = set()
        for (key,) in all_bubble_keys:
            parts = key.split(':')
            if len(parts) >= 2:
                conv_id = parts[1]
                conv_ids_from_bubbles.add(conv_id)
        
        print(f"📊 Found {len(conv_ids_from_bubbles)} unique conversation IDs from bubbles.")

        # 2. Get existing composerData
        cursor.execute("SELECT value FROM ItemTable WHERE key = 'composer.composerData';")
        composer_row = cursor.fetchone()
        
        composer_data = {"allComposers": []}
        if composer_row:
            composer_data = json.loads(composer_row[0])
        
        existing_composer_map = {c.get('id') or c.get('composerId'): c for c in composer_data.get('allComposers', [])}
        
        # 3. Check which conversations are agents vs regular
        agent_convs_with_bubbles = []
        regular_convs_with_bubbles = []
        
        for conv_id in conv_ids_from_bubbles:
            # Check if it's already in composerData
            if conv_id in existing_composer_map:
                composer_entry = existing_composer_map[conv_id]
                unified_mode = composer_entry.get('unifiedMode', '')
                if unified_mode == 'agent':
                    agent_convs_with_bubbles.append(conv_id)
                else:
                    regular_convs_with_bubbles.append(conv_id)
            else:
                # Not in composerData - need to determine if it's an agent
                # Check bubble content for agent patterns
                cursor.execute("SELECT value FROM cursorDiskKV WHERE key LIKE ? LIMIT 1;", (f'bubbleId:{conv_id}%',))
                bubble_row = cursor.fetchone()
                if bubble_row:
                    try:
                        bubble_data = json.loads(bubble_row[0])
                        # Heuristic: if bubble has agent-related fields, it's an agent
                        if isinstance(bubble_data, dict):
                            # Check for agent indicators
                            is_agent = any([
                                'agent' in str(bubble_data).lower(),
                                bubble_data.get('agentType'),
                                bubble_data.get('unifiedMode') == 'agent'
                            ])
                            if is_agent:
                                agent_convs_with_bubbles.append(conv_id)
                            else:
                                regular_convs_with_bubbles.append(conv_id)
                    except:
                        # Default to regular if we can't determine
                        regular_convs_with_bubbles.append(conv_id)
        
        print(f"\n📊 Conversation Type Breakdown:")
        print(f"   🤖 Agent conversations with bubbles: {len(agent_convs_with_bubbles)}")
        print(f"   💬 Regular conversations with bubbles: {len(regular_convs_with_bubbles)}")
        
        # 4. Build new composerData with proper agent/regular distinction
        new_composers = []
        kept_agents = 0
        kept_regular = 0
        added_agents = 0
        added_regular = 0
        
        # Keep existing valid entries
        for composer_entry in composer_data.get('allComposers', []):
            comp_id = composer_entry.get('id') or composer_entry.get('composerId')
            if comp_id in conv_ids_from_bubbles:
                new_composers.append(composer_entry)
                if composer_entry.get('unifiedMode') == 'agent':
                    kept_agents += 1
                else:
                    kept_regular += 1
        
        # Add missing agent conversations
        for conv_id in agent_convs_with_bubbles:
            if conv_id not in existing_composer_map:
                # Create agent composerData entry
                new_composer_entry = {
                    "id": conv_id,
                    "composerId": conv_id,
                    "unifiedMode": "agent",  # CRITICAL: Mark as agent
                    "type": "head",
                    "title": f"Agent Conversation {conv_id[:8]}...",
                    "lastUpdated": datetime.now().isoformat(),
                    "created": datetime.now().isoformat(),
                    "pinned": False,
                    "archived": False,
                    "isDeleted": False,
                    "isHidden": False,
                    "model": "unknown",
                    "context": []
                }
                new_composers.append(new_composer_entry)
                added_agents += 1
                print(f"   ➕ Added agent conversation: {conv_id[:8]}...")
        
        # Add missing regular conversations
        for conv_id in regular_convs_with_bubbles:
            if conv_id not in existing_composer_map:
                # Create regular composerData entry
                new_composer_entry = {
                    "id": conv_id,
                    "composerId": conv_id,
                    "unifiedMode": "chat",  # Regular chat
                    "type": "head",
                    "title": f"Chat {conv_id[:8]}...",
                    "lastUpdated": datetime.now().isoformat(),
                    "created": datetime.now().isoformat(),
                    "pinned": False,
                    "archived": False,
                    "isDeleted": False,
                    "isHidden": False,
                    "model": "unknown",
                    "context": []
                }
                new_composers.append(new_composer_entry)
                added_regular += 1
        
        composer_data['allComposers'] = new_composers
        
        # 5. Update ItemTable with the new composerData
        cursor.execute(
            "INSERT OR REPLACE INTO ItemTable (key, value) VALUES (?, ?);",
            ('composer.composerData', json.dumps(composer_data))
        )
        conn.commit()
        
        print("\n🎉 Agent conversation fix complete!")
        print(f"\n📊 Statistics:")
        print(f"   🤖 Agent conversations:")
        print(f"      ✅ Kept existing: {kept_agents}")
        print(f"      ➕ Added missing: {added_agents}")
        print(f"   💬 Regular conversations:")
        print(f"      ✅ Kept existing: {kept_regular}")
        print(f"      ➕ Added missing: {added_regular}")
        print(f"   ➡️ Total conversations in composerData: {len(new_composers)}")
        print(f"\n💡 IMPORTANT:")
        print(f"   Agent conversations may appear in a separate 'Agents' panel in Cursor.")
        print(f"   Look for an 'Agents' tab or panel in the sidebar.")
        print(f"   Regular chats appear in the main chat history.")
        print("\n🔄 Restart Cursor to see all conversations.\n")
        
        return True

    except sqlite3.Error as e:
        print(f"❌ SQLite error: {e}")
        return False
    except json.JSONDecodeError as e:
        print(f"❌ JSON decoding error: {e}")
        return False
    except Exception as e:
        print(f"❌ An unexpected error occurred: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if conn:
            conn.close()

if __name__ == "__main__":
    fix_agent_visibility()



