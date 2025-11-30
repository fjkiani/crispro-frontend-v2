#!/usr/bin/env python3
"""
Restore conversations from workspace-specific database to global database.
These are conversations that exist only in workspaceStorage, not in the global backup.
"""

import json
import sqlite3
import sys
from pathlib import Path
from datetime import datetime

def restore_workspace_conversations():
    """
    Restore conversations from workspace database to global database.
    """
    workspace_db = "/Users/fahadkiani/Library/Application Support/Cursor/User/workspaceStorage/8eae1effa861c9269e5cc4072b2a492b/state.vscdb"
    global_db = "/Users/fahadkiani/Library/Application Support/Cursor/User/globalStorage/state.vscdb"
    
    print(f"🔄 Restoring workspace-specific conversations...\n")
    print(f"📁 Workspace DB: {workspace_db}")
    print(f"📁 Global DB: {global_db}\n")
    
    if not Path(workspace_db).exists():
        print(f"❌ Workspace DB not found: {workspace_db}")
        return False
    if not Path(global_db).exists():
        print(f"❌ Global DB not found: {global_db}")
        return False

    # Safety check
    response = input("⚠️  Have you closed Cursor completely? (yes/no): ").strip().lower()
    if response != 'yes':
        print("❌ Please close Cursor first using: pkill -f 'Cursor'")
        return False

    try:
        # Get workspace conversations
        conn_workspace = sqlite3.connect(workspace_db)
        cursor_workspace = conn_workspace.cursor()
        
        cursor_workspace.execute("SELECT value FROM ItemTable WHERE key = 'composer.composerData';")
        workspace_row = cursor_workspace.fetchone()
        
        if not workspace_row:
            print("⚠️  No composerData in workspace database")
            conn_workspace.close()
            return False
        
        workspace_composer_data = json.loads(workspace_row[0])
        workspace_composers = workspace_composer_data.get('allComposers', [])
        workspace_conv_ids = {c.get('id') or c.get('composerId') for c in workspace_composers if c.get('id') or c.get('composerId')}
        
        print(f"📊 Found {len(workspace_conv_ids)} conversations in workspace database\n")
        
        # Get global conversations
        conn_global = sqlite3.connect(global_db)
        cursor_global = conn_global.cursor()
        
        cursor_global.execute("SELECT value FROM ItemTable WHERE key = 'composer.composerData';")
        global_row = cursor_global.fetchone()
        
        global_composer_data = {}
        if global_row:
            global_composer_data = json.loads(global_row[0])
        
        if 'allComposers' not in global_composer_data:
            global_composer_data['allComposers'] = []
        
        global_conv_ids = {c.get('id') or c.get('composerId') for c in global_composer_data['allComposers'] if c.get('id') or c.get('composerId')}
        
        # Find unique conversations
        unique_conv_ids = workspace_conv_ids - global_conv_ids
        print(f"🔍 Found {len(unique_conv_ids)} unique conversations to restore\n")
        
        if not unique_conv_ids:
            print("✅ All workspace conversations already in global database!")
            conn_workspace.close()
            conn_global.close()
            return True
        
        # Restore each unique conversation
        total_bubbles = 0
        restored_count = 0
        
        print("🚀 Starting restoration...\n")
        print("=" * 80)
        
        for i, conv_id in enumerate(unique_conv_ids, 1):
            try:
                # Get bubbles from workspace
                cursor_workspace.execute(f"SELECT key, value FROM cursorDiskKV WHERE key LIKE 'bubbleId:{conv_id}:%';")
                bubbles = cursor_workspace.fetchall()
                
                if not bubbles:
                    print(f"[{i}/{len(unique_conv_ids)}] ⏭️  {conv_id[:8]}... (no bubbles)")
                    continue
                
                # Copy bubbles to global
                for key, value in bubbles:
                    cursor_global.execute(
                        "INSERT OR REPLACE INTO cursorDiskKV (key, value) VALUES (?, ?);",
                        (key, value)
                    )
                
                # Get composer metadata
                workspace_composer = next((c for c in workspace_composers if (c.get('id') or c.get('composerId')) == conv_id), None)
                
                if workspace_composer:
                    # Add to global composerData
                    global_composer_data['allComposers'].append(workspace_composer)
                
                total_bubbles += len(bubbles)
                restored_count += 1
                print(f"[{i}/{len(unique_conv_ids)}] ✅ {conv_id[:8]}... ({len(bubbles)} bubbles)")
                
            except Exception as e:
                print(f"[{i}/{len(unique_conv_ids)}] ❌ {conv_id[:8]}... - Error: {e}")
        
        # Update global composerData
        cursor_global.execute(
            "INSERT OR REPLACE INTO ItemTable (key, value) VALUES (?, ?);",
            ('composer.composerData', json.dumps(global_composer_data))
        )
        
        conn_global.commit()
        conn_workspace.close()
        conn_global.close()
        
        print("=" * 80)
        print(f"\n🎉 Restoration complete!\n")
        print(f"📊 Statistics:")
        print(f"   ✅ Restored: {restored_count} conversations")
        print(f"   💬 Total bubbles: {total_bubbles:,}")
        print(f"\n🔄 Restart Cursor to see restored conversations.\n")
        
        return True
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = restore_workspace_conversations()
    sys.exit(0 if success else 1)









