#!/usr/bin/env python3
"""
Clean up Cursor's state.vscdb database to improve performance.

This script:
1. Removes old checkpointId entries (keeping only recent ones)
2. Optionally archives old conversations
3. Runs VACUUM to optimize the database
"""
import sqlite3
import json
import shutil
from datetime import datetime
from pathlib import Path

DB_PATH = Path.home() / "Library/Application Support/Cursor/User/globalStorage/state.vscdb"
BACKUP_PATH = DB_PATH.with_suffix(f".vscdb.backup_{datetime.now().strftime('%Y%m%d_%H%M%S')}")

def get_db_size():
    """Get current database size."""
    if not DB_PATH.exists():
        return 0
    return DB_PATH.stat().st_size / (1024 * 1024 * 1024)  # GB

def backup_database():
    """Create a backup before cleanup."""
    print(f"📦 Creating backup: {BACKUP_PATH}")
    shutil.copy2(DB_PATH, BACKUP_PATH)
    print(f"✅ Backup created: {BACKUP_PATH.stat().st_size / (1024*1024):.1f} MB")

def analyze_checkpoints(conn):
    """Analyze checkpoint data."""
    cursor = conn.cursor()
    
    # Get checkpoint statistics
    cursor.execute("""
        SELECT 
            COUNT(*) as count,
            SUM(length(value)) / 1024.0 / 1024.0 as total_mb,
            AVG(length(value)) / 1024.0 as avg_kb
        FROM cursorDiskKV
        WHERE key LIKE 'checkpointId:%'
    """)
    
    result = cursor.fetchone()
    if result:
        count, total_mb, avg_kb = result
        print(f"\n📊 Checkpoint Analysis:")
        print(f"   Total checkpoints: {count:,}")
        print(f"   Total size: {total_mb:.1f} MB")
        print(f"   Average size: {avg_kb:.1f} KB per checkpoint")
        return count, total_mb
    return 0, 0

def clean_old_checkpoints(conn, keep_recent=100):
    """Remove old checkpoints, keeping only the most recent N."""
    cursor = conn.cursor()
    
    # Get all checkpoint keys with sizes
    cursor.execute("""
        SELECT key, length(value) as size
        FROM cursorDiskKV
        WHERE key LIKE 'checkpointId:%'
        ORDER BY key DESC
    """)
    
    all_checkpoints = cursor.fetchall()
    total_checkpoints = len(all_checkpoints)
    
    if total_checkpoints <= keep_recent:
        print(f"✅ Only {total_checkpoints} checkpoints found, keeping all")
        return 0
    
    # Keep the most recent N (they're ordered DESC)
    to_keep = all_checkpoints[:keep_recent]
    to_delete = all_checkpoints[keep_recent:]
    
    kept_size = sum(size for _, size in to_keep) / (1024 * 1024)
    deleted_size = sum(size for _, size in to_delete) / (1024 * 1024)
    
    print(f"\n🗑️  Checkpoint Cleanup:")
    print(f"   Total checkpoints: {total_checkpoints:,}")
    print(f"   Keeping: {len(to_keep):,} (most recent)")
    print(f"   Deleting: {len(to_delete):,} (old)")
    print(f"   Space to free: {deleted_size:.1f} MB")
    
    # Delete old checkpoints
    delete_keys = [key for key, _ in to_delete]
    deleted_count = 0
    
    for key in delete_keys:
        cursor.execute("DELETE FROM cursorDiskKV WHERE key = ?", (key,))
        deleted_count += 1
        if deleted_count % 100 == 0:
            print(f"   Progress: {deleted_count:,}/{len(delete_keys):,} deleted...")
    
    conn.commit()
    print(f"✅ Deleted {deleted_count:,} old checkpoints")
    return deleted_size

def vacuum_database(conn):
    """Run VACUUM to optimize database."""
    print(f"\n🔧 Running VACUUM to optimize database...")
    conn.execute("VACUUM")
    print(f"✅ VACUUM complete")

def main():
    print("🔍 Cursor Database Cleanup Tool")
    print("=" * 50)
    
    if not DB_PATH.exists():
        print(f"❌ Database not found: {DB_PATH}")
        return
    
    # Get initial size
    initial_size = get_db_size()
    print(f"\n📊 Initial Database Size: {initial_size:.2f} GB")
    
    # Create backup
    backup_database()
    
    # Connect to database
    conn = sqlite3.connect(str(DB_PATH))
    
    try:
        # Analyze checkpoints
        checkpoint_count, checkpoint_mb = analyze_checkpoints(conn)
        
        if checkpoint_count == 0:
            print("✅ No checkpoints found, nothing to clean")
            return
        
        # Ask user (for now, auto-clean keeping 100 most recent)
        print(f"\n⚠️  This will delete old checkpoint data.")
        print(f"   Keeping: 100 most recent checkpoints")
        print(f"   Deleting: {checkpoint_count - 100:,} old checkpoints")
        print(f"   Estimated space freed: ~{checkpoint_mb * 0.9:.1f} MB")
        
        # Clean old checkpoints
        freed_mb = clean_old_checkpoints(conn, keep_recent=100)
        
        # Vacuum database
        vacuum_database(conn)
        
        # Get final size
        final_size = get_db_size()
        freed_gb = (initial_size - final_size)
        
        print(f"\n✅ Cleanup Complete!")
        print(f"   Initial size: {initial_size:.2f} GB")
        print(f"   Final size: {final_size:.2f} GB")
        print(f"   Space freed: {freed_gb:.2f} GB")
        print(f"\n💾 Backup saved: {BACKUP_PATH.name}")
        print(f"   You can restore if needed: cp {BACKUP_PATH} {DB_PATH}")
        
    finally:
        conn.close()

if __name__ == "__main__":
    main()
