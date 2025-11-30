#!/bin/bash
# Clear Cursor's cache directories to force UI refresh

echo "🧹 Clearing Cursor cache directories..."

CURSOR_SUPPORT="$HOME/Library/Application Support/Cursor"

# Close Cursor first
echo "🛑 Closing Cursor..."
pkill -f "Cursor" 2>/dev/null
sleep 2

# Clear cache directories (safe to delete - they'll regenerate)
echo "🗑️  Clearing caches..."

# Main cache
rm -rf "$CURSOR_SUPPORT/Cache" 2>/dev/null
echo "   ✅ Cleared Cache/"

# GPUCache
rm -rf "$CURSOR_SUPPORT/GPUCache" 2>/dev/null
echo "   ✅ Cleared GPUCache/"

# CachedData
rm -rf "$CURSOR_SUPPORT/CachedData" 2>/dev/null
echo "   ✅ Cleared CachedData/"

# WebStorage cache
rm -rf "$CURSOR_SUPPORT/WebStorage" 2>/dev/null
echo "   ✅ Cleared WebStorage/"

# Partition caches (workspace-specific)
rm -rf "$CURSOR_SUPPORT/Partitions"/*/Cache 2>/dev/null
rm -rf "$CURSOR_SUPPORT/Partitions"/*/GPUCache 2>/dev/null
echo "   ✅ Cleared Partition caches/"

echo ""
echo "✅ Cache cleared!"
echo ""
echo "🔄 Now restart Cursor - conversations should appear!"













