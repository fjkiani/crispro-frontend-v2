-- Main Agent Chat Export Script
-- Exports all chat and agent-related data from Cursor's state database
-- Generated: $(date)

-- Create export file
.output main_agent_chat_export.sql

-- Header
SELECT '-- Main Agent Chat Export' AS '';
SELECT '-- Generated from: /Users/fahadkiani/Library/Application Support/Cursor/User/globalStorage/state.vscdb' AS '';
SELECT '-- Export Date: ' || datetime('now') AS '';
SELECT '' AS '';

-- Export chat-related keys
SELECT '-- ============================================' AS '';
SELECT '-- CHAT-RELATED DATA' AS '';
SELECT '-- ============================================' AS '';
SELECT '' AS '';

SELECT 'INSERT INTO ItemTable (key, value) VALUES (' || 
       quote(key) || ', ' || 
       quote(value) || ');'
FROM ItemTable 
WHERE key LIKE '%chat%' 
   OR key LIKE '%agent%'
   OR key LIKE '%composer%'
   OR key LIKE '%aichat%'
ORDER BY key;

-- Export agent layout data
SELECT '' AS '';
SELECT '-- ============================================' AS '';
SELECT '-- AGENT LAYOUT DATA' AS '';
SELECT '-- ============================================' AS '';
SELECT '' AS '';

SELECT 'INSERT INTO ItemTable (key, value) VALUES (' || 
       quote(key) || ', ' || 
       quote(value) || ');'
FROM ItemTable 
WHERE key LIKE '%agentLayout%'
   OR key LIKE '%workbench.view.agents%'
ORDER BY key;

-- Export AI code tracking (composer-generated code)
SELECT '' AS '';
SELECT '-- ============================================' AS '';
SELECT '-- AI CODE TRACKING (Composer History)' AS '';
SELECT '-- ============================================' AS '';
SELECT '' AS '';

SELECT 'INSERT INTO ItemTable (key, value) VALUES (' || 
       quote(key) || ', ' || 
       quote(value) || ');'
FROM ItemTable 
WHERE key = 'aiCodeTrackingLines';

-- Summary statistics
SELECT '' AS '';
SELECT '-- ============================================' AS '';
SELECT '-- EXPORT SUMMARY' AS '';
SELECT '-- ============================================' AS '';
SELECT '' AS '';

SELECT '-- Total chat-related keys: ' || COUNT(*) AS ''
FROM ItemTable 
WHERE key LIKE '%chat%' OR key LIKE '%agent%' OR key LIKE '%composer%' OR key LIKE '%aichat%';

SELECT '-- Export completed: ' || datetime('now') AS '';

.output stdout





