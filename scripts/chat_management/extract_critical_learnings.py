#!/usr/bin/env python3
"""
Extract critical learnings from chat history for .cursorrules integration.
Focuses on: bugs, fixes, patterns, decisions, manager guidance.
"""

import re
import json
from pathlib import Path
from typing import List, Dict
from collections import defaultdict

def extract_lessons_learned_section(content: str) -> List[Dict]:
    """Extract explicit "Lessons Learned" sections."""
    lessons = []
    
    # Pattern for lessons learned sections
    patterns = [
        r'##\s*[🎖️⚔️]*\s*LESSONS?\s*LEARNED[^#]*',
        r'###\s*[Ww]hat\s+[Ww]ent\s+[Ww]rong[^#]*',
        r'[Ll]esson[^:]*:\s*[^\n]+',
        r'[Nn]ever\s+[Aa]gain[^#]*',
        r'[Aa]lways\s+[Dd]o[^#]*',
        r'[Cc]ritical\s+[Ii]nsight[^#]*',
    ]
    
    for pattern in patterns:
        matches = re.finditer(pattern, content, re.IGNORECASE | re.MULTILINE | re.DOTALL)
        for match in matches:
            lesson_text = match.group(0)
            # Extract the actual lesson content
            if len(lesson_text) > 50 and len(lesson_text) < 2000:
                lessons.append({
                    'text': lesson_text.strip(),
                    'position': match.start()
                })
    
    return lessons

def extract_bug_fixes(content: str) -> List[Dict]:
    """Extract bug fixes with context."""
    fixes = []
    
    # Look for bug fix patterns
    bug_patterns = [
        r'BUG\s+FIXED[^#]*?(?=\n\n|\n===|$)',
        r'Fixed.*?issue[^#]*?(?=\n\n|\n===|$)',
        r'Error.*?found[^#]*?(?=\n\n|\n===|$)',
        r'TypeError[^#]*?(?=\n\n|\n===|$)',
        r'ImportError[^#]*?(?=\n\n|\n===|$)',
        r'AttributeError[^#]*?(?=\n\n|\n===|$)',
    ]
    
    for pattern in bug_patterns:
        matches = re.finditer(pattern, content, re.IGNORECASE | re.MULTILINE | re.DOTALL)
        for match in matches:
            fix_text = match.group(0)
            if len(fix_text) > 100 and len(fix_text) < 3000:
                # Extract file references
                file_refs = re.findall(r'`([^`]+\.(py|jsx|tsx|js|ts))`', fix_text)
                fixes.append({
                    'text': fix_text.strip(),
                    'files': [f[0] for f in file_refs],
                    'position': match.start()
                })
    
    return fixes

def extract_manager_policies(content: str) -> List[Dict]:
    """Extract manager guidance and policies."""
    policies = []
    
    # Manager answer patterns
    patterns = [
        r'[Mm]anager.*?[Aa]nswer[^#]*?(?=\n\n|\n===|$)',
        r'[Mm]anager.*?[Pp]olicy[^#]*?(?=\n\n|\n===|$)',
        r'[Cc]ommander.*?[Oo]rder[^#]*?(?=\n\n|\n===|$)',
        r'Q\d+[^#]*?[Aa]nswer[^#]*?(?=\n\n|\n===|$)',
        r'[Pp]olicy[^:]*:\s*[^\n]+',
    ]
    
    for pattern in patterns:
        matches = re.finditer(pattern, content, re.IGNORECASE | re.MULTILINE | re.DOTALL)
        for match in matches:
            policy_text = match.group(0)
            if len(policy_text) > 50 and len(policy_text) < 2000:
                policies.append({
                    'text': policy_text.strip(),
                    'position': match.start()
                })
    
    return policies

def extract_architectural_decisions(content: str) -> List[Dict]:
    """Extract architectural decisions with rationale."""
    decisions = []
    
    # Decision patterns with rationale
    patterns = [
        r'[Aa]rchitectural.*?[Dd]ecision[^#]*?(?=\n\n|\n===|$)',
        r'[Dd]esign.*?[Pp]attern[^#]*?(?=\n\n|\n===|$)',
        r'[Ss]trategy[^:]*:\s*[^#]+?(?=\n\n|\n===|$)',
        r'[Ww]hy.*?[Ww]e.*?[Bb]uilt[^#]*?(?=\n\n|\n===|$)',
        r'[Rr]ationale[^:]*:\s*[^#]+?(?=\n\n|\n===|$)',
    ]
    
    for pattern in patterns:
        matches = re.finditer(pattern, content, re.IGNORECASE | re.MULTILINE | re.DOTALL)
        for match in matches:
            decision_text = match.group(0)
            if len(decision_text) > 100 and len(decision_text) < 3000:
                decisions.append({
                    'text': decision_text.strip(),
                    'position': match.start()
                })
    
    return decisions

def extract_code_patterns(content: str) -> List[Dict]:
    """Extract reusable code patterns and best practices."""
    patterns = []
    
    # Look for code blocks with explanations
    code_block_pattern = r'```(?:python|javascript|typescript|bash|json)[\s\S]*?```'
    code_blocks = re.finditer(code_block_pattern, content)
    
    for match in code_blocks:
        code_block = match.group(0)
        # Look for context before the code block
        start = max(0, match.start() - 200)
        context_before = content[start:match.start()].strip()
        
        # Check if there's explanatory text
        if any(keyword in context_before.lower() for keyword in ['pattern', 'best practice', 'always', 'never', 'should', 'must']):
            patterns.append({
                'code': code_block,
                'context': context_before[-200:],
                'position': match.start()
            })
    
    return patterns

def extract_doctrine_statements(content: str) -> List[Dict]:
    """Extract doctrine statements."""
    doctrines = []
    
    # Doctrine markers
    patterns = [
        r'DOCTRINE[^:]*:\s*[^#]+?(?=\n\n|\n===|$)',
        r'⚔️.*?DOCTRINE[^#]*?(?=\n\n|\n===|$)',
        r'[Tt]he\s+[A-Z][a-z]+\s+[Dd]octrine[^#]*?(?=\n\n|\n===|$)',
        r'[Dd]octrine.*?STATUS[^#]*?(?=\n\n|\n===|$)',
    ]
    
    for pattern in patterns:
        matches = re.finditer(pattern, content, re.IGNORECASE | re.MULTILINE | re.DOTALL)
        for match in matches:
            doctrine_text = match.group(0)
            if len(doctrine_text) > 50 and len(doctrine_text) < 2000:
                doctrines.append({
                    'text': doctrine_text.strip(),
                    'position': match.start()
                })
    
    return doctrines

def main():
    """Extract critical learnings."""
    chat_file = Path('/Users/fahadkiani/Desktop/development/crispr-assistant-main/ALL_CHAT_HISTORY_FROM_REAL_BACKUP.txt')
    
    print("📖 Reading chat history...")
    with open(chat_file, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    print("🔍 Extracting critical learnings...")
    
    # Extract different types of learnings
    lessons = extract_lessons_learned_section(content)
    bug_fixes = extract_bug_fixes(content)
    manager_policies = extract_manager_policies(content)
    arch_decisions = extract_architectural_decisions(content)
    code_patterns_list = extract_code_patterns(content)
    doctrines = extract_doctrine_statements(content)
    
    # Save results
    output_dir = Path('/Users/fahadkiani/Desktop/development/crispr-assistant-main/.cursor/chat_history_analysis')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create comprehensive report
    report_file = output_dir / 'critical_learnings.md'
    with open(report_file, 'w') as f:
        f.write("# Critical Learnings from Chat History\n\n")
        f.write(f"**Source:** {chat_file.name}\n")
        f.write(f"**Total Conversations:** 177\n")
        f.write(f"**Total Messages:** 66,858\n\n")
        
        f.write("## 📚 Lessons Learned\n\n")
        for i, lesson in enumerate(lessons[:50], 1):  # Top 50
            f.write(f"### Lesson {i}\n\n")
            f.write(f"{lesson['text'][:500]}...\n\n")
            f.write("---\n\n")
        
        f.write("## 🐛 Bug Fixes\n\n")
        for i, fix in enumerate(bug_fixes[:50], 1):  # Top 50
            f.write(f"### Fix {i}\n\n")
            if fix['files']:
                f.write(f"**Files:** {', '.join(fix['files'])}\n\n")
            f.write(f"{fix['text'][:500]}...\n\n")
            f.write("---\n\n")
        
        f.write("## 👔 Manager Policies\n\n")
        for i, policy in enumerate(manager_policies[:30], 1):  # Top 30
            f.write(f"### Policy {i}\n\n")
            f.write(f"{policy['text'][:500]}...\n\n")
            f.write("---\n\n")
        
        f.write("## 🏗️ Architectural Decisions\n\n")
        for i, decision in enumerate(arch_decisions[:30], 1):  # Top 30
            f.write(f"### Decision {i}\n\n")
            f.write(f"{decision['text'][:500]}...\n\n")
            f.write("---\n\n")
        
        f.write("## ⚔️ Doctrines\n\n")
        for i, doctrine in enumerate(doctrines[:30], 1):  # Top 30
            f.write(f"### Doctrine {i}\n\n")
            f.write(f"{doctrine['text'][:500]}...\n\n")
            f.write("---\n\n")
    
    # Save JSON for programmatic access
    json_file = output_dir / 'critical_learnings.json'
    with open(json_file, 'w') as f:
        json.dump({
            'lessons_learned': lessons[:100],
            'bug_fixes': bug_fixes[:100],
            'manager_policies': manager_policies[:50],
            'architectural_decisions': arch_decisions[:50],
            'code_patterns': code_patterns_list[:50],
            'doctrines': doctrines[:50]
        }, f, indent=2)
    
    print(f"\n✅ Critical learnings extracted:")
    print(f"   - Lessons Learned: {len(lessons)}")
    print(f"   - Bug Fixes: {len(bug_fixes)}")
    print(f"   - Manager Policies: {len(manager_policies)}")
    print(f"   - Architectural Decisions: {len(arch_decisions)}")
    print(f"   - Code Patterns: {len(code_patterns_list)}")
    print(f"   - Doctrines: {len(doctrines)}")
    print(f"\n📄 Reports saved to:")
    print(f"   - {report_file}")
    print(f"   - {json_file}")

if __name__ == '__main__':
    main()


























