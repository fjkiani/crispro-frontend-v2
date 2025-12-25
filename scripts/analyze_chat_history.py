#!/usr/bin/env python3
"""
Systematic analysis of chat history to extract key learnings.
Focuses on: decisions, patterns, lessons, architectural choices, code patterns.
"""

import re
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set
import sys

def extract_conversation_blocks(content: str) -> List[Dict]:
    """Extract individual conversations from the backup file."""
    conversations = []
    current_conv = None
    current_messages = []
    
    lines = content.split('\n')
    for i, line in enumerate(lines):
        # Conversation header
        if line.startswith('CONVERSATION:'):
            if current_conv:
                current_conv['messages'] = current_messages
                conversations.append(current_conv)
            conv_id = line.split(':')[1].strip()
            current_conv = {
                'id': conv_id,
                'messages': []
            }
            current_messages = []
        
        # Message markers
        elif line.startswith('[') and ('[USER]' in line or '[ASSISTANT]' in line):
            if current_messages:
                current_messages[-1]['content'] += '\n'
            current_messages.append({
                'type': 'USER' if '[USER]' in line else 'ASSISTANT',
                'content': ''
            })
        
        # Content
        elif current_messages and not line.startswith('===') and not line.startswith('---'):
            if current_messages[-1]['content']:
                current_messages[-1]['content'] += '\n' + line
            else:
                current_messages[-1]['content'] = line
    
    if current_conv:
        current_conv['messages'] = current_messages
        conversations.append(current_conv)
    
    return conversations

def extract_key_patterns(content: str) -> Dict:
    """Extract key patterns: decisions, lessons, code patterns, etc."""
    patterns = {
        'decisions': [],
        'lessons_learned': [],
        'code_patterns': [],
        'architectural_choices': [],
        'bug_fixes': [],
        'manager_guidance': [],
        'doctrines': [],
        'test_results': [],
        'deployment_notes': [],
        'api_changes': []
    }
    
    # Decision markers
    decision_patterns = [
        r'✅.*[Dd]ecision',
        r'⚔️.*[Cc]ommander',
        r'[Mm]anager.*[Aa]nswer',
        r'[Pp]roceed with',
        r'[Aa]pproved',
        r'[Ff]inal.*[Dd]ecision'
    ]
    
    # Lesson learned markers
    lesson_patterns = [
        r'[Ll]esson.*[Ll]earned',
        r'[Ww]hat.*[Ww]ent.*[Ww]rong',
        r'[Bb]ug.*[Ff]ound',
        r'[Ff]ix.*[Aa]pplied',
        r'[Dd]octrine',
        r'[Nn]ever.*[Aa]gain',
        r'[Aa]lways.*[Dd]o',
        r'[Cc]ritical.*[Ii]nsight'
    ]
    
    # Code pattern markers
    code_patterns = [
        r'def\s+\w+',
        r'class\s+\w+',
        r'@\w+',
        r'async def',
        r'POST /api/',
        r'GET /api/',
        r'File:.*\.py',
        r'File:.*\.jsx',
        r'Location:.*\.py'
    ]
    
    # Architectural markers
    arch_patterns = [
        r'[Aa]rchitecture',
        r'[Dd]esign.*[Pp]attern',
        r'[Ss]ervice.*[Ss]eparation',
        r'[Mm]odular',
        r'[Oo]rchestrator',
        r'[Ff]ramework',
        r'[Ss]trategy'
    ]
    
    # Bug fix markers
    bug_patterns = [
        r'[Bb]ug.*[Ff]ixed',
        r'[Ff]ixed.*[Ii]ssue',
        r'[Ee]rror.*[Ff]ound',
        r'[Tt]ypeError',
        r'[Ii]mportError',
        r'[Aa]ttributeError',
        r'[Ff]ix.*[Aa]pplied'
    ]
    
    # Manager guidance markers
    manager_patterns = [
        r'[Mm]anager.*[Qq]uestion',
        r'[Mm]anager.*[Aa]nswer',
        r'[Cc]ommander.*[Oo]rder',
        r'[Pp]olicy',
        r'[Dd]octrine.*[Mm]anager'
    ]
    
    # Test result markers
    test_patterns = [
        r'\d+/\d+.*[Tt]est.*[Pp]ass',
        r'[Aa]ll.*[Tt]est.*[Pp]ass',
        r'[Tt]est.*[Cc]overage',
        r'[Aa]UROC',
        r'[Aa]UPRC',
        r'[Pp]recision',
        r'[Rr]ecall'
    ]
    
    # Deployment markers
    deploy_patterns = [
        r'[Dd]eploy',
        r'[Mm]odal.*[Ss]ervice',
        r'[Pp]roduction',
        r'[Rr]elease',
        r'[Vv]ersion'
    ]
    
    # API change markers
    api_patterns = [
        r'POST /api/',
        r'GET /api/',
        r'[Ee]ndpoint.*[Cc]hange',
        r'[Aa]PI.*[Cc]ontract',
        r'[Bb]reaking.*[Cc]hange'
    ]
    
    # Extract patterns
    for pattern_list, key in [
        (decision_patterns, 'decisions'),
        (lesson_patterns, 'lessons_learned'),
        (code_patterns, 'code_patterns'),
        (arch_patterns, 'architectural_choices'),
        (bug_patterns, 'bug_fixes'),
        (manager_patterns, 'manager_guidance'),
        (test_patterns, 'test_results'),
        (deploy_patterns, 'deployment_notes'),
        (api_patterns, 'api_changes')
    ]:
        for pattern in pattern_list:
            matches = re.finditer(pattern, content, re.IGNORECASE | re.MULTILINE)
            for match in matches:
                # Extract context (50 chars before and after)
                start = max(0, match.start() - 50)
                end = min(len(content), match.end() + 50)
                context = content[start:end].strip()
                patterns[key].append({
                    'match': match.group(),
                    'context': context,
                    'position': match.start()
                })
    
    # Extract doctrines (look for doctrine markers)
    doctrine_matches = re.finditer(
        r'[Dd]octrine.*?[:\-]|⚔️.*[Dd]octrine|DOCTRINE.*STATUS',
        content,
        re.MULTILINE
    )
    for match in doctrine_matches:
        start = max(0, match.start() - 100)
        end = min(len(content), match.end() + 200)
        context = content[start:end].strip()
        patterns['doctrines'].append({
            'match': match.group(),
            'context': context,
            'position': match.start()
        })
    
    return patterns

def extract_file_references(content: str) -> Dict[str, List[str]]:
    """Extract file paths and code references."""
    file_refs = defaultdict(list)
    
    # File path patterns
    file_patterns = [
        r'`([^`]+\.(py|jsx|tsx|ts|js|json|md|mdc))`',
        r'File:\s*([^\s]+\.(py|jsx|tsx|ts|js|json|md|mdc))',
        r'Location:\s*([^\s]+\.(py|jsx|tsx|ts|js|json|md|mdc))',
        r'([a-zA-Z0-9_/]+\.(py|jsx|tsx|ts|js|json|md|mdc))',
    ]
    
    for pattern in file_patterns:
        matches = re.finditer(pattern, content)
        for match in matches:
            filepath = match.group(1) if match.lastindex >= 1 else match.group(0)
            if '/' in filepath or filepath.startswith('api/') or filepath.startswith('src/'):
                file_refs[filepath].append(match.start())
    
    return dict(file_refs)

def extract_technical_decisions(content: str) -> List[Dict]:
    """Extract technical decisions and rationale."""
    decisions = []
    
    # Look for decision patterns with rationale
    decision_blocks = re.finditer(
        r'(?:Decision|Choice|Approach|Strategy|Solution).*?(?=\n\n|\n===|\n---|$)',
        content,
        re.IGNORECASE | re.DOTALL
    )
    
    for match in decision_blocks:
        block = match.group(0)
        # Extract rationale
        rationale_match = re.search(r'(?:because|since|due to|rationale|reason|why).*?[.!?]', block, re.IGNORECASE)
        rationale = rationale_match.group(0) if rationale_match else None
        
        decisions.append({
            'decision': block[:200],  # First 200 chars
            'rationale': rationale,
            'position': match.start()
        })
    
    return decisions

def main():
    """Main analysis function."""
    chat_file = Path('/Users/fahadkiani/Desktop/development/crispr-assistant-main/ALL_CHAT_HISTORY_FROM_REAL_BACKUP.txt')
    
    print(f"📊 Analyzing chat history: {chat_file}")
    print(f"📏 File size: {chat_file.stat().st_size / 1024 / 1024:.2f} MB")
    
    # Read in chunks to avoid memory issues
    print("📖 Reading file...")
    with open(chat_file, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    print(f"✅ Read {len(content)} characters")
    print("🔍 Extracting patterns...")
    
    # Extract patterns
    patterns = extract_key_patterns(content)
    
    # Extract file references
    file_refs = extract_file_references(content)
    
    # Extract technical decisions
    decisions = extract_technical_decisions(content)
    
    # Summary statistics
    print("\n" + "="*80)
    print("📊 EXTRACTION SUMMARY")
    print("="*80)
    for key, items in patterns.items():
        print(f"  {key}: {len(items)} items")
    print(f"  File references: {len(file_refs)} unique files")
    print(f"  Technical decisions: {len(decisions)} items")
    
    # Save results
    output_dir = Path('/Users/fahadkiani/Desktop/development/crispr-assistant-main/.cursor/chat_history_analysis')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save patterns
    patterns_file = output_dir / 'extracted_patterns.json'
    with open(patterns_file, 'w') as f:
        json.dump({
            'patterns': {k: v[:100] for k, v in patterns.items()},  # Limit to first 100 per category
            'file_references': {k: len(v) for k, v in list(file_refs.items())[:200]},
            'decisions': decisions[:100]
        }, f, indent=2)
    
    print(f"\n✅ Results saved to: {output_dir}")
    print(f"   - Patterns: {patterns_file}")
    
    # Generate summary report
    summary_file = output_dir / 'analysis_summary.md'
    with open(summary_file, 'w') as f:
        f.write("# Chat History Analysis Summary\n\n")
        f.write(f"**Total Conversations:** 177\n")
        f.write(f"**Total Messages:** 66,858\n")
        f.write(f"**File Size:** {chat_file.stat().st_size / 1024 / 1024:.2f} MB\n\n")
        
        f.write("## Key Patterns Extracted\n\n")
        for key, items in patterns.items():
            f.write(f"### {key.replace('_', ' ').title()}\n")
            f.write(f"**Count:** {len(items)}\n\n")
            if items:
                # Show top 5 examples
                for i, item in enumerate(items[:5], 1):
                    f.write(f"{i}. `{item.get('match', 'N/A')[:50]}`\n")
                    f.write(f"   Context: {item.get('context', '')[:100]}...\n\n")
            f.write("\n")
        
        f.write("## Most Referenced Files\n\n")
        sorted_files = sorted(file_refs.items(), key=lambda x: len(x[1]), reverse=True)
        for filepath, positions in sorted_files[:20]:
            f.write(f"- `{filepath}` ({len(positions)} references)\n")
    
    print(f"   - Summary: {summary_file}")
    print("\n✅ Analysis complete!")

if __name__ == '__main__':
    main()


























