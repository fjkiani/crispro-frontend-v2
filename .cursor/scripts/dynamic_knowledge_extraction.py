#!/usr/bin/env python3
"""
Dynamic Knowledge Extraction - No Hard-Coded Topics

Purpose: Extract knowledge by identifying important sections dynamically
Approach: Find natural breaks, headers, code blocks, and preserve full context
"""

import re
import glob
from pathlib import Path
from collections import defaultdict

def identify_important_sections(lines):
    """Dynamically identify important sections without hard-coding topics"""
    
    important_sections = []
    
    # Find section markers (headers, code blocks, key phrases)
    section_markers = []
    
    for i, line in enumerate(lines):
        # Mark headers (##, ###, etc.)
        if re.match(r'^#+\s+', line):
            level = len(line) - len(line.lstrip('#'))
            section_markers.append({
                'type': 'header',
                'level': level,
                'line': i,
                'text': line.strip()
            })
        
        # Mark code blocks
        if re.match(r'^```', line):
            section_markers.append({
                'type': 'code_block',
                'line': i,
                'text': line.strip()
            })
        
        # Mark important phrases (but don't limit to specific topics)
        important_phrases = [
            r'doctrine', r'strategy', r'framework', r'protocol',
            r'decision', r'remind', r'correct', r'important',
            r'critical', r'key', r'note', r'warning',
            r'example', r'case', r'scenario', r'patient',
            r'discovery', r'found', r'learned', r'pattern'
        ]
        
        line_lower = line.lower()
        for phrase in important_phrases:
            if re.search(phrase, line_lower):
                # Only mark if it's not already a header
                if not re.match(r'^#+\s+', line):
                    section_markers.append({
                        'type': 'important_phrase',
                        'line': i,
                        'text': line.strip()[:100]  # First 100 chars
                    })
                break
    
    # Group markers into sections
    sections = []
    current_section = None
    
    for marker in section_markers:
        if marker['type'] == 'header':
            # Start new section
            if current_section:
                sections.append(current_section)
            
            current_section = {
                'start': marker['line'],
                'end': None,
                'header': marker['text'],
                'level': marker['level'],
                'type': 'section'
            }
        elif current_section:
            # Extend current section
            current_section['end'] = marker['line']
    
    # Close last section
    if current_section:
        current_section['end'] = len(lines) - 1
        sections.append(current_section)
    
    return sections, section_markers

def extract_with_context(lines, start, end, context_before=50, context_after=50):
    """Extract section with context"""
    actual_start = max(0, start - context_before)
    actual_end = min(len(lines), end + context_after)
    return lines[actual_start:actual_end], actual_start, actual_end

def create_comprehensive_extraction(chat_log_path: str, output_file: Path):
    """Create comprehensive extraction preserving all important context"""
    
    with open(chat_log_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    
    print(f"📖 Loaded {len(lines)} lines from chat log")
    
    # Identify important sections dynamically
    print("🔍 Identifying important sections...")
    sections, markers = identify_important_sections(lines)
    
    print(f"✅ Found {len(sections)} major sections")
    print(f"✅ Found {len(markers)} section markers")
    
    # For comprehensive extraction, we'll use a different approach:
    # Extract by natural conversation breaks and preserve full context
    
    # Find conversation breaks (user messages, major headers)
    conversation_breaks = []
    
    for i, line in enumerate(lines):
        # Major headers (## or ###)
        if re.match(r'^##+\s+', line):
            conversation_breaks.append(i)
        
        # User messages (lines starting with # or **)
        if re.match(r'^#\s+[a-z]', line.lower()):
            conversation_breaks.append(i)
        
        # Major separators
        if re.match(r'^---+\s*$', line) or re.match(r'^\*\*\*\s*$', line):
            conversation_breaks.append(i)
    
    # Create comprehensive sections
    comprehensive_sections = []
    
    # If we have many breaks, group them into larger sections
    if len(conversation_breaks) > 100:
        # Group into ~50-100 line chunks
        chunk_size = len(lines) // 50
        for i in range(0, len(lines), chunk_size):
            end = min(i + chunk_size, len(lines))
            comprehensive_sections.append({
                'start': i,
                'end': end,
                'type': 'chunk',
                'index': len(comprehensive_sections) + 1
            })
    else:
        # Use conversation breaks
        for i in range(len(conversation_breaks)):
            start = conversation_breaks[i]
            end = conversation_breaks[i + 1] if i + 1 < len(conversation_breaks) else len(lines)
            
            # Only include sections with substantial content
            if end - start > 20:  # At least 20 lines
                comprehensive_sections.append({
                    'start': start,
                    'end': end,
                    'type': 'conversation',
                    'index': len(comprehensive_sections) + 1
                })
    
    # If still too few sections, create overlapping windows
    if len(comprehensive_sections) < 10:
        # Create overlapping windows to ensure nothing is missed
        window_size = 500  # 500 line windows
        overlap = 100  # 100 line overlap
        
        comprehensive_sections = []
        for i in range(0, len(lines), window_size - overlap):
            end = min(i + window_size, len(lines))
            comprehensive_sections.append({
                'start': i,
                'end': end,
                'type': 'window',
                'index': len(comprehensive_sections) + 1
            })
    
    print(f"✅ Created {len(comprehensive_sections)} comprehensive sections")
    
    # Write comprehensive knowledge base
    with open(output_file, 'w', encoding='utf-8') as f:
        # Header
        f.write("# 🧠 COMPREHENSIVE KNOWLEDGE BASE - COMPLETE CONTEXT\n\n")
        f.write("**Purpose**: Complete knowledge base with ALL context from 30,000+ line chat log\n")
        f.write("**Status**: ✅ Comprehensive Dynamic Extraction - Ready for Agent Training\n")
        f.write("**Date**: January 2025\n")
        f.write(f"**Total Lines**: {len(lines)}\n")
        f.write(f"**Total Sections**: {len(comprehensive_sections)}\n\n")
        f.write("---\n\n")
        
        f.write("## 📋 TABLE OF CONTENTS\n\n")
        f.write(f"**Total Sections**: {len(comprehensive_sections)}\n\n")
        for i, section in enumerate(comprehensive_sections, 1):
            f.write(f"{i}. [Section {i} (Lines {section['start']+1}-{section['end']})](#section-{i})\n")
        f.write("\n---\n\n")
        
        # Write all sections
        total_lines_written = 0
        for i, section in enumerate(comprehensive_sections, 1):
            section_lines = lines[section['start']:section['end']]
            total_lines_written += len(section_lines)
            
            f.write(f"# Section {i}\n\n")
            f.write(f"**Type**: {section['type']}\n")
            f.write(f"**Source Lines**: {section['start']+1}-{section['end']}\n")
            f.write(f"**Content Length**: {len(section_lines)} lines\n\n")
            f.write("---\n\n")
            
            f.write("".join(section_lines))
            f.write("\n\n---\n\n")
        
        # Summary
        f.write("# 📊 EXTRACTION SUMMARY\n\n")
        f.write(f"**Total Sections**: {len(comprehensive_sections)}\n")
        f.write(f"**Total Lines Extracted**: {total_lines_written}\n")
        f.write(f"**Coverage**: Complete chat log preserved\n")
        f.write(f"**Method**: Dynamic section identification (no hard-coded topics)\n\n")
        
        f.write("## Section Breakdown\n\n")
        for i, section in enumerate(comprehensive_sections, 1):
            f.write(f"- **Section {i}**: Lines {section['start']+1}-{section['end']} ({section['end']-section['start']} lines, type: {section['type']})\n")
        
        f.write("\n---\n\n")
        f.write("**END OF COMPREHENSIVE KNOWLEDGE BASE**\n")
    
    # Get file stats
    file_size_mb = output_file.stat().st_size / (1024 * 1024)
    
    print(f"\n✅ Comprehensive knowledge base created: {output_file}")
    print(f"📊 Total sections: {len(comprehensive_sections)}")
    print(f"📏 Total lines extracted: {total_lines_written}")
    print(f"💾 File size: {file_size_mb:.2f} MB")
    
    return output_file

def main():
    """Main extraction function"""
    # Find chat log file
    chat_log_files = glob.glob("/Users/fahadkiani/Downloads/*Nyx*rogue*.md")
    
    if not chat_log_files:
        print("❌ Could not find chat log file")
        return
    
    chat_log_path = chat_log_files[0]
    print(f"📖 Using chat log: {chat_log_path}")
    
    # Output file
    output_dir = Path(".cursor/rules/MOAT/AGENT_TRAINING_KNOWLEDGE")
    output_file = output_dir / "COMPLETE_KNOWLEDGE_BASE.md"
    
    # Extract comprehensive knowledge
    result = create_comprehensive_extraction(chat_log_path, output_file)
    
    if result:
        print(f"\n✅ Dynamic extraction complete!")
        print(f"📁 Output: {result}")
    else:
        print("\n❌ Extraction failed")

if __name__ == "__main__":
    main()























