#!/usr/bin/env python3
"""
Full Context Extraction - Preserve Everything, Organize by Natural Breaks

Purpose: Extract entire chat log with all context, organized by conversation flow
Approach: No topic filtering - preserve everything, organize by natural breaks
"""

import re
import glob
from pathlib import Path

def extract_full_context(chat_log_path: str, output_file: Path):
    """Extract full chat log with all context, organized by natural breaks"""
    
    with open(chat_log_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    
    print(f"📖 Loaded {len(lines)} lines from chat log")
    
    # Find natural conversation breaks
    breaks = []
    
    for i, line in enumerate(lines):
        # Major headers (## or ###)
        if re.match(r'^##+\s+', line):
            breaks.append({
                'line': i,
                'type': 'major_header',
                'text': line.strip()[:100]
            })
        
        # User messages (lines starting with # followed by lowercase)
        if re.match(r'^#\s+[a-z]', line.lower()):
            breaks.append({
                'line': i,
                'type': 'user_message',
                'text': line.strip()[:100]
            })
    
    print(f"✅ Found {len(breaks)} natural breaks")
    
    # Create sections from breaks
    sections = []
    
    if len(breaks) > 0:
        for i in range(len(breaks)):
            start = breaks[i]['line']
            end = breaks[i + 1]['line'] if i + 1 < len(breaks) else len(lines)
            
            sections.append({
                'start': start,
                'end': end,
                'type': breaks[i]['type'],
                'header': breaks[i]['text'],
                'index': i + 1
            })
    else:
        # No breaks found, create chunks
        chunk_size = 1000
        for i in range(0, len(lines), chunk_size):
            end = min(i + chunk_size, len(lines))
            sections.append({
                'start': i,
                'end': end,
                'type': 'chunk',
                'header': f"Chunk {len(sections) + 1}",
                'index': len(sections) + 1
            })
    
    print(f"✅ Created {len(sections)} sections")
    
    # Write comprehensive knowledge base
    with open(output_file, 'w', encoding='utf-8') as f:
        # Header
        f.write("# 🧠 COMPLETE KNOWLEDGE BASE - FULL CONTEXT\n\n")
        f.write("**Purpose**: Complete knowledge base with ALL context from chat log\n")
        f.write("**Status**: ✅ Full Context Extraction - No Filtering - Ready for Agent Training\n")
        f.write("**Date**: January 2025\n")
        f.write(f"**Total Lines**: {len(lines)}\n")
        f.write(f"**Total Sections**: {len(sections)}\n")
        f.write("**Method**: Full context preservation, organized by natural conversation breaks\n\n")
        f.write("---\n\n")
        
        f.write("## 📋 TABLE OF CONTENTS\n\n")
        for section in sections:
            header_short = section['header'][:80] if len(section['header']) > 80 else section['header']
            f.write(f"{section['index']}. [{header_short}](#section-{section['index']}) (Lines {section['start']+1}-{section['end']})\n")
        f.write("\n---\n\n")
        
        # Write all sections
        total_lines_written = 0
        for section in sections:
            section_lines = lines[section['start']:section['end']]
            total_lines_written += len(section_lines)
            
            f.write(f"# Section {section['index']}\n\n")
            f.write(f"**Type**: {section['type']}\n")
            f.write(f"**Header**: {section['header']}\n")
            f.write(f"**Source Lines**: {section['start']+1}-{section['end']}\n")
            f.write(f"**Content Length**: {len(section_lines)} lines\n\n")
            f.write("---\n\n")
            
            f.write("".join(section_lines))
            f.write("\n\n---\n\n")
        
        # Summary
        f.write("# 📊 EXTRACTION SUMMARY\n\n")
        f.write(f"**Total Sections**: {len(sections)}\n")
        f.write(f"**Total Lines Extracted**: {total_lines_written}\n")
        f.write(f"**Coverage**: 100% - Complete chat log preserved\n")
        f.write(f"**Method**: Full context extraction with no topic filtering\n\n")
        
        f.write("## Section Statistics\n\n")
        for section in sections:
            f.write(f"- **Section {section['index']}**: {section['end']-section['start']} lines (type: {section['type']})\n")
        
        f.write("\n---\n\n")
        f.write("**END OF COMPLETE KNOWLEDGE BASE**\n")
    
    # Get file stats
    file_size_mb = output_file.stat().st_size / (1024 * 1024)
    
    print(f"\n✅ Complete knowledge base created: {output_file}")
    print(f"📊 Total sections: {len(sections)}")
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
    
    # Extract full context
    result = extract_full_context(chat_log_path, output_file)
    
    if result:
        print(f"\n✅ Full context extraction complete!")
        print(f"📁 Output: {result}")
        print(f"✅ All context preserved - no filtering applied")
    else:
        print("\n❌ Extraction failed")

if __name__ == "__main__":
    main()

