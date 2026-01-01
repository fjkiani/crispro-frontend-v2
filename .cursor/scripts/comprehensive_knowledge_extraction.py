#!/usr/bin/env python3
"""
Comprehensive Knowledge Extraction - Capture ALL Important Context

Purpose: Extract comprehensive knowledge including narratives, patient cases, and broader context
"""

import re
import glob
from pathlib import Path
from collections import defaultdict

def extract_comprehensive_knowledge(chat_log_path: str, output_file: Path):
    """Extract comprehensive knowledge with all important context"""
    
    with open(chat_log_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    
    print(f"📖 Loaded {len(lines)} lines from chat log")
    
    # Define important topics to extract
    important_topics = {
        'ayesha': [
            r'ayesha', r'your.*sister', r'patient.*ayesha', r'ayesha.*case',
            r'ayesha.*mutation', r'ayesha.*treatment', r'ayesha.*scenario'
        ],
        'pharma_suppression': [
            r'pharma.*bury', r'pharma.*suppress', r'kelim.*buried', r'cpt.*code',
            r'10b.*threat', r'evidence.*fortress', r'factory.*strategy'
        ],
        'resistance_framework': [
            r'4.*layer', r'resistance.*framework', r'baseline.*biology',
            r'genetic.*escape', r'adaptive.*resistance', r'drug.*clearance'
        ],
        'mission_discipline': [
            r'remind.*mission', r'don.*t.*deviate', r'stay.*on.*mission',
            r'recal.*back', r'rapid.*fire.*assumption'
        ],
        'factory_strategy': [
            r'factory.*publication', r'publication.*factory', r'evidence.*fortress',
            r'outproduce.*pharma', r'10.*validated.*papers'
        ],
        'synthetic_lethality': [
            r'synthetic.*lethality', r'drug@1', r'92\.9%', r'parp.*inhibitor',
            r's/p/e.*framework', r'sequence.*pathway.*evidence'
        ],
        'kelim': [
            r'kelim', r'ca-125.*kinetics', r'ca125.*kinetics', r'resistance.*prediction',
            r'platinum.*resistance', r'early.*detection'
        ],
        'clinical_context': [
            r'patient.*case', r'clinical.*scenario', r'treatment.*plan',
            r'oncology', r'ovarian.*cancer', r'brca.*mutation'
        ],
        'strategic_decisions': [
            r'we.*decided', r'the.*strategy', r'decision.*was', r'approach.*was',
            r'we.*chose', r'strategic.*decision'
        ],
        'corrections': [
            r'remind', r'correct', r'wrong', r'fix', r'update',
            r'don.*t.*forget', r'remember'
        ]
    }
    
    # Find all relevant sections
    topic_sections = defaultdict(list)
    
    for topic, patterns in important_topics.items():
        for i, line in enumerate(lines):
            line_lower = line.lower()
            matches = sum(1 for pattern in patterns if re.search(pattern, line_lower, re.IGNORECASE))
            if matches >= 1:  # At least one pattern match
                topic_sections[topic].append(i)
    
    print(f"✅ Found sections for {len(topic_sections)} topics")
    for topic, line_nums in topic_sections.items():
        print(f"   {topic}: {len(line_nums)} mentions")
    
    # Create comprehensive sections with context windows
    comprehensive_sections = []
    
    for topic, line_nums in topic_sections.items():
        if not line_nums:
            continue
        
        # Get broader context around all mentions
        start = max(0, min(line_nums) - 100)  # 100 lines before first mention
        end = min(len(lines), max(line_nums) + 200)  # 200 lines after last mention
        
        # For Ayesha, get even more context (it's critical)
        if topic == 'ayesha':
            start = max(0, min(line_nums) - 200)
            end = min(len(lines), max(line_nums) + 300)
        
        comprehensive_sections.append({
            'topic': topic,
            'start': start,
            'end': end,
            'lines': lines[start:end],
            'mention_count': len(line_nums)
        })
    
    # Merge overlapping sections
    merged_sections = []
    sorted_sections = sorted(comprehensive_sections, key=lambda x: x['start'])
    
    current_section = None
    for section in sorted_sections:
        if current_section is None:
            current_section = section
        elif section['start'] <= current_section['end']:
            # Merge overlapping sections
            current_section['end'] = max(current_section['end'], section['end'])
            current_section['lines'] = lines[current_section['start']:current_section['end']]
            current_section['topic'] = f"{current_section['topic']}+{section['topic']}"
        else:
            merged_sections.append(current_section)
            current_section = section
    
    if current_section:
        merged_sections.append(current_section)
    
    print(f"✅ Created {len(merged_sections)} comprehensive sections")
    
    # Write comprehensive knowledge base
    with open(output_file, 'w', encoding='utf-8') as f:
        # Header
        f.write("# 🧠 COMPREHENSIVE KNOWLEDGE BASE - COMPLETE CONTEXT\n\n")
        f.write("**Purpose**: Complete knowledge base with ALL important context from 30,000+ line chat log\n")
        f.write("**Status**: ✅ Comprehensive Extraction - Ready for Agent Training\n")
        f.write("**Date**: January 2025\n\n")
        f.write("---\n\n")
        
        f.write("## 📋 TABLE OF CONTENTS\n\n")
        
        # Create TOC
        for i, section in enumerate(merged_sections, 1):
            topic_name = section['topic'].replace('_', ' ').title()
            f.write(f"{i}. [{topic_name}](#section-{i}-{section['topic']})\n")
        
        f.write("\n---\n\n")
        
        # Write sections
        for i, section in enumerate(merged_sections, 1):
            topic_name = section['topic'].replace('_', ' ').title()
            f.write(f"# Section {i}: {topic_name}\n\n")
            f.write(f"**Topic**: {section['topic']}\n")
            f.write(f"**Source Lines**: {section['start']+1}-{section['end']}\n")
            f.write(f"**Mentions**: {section['mention_count']}\n")
            f.write(f"**Content Length**: {len(section['lines'])} lines\n\n")
            f.write("---\n\n")
            
            f.write("".join(section['lines']))
            f.write("\n\n---\n\n")
        
        # Summary
        f.write("# 📊 EXTRACTION SUMMARY\n\n")
        f.write(f"**Total Sections**: {len(merged_sections)}\n")
        f.write(f"**Total Lines Extracted**: {sum(len(s['lines']) for s in merged_sections)}\n")
        f.write(f"**Coverage**: Comprehensive context extraction\n\n")
        
        f.write("## Topics Covered\n\n")
        for section in merged_sections:
            f.write(f"- **{section['topic']}**: Lines {section['start']+1}-{section['end']} ({len(section['lines'])} lines)\n")
        
        f.write("\n---\n\n")
        f.write("**END OF COMPREHENSIVE KNOWLEDGE BASE**\n")
    
    # Get file stats
    file_size_mb = output_file.stat().st_size / (1024 * 1024)
    line_count = sum(len(s['lines']) for s in merged_sections)
    
    print(f"\n✅ Comprehensive knowledge base created: {output_file}")
    print(f"📊 Total sections: {len(merged_sections)}")
    print(f"📏 Total lines extracted: {line_count}")
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
    output_file = output_dir / "COMPREHENSIVE_KNOWLEDGE_BASE.md"
    
    # Extract comprehensive knowledge
    result = extract_comprehensive_knowledge(chat_log_path, output_file)
    
    if result:
        print(f"\n✅ Comprehensive extraction complete!")
        print(f"📁 Output: {result}")
    else:
        print("\n❌ Extraction failed")

if __name__ == "__main__":
    main()



