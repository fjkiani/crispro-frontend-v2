#!/usr/bin/env python3
"""
Enhanced Knowledge Extraction Script - Extract all major knowledge domains

Purpose: Comprehensive extraction of doctrines, frameworks, patterns, and context
Output: Organized .mdc files in hierarchical structure
"""

import re
import os
import glob
from pathlib import Path
from collections import defaultdict
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass, field

@dataclass
class KnowledgeSection:
    """Represents an extracted knowledge section"""
    title: str
    category: str  # 'doctrine', 'framework', 'pattern', 'context'
    priority: int  # P0, P1, P2, P3
    line_start: int
    line_end: int
    content: List[str]
    related_sections: List[str] = field(default_factory=list)
    key_quotes: List[str] = field(default_factory=list)

class EnhancedKnowledgeExtractor:
    """Enhanced extractor for comprehensive knowledge extraction"""
    
    def __init__(self, chat_log_path: str):
        self.chat_log_path = chat_log_path
        self.lines = []
        self.sections = []
        self.output_dir = Path(".cursor/rules/MOAT/AGENT_TRAINING_KNOWLEDGE")
        
    def load_chat_log(self):
        """Load chat log file"""
        with open(self.chat_log_path, 'r', encoding='utf-8', errors='ignore') as f:
            self.lines = f.readlines()
        print(f"✅ Loaded {len(self.lines)} lines from chat log")
    
    def find_section_by_patterns(self, patterns: List[str], context_before: int = 50, 
                                context_after: int = 100, min_matches: int = 2) -> Optional[Tuple[int, int]]:
        """Find section by multiple patterns"""
        matching_lines = []
        
        for i, line in enumerate(self.lines):
            line_lower = line.lower()
            matches = sum(1 for pattern in patterns if re.search(pattern, line_lower, re.IGNORECASE))
            if matches >= min_matches:
                matching_lines.append(i)
        
        if not matching_lines:
            return None
        
        # Get context window around all matches
        start = max(0, min(matching_lines) - context_before)
        end = min(len(self.lines), max(matching_lines) + context_after)
        
        return (start, end)
    
    def extract_resistance_framework(self) -> Optional[KnowledgeSection]:
        """Extract 4-layer Resistance Prediction Framework"""
        patterns = [
            r'4.*layer.*resistance',
            r'resistance.*prediction.*framework',
            r'baseline.*biology',
            r'genetic.*escape',
            r'additional.*resistance',
            r'drug.*clearance',
            r'resistancepredictionengine',
            r'layer\s*1.*layer\s*2.*layer\s*3.*layer\s*4'
        ]
        
        result = self.find_section_by_patterns(patterns, context_before=100, context_after=200, min_matches=3)
        if not result:
            return None
        
        start, end = result
        content = self.lines[start:end]
        
        # Extract key quotes
        key_quotes = []
        for line in content:
            if any(term in line.lower() for term in ['4-layer', '4 layer', 'resistance framework', 'baseline biology']):
                if len(line.strip()) > 20:  # Meaningful line
                    key_quotes.append(line.strip()[:200])  # First 200 chars
        
        return KnowledgeSection(
            title="Resistance Prediction Framework",
            category="doctrine",
            priority=0,
            line_start=start + 1,
            line_end=end,
            content=content,
            related_sections=["Pharma Suppression Doctrine", "Mission Discipline Doctrine"],
            key_quotes=key_quotes[:5]  # Top 5 quotes
        )
    
    def extract_code_blocks(self) -> List[KnowledgeSection]:
        """Extract code blocks (Python, JavaScript, etc.)"""
        code_sections = []
        in_code_block = False
        code_start = None
        code_language = None
        code_lines = []
        
        for i, line in enumerate(self.lines):
            # Detect code block start
            if re.match(r'```(\w+)', line):
                if in_code_block:
                    # End previous block
                    if code_start is not None and len(code_lines) > 5:  # At least 5 lines
                        code_sections.append({
                            'start': code_start,
                            'end': i,
                            'language': code_language,
                            'lines': code_lines
                        })
                
                # Start new block
                in_code_block = True
                code_start = i
                code_language = re.match(r'```(\w+)', line).group(1)
                code_lines = []
            elif line.strip() == '```' and in_code_block:
                # End code block
                if code_start is not None and len(code_lines) > 5:
                    code_sections.append({
                        'start': code_start,
                        'end': i,
                        'language': code_language,
                        'lines': code_lines
                    })
                in_code_block = False
                code_start = None
                code_lines = []
            elif in_code_block:
                code_lines.append(line)
        
        # Convert to KnowledgeSection objects
        sections = []
        for idx, code_block in enumerate(code_sections[:20]):  # Limit to top 20
            # Get context before code block
            context_start = max(0, code_block['start'] - 30)
            context_end = min(len(self.lines), code_block['end'] + 10)
            
            # Combine context + code
            full_content = self.lines[context_start:code_block['start']] + \
                          code_block['lines'] + \
                          self.lines[code_block['end']:context_end]
            
            sections.append(KnowledgeSection(
                title=f"Code Example {idx+1} ({code_block['language']})",
                category="pattern",
                priority=2,
                line_start=context_start + 1,
                line_end=context_end,
                content=full_content,
                related_sections=[]
            ))
        
        return sections
    
    def extract_strategic_decisions(self) -> List[KnowledgeSection]:
        """Extract strategic decision points"""
        decision_patterns = [
            r'we\s+decided',
            r'the\s+strategy\s+was',
            r'the\s+approach',
            r'decision\s+was',
            r'we\s+chose',
            r'strategic\s+decision'
        ]
        
        decision_lines = []
        for i, line in enumerate(self.lines):
            line_lower = line.lower()
            if any(re.search(pattern, line_lower) for pattern in decision_patterns):
                decision_lines.append(i)
        
        # Group nearby decisions
        sections = []
        if decision_lines:
            # Get context around first few decisions
            for i, line_num in enumerate(decision_lines[:10]):  # Top 10 decisions
                context_start = max(0, line_num - 50)
                context_end = min(len(self.lines), line_num + 50)
                
                sections.append(KnowledgeSection(
                    title=f"Strategic Decision {i+1}",
                    category="context",
                    priority=3,
                    line_start=context_start + 1,
                    line_end=context_end,
                    content=self.lines[context_start:context_end],
                    related_sections=[]
                ))
        
        return sections
    
    def extract_corrections_reminders(self) -> List[KnowledgeSection]:
        """Extract correction and reminder moments"""
        correction_patterns = [
            r'remind',
            r'correct',
            r'don.*t.*deviate',
            r'stay.*on.*mission',
            r'recal.*back',
            r'remember',
            r'focus'
        ]
        
        correction_lines = []
        for i, line in enumerate(self.lines):
            line_lower = line.lower()
            if any(re.search(pattern, line_lower) for pattern in correction_patterns):
                # Check if nearby context mentions pharma/mission/factory
                context_window = ' '.join([l.lower() for l in self.lines[max(0, i-10):min(len(self.lines), i+10)]])
                if any(term in context_window for term in ['pharma', 'mission', 'factory', 'kelim', 'strategy']):
                    correction_lines.append(i)
        
        sections = []
        for i, line_num in enumerate(correction_lines[:15]):  # Top 15 corrections
            context_start = max(0, line_num - 30)
            context_end = min(len(self.lines), line_num + 50)
            
            sections.append(KnowledgeSection(
                title=f"Correction/Reminder {i+1}",
                category="context",
                priority=3,
                line_start=context_start + 1,
                line_end=context_end,
                content=self.lines[context_start:context_end],
                related_sections=["Mission Discipline Doctrine"]
            ))
        
        return sections
    
    def extract_all_sections(self) -> List[KnowledgeSection]:
        """Extract all knowledge sections"""
        sections = []
        
        print("\n🔍 Extracting knowledge sections...")
        
        # P0: Critical doctrines (already extracted, but re-extract to ensure completeness)
        print("  📋 Extracting Resistance Framework...")
        resistance = self.extract_resistance_framework()
        if resistance:
            sections.append(resistance)
            print(f"    ✅ Extracted: {resistance.title} (lines {resistance.line_start}-{resistance.line_end})")
        
        # P2: Technical patterns
        print("  💻 Extracting code blocks...")
        code_blocks = self.extract_code_blocks()
        sections.extend(code_blocks)
        print(f"    ✅ Extracted: {len(code_blocks)} code blocks")
        
        # P3: Historical context
        print("  📜 Extracting strategic decisions...")
        decisions = self.extract_strategic_decisions()
        sections.extend(decisions)
        print(f"    ✅ Extracted: {len(decisions)} strategic decisions")
        
        print("  🔔 Extracting corrections/reminders...")
        corrections = self.extract_corrections_reminders()
        sections.extend(corrections)
        print(f"    ✅ Extracted: {len(corrections)} corrections/reminders")
        
        return sections
    
    def write_section_to_file(self, section: KnowledgeSection):
        """Write extracted section to .mdc file"""
        # Determine subdirectory based on category
        if section.category == "doctrine":
            subdir = self.output_dir / "CORE_DOCTRINES"
        elif section.category == "framework":
            subdir = self.output_dir / "STRATEGIC_FRAMEWORKS"
        elif section.category == "pattern":
            subdir = self.output_dir / "TECHNICAL_PATTERNS"
        else:
            subdir = self.output_dir / "HISTORICAL_CONTEXT"
        
        subdir.mkdir(parents=True, exist_ok=True)
        
        # Create filename from title
        filename = section.title.lower().replace(" ", "_").replace("/", "_") + ".mdc"
        # Remove special characters
        filename = re.sub(r'[^\w\-_\.]', '', filename)
        filepath = subdir / filename
        
        # Write content
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(f"# {section.title}\n\n")
            f.write(f"**Category**: {section.category}\n")
            f.write(f"**Priority**: P{section.priority}\n")
            f.write(f"**Source**: Chat log lines {section.line_start}-{section.line_end}\n")
            if section.related_sections:
                f.write(f"**Related**: {', '.join(section.related_sections)}\n")
            if section.key_quotes:
                f.write(f"\n**Key Quotes**:\n")
                for quote in section.key_quotes:
                    f.write(f"- {quote}\n")
            f.write("\n---\n\n")
            f.write("".join(section.content))
        
        return filepath
    
    def create_master_index(self, sections: List[KnowledgeSection]):
        """Create master index file"""
        index_path = self.output_dir / "KNOWLEDGE_INDEX.md"
        
        with open(index_path, 'w', encoding='utf-8') as f:
            f.write("# 🧠 KNOWLEDGE BASE INDEX\n\n")
            f.write("**Purpose**: Master index of all extracted knowledge\n")
            f.write(f"**Total Sections**: {len(sections)}\n\n")
            f.write("---\n\n")
            
            # Organize by priority
            by_priority = defaultdict(list)
            for section in sections:
                by_priority[section.priority].append(section)
            
            for priority in sorted(by_priority.keys()):
                f.write(f"## P{priority} Sections ({len(by_priority[priority])} sections)\n\n")
                for section in by_priority[priority]:
                    f.write(f"- **{section.title}** ({section.category})\n")
                    f.write(f"  - Lines: {section.line_start}-{section.line_end}\n")
                    if section.related_sections:
                        f.write(f"  - Related: {', '.join(section.related_sections)}\n")
                    f.write(f"  - File: `{section.category.upper()}/{section.title.lower().replace(' ', '_')}.mdc`\n\n")
        
        print(f"✅ Created master index: {index_path}")

def main():
    """Main extraction function"""
    # Find chat log file
    chat_log_files = glob.glob("/Users/fahadkiani/Downloads/*Nyx*rogue*.md")
    
    if not chat_log_files:
        print("❌ Could not find chat log file")
        return
    
    chat_log_path = chat_log_files[0]
    print(f"📖 Using chat log: {chat_log_path}")
    
    # Initialize extractor
    extractor = EnhancedKnowledgeExtractor(chat_log_path)
    extractor.load_chat_log()
    
    # Extract all sections
    sections = extractor.extract_all_sections()
    print(f"\n✅ Extracted {len(sections)} total sections")
    
    # Write sections to files
    print("\n📝 Writing sections to files...")
    written_files = []
    for section in sections:
        filepath = extractor.write_section_to_file(section)
        written_files.append(filepath)
        print(f"  ✅ Wrote: {filepath.name}")
    
    # Create master index
    print("\n📋 Creating master index...")
    extractor.create_master_index(sections)
    
    print(f"\n✅ Knowledge extraction complete!")
    print(f"📁 Output directory: {extractor.output_dir}")
    print(f"📊 Total files created: {len(written_files)}")

if __name__ == "__main__":
    main()











