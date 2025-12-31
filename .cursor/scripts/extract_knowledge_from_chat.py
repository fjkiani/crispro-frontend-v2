#!/usr/bin/env python3
"""
Knowledge Extraction Script - Organize 30,000+ line chat log for agent training

Purpose: Extract key doctrines, frameworks, and patterns from chat log
Output: Organized .mdc files in hierarchical structure
"""

import re
import os
from pathlib import Path
from collections import defaultdict
from typing import List, Tuple, Dict
from dataclasses import dataclass

@dataclass
class KnowledgeSection:
    """Represents an extracted knowledge section"""
    title: str
    category: str  # 'doctrine', 'framework', 'pattern', 'context'
    priority: int  # P0, P1, P2, P3
    line_start: int
    line_end: int
    content: List[str]
    related_sections: List[str] = None

class KnowledgeExtractor:
    """Extracts and organizes knowledge from chat log"""
    
    def __init__(self, chat_log_path: str):
        self.chat_log_path = chat_log_path
        self.lines = []
        self.sections = []
        
    def load_chat_log(self):
        """Load chat log file"""
        with open(self.chat_log_path, 'r', encoding='utf-8', errors='ignore') as f:
            self.lines = f.readlines()
        print(f"Loaded {len(self.lines)} lines from chat log")
    
    def extract_pharma_suppression_doctrine(self) -> KnowledgeSection:
        """Extract pharma suppression doctrine (P0)"""
        # Find key sections about pharma suppression
        key_lines = []
        
        for i, line in enumerate(self.lines):
            line_lower = line.lower()
            if any(term in line_lower for term in [
                'pharma.*bury', 'kelim.*buried', 'cpt.*code', 
                'suppress.*kelim', '10b.*threat', 'evidence fortress',
                'factory.*strategy', 'outproduce.*pharma'
            ]):
                key_lines.append(i)
        
        if not key_lines:
            return None
        
        # Extract context windows around key lines
        context_start = max(0, min(key_lines) - 100)
        context_end = min(len(self.lines), max(key_lines) + 100)
        
        content = self.lines[context_start:context_end]
        
        return KnowledgeSection(
            title="Pharma Suppression Doctrine",
            category="doctrine",
            priority=0,
            line_start=context_start + 1,
            line_end=context_end,
            content=content,
            related_sections=["Factory Strategy", "Mission Discipline"]
        )
    
    def extract_mission_reminder(self) -> KnowledgeSection:
        """Extract critical mission reminder (P0)"""
        # Find the exact reminder line (26773)
        reminder_line = None
        
        for i, line in enumerate(self.lines):
            if 'recal back' in line.lower() and 'pharma' in line.lower():
                reminder_line = i
                break
        
        if reminder_line is None:
            return None
        
        # Extract context window
        context_start = max(0, reminder_line - 50)
        context_end = min(len(self.lines), reminder_line + 100)
        
        content = self.lines[context_start:context_end]
        
        return KnowledgeSection(
            title="Mission Discipline Reminder",
            category="doctrine",
            priority=0,
            line_start=context_start + 1,
            line_end=context_end,
            content=content,
            related_sections=["Pharma Suppression Doctrine"]
        )
    
    def extract_factory_strategy(self) -> KnowledgeSection:
        """Extract Factory strategy (P1)"""
        # Find Factory strategy sections
        key_lines = []
        
        for i, line in enumerate(self.lines):
            line_lower = line.lower()
            if any(term in line_lower for term in [
                'factory.*strategy', 'publication.*factory', 
                'evidence fortress', 'outproduce.*pharma',
                '10.*validated.*papers'
            ]):
                key_lines.append(i)
        
        if not key_lines:
            return None
        
        # Extract context windows
        context_start = max(0, min(key_lines) - 50)
        context_end = min(len(self.lines), max(key_lines) + 50)
        
        content = self.lines[context_start:context_end]
        
        return KnowledgeSection(
            title="Factory Publication Strategy",
            category="framework",
            priority=1,
            line_start=context_start + 1,
            line_end=context_end,
            content=content,
            related_sections=["Pharma Suppression Doctrine"]
        )
    
    def extract_resistance_framework(self) -> KnowledgeSection:
        """Extract Resistance Prediction Framework (P0)"""
        # Find 4-layer framework sections
        key_lines = []
        
        for i, line in enumerate(self.lines):
            line_lower = line.lower()
            if any(term in line_lower for term in [
                'resistance.*framework', '4.*layer', 
                'baseline.*biology', 'genetic.*escape',
                'adaptive.*resistance', 'drug.*clearance'
            ]):
                key_lines.append(i)
        
        if not key_lines:
            return None
        
        # Extract context windows
        context_start = max(0, min(key_lines) - 50)
        context_end = min(len(self.lines), max(key_lines) + 200)  # Framework might be longer
        
        content = self.lines[context_start:context_end]
        
        return KnowledgeSection(
            title="Resistance Prediction Framework",
            category="doctrine",
            priority=0,
            line_start=context_start + 1,
            line_end=context_end,
            content=content,
            related_sections=[]
        )
    
    def extract_all_sections(self) -> List[KnowledgeSection]:
        """Extract all knowledge sections"""
        sections = []
        
        # P0: Critical doctrines
        pharma_doctrine = self.extract_pharma_suppression_doctrine()
        if pharma_doctrine:
            sections.append(pharma_doctrine)
        
        mission_reminder = self.extract_mission_reminder()
        if mission_reminder:
            sections.append(mission_reminder)
        
        resistance_framework = self.extract_resistance_framework()
        if resistance_framework:
            sections.append(resistance_framework)
        
        # P1: Strategic frameworks
        factory_strategy = self.extract_factory_strategy()
        if factory_strategy:
            sections.append(factory_strategy)
        
        return sections
    
    def write_section_to_file(self, section: KnowledgeSection, output_dir: Path):
        """Write extracted section to .mdc file"""
        # Determine subdirectory based on category
        if section.category == "doctrine":
            subdir = output_dir / "CORE_DOCTRINES"
        elif section.category == "framework":
            subdir = output_dir / "STRATEGIC_FRAMEWORKS"
        elif section.category == "pattern":
            subdir = output_dir / "TECHNICAL_PATTERNS"
        else:
            subdir = output_dir / "HISTORICAL_CONTEXT"
        
        subdir.mkdir(parents=True, exist_ok=True)
        
        # Create filename from title
        filename = section.title.lower().replace(" ", "_") + ".mdc"
        filepath = subdir / filename
        
        # Write content
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(f"# {section.title}\n\n")
            f.write(f"**Category**: {section.category}\n")
            f.write(f"**Priority**: P{section.priority}\n")
            f.write(f"**Source**: Chat log lines {section.line_start}-{section.line_end}\n")
            if section.related_sections:
                f.write(f"**Related**: {', '.join(section.related_sections)}\n")
            f.write("\n---\n\n")
            f.write("".join(section.content))
        
        print(f"✅ Wrote {filepath}")
        return filepath
    
    def create_master_index(self, sections: List[KnowledgeSection], output_dir: Path):
        """Create master index file"""
        index_path = output_dir / "KNOWLEDGE_INDEX.md"
        
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
                f.write(f"## P{priority} Sections\n\n")
                for section in by_priority[priority]:
                    f.write(f"- **{section.title}** ({section.category})\n")
                    f.write(f"  - Lines: {section.line_start}-{section.line_end}\n")
                    if section.related_sections:
                        f.write(f"  - Related: {', '.join(section.related_sections)}\n")
                    f.write(f"  - File: `{section.category}/{section.title.lower().replace(' ', '_')}.mdc`\n\n")
        
        print(f"✅ Created master index: {index_path}")

def main():
    """Main extraction function"""
    # Find chat log file
    import glob
    chat_log_files = glob.glob("/Users/fahadkiani/Downloads/*Nyx*rogue*.md")
    
    if not chat_log_files:
        print("❌ Could not find chat log file")
        return
    
    chat_log_path = chat_log_files[0]
    print(f"📖 Using chat log: {chat_log_path}")
    
    # Initialize extractor
    extractor = KnowledgeExtractor(chat_log_path)
    extractor.load_chat_log()
    
    # Extract all sections
    print("\n🔍 Extracting knowledge sections...")
    sections = extractor.extract_all_sections()
    print(f"✅ Extracted {len(sections)} sections")
    
    # Create output directory
    output_dir = Path(".cursor/rules/MOAT/EXTRACTED_KNOWLEDGE")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Write sections to files
    print("\n📝 Writing sections to files...")
    for section in sections:
        extractor.write_section_to_file(section, output_dir)
    
    # Create master index
    print("\n📋 Creating master index...")
    extractor.create_master_index(sections, output_dir)
    
    print("\n✅ Knowledge extraction complete!")
    print(f"📁 Output directory: {output_dir}")

if __name__ == "__main__":
    main()


