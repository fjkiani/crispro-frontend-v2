#!/usr/bin/env python3
"""
Targeted Extraction: Resistance Prediction Framework (4-Layer Architecture)

Purpose: Extract the complete Resistance Framework doctrine
Priority: P0 (Critical)
"""

import re
import glob
from pathlib import Path

def extract_resistance_framework(chat_log_path: str, output_dir: Path):
    """Extract Resistance Prediction Framework"""
    
    with open(chat_log_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    
    print(f"📖 Loaded {len(lines)} lines from chat log")
    
    # Find all lines mentioning resistance framework or 4-layer architecture
    framework_lines = []
    
    for i, line in enumerate(lines):
        line_lower = line.lower()
        # Look for framework mentions
        if any(term in line_lower for term in [
            '4-layer', '4 layer', 'four layer',
            'resistance.*framework', 'resistance.*prediction.*framework',
            'baseline.*biology', 'genetic.*escape', 'adaptive.*resistance', 'drug.*clearance',
            'layer 1', 'layer 2', 'layer 3', 'layer 4',
            'resistancepredictionengine', 'resistance.*prediction.*engine'
        ]):
            framework_lines.append(i)
    
    if not framework_lines:
        print("❌ No resistance framework sections found")
        return None
    
    print(f"✅ Found {len(framework_lines)} framework-related lines")
    
    # Find the main framework section (look for comprehensive sections)
    # Start from the first mention and expand context
    start_line = max(0, min(framework_lines) - 100)
    end_line = min(len(lines), max(framework_lines) + 200)
    
    # Try to find a more specific section by looking for "THE RESISTANCE PREDICTION FRAMEWORK" or similar headers
    for i in range(len(lines)):
        line_lower = lines[i].lower()
        if any(header in line_lower for header in [
            'resistance.*prediction.*framework',
            '4.*layer.*resistance',
            '⚔️.*resistance.*prediction.*framework',
            'the.*resistance.*prediction.*framework'
        ]):
            # Found a header, extract from here
            start_line = max(0, i - 50)
            # Look for next major section or end of framework
            for j in range(i + 1, min(len(lines), i + 2000)):
                if re.search(r'^#+\s+[A-Z]', lines[j]) and 'framework' not in lines[j].lower():
                    end_line = j
                    break
            else:
                end_line = min(len(lines), i + 2000)
            break
    
    content = lines[start_line:end_line]
    
    print(f"✅ Extracted framework section: lines {start_line+1}-{end_line}")
    print(f"   Content length: {len(content)} lines")
    
    # Write to file
    output_file = output_dir / "CORE_DOCTRINES" / "resistance_framework_doctrine.mdc"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("# Resistance Prediction Framework (4-Layer Architecture)\n\n")
        f.write("**Category**: doctrine\n")
        f.write("**Priority**: P0 (Critical)\n")
        f.write(f"**Source**: Chat log lines {start_line+1}-{end_line}\n")
        f.write("**Related**: Pharma Suppression Doctrine, Mission Discipline Doctrine\n")
        f.write("\n---\n\n")
        f.write("".join(content))
    
    print(f"✅ Wrote: {output_file}")
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
    
    # Output directory
    output_dir = Path(".cursor/rules/MOAT/AGENT_TRAINING_KNOWLEDGE")
    
    # Extract framework
    result = extract_resistance_framework(chat_log_path, output_dir)
    
    if result:
        print(f"\n✅ Resistance Framework extraction complete!")
        print(f"📁 Output: {result}")
    else:
        print("\n❌ Extraction failed")

if __name__ == "__main__":
    main()










