#!/usr/bin/env python3
"""
Extract PFI data from Science Advances supplementary PDF
"""

import sys
from pathlib import Path

try:
    import pdfplumber
    HAS_PDFPLUMBER = True
except ImportError:
    HAS_PDFPLUMBER = False
    print("⚠️  pdfplumber not installed. Installing...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pdfplumber", "-q"])
    import pdfplumber
    HAS_PDFPLUMBER = True

PDF_PATH = Path.home() / "Downloads" / "sciadv.abm1831_sm.pdf"
OUTPUT_PATH = Path("data/serial_sae/gse165897/results/pfi_data_extracted.json")

# GSE165897 patient IDs (confirmed)
GSE165897_PATIENTS = [
    'EOC1005', 'EOC136', 'EOC153', 'EOC227', 'EOC349', 
    'EOC372', 'EOC3', 'EOC443', 'EOC540', 'EOC733', 'EOC87'
]

def extract_pfi_data():
    """Extract PFI data from PDF tables."""
    print(f"📄 Extracting PFI data from: {PDF_PATH}")
    
    if not PDF_PATH.exists():
        print(f"❌ PDF not found: {PDF_PATH}")
        return None
    
    pfi_data = {}
    
    with pdfplumber.open(PDF_PATH) as pdf:
        print(f"✅ Opened PDF ({len(pdf.pages)} pages)")
        
        # Search through all pages for tables
        for page_num, page in enumerate(pdf.pages, 1):
            tables = page.extract_tables()
            
            if not tables:
                continue
            
            print(f"📊 Page {page_num}: Found {len(tables)} table(s)")
            
            for table_idx, table in enumerate(tables):
                if not table or len(table) < 2:
                    continue
                
                # Look for patient ID patterns and PFI columns
                header = table[0] if table else []
                
                # Check if this looks like a patient data table
                header_str = ' '.join([str(c) for c in header if c]).lower()
                
                if 'patient' in header_str or 'eoc' in header_str.lower():
                    print(f"  Table {table_idx + 1}: Potential patient table")
                    print(f"    Header: {header[:5]}...")
                    
                    # Look for PFI column
                    pfi_col_idx = None
                    patient_col_idx = None
                    
                    for idx, col in enumerate(header):
                        if col:
                            col_lower = str(col).lower()
                            if 'pfi' in col_lower or 'platinum' in col_lower:
                                pfi_col_idx = idx
                            if 'patient' in col_lower or 'id' in col_lower:
                                patient_col_idx = idx
                    
                    if pfi_col_idx is not None:
                        print(f"    Found PFI column at index {pfi_col_idx}")
                        
                        # Extract data rows
                        for row_idx, row in enumerate(table[1:], 1):
                            if not row or len(row) <= max(pfi_col_idx, patient_col_idx or 0):
                                continue
                            
                            # Try to find patient ID
                            patient_id = None
                            if patient_col_idx is not None:
                                patient_id = str(row[patient_col_idx]).strip() if row[patient_col_idx] else None
                            else:
                                # Try first column
                                patient_id = str(row[0]).strip() if row[0] else None
                            
                            # Check if it matches GSE165897 patients
                            if patient_id and any(pid in patient_id for pid in GSE165897_PATIENTS):
                                pfi_val = row[pfi_col_idx] if pfi_col_idx < len(row) else None
                                
                                # Extract numeric PFI value
                                pfi_months = None
                                if pfi_val:
                                    pfi_str = str(pfi_val).strip()
                                    # Try to extract number
                                    import re
                                    numbers = re.findall(r'\d+\.?\d*', pfi_str)
                                    if numbers:
                                        pfi_months = float(numbers[0])
                                
                                if pfi_months is not None:
                                    # Find matching patient ID
                                    for pid in GSE165897_PATIENTS:
                                        if pid in patient_id:
                                            pfi_data[pid] = {
                                                'pfi_months': pfi_months,
                                                'pfi_days': pfi_months * 30.44,  # Approximate
                                                'resistant': 'R' if pfi_months < 6 else 'S',
                                                'source': f'Page {page_num}, Table {table_idx + 1}'
                                            }
                                            print(f"    ✅ {pid}: {pfi_months} months ({'R' if pfi_months < 6 else 'S'})")
                                            break
    
    return pfi_data

def main():
    print("="*80)
    print("PFI DATA EXTRACTION FROM PDF")
    print("="*80)
    print()
    
    pfi_data = extract_pfi_data()
    
    if not pfi_data:
        print("❌ No PFI data extracted. Trying alternative extraction method...")
        # Try text extraction and pattern matching
        with pdfplumber.open(PDF_PATH) as pdf:
            print("\n📝 Extracting all text and searching for patient IDs...")
            full_text = ""
            for page in pdf.pages:
                full_text += page.extract_text() + "\n"
            
            # Search for patient IDs with PFI values
            import re
            for patient_id in GSE165897_PATIENTS:
                # Look for patterns like "EOC1005" followed by numbers
                pattern = rf'{patient_id}.*?(\d+\.?\d*)\s*(?:months?|mo|days?|d)'
                matches = re.findall(pattern, full_text, re.IGNORECASE)
                if matches:
                    pfi_val = float(matches[0])
                    if pfi_val < 100:  # Likely months
                        pfi_data[patient_id] = {
                            'pfi_months': pfi_val,
                            'pfi_days': pfi_val * 30.44,
                            'resistant': 'R' if pfi_val < 6 else 'S',
                            'source': 'Text extraction'
                        }
                        print(f"  ✅ {patient_id}: {pfi_val} months")
    
    # Save results
    if pfi_data:
        import json
        OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
        
        output = {
            '_comment': 'PFI data extracted from sciadv.abm1831_sm.pdf',
            '_total_patients': len(pfi_data),
            '_gse165897_patients': GSE165897_PATIENTS,
            **{pid: pfi_data.get(pid, None) for pid in GSE165897_PATIENTS}
        }
        
        with open(OUTPUT_PATH, 'w') as f:
            json.dump(output, f, indent=2)
        
        print(f"\n💾 Saved PFI data: {OUTPUT_PATH}")
        print(f"✅ Extracted PFI for {len(pfi_data)}/{len(GSE165897_PATIENTS)} patients")
        
        # Show summary
        resistant = sum(1 for v in pfi_data.values() if v.get('resistant') == 'R')
        sensitive = sum(1 for v in pfi_data.values() if v.get('resistant') == 'S')
        print(f"   Resistant: {resistant}, Sensitive: {sensitive}")
    else:
        print("❌ No PFI data found. Please check the PDF manually.")
        print(f"   PDF location: {PDF_PATH}")

if __name__ == "__main__":
    main()
