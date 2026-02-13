
import csv

series_matrix_path = 'oncology-coPilot/oncology-backend-minimal/data/serial_sae/gse241908/series_matrix_temp.txt'

def parse_metadata():
    title_line = []
    geo_line = []
    characteristics_line = []
    description_line = []
    
    with open(series_matrix_path, 'r') as f:
        for line in f:
            if line.startswith('!Sample_title'):
                title_line = [x.strip('"').strip() for x in line.strip().split('\t')[1:]]
            if line.startswith('!Sample_geo_accession'):
                geo_line = [x.strip('"').strip() for x in line.strip().split('\t')[1:]]
            if line.startswith('!Sample_characteristics_ch1'):
                # Heuristic to find the meaningful characteristic line
                if "chemotherapy" in line or "time" in line:
                    characteristics_line = [x.strip('"').strip() for x in line.strip().split('\t')[1:]]
            if line.startswith('!Sample_description'):  # Might contain ShV
                 description_line = [x.strip('"').strip() for x in line.strip().split('\t')[1:]]
            
            if title_line and geo_line and characteristics_line and description_line:
                break
    
    # Check if ShV is in title or description
    print(f"{'Index':<5} | {'GSM':<10} | {'Title':<40} | {'ShV Match?':<15} | {'Characteristic'}")
    print("-" * 120)
    
    # If description is empty, handle it
    if not description_line:
        description_line = [""] * len(title_line)
        
    for i, (t, g, c, d) in enumerate(zip(title_line, geo_line, characteristics_line, description_line)):
        # Try to find ShV in title or description
        shv = "???"
        if "ShV" in t:
             shv = t
        elif "ShV" in d:
             shv = d
        # Or maybe it matches the order?
        print(f"{i:<5} | {g:<10} | {t:<40} | {shv:<15} | {c}")

if __name__ == "__main__":
    parse_metadata()
