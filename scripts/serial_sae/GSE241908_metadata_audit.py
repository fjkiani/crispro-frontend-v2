
import csv
import re

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
                if "chemotherapy" in line or "time" in line:
                    characteristics_line = [x.strip('"').strip() for x in line.strip().split('\t')[1:]]
            if line.startswith('!Sample_description'):  
                 description_line = [x.strip('"').strip() for x in line.strip().split('\t')[1:]]
            
            if title_line and geo_line and characteristics_line and description_line:
                break
    
    print(f"{'Index':<5} | {'GSM':<10} | {'ShV (Extracted)':<15} | {'Title':<40}")
    print("-" * 100)
    
    if not description_line:
         description_line = [""] * len(title_line)

    for i, (t, g, c, d) in enumerate(zip(title_line, geo_line, characteristics_line, description_line)):
        # Extract ShV using regex from description or title
        match = re.search(r"ShV[-\s]?(\d+)", d)
        if not match:
             match = re.search(r"ShV[-\s]?(\d+)", t)
        
        shv = f"ShV-{match.group(1)}" if match else "???"
        
        print(f"{i:<5} | {g:<10} | {shv:<15} | {t:<40}")

if __name__ == "__main__":
    parse_metadata()
