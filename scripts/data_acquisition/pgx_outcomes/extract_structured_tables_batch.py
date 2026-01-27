import os
import re
import json
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from dotenv import load_dotenv
from Bio import Entrez

REPO_ROOT = Path(__file__).resolve().parents[3]
BIOMED_ENV = REPO_ROOT / 'scripts' / 'data_acquisition' / 'mcp_servers' / 'BioMed-MCP' / '.env'
OUT_DIR = REPO_ROOT / 'scripts' / 'data_acquisition' / 'pgx_outcomes' / 'receipts' / 'pmc_structured'
OUT_DIR.mkdir(parents=True, exist_ok=True)


def setup_entrez() -> None:
    load_dotenv(BIOMED_ENV)
    email = os.getenv('PUBMED_EMAIL')
    if not email:
        raise RuntimeError('PUBMED_EMAIL not set')
    Entrez.email = email
    Entrez.tool = 'crispr-assistant-pmc-structured-batch'


def pmc_link(pmid: str) -> Optional[str]:
    h = Entrez.elink(dbfrom='pubmed', db='pmc', id=pmid)
    xml = h.read(); h.close()
    if isinstance(xml, (bytes, bytearray)):
        xml = xml.decode('utf-8', errors='replace')
    m = re.search(r'<LinkSetDb>.*?<Id>(\d+)</Id>.*?</LinkSetDb>', xml, re.S)
    return m.group(1) if m else None


def fetch_pmc_xml(pmc_id: str) -> str:
    h = Entrez.efetch(db='pmc', id=pmc_id, rettype='xml')
    xml = h.read(); h.close()
    if isinstance(xml, (bytes, bytearray)):
        xml = xml.decode('utf-8', errors='replace')
    return xml


def _clean_xml_for_et(xml: str) -> str:
    cleaned = re.sub(r'xmlns(:\w+)?="[^"]+"', '', xml)
    cleaned = re.sub(r'<(/?)\w+:', r'<\1', cleaned)
    cleaned = re.sub(r'(\s)\w+:(\w+)=', r'\1\2=', cleaned)
    return cleaned


def _cell_text(el: ET.Element) -> str:
    parts = []
    for t in el.itertext():
        t = (t or '').strip()
        if t:
            parts.append(t)
    return re.sub(r'\s+', ' ', ' '.join(parts)).strip()


def extract_table(root: ET.Element, want_label_prefix: str) -> Dict[str, Any]:
    for tw in root.findall('.//table-wrap'):
        label_el = tw.find('label')
        label = _cell_text(label_el) if label_el is not None else ''
        if label.strip().lower().startswith(want_label_prefix.lower()):
            cap = tw.find('caption')
            title_el = cap.find('title') if cap is not None else None
            title = _cell_text(title_el) if title_el is not None else (_cell_text(cap) if cap is not None else None)

            rows = []
            table_el = tw.find('.//table')
            if table_el is not None:
                for tr in table_el.findall('.//tr'):
                    cells = []
                    for td in list(tr):
                        if (td.tag or '').lower() in ('td','th'):
                            cells.append(_cell_text(td))
                    if any(cells):
                        rows.append(cells)

            return {'label': label, 'title': title, 'n_rows': len(rows), 'rows': rows}

    raise RuntimeError(f'No table found with label starting {want_label_prefix!r}')


def extract(pmid: str, label_prefixes: List[str]) -> Dict[str, Any]:
    pmc_id = pmc_link(pmid)
    if not pmc_id:
        raise RuntimeError(f'No PMC id for {pmid}')
    xml = fetch_pmc_xml(pmc_id)
    root = ET.fromstring(_clean_xml_for_et(xml))

    tables = {}
    for pref in label_prefixes:
        tables[pref] = extract_table(root, pref)

    out = {
        'generated_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'pmid': pmid,
        'pmc_id': pmc_id,
        'tables': tables,
    }

    out_path = OUT_DIR / f'pmid_{pmid}_tables_{"_".join([p.replace(" ","") for p in label_prefixes])}.json'
    out_path.write_text(json.dumps(out, indent=2))
    return out


def main():
    setup_entrez()

    # Two key breakthroughs:
    # 1) PREPARE secondary analysis has true negative controls for toxicity endpoints.
    # 2) CYP2C19 clopidogrel outcomes cohort has intermediate metabolizer outcomes tables.
    tasks = [
        ('39641926', ['Table 1', 'Table 2']),
        ('40944685', ['Table 2', 'Table 4']),
    ]

    for pmid, prefs in tasks:
        print('extract', pmid, prefs)
        out = extract(pmid, prefs)
        for k, t in out['tables'].items():
            print(' ', k, '->', t.get('label'), '|', (t.get('title') or '')[:80], '| rows', t.get('n_rows'))
        time.sleep(0.4)


if __name__ == '__main__':
    main()
