import os
import re
import json
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Optional

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
    Entrez.tool = 'crispr-assistant-pmc-structured-table'


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
    # remove xmlns and namespace prefixes in tags + attributes
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
    s = ' '.join(parts)
    s = re.sub(r'\s+', ' ', s).strip()
    return s


def extract_table(root: ET.Element, want_label_prefix: str) -> Dict[str, Any]:
    for tw in root.findall('.//table-wrap'):
        label_el = tw.find('label')
        label = _cell_text(label_el) if label_el is not None else ''
        if label.strip().lower().startswith(want_label_prefix.lower()):
            title_el = None
            cap = tw.find('caption')
            if cap is not None:
                title_el = cap.find('title')
            title = _cell_text(title_el) if title_el is not None else (_cell_text(cap) if cap is not None else None)

            rows = []
            table_el = tw.find('.//table')
            if table_el is not None:
                for tr in table_el.findall('.//tr'):
                    cells = []
                    for td in list(tr):
                        tag = (td.tag or '').lower()
                        if tag not in ('td','th'):
                            continue
                        cells.append(_cell_text(td))
                    if any(cells):
                        rows.append(cells)

            return {
                'label': label,
                'title': title,
                'n_rows': len(rows),
                'rows': rows,
            }

    raise RuntimeError(f'No table found with label starting {want_label_prefix!r}')


def main():
    setup_entrez()

    # PREPARE secondary analysis paper from receipts
    pmid = '39641926'
    want_label_prefix = 'Table 2'

    pmc_id = pmc_link(pmid)
    if not pmc_id:
        raise RuntimeError('No PMC id')

    xml = fetch_pmc_xml(pmc_id)
    cleaned = _clean_xml_for_et(xml)
    root = ET.fromstring(cleaned)

    tbl = extract_table(root, want_label_prefix)

    out = {
        'generated_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'pmid': pmid,
        'pmc_id': pmc_id,
        'requested_label_prefix': want_label_prefix,
        'table': tbl,
    }

    out_path = OUT_DIR / f'pmid_{pmid}_{want_label_prefix.replace(" ", "_")}.json'
    out_path.write_text(json.dumps(out, indent=2))
    print('wrote', out_path)
    print('label:', tbl.get('label'))
    print('title:', (tbl.get('title') or '')[:160])
    print('rows:', tbl.get('n_rows'))
    # print first 8 rows for inspection
    for r in tbl.get('rows', [])[:8]:
        print(r)


if __name__ == '__main__':
    main()
