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
OUT_DIR = REPO_ROOT / 'scripts' / 'data_acquisition' / 'pgx_outcomes' / 'receipts' / 'pmc_tables'
OUT_DIR.mkdir(parents=True, exist_ok=True)


def setup_entrez() -> None:
    load_dotenv(BIOMED_ENV)
    email = os.getenv('PUBMED_EMAIL')
    if not email:
        raise RuntimeError('PUBMED_EMAIL not set')
    Entrez.email = email
    Entrez.tool = 'crispr-assistant-pmc-table-extractor'
    api_key = os.getenv('PUBMED_API_KEY')
    if api_key:
        Entrez.api_key = api_key


def pmc_link(pmid: str) -> Optional[str]:
    h = Entrez.elink(dbfrom='pubmed', db='pmc', id=pmid)
    xml = h.read()
    h.close()
    if isinstance(xml, (bytes, bytearray)):
        xml = xml.decode('utf-8', errors='replace')
    m = re.search(r'<LinkSetDb>.*?<Id>(\d+)</Id>.*?</LinkSetDb>', xml, re.S)
    return m.group(1) if m else None


def fetch_pmc_xml(pmc_id: str) -> str:
    h = Entrez.efetch(db='pmc', id=pmc_id, rettype='xml')
    xml = h.read()
    h.close()
    if isinstance(xml, (bytes, bytearray)):
        xml = xml.decode('utf-8', errors='replace')
    return xml


def text_content(elem: ET.Element) -> str:
    # join all text nodes
    parts: List[str] = []
    for t in elem.itertext():
        t = (t or '').strip()
        if t:
            parts.append(t)
    s = ' '.join(parts)
    s = re.sub(r'\s+', ' ', s).strip()
    return s


def extract_tables(pmc_xml: str) -> List[Dict[str, Any]]:
    # PMC XML can have namespaces; strip them by regex replacement for ElementTree simplicity
    cleaned = re.sub(r'xmlns(:\w+)?="[^"]+"', '', pmc_xml)
    cleaned = re.sub(r'<(/?)\w+:', r'<\1', cleaned)  # drop prefix
    # Drop namespace prefixes in attribute names too (e.g., xlink:href -> href)
    cleaned = re.sub(r'(\s)\w+:(\w+)=', r'\1\2=', cleaned)

    root = ET.fromstring(cleaned)
    tables: List[Dict[str, Any]] = []

    for tw in root.findall('.//table-wrap'):
        caption = tw.find('caption')
        label = tw.find('label')
        title = None
        if caption is not None:
            # often caption has <title>
            t = caption.find('title')
            title = text_content(t) if t is not None else text_content(caption)
        tables.append({
            'label': text_content(label) if label is not None else None,
            'title': title,
            'text': text_content(tw)[:20000],
        })

    return tables


def filter_tables(tables: List[Dict[str, Any]], keywords: List[str]) -> List[Dict[str, Any]]:
    kws = [k.lower() for k in keywords]
    out = []
    for t in tables:
        blob = ((t.get('title') or '') + ' ' + (t.get('text') or '')).lower()
        if any(k in blob for k in kws):
            out.append(t)
    return out


def run(pmid: str, keywords: List[str]) -> Dict[str, Any]:
    pmc_id = pmc_link(pmid)
    if not pmc_id:
        raise RuntimeError(f'No PMC ID for PMID {pmid}')

    xml = fetch_pmc_xml(pmc_id)
    tables = extract_tables(xml)
    hits = filter_tables(tables, keywords)

    out = {
        'pmid': pmid,
        'pmc_id': pmc_id,
        'n_tables': len(tables),
        'keywords': keywords,
        'hit_tables': hits,
        'generated_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
    }

    (OUT_DIR / f'pmid_{pmid}_tables.json').write_text(json.dumps(out, indent=2))
    return out


def main():
    setup_entrez()

    targets = [
        ('39641926', ['dpyd', 'ugt1a1', 'fluoropyrimidine', 'irinotecan', 'toxicity', 'grade']),
        ('40773711', ['dpyd', 'ugt1a1', 'prospective', 'trial', 'toxicity', 'grade']),
        ('39516829', ['dpyd', 'genotype-guided', 'toxicity', 'grade']),
    ]

    for pmid, kws in targets:
        print('extracting tables for', pmid)
        out = run(pmid, kws)
        print('  pmc', out['pmc_id'], 'tables', out['n_tables'], 'hits', len(out['hit_tables']))
        time.sleep(0.34)


if __name__ == '__main__':
    main()
