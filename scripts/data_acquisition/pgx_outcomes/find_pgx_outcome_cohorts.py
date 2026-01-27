import os
import re
import json
import time
from pathlib import Path
from typing import Dict, Any, List, Optional

from dotenv import load_dotenv
from Bio import Entrez

REPO_ROOT = Path(__file__).resolve().parents[3]
BIOMED_ENV = REPO_ROOT / 'scripts' / 'data_acquisition' / 'mcp_servers' / 'BioMed-MCP' / '.env'
OUT_DIR = REPO_ROOT / 'scripts' / 'data_acquisition' / 'pgx_outcomes' / 'receipts'
OUT_DIR.mkdir(parents=True, exist_ok=True)


def _setup_entrez() -> None:
    load_dotenv(BIOMED_ENV)
    email = os.getenv('PUBMED_EMAIL')
    if not email:
        raise RuntimeError('PUBMED_EMAIL not set (expected in BioMed-MCP/.env)')
    Entrez.email = email
    Entrez.tool = 'crispr-assistant-pgx-outcomes'
    api_key = os.getenv('PUBMED_API_KEY')
    if api_key:
        Entrez.api_key = api_key


def esearch(query: str, retmax: int = 25) -> List[str]:
    h = Entrez.esearch(db='pubmed', term=query, retmax=str(retmax))
    xml = h.read()
    h.close()
    if isinstance(xml, (bytes, bytearray)):
        xml = xml.decode('utf-8', errors='replace')
    ids = re.findall(r'<Id>(\d+)</Id>', xml)
    return ids


def efetch_details(pmids: List[str]) -> Dict[str, Dict[str, Any]]:
    # Fetch XML for batch
    if not pmids:
        return {}
    h = Entrez.efetch(db='pubmed', id=','.join(pmids), rettype='xml')
    xml = h.read()
    h.close()
    if isinstance(xml, (bytes, bytearray)):
        xml = xml.decode('utf-8', errors='replace')

    # Extremely lightweight parsing using regex (good enough for discovery)
    # For robust parsing we'd use ElementTree, but discovery is fine.
    out: Dict[str, Dict[str, Any]] = {}
    for pmid in pmids:
        out[pmid] = {'pmid': pmid}

    # Split articles
    articles = xml.split('<PubmedArticle>')
    for block in articles[1:]:
        pmid_m = re.search(r'<PMID[^>]*>(\d+)</PMID>', block)
        if not pmid_m:
            continue
        pmid = pmid_m.group(1)
        title_m = re.search(r'<ArticleTitle>(.*?)</ArticleTitle>', block, re.S)
        abstract_m = re.search(r'<Abstract>(.*?)</Abstract>', block, re.S)
        doi_m = re.search(r'<ArticleId IdType="doi">(.*?)</ArticleId>', block, re.S)
        year_m = re.search(r'<PubDate>.*?<Year>(\d+)</Year>.*?</PubDate>', block, re.S)

        title = _strip_xml(title_m.group(1)) if title_m else None
        abstract = _strip_xml(abstract_m.group(1)) if abstract_m else None
        doi = _strip_xml(doi_m.group(1)) if doi_m else None
        year = year_m.group(1) if year_m else None

        out[pmid].update({'title': title, 'abstract': abstract, 'doi': doi, 'year': year})

    return out


def _strip_xml(s: str) -> str:
    s = re.sub(r'<[^>]+>', ' ', s or '')
    s = re.sub(r'\s+', ' ', s).strip()
    return s


def pmc_link(pmid: str) -> Optional[str]:
    h = Entrez.elink(dbfrom='pubmed', db='pmc', id=pmid)
    xml = h.read()
    h.close()
    if isinstance(xml, (bytes, bytearray)):
        xml = xml.decode('utf-8', errors='replace')
    m = re.search(r'<LinkSetDb>.*?<Id>(\d+)</Id>.*?</LinkSetDb>', xml, re.S)
    return m.group(1) if m else None


def score_candidate(article: Dict[str, Any]) -> Dict[str, Any]:
    text = (article.get('title') or '') + ' ' + (article.get('abstract') or '')
    tl = text.lower()

    flags = {
        'mentions_genotype_guided': 'genotype-guided' in tl or 'genotype guided' in tl,
        'mentions_prospective': 'prospective' in tl,
        'mentions_randomized': 'randomized' in tl or 'randomised' in tl,
        'mentions_toxicity': 'toxicity' in tl or 'adverse' in tl,
        'mentions_grade': 'grade 3' in tl or 'grade iii' in tl or 'ctcae' in tl,
        'mentions_outcome': 'outcome' in tl or 'survival' in tl or 'event' in tl,
        'mentions_genotype': 'genotype' in tl or 'dypd' in tl or 'tpmt' in tl or 'ugt1a1' in tl or 'cyp2c19' in tl or 'cyp2d6' in tl,
        'mentions_clopidogrel': 'clopidogrel' in tl,
        'mentions_fluoropyrimidine': '5-fu' in tl or 'fluoropyrimidine' in tl or 'capecitabine' in tl,
        'mentions_thiopurine': 'thiopurine' in tl or 'azathioprine' in tl or 'mercaptopurine' in tl,
    }

    # crude count of numeric patterns suggesting cohort reporting
    n_like = len(re.findall(r'\bn\s*=\s*\d+\b', tl)) + len(re.findall(r'\b\d+\s*patients\b', tl))

    score = 0
    score += 2 if flags['mentions_genotype_guided'] else 0
    score += 1 if flags['mentions_prospective'] else 0
    score += 1 if flags['mentions_randomized'] else 0
    score += 2 if flags['mentions_toxicity'] else 0
    score += 1 if flags['mentions_grade'] else 0
    score += 1 if flags['mentions_outcome'] else 0
    score += min(2, n_like)

    return {'score': score, 'flags': flags, 'n_like': n_like}


def run_query(name: str, query: str, retmax: int = 25) -> Dict[str, Any]:
    pmids = esearch(query, retmax=retmax)
    details = efetch_details(pmids)

    rows = []
    for pmid in pmids:
        a = details.get(pmid, {'pmid': pmid})
        pmc = pmc_link(pmid)
        scored = score_candidate(a)
        rows.append({
            **a,
            'pmc_id': pmc,
            **scored,
        })
        time.sleep(0.34)

    # sort high-signal first
    rows.sort(key=lambda r: (r.get('pmc_id') is not None, r.get('score', 0)), reverse=True)

    out = {
        'query_name': name,
        'query': query,
        'generated_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'n_results': len(rows),
        'top10': rows[:10],
        'rows': rows,
    }

    (OUT_DIR / f'{name}.json').write_text(json.dumps(out, indent=2))
    return out


def main() -> None:
    _setup_entrez()

    queries = {
        # DPYD: need cohorts with BOTH toxicity and no toxicity under exposure
        'dpyd_prospective_toxicity': '(DPYD OR dihydropyrimidine dehydrogenase) AND (fluoropyrimidine OR capecitabine OR 5-fluorouracil OR 5-FU) AND (genotype-guided OR genotype guided OR dose individualisation OR dose individualization) AND (toxicity OR adverse)',
        # TPMT: thiopurine myelosuppression studies
        'tpmt_thiopurine_outcomes': '(TPMT OR thiopurine methyltransferase) AND (azathioprine OR mercaptopurine OR thiopurine) AND (toxicity OR myelosuppression) AND (genotype OR genotyping)',
        # UGT1A1: irinotecan neutropenia/diarrhea
        'ugt1a1_irinotecan_outcomes': '(UGT1A1 OR UDP glucuronosyltransferase) AND irinotecan AND (neutropenia OR diarrhea OR toxicity) AND (genotype OR genotyping)',
        # CYP2C19 clopidogrel intermediate metabolizer outcomes
        'cyp2c19_clopidogrel_im': '(CYP2C19) AND clopidogrel AND (intermediate metabolizer OR loss-of-function OR *1/*2 OR heterozygous) AND (outcome OR event OR thrombosis OR stent)',
        # CYP2D6 codeine: reduced function & toxicity/efficacy outcomes
        'cyp2d6_codeine_borderline': '(CYP2D6) AND codeine AND (intermediate metabolizer OR reduced function OR *10 OR *41) AND (pain OR toxicity OR adverse)',
    }

    receipts = []
    for name, q in queries.items():
        print('Running query:', name)
        receipts.append(run_query(name, q, retmax=25))

    summary = {
        'generated_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'receipts': [r['query_name'] for r in receipts],
        'output_dir': str(OUT_DIR),
    }
    (OUT_DIR / 'SUMMARY.json').write_text(json.dumps(summary, indent=2))
    print('Wrote receipts to', OUT_DIR)


if __name__ == '__main__':
    main()
