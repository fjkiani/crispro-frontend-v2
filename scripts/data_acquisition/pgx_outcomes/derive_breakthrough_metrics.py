import json
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
IN_DIR = REPO_ROOT / 'scripts' / 'data_acquisition' / 'pgx_outcomes' / 'receipts' / 'pmc_structured'
OUT_DIR = REPO_ROOT / 'scripts' / 'data_acquisition' / 'pgx_outcomes' / 'receipts'
OUT_DIR.mkdir(parents=True, exist_ok=True)


def parse_int(s: str) -> int:
    return int(str(s).strip())


def parse_events_cell(s: str) -> int:
    # formats like "54 (17.4)" or "0"
    s = (s or '').strip()
    if s == '0':
        return 0
    num = s.split('(')[0].strip()
    return int(num)


def prepare_toxicity_from_table2(path: Path):
    d = json.loads(path.read_text())
    tbl = d['table'] if 'table' in d else d['tables']['Table 2']
    rows = tbl['rows']

    # Build lookup by section/arm
    # Expected layout:
    # [header]
    # ["Any actionable genotype carriers"]
    # ["Control arm", N, tox]
    # ["Intervention arm", N, tox]
    # ["All patients"]
    # ["Control arm", N, tox]
    # ["Intervention arm", N, tox]
    sections = {}
    current = None
    for r in rows[1:]:
        if len(r) >= 1 and r[1:] == ['', '', '', '']:
            current = r[0]
            sections[current] = {}
            continue
        if current and r[0] in ('Control arm', 'Intervention arm'):
            sections[current][r[0]] = {
                'n': parse_int(r[1]),
                'tox_events': parse_events_cell(r[2]),
                'tox_cell': r[2],
            }

    # Derive nonactionable counts by subtraction using "All patients" - "Any actionable genotype carriers"
    all_ctrl = sections['All patients']['Control arm']
    all_int = sections['All patients']['Intervention arm']
    act_ctrl = sections['Any actionable genotype carriers']['Control arm']
    act_int = sections['Any actionable genotype carriers']['Intervention arm']

    non_ctrl = {
        'n': all_ctrl['n'] - act_ctrl['n'],
        'tox_events': all_ctrl['tox_events'] - act_ctrl['tox_events'],
    }
    non_int = {
        'n': all_int['n'] - act_int['n'],
        'tox_events': all_int['tox_events'] - act_int['tox_events'],
    }

    return {
        'pmid': d['pmid'],
        'pmc_id': d['pmc_id'],
        'endpoint': 'clinically relevant toxic effects',
        'arms': {
            'control': {
                'all_patients': all_ctrl,
                'actionable_carriers': act_ctrl,
                'nonactionable': non_ctrl,
            },
            'intervention': {
                'all_patients': all_int,
                'actionable_carriers': act_int,
                'nonactionable': non_int,
            },
        },
        'note': 'Derived from Table 2 by subtraction; provides true negative controls (nonactionable) and carrier vs noncarrier event counts.'
    }


def cyp2c19_outcomes_from_table(path: Path, table_key: str):
    d = json.loads(path.read_text())
    tbl = d['tables'][table_key]
    rows = tbl['rows']

    # header row gives N
    header = rows[0]
    em_n = int(header[1].split('n = ')[1].split(')')[0])
    pm_im_n = int(header[2].split('n = ')[1].split(')')[0])

    # find the block for symptomatic ischemic stroke/TIA
    # pattern: ["Symptomatic ischemic stroke/TIA",...]
    # then ["Events", em, pm_im]
    ev_em = None
    ev_pmim = None
    for i,r in enumerate(rows):
        if r and r[0].strip().lower() == 'symptomatic ischemic stroke/tia':
            ev_row = rows[i+1]
            if ev_row[0].strip().lower() == 'events':
                ev_em = parse_events_cell(ev_row[1])
                ev_pmim = parse_events_cell(ev_row[2])
            break

    if ev_em is None:
        raise RuntimeError('Could not find events row')

    return {
        'pmid': d['pmid'],
        'pmc_id': d['pmc_id'],
        'population': 'clopidogrel-treated' if table_key == 'Table 4' else 'overall',
        'endpoint': 'symptomatic ischemic stroke/TIA',
        'groups': {
            'extensive_metabolizer': {'n': em_n, 'events': ev_em},
            'poor_or_intermediate_metabolizer': {'n': pm_im_n, 'events': ev_pmim},
        },
        'note': 'Counts extracted directly from outcome table; provides borderline IM cases and outcome-linked efficacy signal.'
    }


def main():
    generated_at = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())

    # PREPARE: use the dedicated structured table 2 file
    prepare_t2 = IN_DIR / 'pmid_39641926_Table_2.json'
    prep = prepare_toxicity_from_table2(prepare_t2)

    # CYP2C19 cohort: use batch extraction file
    cyp_batch = IN_DIR / 'pmid_40944685_tables_Table2_Table4.json'
    cyp_overall = cyp2c19_outcomes_from_table(cyp_batch, 'Table 2')
    cyp_clop = cyp2c19_outcomes_from_table(cyp_batch, 'Table 4')

    out = {
        'generated_at': generated_at,
        'breakthroughs': {
            'dpyd_ugt1a1_prepare_rct': prep,
            'cyp2c19_clopidogrel_outcomes': {
                'overall': cyp_overall,
                'clopidogrel_treated': cyp_clop,
            }
        },
        'why_this_matters': [
            'Adds true outcome-linked negatives (nonactionable genotypes with observed toxicity vs no toxicity).',
            'Adds borderline CYP2C19 intermediate metabolizer efficacy outcomes on clopidogrel.',
            'Supports moving beyond DPYD-only cohorts toward more representative PGx validation.'
        ]
    }

    (OUT_DIR / 'BREAKTHROUGH_OUTCOME_LINKED_METRICS.json').write_text(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
