#!/usr/bin/env python3
import json
import os
import sys
import math
import urllib.request
from typing import List, Tuple

BACKEND = os.environ.get('BACKEND_URL', 'https://crispro-oncology-backend-minimal.vercel.app')
PANEL = os.environ.get('PANEL_PATH', 'docs/benchmarks/myeloma_panel.json')
MODEL = os.environ.get('MODEL_ID', 'evo2_7b')
OUT_SUMMARY = os.environ.get('OUT_SUMMARY', 'benchmark_output.json')
OUT_RAW = os.environ.get('OUT_RAW', 'benchmark_raw.json')
OUT_METRICS = os.environ.get('OUT_METRICS', 'benchmark_metrics.json')
CHUNK_SIZE = int(os.environ.get('BENCHMARK_CHUNK_SIZE', '20'))


def post_json(url: str, payload: dict, timeout: int = 600) -> dict:
    data = json.dumps(payload).encode('utf-8')
    req = urllib.request.Request(url, data=data, headers={'Content-Type': 'application/json'})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        body = resp.read().decode('utf-8')
        return json.loads(body)


def _safe_get(d: dict, path: List[str], default=None):
    cur = d
    for k in path:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


def _auc_roc(y_true: List[int], y_score: List[float]) -> float:
    paired = sorted(zip(y_score, y_true), key=lambda x: x[0])
    n_pos = sum(y_true)
    n_neg = len(y_true) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float('nan')
    rank_sum = 0.0
    for idx, (_, y) in enumerate(paired, start=1):
        if y == 1:
            rank_sum += idx
    auc = (rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    return float(auc)


def _auc_pr(y_true: List[int], y_score: List[float]) -> float:
    paired = sorted(zip(y_score, y_true), key=lambda x: -x[0])
    tp = 0
    fp = 0
    prev_recall = 0.0
    area = 0.0
    pos_total = sum(y_true)
    if pos_total == 0:
        return float('nan')
    for score, y in paired:
        if y == 1:
            tp += 1
        else:
            fp += 1
        recall = tp / pos_total
        precision = tp / max(1, (tp + fp))
        area += precision * (recall - prev_recall)
        prev_recall = recall
    return float(area)


def _sigmoid(x: float) -> float:
    try:
        return 1.0 / (1.0 + math.exp(-x))
    except OverflowError:
        return 0.0 if x < 0 else 1.0


def compute_metrics(panel: List[dict], result: dict) -> dict:
    detailed = result.get('detailed_analysis') or []
    if not isinstance(detailed, list):
        detailed = []

    labels: List[int] = []
    scores: List[float] = []
    confs: List[float] = []

    for idx, entry in enumerate(panel):
        label = entry.get('label_is_pathogenic', None)
        if label is None:
            continue
        if idx >= len(detailed):
            continue
        evo = detailed[idx].get('evo2_result') or {}
        min_delta = evo.get('min_delta')
        if isinstance(min_delta, (int, float)):
            score = -float(min_delta)
        else:
            z = evo.get('zeta_score')
            score = -float(z) if isinstance(z, (int, float)) else 0.0
        prob_like = _sigmoid(score)
        labels.append(1 if bool(label) else 0)
        scores.append(prob_like)
        conf = evo.get('confidence_score')
        confs.append(float(conf) if isinstance(conf, (int, float)) else None)

    metrics = {
        'num_total': len(panel),
        'num_labeled': len(labels),
        'has_both_classes': (sum(labels) > 0 and sum(labels) < len(labels)),
        'auroc': None,
        'auprc': None,
        'calibration': None,
    }

    if len(labels) >= 2 and metrics['has_both_classes']:
        metrics['auroc'] = _auc_roc(labels, scores)
        metrics['auprc'] = _auc_pr(labels, scores)

    cal = []
    if confs and len(labels) == len(confs):
        bins = [0.0, 0.25, 0.5, 0.75, 1.01]
        for i in range(len(bins) - 1):
            lo, hi = bins[i], bins[i + 1]
            items = [(labels[j], confs[j]) for j in range(len(labels)) if confs[j] is not None and lo <= confs[j] < hi]
            if items:
                avg_conf = sum(c for _, c in items) / len(items)
                avg_label = sum(y for y, _ in items) / len(items)
                cal.append({'bin_lo': lo, 'bin_hi': hi, 'avg_confidence': avg_conf, 'empirical_positive_rate': avg_label, 'n': len(items)})
    metrics['calibration'] = cal if cal else None

    return metrics


def chunks(lst: List[dict], n: int) -> List[List[dict]]:
    return [lst[i:i+n] for i in range(0, len(lst), n)]


def main():
    with open(PANEL, 'r') as f:
        panel = json.load(f)
    url = f"{BACKEND.rstrip('/')}/api/predict/myeloma_drug_response"

    # Chunked execution to avoid gateway timeouts
    all_details: List[dict] = []
    pathway_scores_agg = {'summed_impact_ras_pathway': 0.0, 'summed_impact_tp53': 0.0}
    last_resp = None
    for chunk in chunks(panel, CHUNK_SIZE):
        payload = {'model_id': MODEL, 'mutations': chunk}
        resp = post_json(url, payload)
        last_resp = resp
        details = resp.get('detailed_analysis') or []
        all_details.extend(details)
        ps = resp.get('pathway_scores') or {}
        # Aggregate naïvely; final label is re-derived from sum
        pathway_scores_agg['summed_impact_ras_pathway'] += float(ps.get('summed_impact_ras_pathway', 0.0) or 0.0)
        pathway_scores_agg['summed_impact_tp53'] += float(ps.get('summed_impact_tp53', 0.0) or 0.0)

    prediction_label = "Likely Resistant" if pathway_scores_agg['summed_impact_ras_pathway'] >= 2.0 else "Likely Sensitive"

    result_combined = {
        'prediction': prediction_label,
        'pathway_scores': pathway_scores_agg,
        'mode': last_resp.get('mode') if last_resp else 'live',
        'selected_model': last_resp.get('selected_model') if last_resp else MODEL,
        'upstream_service': last_resp.get('upstream_service') if last_resp else None,
        'run_signature': last_resp.get('run_signature') if last_resp else None,
        'detailed_analysis': all_details,
    }

    with open(OUT_RAW, 'w') as f:
        json.dump(result_combined, f)

    metrics = compute_metrics(panel, result_combined)
    with open(OUT_METRICS, 'w') as f:
        json.dump(metrics, f)

    summary = {
        'prediction': result_combined.get('prediction'),
        'pathway_scores': result_combined.get('pathway_scores'),
        'mode': result_combined.get('mode'),
        'selected_model': result_combined.get('selected_model'),
        'upstream_service': result_combined.get('upstream_service'),
        'run_signature': result_combined.get('run_signature'),
        'metrics': metrics,
    }
    print(json.dumps(summary, indent=2))
    with open(OUT_SUMMARY, 'w') as f:
        json.dump(summary, f)


if __name__ == '__main__':
    main()