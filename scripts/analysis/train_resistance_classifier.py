"""
Train a resistance classifier on the merged TCGA-OV dataset (166 patients).

Goal: Achieve AUROC >= 0.70 using engineered features + class weighting.

Features:
- brca_mut (any BRCA1/2 mutation)
- ddr_count (DDR genes mutated)
- mapk_any (MAPK pathway gene mutated)
- pi3k_any (PI3K/AKT/mTOR gene mutated)
- nf1_any (NF1 mutated)
- tp53_any (TP53 mutated)
- stage_iv (IV/IVA/IVB)
- tmb (clinical_raw.TMB_NONSYNONYMOUS)
- mutation_count (clinical_raw.MUTATION_COUNT)

Author: Zo
Date: Jan 14, 2025
"""

import json
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

DDR = {"BRCA1","BRCA2","RAD51","RAD51C","RAD51D","RAD50","PALB2","ATM","ATR","CHEK1","CHEK2"}
MAPK = {"KRAS","NRAS","BRAF","MAP2K1","MAP2K2","MAPK1","MAPK3"}
PI3K = {"PIK3CA","PIK3CB","PIK3R1","AKT1","AKT2","MTOR","PTEN"}


def extract_features(rec: dict) -> dict:
    muts = rec.get('mutations', [])
    genes = {m.get('gene') for m in muts if m.get('gene')}
    clinical = rec.get('clinical_raw', {}) or {}

    brca_mut = int(bool(genes & {"BRCA1","BRCA2"}))
    ddr_count = sum(1 for g in genes if g in DDR)
    mapk_any = int(bool(genes & MAPK))
    pi3k_any = int(bool(genes & PI3K))
    nf1_any = int('NF1' in genes)
    tp53_any = int('TP53' in genes)

    stage = (rec.get('stage') or '').upper()
    stage_iv = int(stage.startswith('IV'))

    def _to_float(x):
        try:
            return float(x)
        except Exception:
            return np.nan

    tmb = _to_float(clinical.get('TMB_NONSYNONYMOUS'))
    mutation_count = _to_float(clinical.get('MUTATION_COUNT'))

    return {
        'brca_mut': brca_mut,
        'ddr_count': ddr_count,
        'mapk_any': mapk_any,
        'pi3k_any': pi3k_any,
        'nf1_any': nf1_any,
        'tp53_any': tp53_any,
        'stage_iv': stage_iv,
        'tmb': tmb if not np.isnan(tmb) else 0.0,
        'mutation_count': mutation_count if not np.isnan(mutation_count) else 0.0,
    }


def load_dataset(path: Path):
    js = json.load(open(path))
    pts = js['patients']
    X_rows, y = [], []
    for r in pts:
        feats = extract_features(r)
        X_rows.append(feats)
        y.append(int(r.get('resistant_binary', 1 if r.get('platinum_response') != 'sensitive' else 0)))
    X = pd.DataFrame(X_rows)
    y = np.array(y)
    return X, y


def cross_validate(X: pd.DataFrame, y: np.ndarray, n_splits: int = 5) -> float:
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    aucs = []
    for train_idx, test_idx in skf.split(X, y):
        Xtr, Xte = X.iloc[train_idx], X.iloc[test_idx]
        ytr, yte = y[train_idx], y[test_idx]
        clf = LogisticRegression(max_iter=500, class_weight='balanced', solver='liblinear')
        clf.fit(Xtr, ytr)
        prob = clf.predict_proba(Xte)[:,1]
        auc = roc_auc_score(yte, prob)
        aucs.append(auc)
    return float(np.mean(aucs)), float(np.std(aucs))


def main():
    path = Path('tools/benchmarks/tcga_ov_platinum_response_with_genomics.json')
    if not path.exists():
        logger.error(f"Dataset not found: {path}")
        return
    X, y = load_dataset(path)
    logger.info(f"Dataset: N={len(y)}, positives={int(y.sum())} ({y.mean():.1%})")

    mean_auc, std_auc = cross_validate(X, y, n_splits=5)
    logger.info(f"5-fold AUROC: {mean_auc:.3f} ± {std_auc:.3f}")

    out = Path('results/resistance_prophet_training')
    out.mkdir(parents=True, exist_ok=True)
    report = {
        'n': int(len(y)),
        'prevalence': float(y.mean()),
        'auroc_cv_mean': mean_auc,
        'auroc_cv_std': std_auc,
        'features': list(X.columns),
    }
    json.dump(report, open(out / 'training_report.json', 'w'), indent=2)
    logger.info(f"Report saved: {out / 'training_report.json'}")


if __name__ == '__main__':
    main()

