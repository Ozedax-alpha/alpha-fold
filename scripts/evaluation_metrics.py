"""
Held-out evaluation: pathogenic vs benign on numeric structure columns.

Uses ClinVar ``germline_date_last_evaluated`` for a time-based split when enough
dated rows exist; otherwise a stratified pseudo-random split (fixed seed).
Pure NumPy (no sklearn).
"""

from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


def _parse_eval_date(s: str) -> tuple[int, int, int] | None:
    if not isinstance(s, str) or not s.strip():
        return None
    m = re.match(r"^(\d{4})-(\d{2})-(\d{2})", s.strip())
    if not m:
        return None
    return int(m.group(1)), int(m.group(2)), int(m.group(3))


def _date_key(t: tuple[int, int, int] | None) -> int:
    if t is None:
        return 0
    y, mo, d = t
    return y * 10000 + mo * 100 + d


def roc_auc_binary(y: np.ndarray, scores: np.ndarray) -> float:
    """y in {0,1}; higher score => class 1. Mann–Whitney rank equivalent."""
    y = np.asarray(y, dtype=int)
    scores = np.asarray(scores, dtype=float)
    mask = np.isfinite(scores)
    y, scores = y[mask], scores[mask]
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(scores)
    y_sorted = y[order]
    ranks = np.empty_like(scores, dtype=float)
    ranks[order] = np.arange(1, len(scores) + 1, dtype=float)
    sum_ranks_pos = ranks[y == 1].sum()
    auc = (sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    return float(auc)


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    v = np.concatenate([a, b])
    sd = float(np.std(v, ddof=1))
    if sd == 0:
        return float("nan")
    return (float(np.mean(a)) - float(np.mean(b))) / sd


def stratified_split_indices(
    y: np.ndarray, *, test_fraction: float, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    train_parts: list[int] = []
    test_parts: list[int] = []
    for cls in (0, 1):
        idx = np.where(y == cls)[0]
        rng.shuffle(idx)
        n_test = max(1, int(round(len(idx) * test_fraction)))
        test_parts.extend(idx[:n_test].tolist())
        train_parts.extend(idx[n_test:].tolist())
    return np.array(sorted(train_parts)), np.array(sorted(test_parts))


def time_split_indices(
    dates: list[str | float | None], n_rows: int, *, test_fraction: float
) -> tuple[np.ndarray, np.ndarray] | None:
    keys = np.array([_date_key(_parse_eval_date(str(x))) if pd.notna(x) and str(x).strip() else 0 for x in dates])
    valid = keys > 0
    if int(valid.sum()) < 50:
        return None
    order = np.argsort(keys)
    n = n_rows
    n_test = max(20, int(round(n * test_fraction)))
    test_rows = set(order[-n_test:].tolist())
    train_idx = np.array([i for i in range(n) if i not in test_rows], dtype=int)
    test_idx = np.array(sorted(test_rows), dtype=int)
    return train_idx, test_idx


def main() -> None:
    ap = argparse.ArgumentParser(description="Simple held-out metrics for missense CSV.")
    ap.add_argument(
        "--run-dir",
        type=str,
        default=None,
        help="Run directory containing pipeline_config.json and data/processed/.",
    )
    ap.add_argument(
        "--missense-csv",
        type=str,
        default=None,
        help="Override path to *_missense_mappable.csv.",
    )
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--test-fraction", type=float, default=0.2)
    ap.add_argument(
        "--out-json",
        type=str,
        default=None,
        help="Write metrics JSON (default: <processed>/evaluation_metrics.json).",
    )
    args = ap.parse_args()

    repo = Path(__file__).resolve().parents[1]
    if args.run_dir:
        rd = Path(args.run_dir)
        if not rd.is_absolute():
            rd = (repo / rd).resolve()
        os.environ["PIPELINE_CONFIG_PATH"] = str(rd / "pipeline_config.json")
        os.environ["PIPELINE_OUTPUT_DIR"] = str(rd / "data")

    from pipeline_runtime import PROCESSED_DIR, load_config

    cfg = load_config()
    base = str(cfg["output_basename"])
    missense_path = (
        Path(args.missense_csv)
        if args.missense_csv
        else (PROCESSED_DIR / f"{base}_missense_mappable.csv")
    )
    if not missense_path.is_file():
        raise SystemExit(f"Missing missense CSV: {missense_path}")

    df = pd.read_csv(missense_path)
    b = df.get("clinical_significance_bucket")
    if b is None:
        raise SystemExit("Expected clinical_significance_bucket column; re-run process stage.")

    harm_mask = b.isin(["Pathogenic", "Likely pathogenic"])
    ben_mask = b.isin(["Benign", "Likely benign"])
    sub = df[harm_mask | ben_mask].copy()
    if len(sub) < 40:
        raise SystemExit(f"Too few path/benign rows for split: {len(sub)}")

    y = sub["clinical_significance_bucket"].isin(["Pathogenic", "Likely pathogenic"]).astype(int).to_numpy()

    split_mode = "random_stratified"
    train_idx: np.ndarray
    test_idx: np.ndarray

    if "germline_date_last_evaluated" in sub.columns:
        ts = time_split_indices(
            sub["germline_date_last_evaluated"].tolist(),
            len(sub),
            test_fraction=args.test_fraction,
        )
        if ts is not None:
            train_idx, test_idx = ts
            split_mode = "time_based_on_germline_date"

    if split_mode == "random_stratified":
        train_idx, test_idx = stratified_split_indices(
            y, test_fraction=args.test_fraction, seed=args.seed
        )

    feat_cols = [
        c
        for c in (
            "alphafold_confidence_score",
            "plddt_neighbor_window_5_mean",
            "ca_distance_to_centroid_angstrom",
            "residue_sasa_angstrom2",
        )
        if c in sub.columns
    ]
    if not feat_cols:
        raise SystemExit("No numeric feature columns found in missense CSV.")

    report: dict = {
        "split_mode": split_mode,
        "n_train": int(len(train_idx)),
        "n_test": int(len(test_idx)),
        "features": {},
    }

    y_tr, y_te = y[train_idx], y[test_idx]
    pos_tr = y_tr == 1
    neg_tr = y_tr == 0

    for col in feat_cols:
        x = pd.to_numeric(sub[col], errors="coerce").to_numpy(dtype=float)
        x_tr, x_te = x[train_idx], x[test_idx]
        d = cohens_d(x_tr[pos_tr], x_tr[neg_tr])
        auc = roc_auc_binary(y_te, x_te)
        report["features"][col] = {
            "cohens_d_train": d,
            "roc_auc_test_score_higher_is_pathogenic": auc,
        }

    out_path = (
        Path(args.out_json)
        if args.out_json
        else (missense_path.parent / "evaluation_metrics.json")
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Wrote {out_path}")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
