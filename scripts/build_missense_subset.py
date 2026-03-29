"""
Build missense + structure-mappable subset and join AlphaFold per-residue
confidence when a confidence JSON is present (e.g. AF-...-confidence_v6.json).
"""

from __future__ import annotations

import json

import pandas as pd

from pipeline_runtime import PROCESSED_DIR, RAW_DIR, load_config


def _confidence_maps(conf_path: Path) -> tuple[dict[int, float], dict[int, str]]:
    j = json.loads(conf_path.read_text(encoding="utf-8"))
    nums = j.get("residueNumber")
    scores = j.get("confidenceScore")
    cats = j.get("confidenceCategory")
    if not isinstance(nums, list) or not isinstance(scores, list):
        return {}, {}
    smap = {int(a): float(b) for a, b in zip(nums, scores)}
    cmap: dict[int, str] = {}
    if isinstance(cats, list) and len(cats) == len(nums):
        cmap = {int(a): str(c) for a, c in zip(nums, cats)}
    return smap, cmap


def main() -> None:
    cfg = load_config()
    uid = cfg["uniprot_id"]
    frag = cfg["alphafold_fragment"]
    base = cfg["output_basename"]

    basic_path = PROCESSED_DIR / f"{base}_variants_basic.csv"
    if not basic_path.is_file():
        raise FileNotFoundError(
            f"Missing {basic_path}. Run process_tp53_variants.py first."
        )

    df = pd.read_csv(basic_path)
    missense = df["molecular_consequence"].fillna("").str.contains(
        "missense", case=False
    )
    sub = df[missense & df["protein_position"].notna() & df["in_alphafold_structure"]].copy()

    conf_glob = sorted(RAW_DIR.glob(f"AF-{uid}-{frag}-confidence*.json"))
    if conf_glob:
        conf_path = conf_glob[-1]
        smap, cmap = _confidence_maps(conf_path)
        pos = sub["protein_position"].astype(int)
        sub["alphafold_confidence_score"] = pos.map(lambda p: smap.get(int(p)))
        if cmap:
            sub["alphafold_confidence_category"] = pos.map(lambda p: cmap.get(int(p)))
        print(f"Joined confidence from {conf_path.name}")
    else:
        print("No AF confidence JSON found; subset written without pLDDT columns.")

    out = PROCESSED_DIR / f"{base}_missense_mappable.csv"
    sub.to_csv(out, index=False)
    print(f"Wrote {len(sub)} rows to {out}")


if __name__ == "__main__":
    main()
