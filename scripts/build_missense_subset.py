"""
Build missense + structure-mappable subset and join AlphaFold per-residue
confidence when a confidence JSON is present (e.g. AF-...-confidence_v6.json).
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from pipeline_runtime import PROCESSED_DIR, RAW_DIR, load_config
from structure_features import (
    ca_distance_to_centroid_angstrom,
    load_ca_coords_by_resseq,
    neighbor_plddt_mean,
    protein_centroid,
)
from structure_sasa import residue_sasa_by_resseq


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
            f"Missing variants table at {basic_path}.\n"
            "Fix: run the variants stage after download, e.g.\n"
            "  py -3 -m avi run --run-dir <run-dir> --stage variants\n"
            "or: python scripts/process_tp53_variants.py"
        )

    df = pd.read_csv(basic_path)
    missense = df["molecular_consequence"].fillna("").str.contains(
        "missense", case=False
    )
    sub = df[missense & df["protein_position"].notna() & df["in_alphafold_structure"]].copy()
    pos_col = "model_residue_index" if "model_residue_index" in sub.columns else "protein_position"

    smap: dict[int, float] = {}
    conf_glob = sorted(RAW_DIR.glob(f"AF-{uid}-{frag}-confidence*.json"))
    if conf_glob:
        conf_path = conf_glob[-1]
        smap, cmap = _confidence_maps(conf_path)
        pos = sub[pos_col].astype(int)
        sub["alphafold_confidence_score"] = pos.map(lambda p: smap.get(int(p)))
        if cmap:
            sub["alphafold_confidence_category"] = pos.map(lambda p: cmap.get(int(p)))
        print(f"Joined confidence from {conf_path.name}")
    else:
        print("No AF confidence JSON found; subset written without pLDDT columns.")

    pdb_path = RAW_DIR / f"AF-{uid}-{frag}-alphafold.pdb"
    if pdb_path.is_file() and len(sub) > 0:
        try:
            ca_by_res = load_ca_coords_by_resseq(pdb_path)
            cent = protein_centroid(ca_by_res)
            pos_arr = sub[pos_col].astype(int).to_numpy()
            if cent is not None:
                sub["ca_distance_to_centroid_angstrom"] = ca_distance_to_centroid_angstrom(
                    pos_arr, ca_by_res, cent
                )
            if smap:
                sub["plddt_neighbor_window_5_mean"] = neighbor_plddt_mean(
                    pos_arr, smap, window=5
                )
            try:
                sasa_map = residue_sasa_by_resseq(pdb_path)
                sub["residue_sasa_angstrom2"] = [float(sasa_map.get(int(p), float("nan"))) for p in pos_arr]
            except Exception as e2:
                print(f"Warning: SASA computation skipped ({type(e2).__name__}: {e2})")
        except Exception as e:
            print(f"Warning: structure context features skipped ({type(e).__name__}: {e})")
    elif not pdb_path.is_file():
        print(f"Warning: no PDB at {pdb_path}; skipping centroid distance / neighbor pLDDT.")

    out = PROCESSED_DIR / f"{base}_missense_mappable.csv"
    sub.to_csv(out, index=False)
    print(f"Wrote {len(sub)} rows to {out}")


if __name__ == "__main__":
    main()
