"""
Parse ClinVar JSON from download_tp53_data.py, extract variant features, and map
protein positions to AlphaFold residue indices (UniProt numbering + optional offset).
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from Bio.PDB import PDBParser

from pipeline_runtime import PROCESSED_DIR, RAW_DIR, load_config
from variant_metadata import clinical_significance_bucket, hgvs_protein_short
from variant_parse import germline_date_last_evaluated, resolve_missense_position

PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

_CFG = load_config()
UNIPROT_ID = str(_CFG["uniprot_id"])
_FRAG = str(_CFG["alphafold_fragment"])
_OUT_BASE = str(_CFG["output_basename"])
PDB_PATH = RAW_DIR / f"AF-{UNIPROT_ID}-{_FRAG}-alphafold.pdb"
CLINVAR_JSON = RAW_DIR / f"clinvar_{_OUT_BASE}_variants.json"
OUT_CSV = PROCESSED_DIR / f"{_OUT_BASE}_variants_basic.csv"


def _offset(cfg: dict) -> int:
    v = cfg.get("uniprot_residue_offset")
    try:
        return int(v) if v is not None and str(v).strip() != "" else 0
    except (TypeError, ValueError):
        return 0


def _preferred_transcript(cfg: dict) -> str | None:
    p = cfg.get("preferred_transcript_prefix")
    if p is None or (isinstance(p, str) and not p.strip()):
        return None
    return str(p).strip()


def load_alphafold_residue_ids(pdb_path: Path) -> set[int]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("af", str(pdb_path))
    ids: set[int] = set()
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                resseq = residue.id[1]
                if isinstance(resseq, int):
                    ids.add(resseq)
    return ids


def _grch38_loc(variation_set: list) -> tuple[str | None, str | None, str | None]:
    if not variation_set or not isinstance(variation_set, list):
        return None, None, None
    locs = variation_set[0].get("variation_loc") or []
    for loc in locs:
        if loc.get("status") == "current" and loc.get("assembly_name") == "GRCh38":
            return loc.get("chr"), loc.get("start"), loc.get("stop")
    return None, None, None


def _clinical_significance(rec: dict) -> str | None:
    g = rec.get("germline_classification") or {}
    desc = g.get("description")
    if desc:
        return str(desc)
    o = rec.get("oncogenicity_classification") or {}
    if o.get("description"):
        return str(o["description"])
    return None


def iter_variant_rows(
    summaries: dict[str, dict], *, cfg: dict, af_residues: set[int]
) -> list[dict]:
    rows: list[dict] = []
    pref = _preferred_transcript(cfg)
    off = _offset(cfg)

    for uid, rec in summaries.items():
        if not isinstance(rec, dict):
            continue
        vset = rec.get("variation_set") or []
        if not vset:
            continue
        first = vset[0]
        variation_name = first.get("variation_name") or rec.get("title") or ""
        protein_change = (rec.get("protein_change") or "").strip()
        chr38, start38, stop38 = _grch38_loc(vset)
        sig = _clinical_significance(rec)
        gdate = germline_date_last_evaluated(rec)

        ref, pos, alt = resolve_missense_position(
            protein_change_field=protein_change,
            variation_name=str(variation_name),
            preferred_transcript_prefix=pref,
        )

        mol_cons = rec.get("molecular_consequence_list") or []
        mol_cons_str = ";".join(mol_cons) if isinstance(mol_cons, list) else ""

        hgvs_short = hgvs_protein_short(ref, pos, alt)
        model_res = int(pos) + off if pos is not None else None
        in_af = model_res in af_residues if model_res is not None else False
        mapping_note = None
        if pos is not None and model_res is not None and not in_af:
            mapping_note = "model_residue_not_in_alphafold_pdb"

        rows.append(
            {
                "clinvar_uid": uid,
                "accession": rec.get("accession"),
                "variation_name": variation_name,
                "grch38_chromosome": chr38,
                "grch38_start": start38,
                "grch38_stop": stop38,
                "clinical_significance": sig,
                "clinical_significance_bucket": clinical_significance_bucket(sig),
                "germline_date_last_evaluated": gdate,
                "variant_type": first.get("variant_type"),
                "molecular_consequence": mol_cons_str,
                "protein_ref": ref,
                "protein_position": pos,
                "protein_alt": alt,
                "protein_change_field": protein_change or None,
                "hgvs_protein_short": hgvs_short,
                "model_residue_index": model_res,
                "mapping_note": mapping_note,
            }
        )
    return rows


def main() -> None:
    if not CLINVAR_JSON.is_file():
        raise FileNotFoundError(
            f"Missing ClinVar JSON at {CLINVAR_JSON}.\n"
            "Fix: run the download stage first, e.g.\n"
            "  py -3 -m avi run --run-dir <run-dir> --stage download\n"
            "or: python scripts/download_tp53_data.py"
        )
    if not PDB_PATH.is_file():
        raise FileNotFoundError(
            f"Missing AlphaFold PDB at {PDB_PATH}.\n"
            "Fix: run the download stage first (same as ClinVar), with network access."
        )

    doc = json.loads(CLINVAR_JSON.read_text(encoding="utf-8"))
    summaries = doc.get("summaries")
    if not isinstance(summaries, dict):
        raise ValueError("Expected top-level 'summaries' object in ClinVar JSON.")

    cfg = load_config()
    af_residues = load_alphafold_residue_ids(PDB_PATH)
    rows = iter_variant_rows(summaries, cfg=cfg, af_residues=af_residues)
    df = pd.DataFrame(rows)

    df["alphafold_residue_index"] = df["model_residue_index"]
    df["in_alphafold_structure"] = df["model_residue_index"].apply(
        lambda p: int(p) in af_residues if pd.notna(p) else False
    )

    df.to_csv(OUT_CSV, index=False)
    print(f"Wrote {len(df)} rows to {OUT_CSV}")


if __name__ == "__main__":
    main()
