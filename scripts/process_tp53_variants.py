"""
Parse ClinVar JSON from download_tp53_data.py, extract variant features, and map
protein positions to AlphaFold residue indices (UniProt P04637 numbering).
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
from Bio.PDB import PDBParser

from pipeline_runtime import PROCESSED_DIR, RAW_DIR, load_config

PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

_CFG = load_config()
UNIPROT_ID = str(_CFG["uniprot_id"])
_FRAG = str(_CFG["alphafold_fragment"])
_OUT_BASE = str(_CFG["output_basename"])
PDB_PATH = RAW_DIR / f"AF-{UNIPROT_ID}-{_FRAG}-alphafold.pdb"
CLINVAR_JSON = RAW_DIR / f"clinvar_{_OUT_BASE}_variants.json"
OUT_CSV = PROCESSED_DIR / f"{_OUT_BASE}_variants_basic.csv"

THREE_TO_ONE = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
    "Sec": "U",
}


def _parse_protein_change_short(text: str) -> tuple[str | None, int | None, str | None]:
    """
    NCBI esummary often uses protein_change like 'R175H' or 'S15I'.
    """
    if not text or not isinstance(text, str):
        return None, None, None
    text = text.strip()
    m = re.match(r"^([A-Za-z*?])(\d+)([A-Za-z*?])$", text)
    if not m:
        return None, None, None
    ref, pos_s, alt = m.group(1).upper(), m.group(2), m.group(3).upper()
    return ref, int(pos_s), alt


def _parse_protein_paren_hgvs(text: str) -> tuple[str | None, int | None, str | None]:
    """
    Extract p.Leu123Arg from strings like 'NM_000546.6(TP53):c.123A>G (p.Leu123Arg)'.
    """
    if not text:
        return None, None, None
    m = re.search(
        r"\(p\.([A-Za-z]{3})(\d+)([A-Za-z]{3}|\*)\)",
        text,
    )
    if not m:
        return None, None, None
    ref3, pos_s, alt3 = m.group(1), m.group(2), m.group(3)
    ref = THREE_TO_ONE.get(ref3)
    alt = THREE_TO_ONE.get(alt3) if alt3 != "*" else "*"
    if ref is None or alt is None:
        return None, None, None
    return ref, int(pos_s), alt


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


def load_alphafold_residue_ids(pdb_path: Path) -> set[int]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("tp53", str(pdb_path))
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


def iter_variant_rows(summaries: dict[str, dict]) -> list[dict]:
    rows: list[dict] = []
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

        ref, pos, alt = _parse_protein_change_short(protein_change)
        if pos is None:
            ref, pos, alt = _parse_protein_paren_hgvs(variation_name)

        mol_cons = rec.get("molecular_consequence_list") or []
        mol_cons_str = ";".join(mol_cons) if isinstance(mol_cons, list) else ""

        rows.append(
            {
                "clinvar_uid": uid,
                "accession": rec.get("accession"),
                "variation_name": variation_name,
                "grch38_chromosome": chr38,
                "grch38_start": start38,
                "grch38_stop": stop38,
                "clinical_significance": sig,
                "variant_type": first.get("variant_type"),
                "molecular_consequence": mol_cons_str,
                "protein_ref": ref,
                "protein_position": pos,
                "protein_alt": alt,
                "protein_change_field": protein_change or None,
            }
        )
    return rows


def main() -> None:
    if not CLINVAR_JSON.is_file():
        raise FileNotFoundError(
            f"Missing {CLINVAR_JSON}. Run scripts/download_tp53_data.py first."
        )
    if not PDB_PATH.is_file():
        raise FileNotFoundError(
            f"Missing {PDB_PATH}. Run scripts/download_tp53_data.py first."
        )

    doc = json.loads(CLINVAR_JSON.read_text(encoding="utf-8"))
    summaries = doc.get("summaries")
    if not isinstance(summaries, dict):
        raise ValueError("Expected top-level 'summaries' object in ClinVar JSON.")

    af_residues = load_alphafold_residue_ids(PDB_PATH)
    rows = iter_variant_rows(summaries)
    df = pd.DataFrame(rows)

    df["alphafold_residue_index"] = df["protein_position"]
    df["in_alphafold_structure"] = df["protein_position"].apply(
        lambda p: int(p) in af_residues if pd.notna(p) else False
    )

    df.to_csv(OUT_CSV, index=False)
    print(f"Wrote {len(df)} rows to {OUT_CSV}")


if __name__ == "__main__":
    main()
