"""Per-residue SASA (Å²) via Shrake–Rupley (Biopython; no external DSSP binary)."""

from __future__ import annotations

from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley


def residue_sasa_by_resseq(pdb_path: Path, *, probe_radius: float = 1.4) -> dict[int, float]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("af", str(pdb_path))
    sr = ShrakeRupley(probe_radius=probe_radius)
    sr.compute(structure, level="R")
    out: dict[int, float] = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                resseq = residue.id[1]
                if not isinstance(resseq, int):
                    continue
                sasa = getattr(residue, "sasa", None)
                if sasa is None:
                    continue
                if resseq not in out:
                    out[resseq] = float(sasa)
    return out
