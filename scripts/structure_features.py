"""AlphaFold PDB-derived context for variant positions (CA geometry)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser


def load_ca_coords_by_resseq(pdb_path: Path) -> dict[int, np.ndarray]:
    """Map PDB residue sequence number -> CA coordinate (first occurrence)."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("af", str(pdb_path))
    out: dict[int, np.ndarray] = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                resseq = residue.id[1]
                if not isinstance(resseq, int):
                    continue
                if "CA" not in residue:
                    continue
                if resseq not in out:
                    out[resseq] = np.asarray(residue["CA"].coord, dtype=float)
    return out


def protein_centroid(ca_by_res: dict[int, np.ndarray]) -> np.ndarray | None:
    if not ca_by_res:
        return None
    pts = np.stack(list(ca_by_res.values()), axis=0)
    return pts.mean(axis=0)


def ca_distance_to_centroid_angstrom(
    positions: np.ndarray,
    ca_by_res: dict[int, np.ndarray],
    centroid: np.ndarray,
) -> np.ndarray:
    """Euclidean distance (Å) from each residue's CA to the global CA centroid."""
    out = np.full(len(positions), np.nan, dtype=float)
    for i, p in enumerate(positions):
        pi = int(p)
        coord = ca_by_res.get(pi)
        if coord is None:
            continue
        out[i] = float(np.linalg.norm(coord - centroid))
    return out


def neighbor_plddt_mean(
    positions: np.ndarray,
    smap: dict[int, float],
    *,
    window: int = 5,
) -> np.ndarray:
    """Rolling mean pLDDT over residue index neighbors [pos-window, pos+window]."""
    out = np.full(len(positions), np.nan, dtype=float)
    for i, p in enumerate(positions):
        pi = int(p)
        vals: list[float] = []
        for j in range(pi - window, pi + window + 1):
            v = smap.get(j)
            if v is not None:
                vals.append(float(v))
        if vals:
            out[i] = float(np.mean(vals))
    return out
