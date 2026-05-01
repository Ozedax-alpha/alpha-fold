from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from structure_features import neighbor_plddt_mean


def test_neighbor_plddt_mean_full_window():
    pos = np.array([10])
    smap = {i: float(i) for i in range(5, 16)}
    out = neighbor_plddt_mean(pos, smap, window=5)
    assert out.shape == (1,)
    assert abs(out[0] - np.mean(list(range(5, 16)))) < 1e-9
