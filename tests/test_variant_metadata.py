from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from variant_metadata import clinical_significance_bucket, hgvs_protein_short


def test_hgvs_protein_short():
    assert hgvs_protein_short("R", 175, "H") == "p.R175H"
    assert hgvs_protein_short(None, 1, "A") is None


def test_clinical_significance_bucket():
    assert clinical_significance_bucket("Pathogenic") == "Pathogenic"
    assert clinical_significance_bucket("Likely pathogenic") == "Likely pathogenic"
    assert clinical_significance_bucket("Benign") == "Benign"
    assert clinical_significance_bucket("Uncertain significance") == "VUS / conflicting"
    assert clinical_significance_bucket(None) == "Uncertain / other"
