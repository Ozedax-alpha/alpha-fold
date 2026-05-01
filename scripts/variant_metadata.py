"""Shared variant display labels and coarse clinical buckets for reporting."""

from __future__ import annotations


def hgvs_protein_short(ref: str | None, pos: int | None, alt: str | None) -> str | None:
    """Canonical short form p.Ref123Alt when all parts are known."""
    if ref is None or pos is None or alt is None:
        return None
    r, a = str(ref).strip().upper(), str(alt).strip().upper()
    if not r or not a:
        return None
    return f"p.{r}{int(pos)}{a}"


def clinical_significance_bucket(description: str | None) -> str:
    """
    Map ClinVar free-text germline classification to a small set for plots/stats.

    Unknown / missing values become 'Uncertain / other'.
    """
    if not description or not isinstance(description, str):
        return "Uncertain / other"
    s = description.strip().lower()
    if not s:
        return "Uncertain / other"
    if "pathogenic" in s and "likely" not in s and "benign" not in s:
        return "Pathogenic"
    if "likely pathogenic" in s or s.startswith("likely pathogenic"):
        return "Likely pathogenic"
    if "benign" in s and "likely" not in s:
        return "Benign"
    if "likely benign" in s or s.startswith("likely benign"):
        return "Likely benign"
    if "uncertain" in s or "conflicting" in s:
        return "VUS / conflicting"
    return "Uncertain / other"
