"""HGVS / ClinVar protein position parsing (no heavy dependencies)."""

from __future__ import annotations

import re
from typing import Any

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


def slice_variation_name_for_transcript(
    variation_name: str, preferred_transcript_prefix: str | None
) -> str:
    """
    When ClinVar packs multiple transcripts into one string, prefer the segment
    for NM_... matching ``preferred_transcript_prefix`` (e.g. NM_000546).
    """
    if not variation_name or not preferred_transcript_prefix:
        return variation_name or ""
    pref = preferred_transcript_prefix.strip()
    if not pref:
        return variation_name
    # Split on " NM_" boundaries (space before alternate RefSeq).
    parts = re.split(r"(?=\sNM_\d)", variation_name)
    for part in parts:
        p = part.strip()
        if p.startswith(pref) or pref in p[: len(pref) + 8]:
            return p
    return variation_name


def parse_protein_change_short(text: str) -> tuple[str | None, int | None, str | None]:
    if not text or not isinstance(text, str):
        return None, None, None
    text = text.strip()
    m = re.match(r"^([A-Za-z*?])(\d+)([A-Za-z*?])$", text)
    if not m:
        return None, None, None
    ref, pos_s, alt = m.group(1).upper(), m.group(2), m.group(3).upper()
    return ref, int(pos_s), alt


def parse_protein_paren_hgvs(text: str) -> tuple[str | None, int | None, str | None]:
    """e.g. (p.Leu123Arg) inside a longer variation name."""
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


def parse_hgvs_bracket_form(text: str) -> tuple[str | None, int | None, str | None]:
    """p.(Leu123Arg) without outer transcript parens."""
    if not text:
        return None, None, None
    m = re.search(r"p\.\(([A-Za-z]{3})(\d+)([A-Za-z]{3}|\*)\)", text)
    if not m:
        return None, None, None
    ref3, pos_s, alt3 = m.group(1), m.group(2), m.group(3)
    ref = THREE_TO_ONE.get(ref3)
    alt = THREE_TO_ONE.get(alt3) if alt3 != "*" else "*"
    if ref is None or alt is None:
        return None, None, None
    return ref, int(pos_s), alt


def parse_hgvs_three_letter_loose(text: str) -> tuple[str | None, int | None, str | None]:
    """p.Leu123Arg (no extra parentheses)."""
    if not text:
        return None, None, None
    m = re.search(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3}|\*)", text)
    if not m:
        return None, None, None
    ref3, pos_s, alt3 = m.group(1), m.group(2), m.group(3)
    ref = THREE_TO_ONE.get(ref3)
    alt = THREE_TO_ONE.get(alt3) if alt3 != "*" else "*"
    if ref is None or alt is None:
        return None, None, None
    return ref, int(pos_s), alt


def parse_hgvs_one_letter(text: str) -> tuple[str | None, int | None, str | None]:
    """p.W53* or p.R175H (one-letter)."""
    if not text:
        return None, None, None
    m = re.search(r"p\.([A-Za-z*?])(\d+)([A-Za-z*?])", text)
    if not m:
        return None, None, None
    ref, pos_s, alt = m.group(1).upper(), m.group(2), m.group(3).upper()
    return ref, int(pos_s), alt


def resolve_missense_position(
    *,
    protein_change_field: str,
    variation_name: str,
    preferred_transcript_prefix: str | None,
) -> tuple[str | None, int | None, str | None]:
    """Apply parsers in a stable order."""
    ref, pos, alt = parse_protein_change_short(protein_change_field)
    if pos is not None:
        return ref, pos, alt

    vslice = slice_variation_name_for_transcript(variation_name, preferred_transcript_prefix)
    for parser in (
        parse_protein_paren_hgvs,
        parse_hgvs_bracket_form,
        parse_hgvs_three_letter_loose,
        parse_hgvs_one_letter,
    ):
        ref, pos, alt = parser(vslice)
        if pos is not None:
            return ref, pos, alt

    # Last resort: full variation_name without transcript slicing
    if vslice != variation_name:
        for parser in (
            parse_protein_paren_hgvs,
            parse_hgvs_bracket_form,
            parse_hgvs_three_letter_loose,
            parse_hgvs_one_letter,
        ):
            ref, pos, alt = parser(variation_name)
            if pos is not None:
                return ref, pos, alt

    return None, None, None


def parse_clinvar_last_evaluated_tokens(s: str) -> tuple[int, int, int] | None:
    """
    Parse ClinVar esummary date strings to a calendar day.

    NCBI often returns ``2026/01/15 00:00``; we also accept ISO ``2026-01-15``.
    """
    if not isinstance(s, str) or not s.strip():
        return None
    t = s.strip()
    m = re.match(r"^(\d{4})-(\d{2})-(\d{2})", t)
    if not m:
        m = re.match(r"^(\d{4})/(\d{2})/(\d{2})", t)
    if not m:
        return None
    return int(m.group(1)), int(m.group(2)), int(m.group(3))


def germline_date_last_evaluated(rec: dict[str, Any]) -> str | None:
    g = rec.get("germline_classification")
    if not isinstance(g, dict):
        return None
    raw: str | None = None
    for key in ("date_last_evaluated", "last_evaluated", "review_date"):
        v = g.get(key)
        if isinstance(v, str) and v.strip():
            raw = v.strip()
            break
    if not raw:
        return None
    tok = parse_clinvar_last_evaluated_tokens(raw)
    if tok is None:
        return None
    y, mo, d = tok
    return f"{y:04d}-{mo:02d}-{d:02d}"
