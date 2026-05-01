from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from variant_parse import (
    germline_date_last_evaluated,
    parse_clinvar_last_evaluated_tokens,
    resolve_missense_position,
    slice_variation_name_for_transcript,
)


def test_slice_variation_name_prefers_transcript():
    name = (
        "NM_000546.5(TP53):c.524G>A (p.Arg175His) "
        "NM_001126112.2(TP53):c.524G>A (p.Arg175His)"
    )
    s = slice_variation_name_for_transcript(name, "NM_000546")
    assert s.startswith("NM_000546")
    ref, pos, alt = resolve_missense_position(
        protein_change_field="",
        variation_name=s,
        preferred_transcript_prefix="NM_000546",
    )
    assert (ref, pos, alt) == ("R", 175, "H")


def test_resolve_bracket_hgvs():
    ref, pos, alt = resolve_missense_position(
        protein_change_field="",
        variation_name="p.(Arg175His)",
        preferred_transcript_prefix=None,
    )
    assert (ref, pos, alt) == ("R", 175, "H")


def test_resolve_one_letter_hgvs():
    ref, pos, alt = resolve_missense_position(
        protein_change_field="",
        variation_name="NM_000546.5(TP53):c.524G>A (p.R175H)",
        preferred_transcript_prefix=None,
    )
    assert (ref, pos, alt) == ("R", 175, "H")


def test_parse_clinvar_date_slash_and_iso():
    assert parse_clinvar_last_evaluated_tokens("2026/01/15 00:00") == (2026, 1, 15)
    assert parse_clinvar_last_evaluated_tokens("2020-12-03") == (2020, 12, 3)


def test_germline_date_normalizes_to_iso():
    rec = {
        "germline_classification": {
            "description": "Pathogenic",
            "last_evaluated": "2024/06/01 12:00",
        }
    }
    assert germline_date_last_evaluated(rec) == "2024-06-01"
