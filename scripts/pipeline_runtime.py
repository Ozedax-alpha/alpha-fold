"""Shared paths and optional pipeline_config.json (project root)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

BASE_DIR = Path(__file__).resolve().parents[1]
RAW_DIR = BASE_DIR / "data" / "raw"
PROCESSED_DIR = BASE_DIR / "data" / "processed"
CONFIG_PATH = BASE_DIR / "pipeline_config.json"

_DEFAULTS: dict[str, Any] = {
    "uniprot_id": "P04637",
    "gene_symbol": "TP53",
    "clinvar_esearch_term": "TP53[gene]",
    "output_basename": "tp53",
    "alphafold_fragment": "F1",
}


def load_config() -> dict[str, Any]:
    cfg = dict(_DEFAULTS)
    if CONFIG_PATH.is_file():
        user = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
        if not isinstance(user, dict):
            raise ValueError(f"{CONFIG_PATH} must contain a JSON object.")
        cfg.update(user)
    return cfg
