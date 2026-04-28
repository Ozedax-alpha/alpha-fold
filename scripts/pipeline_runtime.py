"""Shared paths and optional pipeline_config.json (project root).

Supports per-run overrides via env vars:
- PIPELINE_CONFIG_PATH: JSON config path
- PIPELINE_OUTPUT_DIR: directory containing raw/ and processed/
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

BASE_DIR = Path(__file__).resolve().parents[1]


def _path_from_env(var: str) -> Path | None:
    v = os.environ.get(var)
    if not v:
        return None
    p = Path(v)
    if not p.is_absolute():
        p = BASE_DIR / p
    return p


_cfg_path = _path_from_env("PIPELINE_CONFIG_PATH")
CONFIG_PATH = _cfg_path if _cfg_path is not None else (BASE_DIR / "pipeline_config.json")

_out_dir = _path_from_env("PIPELINE_OUTPUT_DIR")
if _out_dir is None:
    RAW_DIR = BASE_DIR / "data" / "raw"
    PROCESSED_DIR = BASE_DIR / "data" / "processed"
else:
    RAW_DIR = _out_dir / "raw"
    PROCESSED_DIR = _out_dir / "processed"

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
