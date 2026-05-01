"""Resolved artifact paths from pipeline config + data directory (matches scripts/)."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def default_config_path(base_dir: Path) -> Path:
    return base_dir / "pipeline_config.json"


def artifact_paths(cfg: dict[str, Any], data_dir: Path) -> dict[str, Path]:
    uid = str(cfg["uniprot_id"])
    frag = str(cfg["alphafold_fragment"])
    base = str(cfg["output_basename"])
    raw = data_dir / "raw"
    proc = data_dir / "processed"
    safe_base = base.replace("/", "_").replace("\\", "_")
    return {
        "raw_dir": raw,
        "processed_dir": proc,
        "alphafold_pdb": raw / f"AF-{uid}-{frag}-alphafold.pdb",
        "clinvar_json": raw / f"clinvar_{base}_variants.json",
        "variants_basic_csv": proc / f"{base}_variants_basic.csv",
        "missense_mappable_csv": proc / f"{base}_missense_mappable.csv",
        "repro_manifest": proc / f"repro_manifest_{safe_base}.json",
        "evaluation_metrics_json": proc / "evaluation_metrics.json",
    }


def format_paths_summary(cfg: dict[str, Any], data_dir: Path) -> str:
    art = artifact_paths(cfg, data_dir)
    lines = [
        "Effective paths:",
        f"  Data root (raw/, processed/): {data_dir}",
        "",
        "Artifacts:",
        f"  AlphaFold PDB:            {art['alphafold_pdb'].name}",
        f"  ClinVar JSON:             {art['clinvar_json'].name}",
        f"  Variants CSV:             {art['variants_basic_csv'].name}",
        f"  Missense CSV:             {art['missense_mappable_csv'].name}",
        f"  Repro manifest (--trail): {art['repro_manifest'].name}",
        f"  Evaluation metrics:       {art['evaluation_metrics_json'].name} (from `avi evaluate`)",
    ]
    return "\n".join(lines)
