from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def _synthetic_missense_csv(path: Path, *, n: int = 120) -> None:
    rows = []
    for i in range(n):
        pathogenic = i % 2 == 0
        rows.append(
            {
                "clinical_significance_bucket": "Pathogenic" if pathogenic else "Benign",
                "protein_position": 100 + (i % 20),
                "alphafold_confidence_score": 40.0 + i * 0.1 if pathogenic else 80.0 - i * 0.05,
                "plddt_neighbor_window_5_mean": 45.0 + i * 0.05 if pathogenic else 78.0,
                "ca_distance_to_centroid_angstrom": 12.0 + (i % 5) if pathogenic else 8.0,
                "residue_sasa_angstrom2": 30.0 if pathogenic else 90.0,
                "germline_date_last_evaluated": f"2020-01-{(i % 28) + 1:02d}",
            }
        )
    pd.DataFrame(rows).to_csv(path, index=False)


def test_evaluation_metrics_script_writes_json(tmp_path: Path) -> None:
    run_dir = tmp_path / "run1"
    proc = run_dir / "data" / "processed"
    proc.mkdir(parents=True)
    cfg = {
        "uniprot_id": "P00000",
        "gene_symbol": "TEST",
        "clinvar_esearch_term": "TEST[gene]",
        "output_basename": "syn",
        "alphafold_fragment": "F1",
    }
    (run_dir / "pipeline_config.json").write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
    _synthetic_missense_csv(proc / "syn_missense_mappable.csv")

    subprocess.check_call(
        [
            sys.executable,
            str(ROOT / "scripts" / "evaluation_metrics.py"),
            "--run-dir",
            str(run_dir),
            "--seed",
            "0",
        ],
        cwd=str(ROOT),
    )
    out = proc / "evaluation_metrics.json"
    assert out.is_file()
    doc = json.loads(out.read_text(encoding="utf-8"))
    assert doc["split_mode"] in ("time_based_on_germline_date", "random_stratified")
    assert "split_metadata" in doc
    assert doc["split_metadata"]["n_rows_with_parseable_germline_date"] == 120
    assert "alphafold_confidence_score" in doc["features"]


def test_evaluation_metrics_handles_single_class(tmp_path: Path) -> None:
    run_dir = tmp_path / "run2"
    proc = run_dir / "data" / "processed"
    proc.mkdir(parents=True)
    cfg = {
        "uniprot_id": "P00000",
        "gene_symbol": "TEST",
        "clinvar_esearch_term": "TEST[gene]",
        "output_basename": "syn",
        "alphafold_fragment": "F1",
    }
    (run_dir / "pipeline_config.json").write_text(
        json.dumps(cfg, indent=2) + "\n", encoding="utf-8"
    )
    # All pathogenic => single-class labeled subset.
    rows = [
        {
            "clinical_significance_bucket": "Pathogenic",
            "alphafold_confidence_score": 50.0,
            "germline_date_last_evaluated": "2020-01-01",
        }
        for _ in range(60)
    ]
    pd.DataFrame(rows).to_csv(proc / "syn_missense_mappable.csv", index=False)

    subprocess.check_call(
        [
            sys.executable,
            str(ROOT / "scripts" / "evaluation_metrics.py"),
            "--run-dir",
            str(run_dir),
            "--seed",
            "0",
        ],
        cwd=str(ROOT),
    )
    out = proc / "evaluation_metrics.json"
    doc = json.loads(out.read_text(encoding="utf-8"))
    assert doc["status"] == "insufficient_class_balance"
    assert doc["split_mode"] == "not_applicable_single_class"
