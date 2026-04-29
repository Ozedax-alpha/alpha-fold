from __future__ import annotations

import argparse
import json
from pathlib import Path

import avi.cli as cli


def _write_cfg(run_dir: Path, *, base: str = "cftr") -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    cfg = {
        "uniprot_id": "P13569",
        "gene_symbol": "CFTR",
        "clinvar_esearch_term": "CFTR[gene]",
        "output_basename": base,
        "alphafold_fragment": "F1",
    }
    (run_dir / "pipeline_config.json").write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")


def test_report_preflight_fails_when_variants_csv_missing(tmp_path: Path):
    run_dir = tmp_path / "runs" / "cftr" / "20260101_000000"
    _write_cfg(run_dir)

    ns = argparse.Namespace(run_dir=str(run_dir), notebook=None)
    try:
        cli.cmd_report(ns)
    except SystemExit as e:
        msg = str(e)
        assert "Missing required artifact:" in msg
        assert "Run first: py -3 -m avi run --run-dir" in msg
    else:
        raise AssertionError("Expected cmd_report to fail preflight when variants CSV is missing.")


def test_report_uses_custom_notebook_path(tmp_path: Path, monkeypatch):
    run_dir = tmp_path / "runs" / "cftr" / "20260101_000001"
    _write_cfg(run_dir)
    data_proc = run_dir / "data" / "processed"
    data_proc.mkdir(parents=True, exist_ok=True)
    # Satisfy preflight.
    (data_proc / "cftr_variants_basic.csv").write_text("a,b\n1,2\n", encoding="utf-8")

    custom_nb = tmp_path / "custom.ipynb"
    custom_nb.write_text("{}", encoding="utf-8")

    captured: dict[str, object] = {}

    def _fake_run_module(module: str, args: list[str], *, env: dict[str, str] | None):
        captured["module"] = module
        captured["args"] = list(args)
        captured["env"] = dict(env or {})
        return 0

    monkeypatch.setattr(cli, "_run_module", _fake_run_module)

    ns = argparse.Namespace(run_dir=str(run_dir), notebook=str(custom_nb))
    rc = cli.cmd_report(ns)
    assert rc == 0
    assert captured["module"] == "jupyter"
    args = captured["args"]
    assert isinstance(args, list)
    assert str(custom_nb) in args
