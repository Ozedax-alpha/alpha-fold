"""Fast checks that CLIs parse (no network)."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _run_help(args: list[str]) -> None:
    subprocess.check_call(
        [sys.executable, *args],
        cwd=ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def test_run_pipeline_help():
    _run_help(["scripts/run_pipeline.py", "--help"])


def test_download_help():
    _run_help(["scripts/download_tp53_data.py", "--help"])


def test_avi_help():
    _run_help(["-m", "avi", "--help"])

def test_avi_init_help():
    _run_help(["-m", "avi", "init", "--help"])

def test_avi_add_preset_help():
    _run_help(["-m", "avi", "add-preset", "--help"])

def test_avi_run_help():
    _run_help(["-m", "avi", "run", "--help"])

def test_avi_explain_help():
    _run_help(["-m", "avi", "explain", "--help"])

def test_avi_doctor_help():
    _run_help(["-m", "avi", "doctor", "--help"])

def test_avi_repro_help():
    _run_help(["-m", "avi", "repro", "--help"])

def test_avi_clean_help():
    _run_help(["-m", "avi", "clean", "--help"])

def test_avi_batch_help():
    _run_help(["-m", "avi", "batch", "--help"])

def test_avi_notebook_help():
    _run_help(["-m", "avi", "notebook", "--help"])

def test_avi_report_help():
    _run_help(["-m", "avi", "report", "--help"])

def test_target_notebook_file_present():
    nb = ROOT / "notebooks" / "01_target_exploration.ipynb"
    assert nb.is_file()

def test_avi_open_report_help():
    _run_help(["-m", "avi", "open-report", "--help"])

def test_avi_list_runs_help():
    _run_help(["-m", "avi", "list-runs", "--help"])

def test_avi_gc_help():
    _run_help(["-m", "avi", "gc", "--help"])


def test_notebook_file_present():
    nb = ROOT / "notebooks" / "01_tp53_exploration.ipynb"
    assert nb.is_file()


def _write_min_cfg(run_dir: Path, *, uid: str = "P13569", base: str = "cftr", frag: str = "F1") -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    cfg = {
        "uniprot_id": uid,
        "gene_symbol": base.upper(),
        "clinvar_esearch_term": f"{base.upper()}[gene]",
        "output_basename": base,
        "alphafold_fragment": frag,
    }
    (run_dir / "pipeline_config.json").write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
    (run_dir / "run.json").write_text(json.dumps({"created_utc": "test"}, indent=2) + "\n", encoding="utf-8")


def test_avi_run_resume_skips_completed_stages(tmp_path: Path):
    # Create a fake run directory with completed artifacts so --resume can skip without executing scripts.
    run_dir = tmp_path / "runs" / "cftr" / "20260101_000000"
    _write_min_cfg(run_dir)

    data_root = run_dir / "data"
    (data_root / "raw").mkdir(parents=True, exist_ok=True)
    (data_root / "processed").mkdir(parents=True, exist_ok=True)

    # Artifacts for download + variants already done.
    (data_root / "raw" / "AF-P13569-F1-alphafold.pdb").write_text("pdb", encoding="utf-8")
    (data_root / "raw" / "clinvar_cftr_variants.json").write_text("{}", encoding="utf-8")
    (data_root / "processed" / "cftr_variants_basic.csv").write_text("a,b\n1,2\n", encoding="utf-8")

    # Ask to run stages, but --resume should skip them due to existing artifacts.
    subprocess.check_call(
        [sys.executable, "-m", "avi", "run", "--run-dir", str(run_dir), "--resume", "--stage", "download", "--stage", "variants"],
        cwd=ROOT,
    )

    state = run_dir / "run_state.json"
    assert state.is_file()
    doc = json.loads(state.read_text(encoding="utf-8"))
    assert doc["stages"]["download"]["skipped"] is True
    assert doc["stages"]["variants"]["skipped"] is True


def test_avi_clean_accepts_dry_run_flag():
    # Just verify parsing/exit code; actual deletions are guarded by --yes.
    subprocess.check_call(
        [sys.executable, "-m", "avi", "clean", "--keep", "1", "--dry-run"],
        cwd=ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
