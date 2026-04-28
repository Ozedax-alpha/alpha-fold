"""Fast checks that CLIs parse (no network)."""

from __future__ import annotations

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

def test_avi_run_help():
    _run_help(["-m", "avi", "run", "--help"])

def test_avi_explain_help():
    _run_help(["-m", "avi", "explain", "--help"])

def test_avi_doctor_help():
    _run_help(["-m", "avi", "doctor", "--help"])

def test_avi_batch_help():
    _run_help(["-m", "avi", "batch", "--help"])

def test_avi_notebook_help():
    _run_help(["-m", "avi", "notebook", "--help"])

def test_avi_report_help():
    _run_help(["-m", "avi", "report", "--help"])


def test_notebook_file_present():
    nb = ROOT / "notebooks" / "01_tp53_exploration.ipynb"
    assert nb.is_file()
