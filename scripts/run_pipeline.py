"""
Run ingest → variant table → missense subset (or individual stages).

From project root:
  python scripts/run_pipeline.py
  python scripts/run_pipeline.py --download-only
  python scripts/run_pipeline.py --process-only
  python scripts/run_pipeline.py --subset-only
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parents[1]


def _run(script_name: str, extra_args: list[str] | None = None) -> None:
    script = BASE_DIR / "scripts" / script_name
    subprocess.check_call(
        [sys.executable, str(script), *(extra_args or [])],
        cwd=BASE_DIR,
    )


def main() -> None:
    p = argparse.ArgumentParser(description="AlphaFold + ClinVar pipeline")
    p.add_argument(
        "--download-only",
        action="store_true",
        help="Only run download_tp53_data.py",
    )
    p.add_argument(
        "--process-only",
        action="store_true",
        help="Run process_tp53_variants.py + build_missense_subset.py (needs raw data)",
    )
    p.add_argument(
        "--subset-only",
        action="store_true",
        help="Only run build_missense_subset.py (needs variants_basic CSV)",
    )
    p.add_argument(
        "--no-subset",
        action="store_true",
        help="With default run, skip build_missense_subset.py",
    )
    p.add_argument(
        "--force-download",
        action="store_true",
        help="Pass --force to download_tp53_data.py (redownload raw PDB/ClinVar JSON)",
    )
    args = p.parse_args()

    flags = sum(
        1
        for x in (
            args.download_only,
            args.process_only,
            args.subset_only,
        )
        if x
    )
    if flags > 1:
        p.error("Choose at most one of --download-only / --process-only / --subset-only")

    dl_extra = ["--force"] if args.force_download else []

    if args.download_only:
        _run("download_tp53_data.py", dl_extra)
        return
    if args.subset_only:
        _run("build_missense_subset.py")
        return
    if args.process_only:
        _run("process_tp53_variants.py")
        _run("build_missense_subset.py")
        print("Process + subset finished.")
        return

    _run("download_tp53_data.py", dl_extra)
    _run("process_tp53_variants.py")
    if not args.no_subset:
        _run("build_missense_subset.py")
    print(
        "Pipeline finished: data/raw/, data/processed/<basename>_variants_basic.csv, "
        "data/processed/<basename>_missense_mappable.csv"
    )


if __name__ == "__main__":
    main()
