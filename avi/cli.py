"""Single front-door CLI: ``py -m avi`` (wraps existing scripts/)."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from avi.paths import artifact_paths, default_config_path, format_paths_summary

BASE_DIR = Path(__file__).resolve().parents[1]
SCRIPTS = BASE_DIR / "scripts"


def _run_script(name: str, extra: list[str], *, env: dict[str, str] | None = None) -> int:
    cmd = [sys.executable, str(SCRIPTS / name), *extra]
    print("+", " ".join(cmd))
    return subprocess.call(cmd, cwd=BASE_DIR, env=env)


def _load_json(path: Path) -> dict[str, Any]:
    doc = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(doc, dict):
        raise SystemExit(f"{path} must contain a JSON object.")
    return doc


def _preset_config(preset: str) -> dict[str, str]:
    r = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "set_pipeline_config.py"),
            "--preset",
            preset,
            "--dry-run",
        ],
        cwd=BASE_DIR,
        capture_output=True,
        text=True,
        check=False,
    )
    if r.returncode != 0:
        sys.stderr.write(r.stderr or r.stdout or "set_pipeline_config failed.\n")
        raise SystemExit(r.returncode or 1)
    doc = json.loads(r.stdout.strip())
    if not isinstance(doc, dict):
        raise SystemExit("set_pipeline_config.py --dry-run returned non-object JSON.")
    return {str(k): str(v) for k, v in doc.items()}


def _git_sha_short() -> str | None:
    r = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"],
        cwd=BASE_DIR,
        capture_output=True,
        text=True,
    )
    if r.returncode != 0:
        return None
    return r.stdout.strip() or None


def _pipeline_env_for_run_dir(run_dir: Path) -> dict[str, str]:
    cfg = run_dir / "pipeline_config.json"
    if not cfg.is_file():
        raise SystemExit(f"Missing {cfg}. Run: py -m avi init --preset <name>")
    data_root = run_dir / "data"
    env = dict(os.environ)
    env["PIPELINE_CONFIG_PATH"] = str(cfg.resolve())
    env["PIPELINE_OUTPUT_DIR"] = str(data_root.resolve())
    return env


def cmd_init(ns: argparse.Namespace) -> int:
    ts = datetime.now(UTC).strftime("%Y%m%d_%H%M%S")
    cfg_preview = _preset_config(ns.preset)
    basename = cfg_preview.get("output_basename", "").strip() or ns.preset
    safe = basename.replace("/", "_").replace("\\", "_")

    runs_parent = Path(ns.runs_parent)
    if not runs_parent.is_absolute():
        runs_parent = BASE_DIR / runs_parent
    run_dir = runs_parent / safe / ts
    run_dir.mkdir(parents=True, exist_ok=False)

    cfg_out = run_dir / "pipeline_config.json"
    rc = _run_script(
        "set_pipeline_config.py",
        ["--preset", ns.preset, "--output", str(cfg_out)],
        env=None,
    )
    if rc != 0:
        return rc

    manifest = {
        "created_utc": datetime.now(UTC).replace(microsecond=0).isoformat(),
        "command": "avi init",
        "preset": ns.preset,
        "run_dir": str(run_dir.resolve()),
        "pipeline_config": str(cfg_out.resolve()),
        "git_sha": _git_sha_short(),
    }
    (run_dir / "run.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"Initialized run directory:\n  {run_dir.resolve()}")
    print("\nNext:")
    print(f"  py -m avi run --run-dir {run_dir}")
    return 0


def cmd_run(ns: argparse.Namespace) -> int:
    stages = ns.stage or []
    env: dict[str, str] | None = None
    run_dir: Path | None = None
    if ns.run_dir:
        run_dir = Path(ns.run_dir)
        if not run_dir.is_absolute():
            run_dir = (BASE_DIR / run_dir).resolve()
        env = _pipeline_env_for_run_dir(run_dir)

    started = time.time()
    rc = 0

    if not stages:
        extra: list[str] = []
        if ns.force_download:
            extra.append("--force-download")
        if ns.no_subset:
            extra.append("--no-subset")
        rc = _run_script("run_pipeline.py", extra, env=env)
    else:
        order = ["download", "variants", "subset", "report", "repro"]
        wanted: list[str] = []
        for s in stages:
            if s not in order:
                raise SystemExit(
                    f"Unknown stage {s!r}. Choose from: {', '.join(order)}"
                )
            if s not in wanted:
                wanted.append(s)

        for s in order:
            if s not in wanted:
                continue
            if s == "download":
                dl = ["--force"] if ns.force_download else []
                rc = _run_script("download_tp53_data.py", dl, env=env)
            elif s == "variants":
                rc = _run_script("process_tp53_variants.py", [], env=env)
            elif s == "subset":
                rc = _run_script("build_missense_subset.py", [], env=env)
            elif s == "report":
                rc = _run_script("export_report.py", [], env=env)
            elif s == "repro":
                rc = _run_script("repro_manifest.py", ["--trail"], env=env)
            if rc != 0:
                break

    elapsed = round(time.time() - started, 3)
    if run_dir is not None:
        run_json = run_dir / "run.json"
        prior: dict[str, Any] = {}
        if run_json.is_file():
            try:
                prior = json.loads(run_json.read_text(encoding="utf-8-sig"))
            except json.JSONDecodeError:
                prior = {}
        if not isinstance(prior, dict):
            prior = {}
        prior["last_run"] = {
            "finished_utc": datetime.now(UTC).replace(microsecond=0).isoformat(),
            "elapsed_sec": elapsed,
            "return_code": rc,
            "stages": stages or ["full_via_run_pipeline"],
        }
        run_json.write_text(
            json.dumps(prior, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )

    return rc


def cmd_explain(ns: argparse.Namespace) -> int:
    cfg_path: Path
    data_dir: Path

    if ns.run_dir:
        rd = Path(ns.run_dir)
        if not rd.is_absolute():
            rd = BASE_DIR / rd
        cfg_path = rd / "pipeline_config.json"
        data_dir = rd / "data"
    elif ns.config:
        cfg_path = Path(ns.config)
        if not cfg_path.is_absolute():
            cfg_path = BASE_DIR / cfg_path
        data_dir = Path(ns.data_dir) if ns.data_dir else (BASE_DIR / "data")
        if not data_dir.is_absolute():
            data_dir = BASE_DIR / data_dir
    else:
        cfg_path = default_config_path(BASE_DIR)
        data_dir = BASE_DIR / "data"

    if not cfg_path.is_file():
        raise SystemExit(f"Missing config: {cfg_path}")
    cfg = _load_json(cfg_path)

    print(format_paths_summary(cfg, data_dir))
    print()
    art = artifact_paths(cfg, data_dir)
    for key, p in art.items():
        if key.endswith("_dir"):
            continue
        print(f"  {key}: {p}")
    return 0


def cmd_doctor(_: argparse.Namespace) -> int:
    ok = True
    print(f"Project root: {BASE_DIR}")

    for mod in ("requests", "pandas", "Bio", "numpy"):
        try:
            __import__(mod)
            print(f"  import {mod}: ok")
        except ImportError as e:
            ok = False
            print(f"  import {mod}: FAILED ({e})")

    cfg = default_config_path(BASE_DIR)
    if cfg.is_file():
        print(f"  pipeline_config.json: present ({cfg})")
        try:
            _load_json(cfg)
            print("  pipeline_config.json: valid JSON object")
        except Exception as e:
            ok = False
            print(f"  pipeline_config.json: INVALID ({e})")
    else:
        print("  pipeline_config.json: missing (copy pipeline_config.example.json)")

    email = os.environ.get("ENTREZ_EMAIL")
    key = os.environ.get("NCBI_API_KEY")
    print(
        f"  ENTREZ_EMAIL: {'set' if email else 'unset (optional, recommended for E-utilities)'}"
    )
    print(f"  NCBI_API_KEY: {'set' if key else 'unset (optional)'}")

    runs = BASE_DIR / "runs"
    try:
        runs.mkdir(parents=True, exist_ok=True)
        probe = runs / ".avi_write_test"
        probe.write_text("ok", encoding="utf-8")
        probe.unlink(missing_ok=True)
        print(f"  writable: {runs}")
    except OSError as e:
        ok = False
        print(f"  writable runs/: FAILED ({e})")

    return 0 if ok else 1


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="avi",
        description="AlphaFold x ClinVar pipeline - front-door CLI (wraps scripts/).",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    p_init = sub.add_parser(
        "init",
        help="Create a timestamped run dir + pipeline_config.json from a preset.",
    )
    p_init.add_argument(
        "--preset",
        required=True,
        metavar="NAME",
        help="Preset from scripts/workflow_targets.json (see set_pipeline_config.py --list-presets).",
    )
    p_init.add_argument(
        "--runs-parent",
        default="runs",
        help="Directory under project root to store runs (default: runs).",
    )
    p_init.set_defaults(func=cmd_init)

    p_run = sub.add_parser(
        "run",
        help="Run pipeline stages (default: full download->variants->subset like scripts/run_pipeline.py).",
    )
    p_run.add_argument(
        "--run-dir",
        type=str,
        default=None,
        help="Run directory containing pipeline_config.json; outputs go to <run-dir>/data/.",
    )
    p_run.add_argument(
        "--stage",
        action="append",
        choices=("download", "variants", "subset", "report", "repro"),
        help="Run only these stages (repeat flag for multiple). Default: full pipeline.",
    )
    p_run.add_argument(
        "--force-download",
        action="store_true",
        help="Force refresh of raw AlphaFold / ClinVar downloads.",
    )
    p_run.add_argument(
        "--no-subset",
        action="store_true",
        help="When running full pipeline, skip build_missense_subset.py.",
    )
    p_run.set_defaults(func=cmd_run)

    p_batch = sub.add_parser(
        "batch",
        help="Run many targets (forwards args to scripts/run_batch.py).",
    )
    p_batch.add_argument(
        "args",
        nargs=argparse.REMAINDER,
        help="Arguments forwarded to scripts/run_batch.py (use '--' if needed).",
    )
    p_batch.set_defaults(func=None)

    p_explain = sub.add_parser(
        "explain",
        help="Show resolved artifact filenames/paths for a config.",
    )
    g = p_explain.add_mutually_exclusive_group()
    g.add_argument("--run-dir", type=str, default=None, help="Run dir from avi init.")
    g.add_argument("--config", type=str, default=None, help="Path to pipeline_config.json.")
    p_explain.add_argument(
        "--data-dir",
        type=str,
        default=None,
        help="With --config only: directory containing raw/ and processed/ (default: ./data).",
    )
    p_explain.set_defaults(func=cmd_explain)

    p_doc = sub.add_parser(
        "doctor",
        help="Preflight checks (imports, config JSON, writable dirs).",
    )
    p_doc.set_defaults(func=cmd_doctor)

    return p


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    ns = parser.parse_args(argv)

    if ns.cmd == "batch":
        cmd = [sys.executable, str(SCRIPTS / "run_batch.py"), *list(ns.args)]
        print("+", " ".join(cmd))
        return subprocess.call(cmd, cwd=BASE_DIR)

    assert ns.func is not None
    return int(ns.func(ns))


if __name__ == "__main__":
    raise SystemExit(main())

