"""Single front-door CLI: ``py -m avi`` (wraps existing scripts/)."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
import webbrowser
import hashlib
import shutil
import re
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from avi.paths import artifact_paths, default_config_path, format_paths_summary
from avi.pipeline import STAGE_ORDER, StageContext, run_stages
from avi import db as avi_db
from avi import webui as avi_webui

try:
    import requests  # type: ignore
except Exception:  # pragma: no cover
    requests = None  # type: ignore

BASE_DIR = Path(__file__).resolve().parents[1]
SCRIPTS = BASE_DIR / "scripts"
NOTEBOOK_PATH = BASE_DIR / "notebooks" / "01_target_exploration.ipynb"
PRESETS_PATH = Path(__file__).with_name("presets.json")
DEFAULT_DB_PATH = avi_db.default_db_path(BASE_DIR)


def _run_script(name: str, extra: list[str], *, env: dict[str, str] | None = None) -> int:
    cmd = [sys.executable, str(SCRIPTS / name), *extra]
    print("+", " ".join(cmd))
    return subprocess.call(cmd, cwd=BASE_DIR, env=env)


def _load_json(path: Path) -> dict[str, Any]:
    doc = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(doc, dict):
        raise SystemExit(f"{path} must contain a JSON object.")
    return doc


def _load_presets() -> dict[str, dict[str, Any]]:
    if not PRESETS_PATH.is_file():
        return {}
    try:
        doc = json.loads(PRESETS_PATH.read_text(encoding="utf-8-sig"))
    except Exception as e:
        raise SystemExit(f"Failed to read {PRESETS_PATH}: {e}")
    if not isinstance(doc, dict):
        raise SystemExit(f"{PRESETS_PATH} must be a JSON object.")
    out: dict[str, dict[str, Any]] = {}
    for k, v in doc.items():
        if isinstance(k, str) and isinstance(v, dict):
            out[k] = dict(v)
    return out


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)


def _save_presets(presets: dict[str, dict[str, Any]]) -> None:
    blob = json.dumps(presets, indent=2, ensure_ascii=False) + "\n"
    _atomic_write_text(PRESETS_PATH, blob)


def _slugify_basename(s: str) -> str:
    t = (s or "").strip().lower()
    t = re.sub(r"\s+", "_", t)
    t = re.sub(r"[^a-z0-9_]+", "", t)
    return t or "output"


def _fetch_uniprot_gene_symbol(uniprot_id: str) -> str | None:
    """Return primary gene symbol for a UniProt accession, if available."""
    if requests is None:
        raise SystemExit("Missing dependency: requests (install requirements.txt).")
    uid = str(uniprot_id).strip()
    if not uid:
        return None
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.json"
    r = requests.get(url, timeout=30)
    if r.status_code == 404:
        raise SystemExit(f"UniProt accession not found: {uid}")
    r.raise_for_status()
    doc = r.json()
    genes = doc.get("genes") if isinstance(doc, dict) else None
    if not isinstance(genes, list):
        return None
    for g in genes:
        if not isinstance(g, dict):
            continue
        gn = g.get("geneName")
        if isinstance(gn, dict) and isinstance(gn.get("value"), str) and gn["value"].strip():
            return gn["value"].strip()
    return None


def _preset_config(preset: str) -> dict[str, Any]:
    key = str(preset).strip()
    if not key:
        raise SystemExit("Empty preset name.")
    doc = _load_presets()
    entry = doc.get(key)
    if not isinstance(entry, dict):
        avail = ", ".join(sorted(doc))
        raise SystemExit(f"Unknown preset {key!r}. Known: {avail}")
    out: dict[str, Any] = {str(k): v for k, v in entry.items()}
    required = (
        "uniprot_id",
        "gene_symbol",
        "clinvar_esearch_term",
        "output_basename",
        "alphafold_fragment",
    )
    missing = [k for k in required if k not in out or out[k] is None or not str(out[k]).strip()]
    if missing:
        raise SystemExit(f"Preset {key!r} missing keys: {missing}")
    return out


def _write_config(cfg: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cfg, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


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


def _sha256_file(path: Path) -> str | None:
    if not path.is_file():
        return None
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_run_manifest(*, run_dir: Path) -> Path:
    cfg_path = run_dir / "pipeline_config.json"
    if not cfg_path.is_file():
        raise SystemExit(f"Missing config: {cfg_path}")
    cfg = _load_json(cfg_path)
    base = str(cfg.get("output_basename") or "run").replace("/", "_").replace("\\", "_")
    data_root = run_dir / "data"
    raw = data_root / "raw"
    proc = data_root / "processed"
    proc.mkdir(parents=True, exist_ok=True)
    art = artifact_paths(cfg, data_root)

    doc = {
        "generated_at_utc": datetime.now(UTC).replace(microsecond=0).isoformat(),
        "run_dir": str(run_dir.resolve()),
        "pipeline_config": str(cfg_path.resolve()),
        "config_effective": cfg,
        "raw_files": [
            {
                "path": str(art["alphafold_pdb"]),
                "exists": art["alphafold_pdb"].is_file(),
                "sha256": _sha256_file(art["alphafold_pdb"]),
            },
            {
                "path": str(art["clinvar_json"]),
                "exists": art["clinvar_json"].is_file(),
                "sha256": _sha256_file(art["clinvar_json"]),
            },
        ],
    }
    out = proc / f"repro_manifest_{base}.json"
    out.write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return out


def _run_module(module: str, args: list[str], *, env: dict[str, str] | None) -> int:
    cmd = [sys.executable, "-m", module, *args]
    print("+", " ".join(cmd))
    return subprocess.call(cmd, cwd=BASE_DIR, env=env)


def _export_html_report(*, run_dir: Path, env: dict[str, str], notebook_path: Path | None = None) -> int:
    out_dir = run_dir / "reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    nb = notebook_path or NOTEBOOK_PATH
    if not nb.is_file():
        raise SystemExit(f"Missing notebook: {nb}")

    # Preflight: the default notebook expects variants CSV to exist.
    cfg_path = run_dir / "pipeline_config.json"
    if cfg_path.is_file():
        try:
            cfg = _load_json(cfg_path)
            art = artifact_paths(cfg, run_dir / "data")
            vb = art["variants_basic_csv"]
            if not vb.is_file():
                raise SystemExit(
                    f"Missing required artifact: {vb}\nRun first: py -3 -m avi run --run-dir {run_dir}"
                )
        except SystemExit:
            raise
        except Exception:
            # If config parsing fails, fall back to letting nbconvert show an error.
            pass
    return _run_module(
        "jupyter",
        [
            "nbconvert",
            "--to",
            "html",
            "--execute",
            str(nb),
            "--output",
            "exploration_report",
            "--output-dir",
            str(out_dir),
        ],
        env=env,
    )


def _report_html_path(run_dir: Path) -> Path:
    return run_dir / "reports" / "exploration_report.html"


def _runs_root() -> Path:
    return BASE_DIR / "runs"


def _iter_run_dirs(runs_root: Path) -> list[Path]:
    if not runs_root.is_dir():
        return []
    out: list[Path] = []
    for target_dir in runs_root.iterdir():
        if not target_dir.is_dir():
            continue
        for run_dir in target_dir.iterdir():
            if run_dir.is_dir() and (run_dir / "pipeline_config.json").is_file():
                out.append(run_dir)
    out.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return out


def cmd_list_runs(ns: argparse.Namespace) -> int:
    runs_root = Path(ns.runs_root) if ns.runs_root else _runs_root()
    if not runs_root.is_absolute():
        runs_root = (BASE_DIR / runs_root).resolve()
    limit = int(ns.limit)

    runs = _iter_run_dirs(runs_root)[:limit]
    if not runs:
        print(f"No runs found under {runs_root}")
        return 0

    for rd in runs:
        cfg = rd / "pipeline_config.json"
        run_json = rd / "run.json"
        out = rd / "data" / "processed"
        vb = None
        mm = None
        try:
            cdoc = json.loads(cfg.read_text(encoding="utf-8-sig"))
            if isinstance(cdoc, dict):
                vb = cdoc.get("output_basename")
        except Exception:
            vb = None
        mm_path = out / (f"{vb}_missense_mappable.csv" if vb else "UNKNOWN")
        vb_path = out / (f"{vb}_variants_basic.csv" if vb else "UNKNOWN")
        mm = "yes" if vb and mm_path.is_file() else "no"
        vb_ok = "yes" if vb and vb_path.is_file() else "no"
        last_rc = "?"
        if run_json.is_file():
            try:
                rdoc = json.loads(run_json.read_text(encoding="utf-8-sig"))
                if isinstance(rdoc, dict):
                    lr = (rdoc.get("last_run") or {}).get("return_code")
                    if lr is not None:
                        last_rc = str(lr)
            except Exception:
                pass
        rel = rd.relative_to(BASE_DIR) if rd.is_absolute() else rd
        print(f"{rel}  basename={vb or '?'}  processed={{basic:{vb_ok},subset:{mm}}}  last_rc={last_rc}")
    return 0


def cmd_add_preset(ns: argparse.Namespace) -> int:
    presets = _load_presets()

    key = str(ns.key).strip()
    if not key:
        raise SystemExit("Preset key is required.")

    overwrite = bool(ns.overwrite)
    if key in presets and not overwrite:
        raise SystemExit(f"Preset {key!r} already exists. Use --overwrite to replace.")

    def ask(prompt: str, default: str | None = None) -> str:
        if ns.non_interactive:
            if default is None or not str(default).strip():
                raise SystemExit(
                    f"Missing required field for non-interactive mode: {prompt}"
                )
            return str(default)
        if getattr(ns, "auto", False):
            if default is None or not str(default).strip():
                raise SystemExit(f"Missing required field for auto mode: {prompt}")
            return str(default)
        suffix = f" [{default}]" if default else ""
        v = input(f"{prompt}{suffix}: ").strip()
        return v or (default or "")

    from_uid = str(ns.from_uniprot).strip() if ns.from_uniprot else ""
    auto_gene = None
    if from_uid:
        auto_gene = _fetch_uniprot_gene_symbol(from_uid)

    uniprot_id = ns.uniprot_id or from_uid or ask("UniProt accession", None)
    gene_symbol = ns.gene_symbol or auto_gene or ask("Gene symbol", None)
    default_term = f"{gene_symbol}[gene]" if gene_symbol else None
    clinvar_term = ns.clinvar_esearch_term or ask(
        "ClinVar esearch term (e.g. TP53[gene])", default_term
    )
    default_base = _slugify_basename(gene_symbol) if gene_symbol else key
    output_basename = ns.output_basename or ask("Output basename (e.g. tp53)", default_base)
    alphafold_fragment = ns.alphafold_fragment or ask("AlphaFold fragment (usually F1)", "F1")

    entry = {
        "uniprot_id": str(uniprot_id).strip(),
        "gene_symbol": str(gene_symbol).strip(),
        "clinvar_esearch_term": str(clinvar_term).strip(),
        "output_basename": str(output_basename).strip(),
        "alphafold_fragment": str(alphafold_fragment).strip(),
    }
    missing = [k for k, v in entry.items() if not v]
    if missing:
        raise SystemExit(f"Missing fields: {', '.join(missing)}")

    presets[key] = entry
    _save_presets(presets)
    print(f"Wrote preset {key!r} to {PRESETS_PATH}")
    return 0

def cmd_init(ns: argparse.Namespace) -> int:
    ts = datetime.now(UTC).strftime("%Y%m%d_%H%M%S")
    if ns.preset:
        cfg_preview = _preset_config(ns.preset)
        run_key = str(ns.preset)
    else:
        uid = str(ns.from_uniprot or "").strip()
        if not uid:
            raise SystemExit("Specify --preset or --from-uniprot.")
        gene = str(ns.gene_symbol or "").strip() or (_fetch_uniprot_gene_symbol(uid) or "")
        if not gene:
            raise SystemExit(f"Could not resolve gene symbol for {uid}.")
        base = str(ns.output_basename or "").strip() or _slugify_basename(gene)
        frag = str(ns.alphafold_fragment or "F1").strip()
        term = str(ns.clinvar_esearch_term or f"{gene}[gene]").strip()
        cfg_preview = {
            "uniprot_id": uid,
            "gene_symbol": gene,
            "clinvar_esearch_term": term,
            "output_basename": base,
            "alphafold_fragment": frag,
        }
        run_key = str(ns.key or base)
    basename = cfg_preview.get("output_basename", "").strip() or run_key
    safe = basename.replace("/", "_").replace("\\", "_")

    runs_parent = Path(ns.runs_parent)
    if not runs_parent.is_absolute():
        runs_parent = BASE_DIR / runs_parent
    run_dir = runs_parent / safe / ts
    run_dir.mkdir(parents=True, exist_ok=False)

    cfg_out = run_dir / "pipeline_config.json"
    _write_config(cfg_preview, cfg_out)

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
    run_abs = str(run_dir.resolve())
    print("  # PowerShell: copy/paste")
    print(f"  $run = {run_abs!r}")
    print("  py -3 -m avi run --run-dir $run --force-download")
    print("  py -3 -m avi explain --run-dir $run")
    print("  py -3 -m avi notebook --run-dir $run")
    print("  py -3 -m avi report --run-dir $run")

    if getattr(ns, "run", False):
        env = _pipeline_env_for_run_dir(run_dir)
        extra: list[str] = []
        if getattr(ns, "force_download", False):
            extra.append("--force-download")
        rc = _run_script("run_pipeline.py", extra, env=env)
        if rc != 0:
            return rc

        # Optional post-actions (in order).
        if getattr(ns, "report", False):
            rc = _export_html_report(run_dir=run_dir, env=env)
            if rc != 0:
                return rc
        if getattr(ns, "open_report", False):
            # Requires report to exist; if user didn't request --report, try to open anyway.
            cmd_open_report(argparse.Namespace(run_dir=str(run_dir)))
        if getattr(ns, "notebook", False):
            return cmd_notebook(argparse.Namespace(run_dir=str(run_dir)))

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

    # Default behavior stays aligned with scripts/run_pipeline.py unless --resume
    # is requested (resume always runs step-by-step).
    if not stages and not getattr(ns, "resume", False):
        extra: list[str] = []
        if ns.force_download:
            extra.append("--force-download")
        if ns.no_subset:
            extra.append("--no-subset")
        rc = _run_script("run_pipeline.py", extra, env=env)
    else:
        wanted: list[str]
        if stages:
            wanted = list(stages)
        else:
            wanted = ["download", "variants"]
            if not ns.no_subset:
                wanted.append("subset")

        cfg: dict[str, Any] | None = None
        if run_dir is not None:
            cfg_path = run_dir / "pipeline_config.json"
            if cfg_path.is_file():
                cfg = _load_json(cfg_path)

        def _runner(stage: str) -> int:
            nonlocal run_dir, env
            if stage == "download":
                dl = ["--force"] if ns.force_download else []
                return _run_script("download_tp53_data.py", dl, env=env)
            if stage == "variants":
                return _run_script("process_tp53_variants.py", [], env=env)
            if stage == "subset":
                return _run_script("build_missense_subset.py", [], env=env)
            if stage == "report":
                if run_dir is None or env is None:
                    raise SystemExit("Stage 'report' requires --run-dir.")
                return _export_html_report(run_dir=run_dir, env=env)
            if stage == "repro":
                if run_dir is None:
                    raise SystemExit("Stage 'repro' requires --run-dir.")
                try:
                    out = _write_run_manifest(run_dir=run_dir)
                    print(f"Wrote {out}")
                    return 0
                except Exception as e:
                    print(f"Manifest error: {e}", file=sys.stderr)
                    return 1
            raise SystemExit(
                f"Unknown stage {stage!r}. Choose from: {', '.join(STAGE_ORDER)}"
            )

        if run_dir is None:
            raise SystemExit("Stage execution requires --run-dir.")
        ctx = StageContext(
            run_dir=run_dir,
            env=env,
            cfg=cfg,
            resume=bool(getattr(ns, 'resume', False)),
            report_html=_report_html_path(run_dir),
        )
        rc = run_stages(
            ctx=ctx,
            stages=wanted,
            stage_runner=_runner,
            atomic_write_text=_atomic_write_text,
        )

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
        if rc == 0 and not getattr(ns, "no_manifest", False):
            try:
                # If --resume, avoid rewriting a manifest that's already present.
                cfg_path = run_dir / "pipeline_config.json"
                cfg = _load_json(cfg_path) if cfg_path.is_file() else {}
                art = artifact_paths(cfg, run_dir / "data") if isinstance(cfg, dict) else {}
                existing = art.get("repro_manifest") if isinstance(art, dict) else None
                if getattr(ns, "resume", False) and isinstance(existing, Path) and existing.is_file():
                    print(f"= Repro manifest already present: {existing.name}")
                else:
                    out = _write_run_manifest(run_dir=run_dir)
                    print(f"Wrote {out}")
            except Exception as e:
                print(f"Manifest error: {e}", file=sys.stderr)

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
    hc = os.environ.get("AVI_USE_HTTP_CACHE", "")
    print(
        f"  AVI_USE_HTTP_CACHE: {hc!r} (set to 1 + pip install '.[cache]' for on-disk HTTP cache)"
    )

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

    # Jupyter availability
    try:
        __import__("jupyter")
        print("  jupyter: ok (python -m jupyter ...)")
    except ImportError:
        print("  jupyter: missing (install requirements.txt)")

    return 0 if ok else 1


def cmd_repro(ns: argparse.Namespace) -> int:
    run_dir = Path(ns.run_dir)
    if not run_dir.is_absolute():
        run_dir = (BASE_DIR / run_dir).resolve()
    out = _write_run_manifest(run_dir=run_dir)
    print(f"Wrote {out}")
    return 0


def cmd_clean(ns: argparse.Namespace) -> int:
    runs_root = Path(ns.runs_root) if ns.runs_root else _runs_root()
    if not runs_root.is_absolute():
        runs_root = (BASE_DIR / runs_root).resolve()
    keep = int(ns.keep)
    dry = bool(getattr(ns, "dry_run", False)) or not bool(ns.yes)

    runs = _iter_run_dirs(runs_root)
    victims = runs[keep:]
    if not victims:
        print("Nothing to clean.")
        return 0

    for rd in victims:
        rel = rd.relative_to(BASE_DIR) if rd.is_absolute() else rd
        if dry:
            print(f"DRY-RUN delete: {rel}")
        else:
            shutil.rmtree(rd)
            print(f"Deleted: {rel}")
    if dry:
        print("Dry run only. Re-run with --yes to delete.")
    return 0


def cmd_gc(ns: argparse.Namespace) -> int:
    runs_root = Path(ns.runs_root) if ns.runs_root else _runs_root()
    if not runs_root.is_absolute():
        runs_root = (BASE_DIR / runs_root).resolve()

    dry = bool(getattr(ns, "dry_run", False)) or not bool(ns.yes)
    older_days = int(ns.older_than_days) if getattr(ns, "older_than_days", None) is not None else None
    now = time.time()

    targets: list[Path] = []
    for target_dir in runs_root.iterdir() if runs_root.is_dir() else []:
        if not target_dir.is_dir():
            continue
        for run_dir in target_dir.iterdir():
            if not run_dir.is_dir():
                continue

            cfg = run_dir / "pipeline_config.json"
            if getattr(ns, "delete_empty", False) and not cfg.is_file():
                targets.append(run_dir)
                continue

            if older_days is not None:
                age_days = (now - run_dir.stat().st_mtime) / (60 * 60 * 24)
                if age_days < older_days:
                    continue

            if getattr(ns, "delete_partial", False):
                if not cfg.is_file():
                    continue
                try:
                    cdoc = json.loads(cfg.read_text(encoding="utf-8-sig"))
                except Exception:
                    continue
                base = cdoc.get("output_basename") if isinstance(cdoc, dict) else None
                if not base:
                    continue
                vb = run_dir / "data" / "processed" / f"{base}_variants_basic.csv"
                if not vb.is_file():
                    targets.append(run_dir)

    # De-dup and sort oldest first (so output is stable).
    uniq: dict[str, Path] = {str(p.resolve()): p for p in targets}
    victims = sorted(uniq.values(), key=lambda p: p.stat().st_mtime)

    if not victims:
        print("Nothing to garbage-collect.")
        return 0

    for rd in victims:
        rel = rd.relative_to(BASE_DIR) if rd.is_absolute() else rd
        if dry:
            print(f"DRY-RUN delete: {rel}")
        else:
            shutil.rmtree(rd)
            print(f"Deleted: {rel}")
    if dry:
        print("Dry run only. Re-run with --yes to delete.")
    return 0


def cmd_batch(ns: argparse.Namespace) -> int:
    presets = _load_presets()

    from_uids: list[str] | None = None
    if getattr(ns, "from_uniprot", None):
        from_uids = [p.strip() for p in str(ns.from_uniprot).split(",") if p.strip()]
        if not from_uids:
            raise SystemExit("No UniProt accessions provided.")

    wanted: list[str] = []
    if from_uids is None:
        if not presets:
            raise SystemExit(f"No presets found in {PRESETS_PATH}")
        if ns.presets:
            wanted = [p.strip() for p in ns.presets.split(",") if p.strip()]
        else:
            wanted = sorted(presets)
        if ns.limit is not None:
            wanted = wanted[: int(ns.limit)]
        missing = [k for k in wanted if k not in presets]
        if missing:
            raise SystemExit(f"Unknown preset(s): {', '.join(missing)}")
    else:
        wanted = from_uids
        if ns.limit is not None:
            wanted = wanted[: int(ns.limit)]

    batch_ts = datetime.now(UTC).strftime("%Y%m%d_%H%M%S")
    runs_parent = Path(ns.runs_parent)
    if not runs_parent.is_absolute():
        runs_parent = BASE_DIR / runs_parent
    batch_out = runs_parent / f"batch_status_{batch_ts}.json"

    status: dict[str, Any] = {
        "started_utc": datetime.now(UTC).replace(microsecond=0).isoformat(),
        "runs_parent": str(runs_parent.resolve()),
        "presets": wanted,
        "results": [],
    }

    for idx, key in enumerate(wanted):
        run_ts = datetime.now(UTC).strftime("%Y%m%d_%H%M%S")
        if from_uids is None:
            cfg_preview = _preset_config(key)
        else:
            uid = key
            gene = _fetch_uniprot_gene_symbol(uid) or ""
            if not gene:
                raise SystemExit(f"Could not resolve gene symbol for {uid}.")
            base = _slugify_basename(gene)
            cfg_preview = {
                "uniprot_id": uid,
                "gene_symbol": gene,
                "clinvar_esearch_term": f"{gene}[gene]",
                "output_basename": base,
                "alphafold_fragment": "F1",
            }
        basename = str(cfg_preview.get("output_basename", "")).strip() or str(key)
        safe = basename.replace("/", "_").replace("\\", "_")
        # Ensure unique run timestamps (batch can run faster than 1/sec).
        suffix = f"_{idx:03d}"
        run_dir = runs_parent / safe / f"{run_ts}{suffix}"
        run_dir.mkdir(parents=True, exist_ok=True)
        _write_config(cfg_preview, run_dir / "pipeline_config.json")
        (run_dir / "run.json").write_text(
            json.dumps(
                {
                    "created_utc": datetime.now(UTC).replace(microsecond=0).isoformat(),
                    "command": "avi batch",
                    "preset": key,
                    "run_dir": str(run_dir.resolve()),
                },
                indent=2,
                ensure_ascii=False,
            )
            + "\n",
            encoding="utf-8",
        )

        t0 = time.time()
        rc = 0
        err: str | None = None
        try:
            env = _pipeline_env_for_run_dir(run_dir)
            extra: list[str] = []
            if ns.force_download:
                extra.append("--force-download")
            rc = _run_script("run_pipeline.py", extra, env=env)
            if rc == 0:
                _write_run_manifest(run_dir=run_dir)
                if getattr(ns, "report", False):
                    try:
                        _export_html_report(run_dir=run_dir, env=env, notebook_path=None)
                    except Exception as e:
                        rc = 1
                        err = f"ReportError: {e}"
                if rc == 0 and getattr(ns, "evaluate", False):
                    try:
                        ev_args = [
                            "--run-dir",
                            str(run_dir),
                            "--seed",
                            str(getattr(ns, "evaluate_seed", 42)),
                            "--min-label-rows",
                            str(getattr(ns, "evaluate_min_label_rows", 40)),
                            "--min-dated-rows",
                            str(getattr(ns, "evaluate_min_dated_rows", 50)),
                        ]
                        rc = _run_script("evaluation_metrics.py", ev_args, env=None)
                    except Exception as e:
                        rc = 1
                        err = f"EvaluateError: {e}"
        except Exception as e:
            rc = 1
            err = f"{type(e).__name__}: {e}"

        if rc != 0 and not err:
            err = "run_pipeline.py exited with non-zero status (see logged commands above)."

        status["results"].append(
            {
                "preset": key if from_uids is None else None,
                "uniprot_id": key if from_uids is not None else None,
                "run_dir": str(run_dir.resolve()),
                "return_code": int(rc),
                "error": err,
                "elapsed_sec": round(time.time() - t0, 3),
            }
        )
        if rc != 0 and not ns.continue_on_error:
            break

    status["finished_utc"] = datetime.now(UTC).replace(microsecond=0).isoformat()
    batch_out.parent.mkdir(parents=True, exist_ok=True)
    batch_out.write_text(json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {batch_out}")

    failed = [r for r in status["results"] if int(r.get("return_code") or 0) != 0]
    if failed:
        print(f"Batch finished with {len(failed)} failure(s) (see {batch_out.name}):")
        for r in failed:
            label = r.get("preset") or r.get("uniprot_id") or "?"
            err = r.get("error") or ""
            rd = r.get("run_dir") or ""
            print(f"  {label!r} rc={r.get('return_code')} {rd}")
            if err:
                print(f"    {err}")
            if rd:
                print(f"    Re-run: py -3 -m avi run --run-dir {rd!r}")
    if any(r["return_code"] != 0 for r in status["results"]):
        raise SystemExit(1)
    return 0


def cmd_notebook(ns: argparse.Namespace) -> int:
    run_dir = Path(ns.run_dir)
    if not run_dir.is_absolute():
        run_dir = (BASE_DIR / run_dir).resolve()
    env = _pipeline_env_for_run_dir(run_dir)
    # Start server at repo root so users can browse notebooks/; do not pass a
    # notebook path here (avoids root-path confusion/404s).
    return _run_module(
        "jupyter",
        ["notebook", "--notebook-dir", str(BASE_DIR)],
        env=env,
    )


def cmd_report(ns: argparse.Namespace) -> int:
    run_dir = Path(ns.run_dir)
    if not run_dir.is_absolute():
        run_dir = (BASE_DIR / run_dir).resolve()
    env = _pipeline_env_for_run_dir(run_dir)
    nb = Path(ns.notebook) if getattr(ns, "notebook", None) else None
    if nb is not None and not nb.is_absolute():
        nb = (BASE_DIR / nb).resolve()
    return _export_html_report(run_dir=run_dir, env=env, notebook_path=nb)


def cmd_open_report(ns: argparse.Namespace) -> int:
    run_dir = Path(ns.run_dir)
    if not run_dir.is_absolute():
        run_dir = (BASE_DIR / run_dir).resolve()
    html = _report_html_path(run_dir)
    if not html.is_file():
        raise SystemExit(
            f"Missing report: {html}\nRun: py -3 -m avi report --run-dir {run_dir}"
        )
    webbrowser.open(html.resolve().as_uri())
    print(f"Opened: {html}")
    return 0


def cmd_evaluate(ns: argparse.Namespace) -> int:
    run_dir = Path(ns.run_dir)
    if not run_dir.is_absolute():
        run_dir = (BASE_DIR / run_dir).resolve()
    if not (run_dir / "pipeline_config.json").is_file():
        raise SystemExit(f"Missing pipeline_config.json under {run_dir}")
    extra = [
        "--run-dir",
        str(run_dir),
        "--seed",
        str(ns.seed),
        "--min-label-rows",
        str(ns.min_label_rows),
        "--min-dated-rows",
        str(ns.min_dated_rows),
    ]
    return _run_script("evaluation_metrics.py", extra, env=None)


def cmd_db_init(ns: argparse.Namespace) -> int:
    db_path = Path(ns.db)
    if not db_path.is_absolute():
        db_path = (BASE_DIR / db_path).resolve()
    con = avi_db.connect(db_path)
    try:
        avi_db.init_db(con)
    finally:
        con.close()
    print(f"Initialized DB: {db_path}")
    return 0


def cmd_db_import_presets(ns: argparse.Namespace) -> int:
    db_path = Path(ns.db)
    if not db_path.is_absolute():
        db_path = (BASE_DIR / db_path).resolve()
    presets = _load_presets()
    if not presets:
        raise SystemExit(f"No presets found in {PRESETS_PATH}")
    con = avi_db.connect(db_path)
    try:
        avi_db.init_db(con)
        n = 0
        for key, cfg in presets.items():
            avi_db.upsert_protein(
                con,
                preset_key=key,
                uniprot_id=str(cfg.get("uniprot_id") or "").strip(),
                gene_symbol=str(cfg.get("gene_symbol") or "").strip() or None,
                clinvar_esearch_term=str(cfg.get("clinvar_esearch_term") or "").strip() or None,
                output_basename=str(cfg.get("output_basename") or "").strip() or None,
                alphafold_fragment=str(cfg.get("alphafold_fragment") or "").strip() or None,
                synonyms=None,
            )
            n += 1
    finally:
        con.close()
    print(f"Imported {n} protein(s) from presets into {db_path.name}")
    return 0


def cmd_db_import_manifest(ns: argparse.Namespace) -> int:
    import csv

    db_path = Path(ns.db)
    if not db_path.is_absolute():
        db_path = (BASE_DIR / db_path).resolve()
    manifest = Path(ns.manifest)
    if not manifest.is_absolute():
        manifest = (BASE_DIR / manifest).resolve()
    if not manifest.is_file():
        raise SystemExit(f"Manifest not found: {manifest}")

    delim = "\t" if manifest.suffix.lower() in (".tsv",) else ","
    con = avi_db.connect(db_path)
    try:
        avi_db.init_db(con)
        with manifest.open("r", encoding="utf-8-sig", newline="") as f:
            rdr = csv.DictReader(f, delimiter=delim)
            if rdr.fieldnames is None:
                raise SystemExit("Manifest has no header row.")
            n = 0
            for row in rdr:
                uid = str(row.get("uniprot_id") or "").strip()
                if not uid:
                    continue
                preset_key = str(row.get("preset_key") or "").strip() or None
                gene = str(row.get("gene_symbol") or "").strip() or None
                clin = str(row.get("clinvar_esearch_term") or "").strip() or None
                outb = str(row.get("output_basename") or "").strip() or None
                frag = str(row.get("alphafold_fragment") or "").strip() or None
                syn = str(row.get("synonyms") or "").strip()
                synonyms = [s.strip() for s in syn.split(";") if s.strip()] if syn else None
                avi_db.upsert_protein(
                    con,
                    preset_key=preset_key,
                    uniprot_id=uid,
                    gene_symbol=gene,
                    clinvar_esearch_term=clin,
                    output_basename=outb,
                    alphafold_fragment=frag,
                    synonyms=synonyms,
                )
                n += 1
    finally:
        con.close()
    print(f"Imported {n} protein(s) from manifest into {db_path.name}")
    return 0


def cmd_db_index_runs(ns: argparse.Namespace) -> int:
    db_path = Path(ns.db)
    if not db_path.is_absolute():
        db_path = (BASE_DIR / db_path).resolve()
    runs_root = Path(ns.runs_root)
    if not runs_root.is_absolute():
        runs_root = (BASE_DIR / runs_root).resolve()

    con = avi_db.connect(db_path)
    try:
        avi_db.init_db(con)
        indexed = 0
        for run_dir in avi_db.iter_run_dirs(runs_root):
            cfg_path = run_dir / "pipeline_config.json"
            if not cfg_path.is_file():
                continue
            try:
                cfg = _load_json(cfg_path)
            except Exception:
                continue

            # Ensure protein row exists.
            preset_key = None
            # When indexing run dirs (not presets), we don't necessarily know the preset key.
            pid = avi_db.upsert_protein(
                con,
                preset_key=preset_key,
                uniprot_id=str(cfg.get("uniprot_id") or "").strip(),
                gene_symbol=str(cfg.get("gene_symbol") or "").strip() or None,
                clinvar_esearch_term=str(cfg.get("clinvar_esearch_term") or "").strip() or None,
                output_basename=str(cfg.get("output_basename") or "").strip() or None,
                alphafold_fragment=str(cfg.get("alphafold_fragment") or "").strip() or None,
                synonyms=None,
            )

            # Artifacts.
            data_root = run_dir / "data"
            art = artifact_paths(cfg, data_root)
            report_html = run_dir / "reports" / "exploration_report.html"
            eval_json = art.get("evaluation_metrics_json")
            eval_exists = eval_json.is_file() if isinstance(eval_json, Path) else False
            status = "ok"
            err = None

            # Created timestamp (best effort).
            created_utc = None
            run_meta = run_dir / "run.json"
            if run_meta.is_file():
                try:
                    rj = json.loads(run_meta.read_text(encoding="utf-8-sig"))
                    if isinstance(rj, dict) and isinstance(rj.get("created_utc"), str):
                        created_utc = str(rj["created_utc"])
                except Exception:
                    pass

            run_id = avi_db.upsert_run(
                con,
                run_dir=run_dir,
                protein_id=pid,
                created_utc=created_utc,
                status=status,
                report_html_path=report_html if report_html.is_file() else None,
                evaluation_json_path=eval_json if eval_exists else None,
                error=err,
            )
            if eval_exists:
                try:
                    doc = json.loads(Path(eval_json).read_text(encoding="utf-8-sig"))
                    if isinstance(doc, dict):
                        avi_db.upsert_run_metrics(con, run_id=run_id, metrics=doc)
                except Exception:
                    pass
            indexed += 1
    finally:
        con.close()
    print(f"Indexed {indexed} run dir(s) under {runs_root}")
    return 0


def cmd_ui(ns: argparse.Namespace) -> int:
    db_path = Path(ns.db)
    if not db_path.is_absolute():
        db_path = (BASE_DIR / db_path).resolve()
    con = avi_db.connect(db_path)
    try:
        avi_db.init_db(con)
    finally:
        con.close()
    if getattr(ns, "open", False):
        try:
            webbrowser.open(f"http://{ns.host}:{int(ns.port)}/")
        except Exception:
            pass
    avi_webui.serve(base_dir=BASE_DIR, db_path=db_path, host=str(ns.host), port=int(ns.port))
    return 0


def cmd_dataset_run(ns: argparse.Namespace) -> int:
    """
    Batch runner driven by a protein manifest -> runs/ + avi.db.

    This is the "universal" entry point for scaling to many proteins: load targets into SQLite,
    execute pipeline per target, and index artifacts/metrics for the UI.
    """
    db_path = Path(ns.db)
    if not db_path.is_absolute():
        db_path = (BASE_DIR / db_path).resolve()
    runs_parent = Path(ns.runs_parent)
    if not runs_parent.is_absolute():
        runs_parent = (BASE_DIR / runs_parent).resolve()

    # Ensure DB exists and load targets.
    con = avi_db.connect(db_path)
    try:
        avi_db.init_db(con)
    finally:
        con.close()

    # Import targets from manifest if provided.
    if ns.manifest:
        # Reuse the db import-manifest command to keep behavior consistent.
        cmd_db_import_manifest(
            argparse.Namespace(db=str(db_path), manifest=str(ns.manifest))
        )

    # Select proteins to run.
    con = avi_db.connect(db_path)
    try:
        avi_db.init_db(con)
        wanted_keys: list[str] | None = None
        if ns.presets:
            wanted_keys = [p.strip() for p in str(ns.presets).split(",") if p.strip()]
        wanted_ids: set[int] | None = None
        if getattr(ns, "protein_ids", None):
            wanted_ids = set()
            for t in str(ns.protein_ids).split(","):
                t = t.strip()
                if not t:
                    continue
                wanted_ids.add(int(t))

        rows = con.execute("SELECT * FROM proteins ORDER BY updated_utc DESC;").fetchall()
        selected = []
        for r in rows:
            if wanted_ids is not None and int(r["id"]) not in wanted_ids:
                continue
            if wanted_keys is not None:
                pk = (r["preset_key"] or "").strip()
                if pk not in wanted_keys:
                    continue
            selected.append(r)
        if ns.limit is not None:
            selected = selected[: int(ns.limit)]
    finally:
        con.close()

    if not selected:
        raise SystemExit("No proteins selected. Import a manifest or set --presets.")

    failures: list[tuple[str, str]] = []
    for idx, p in enumerate(selected):
        cfg_preview = {
            "uniprot_id": str(p["uniprot_id"]),
            "gene_symbol": str(p["gene_symbol"] or ""),
            "clinvar_esearch_term": str(p["clinvar_esearch_term"] or f"{p['gene_symbol']}[gene]"),
            "output_basename": str(p["output_basename"] or (p["gene_symbol"] or p["uniprot_id"])).strip(),
            "alphafold_fragment": str(p["alphafold_fragment"] or "F1"),
        }
        basename = str(cfg_preview["output_basename"]).strip() or str(p["uniprot_id"])
        safe = basename.replace("/", "_").replace("\\", "_")
        run_ts = datetime.now(UTC).strftime("%Y%m%d_%H%M%S")
        suffix = f"_{idx:03d}"
        run_dir = runs_parent / safe / f"{run_ts}{suffix}"
        run_dir.mkdir(parents=True, exist_ok=True)
        _write_config(cfg_preview, run_dir / "pipeline_config.json")
        (run_dir / "run.json").write_text(
            json.dumps(
                {
                    "created_utc": datetime.now(UTC).replace(microsecond=0).isoformat(),
                    "command": "avi dataset run",
                    "protein_id": int(p["id"]),
                    "run_dir": str(run_dir.resolve()),
                },
                indent=2,
                ensure_ascii=False,
            )
            + "\n",
            encoding="utf-8",
        )

        # Mark running early so UI can show progress.
        try:
            con = avi_db.connect(db_path)
            try:
                avi_db.init_db(con)
                avi_db.upsert_run(
                    con,
                    run_dir=run_dir,
                    protein_id=int(p["id"]),
                    created_utc=None,
                    status="running",
                    report_html_path=None,
                    evaluation_json_path=None,
                    error=None,
                )
            finally:
                con.close()
        except Exception:
            pass

        rc = 0
        err: str | None = None
        try:
            env = _pipeline_env_for_run_dir(run_dir)
            extra: list[str] = []
            if ns.force_download:
                extra.append("--force-download")
            rc = _run_script("run_pipeline.py", extra, env=env)
            if rc == 0:
                _write_run_manifest(run_dir=run_dir)
                if ns.report:
                    _export_html_report(run_dir=run_dir, env=env, notebook_path=None)
                if ns.evaluate:
                    ev_args = [
                        "--run-dir",
                        str(run_dir),
                        "--seed",
                        str(ns.evaluate_seed),
                        "--min-label-rows",
                        str(ns.evaluate_min_label_rows),
                        "--min-dated-rows",
                        str(ns.evaluate_min_dated_rows),
                    ]
                    rc = _run_script("evaluation_metrics.py", ev_args, env=None)
        except Exception as e:
            rc = 1
            err = f"{type(e).__name__}: {e}"

        # Index into DB (best effort).
        try:
            con = avi_db.connect(db_path)
            try:
                avi_db.init_db(con)
                cfg = _load_json(run_dir / "pipeline_config.json")
                art = artifact_paths(cfg, run_dir / "data")
                report_html = run_dir / "reports" / "exploration_report.html"
                eval_json = art.get("evaluation_metrics_json")
                eval_path = eval_json if isinstance(eval_json, Path) and eval_json.is_file() else None
                run_id = avi_db.upsert_run(
                    con,
                    run_dir=run_dir,
                    protein_id=int(p["id"]),
                    created_utc=None,
                    status="ok" if rc == 0 else "failed",
                    report_html_path=report_html if report_html.is_file() else None,
                    evaluation_json_path=eval_path,
                    error=err,
                )
                if eval_path is not None:
                    try:
                        doc = json.loads(eval_path.read_text(encoding="utf-8-sig"))
                        if isinstance(doc, dict):
                            avi_db.upsert_run_metrics(con, run_id=run_id, metrics=doc)
                    except Exception:
                        pass
            finally:
                con.close()
        except Exception:
            pass

        if rc != 0:
            label = str(p["gene_symbol"] or p["preset_key"] or p["uniprot_id"])
            failures.append((label, str(run_dir)))
            if not ns.continue_on_error:
                break

    if failures:
        print(f"Dataset run finished with {len(failures)} failure(s):")
        for lab, rd in failures:
            print(f"  {lab}: {rd}")
        raise SystemExit(1)
    return 0
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="avi",
        description="AlphaFold x ClinVar pipeline - front-door CLI (wraps scripts/).",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    p_init = sub.add_parser(
        "init",
        help="Create a timestamped run dir + pipeline_config.json.",
    )
    src = p_init.add_mutually_exclusive_group(required=True)
    src.add_argument(
        "--preset",
        metavar="NAME",
        help="Preset name from avi/presets.json (edit this file to add targets).",
    )
    src.add_argument(
        "--from-uniprot",
        metavar="ACC",
        help="Create a one-off run by fetching gene symbol from UniProt (e.g. P38398).",
    )
    p_init.add_argument(
        "--key",
        default=None,
        help="Optional label used for defaults when using --from-uniprot.",
    )
    p_init.add_argument(
        "--gene-symbol",
        default=None,
        help="With --from-uniprot: override the resolved gene symbol.",
    )
    p_init.add_argument(
        "--output-basename",
        default=None,
        help="With --from-uniprot: override output basename directory/file prefix.",
    )
    p_init.add_argument(
        "--clinvar-esearch-term",
        default=None,
        help="With --from-uniprot: override ClinVar esearch term.",
    )
    p_init.add_argument(
        "--alphafold-fragment",
        default=None,
        help="With --from-uniprot: override AlphaFold fragment (usually F1).",
    )
    p_init.add_argument(
        "--run",
        action="store_true",
        help="After creating the run directory, immediately run the pipeline.",
    )
    p_init.add_argument(
        "--force-download",
        dest="force_download",
        action="store_true",
        help="With --run: force refresh of raw downloads.",
    )
    p_init.add_argument(
        "--report",
        action="store_true",
        help="With --run: export an executed HTML report into <run-dir>/reports/.",
    )
    p_init.add_argument(
        "--open-report",
        action="store_true",
        help="With --run: open the HTML report in your browser.",
    )
    p_init.add_argument(
        "--notebook",
        action="store_true",
        help="With --run: launch Jupyter Notebook for this run (blocks).",
    )
    p_init.add_argument(
        "--runs-parent",
        default="runs",
        help="Directory under project root to store runs (default: runs).",
    )
    p_init.set_defaults(func=cmd_init)

    p_add = sub.add_parser(
        "add-preset",
        help="Add or update a preset in avi/presets.json.",
    )
    p_add.add_argument("key", help="Preset key (e.g. brca1, mygene).")
    p_add.add_argument(
        "--from-uniprot",
        default=None,
        help="Fetch gene symbol from UniProt and use it to prefill fields.",
    )
    p_add.add_argument("--uniprot-id", dest="uniprot_id", default=None)
    p_add.add_argument("--gene-symbol", dest="gene_symbol", default=None)
    p_add.add_argument("--clinvar-esearch-term", dest="clinvar_esearch_term", default=None)
    p_add.add_argument("--output-basename", dest="output_basename", default=None)
    p_add.add_argument("--alphafold-fragment", dest="alphafold_fragment", default=None)
    p_add.add_argument("--overwrite", action="store_true", help="Replace if key exists.")
    p_add.add_argument(
        "--auto",
        action="store_true",
        help="Non-interactive using defaults (e.g. inferred from --from-uniprot).",
    )
    p_add.add_argument(
        "--non-interactive",
        action="store_true",
        help="Do not prompt; require all fields as flags.",
    )
    p_add.set_defaults(func=cmd_add_preset)

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
        "--resume",
        action="store_true",
        help="Resume a run by skipping stages whose outputs already exist (writes run_state.json).",
    )
    p_run.add_argument(
        "--no-subset",
        action="store_true",
        help="When running full pipeline, skip build_missense_subset.py.",
    )
    p_run.add_argument(
        "--no-manifest",
        action="store_true",
        help="Skip writing repro manifest JSON to <run-dir>/data/processed/.",
    )
    p_run.set_defaults(func=cmd_run)

    p_batch = sub.add_parser(
        "batch",
        help="Run the pipeline for multiple presets or UniProt accessions into runs/.",
    )
    srcb = p_batch.add_mutually_exclusive_group()
    srcb.add_argument(
        "--presets",
        default=None,
        help="Comma-separated preset keys. Omit to run all presets in avi/presets.json.",
    )
    srcb.add_argument(
        "--from-uniprot",
        default=None,
        help="Comma-separated UniProt accessions (runs without needing presets).",
    )
    p_batch.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Run only first N selected presets (useful for smoke runs).",
    )
    p_batch.add_argument(
        "--runs-parent",
        default="runs",
        help="Directory under project root to store runs (default: runs).",
    )
    p_batch.add_argument(
        "--force-download",
        action="store_true",
        help="Force refresh of raw downloads for each run.",
    )
    p_batch.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue remaining presets if one fails.",
    )
    p_batch.add_argument(
        "--report",
        action="store_true",
        help="After each successful run, export an executed HTML report into <run-dir>/reports/.",
    )
    p_batch.add_argument(
        "--evaluate",
        action="store_true",
        help="After each successful run, write evaluation_metrics.json into data/processed/.",
    )
    p_batch.add_argument(
        "--evaluate-seed",
        type=int,
        default=42,
        help="With --evaluate: stratified split RNG seed.",
    )
    p_batch.add_argument(
        "--evaluate-min-label-rows",
        type=int,
        default=40,
        metavar="N",
        help="With --evaluate: minimum Path/LP + Ben/LB rows required.",
    )
    p_batch.add_argument(
        "--evaluate-min-dated-rows",
        type=int,
        default=50,
        metavar="N",
        help="With --evaluate: minimum parseable germline_date_last_evaluated rows for time split.",
    )
    p_batch.set_defaults(func=cmd_batch)

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

    p_repro = sub.add_parser(
        "repro",
        help="Write a repro manifest JSON for a run (config + raw file hashes).",
    )
    p_repro.add_argument("--run-dir", required=True)
    p_repro.set_defaults(func=cmd_repro)

    p_clean = sub.add_parser(
        "clean",
        help="Prune old runs under runs/ (dry-run unless --yes).",
    )
    p_clean.add_argument("--runs-root", default=None)
    p_clean.add_argument("--keep", type=int, default=20, help="Keep newest N runs (default: 20).")
    p_clean.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview deletions without removing anything (default unless --yes).",
    )
    p_clean.add_argument("--yes", action="store_true", help="Actually delete (overrides dry-run).")
    p_clean.set_defaults(func=cmd_clean)

    p_gc = sub.add_parser(
        "gc",
        help="Garbage-collect problematic runs (dry-run unless --yes).",
    )
    p_gc.add_argument("--runs-root", default=None)
    p_gc.add_argument(
        "--delete-empty",
        action="store_true",
        help="Delete run dirs missing pipeline_config.json.",
    )
    p_gc.add_argument(
        "--delete-partial",
        action="store_true",
        help="Delete run dirs that have a config but no processed *_variants_basic.csv.",
    )
    p_gc.add_argument(
        "--older-than-days",
        type=int,
        default=None,
        help="Only target runs older than N days (by mtime).",
    )
    p_gc.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview deletions without removing anything (default unless --yes).",
    )
    p_gc.add_argument("--yes", action="store_true", help="Actually delete (overrides dry-run).")
    p_gc.set_defaults(func=cmd_gc)

    p_nb = sub.add_parser(
        "notebook",
        help="Launch Jupyter Notebook with this run's env vars set.",
    )
    p_nb.add_argument(
        "--run-dir",
        required=True,
        help="Run directory containing pipeline_config.json (sets PIPELINE_* env vars).",
    )
    p_nb.set_defaults(func=cmd_notebook)

    p_rep = sub.add_parser(
        "report",
        help="Execute the notebook and write an HTML report under <run-dir>/reports/.",
    )
    p_rep.add_argument(
        "--run-dir",
        required=True,
        help="Run directory containing pipeline_config.json (sets PIPELINE_* env vars).",
    )
    p_rep.add_argument(
        "--notebook",
        default=None,
        help="Override notebook path to execute (default: notebooks/01_target_exploration.ipynb).",
    )
    p_rep.set_defaults(func=cmd_report)

    p_open = sub.add_parser(
        "open-report",
        help="Open the HTML report for a run in your browser.",
    )
    p_open.add_argument(
        "--run-dir",
        required=True,
        help="Run directory containing reports/exploration_report.html.",
    )
    p_open.set_defaults(func=cmd_open_report)

    p_ev = sub.add_parser(
        "evaluate",
        help="Write held-out evaluation_metrics.json (path vs benign on structure features).",
    )
    p_ev.add_argument(
        "--run-dir",
        required=True,
        help="Run directory with data/processed/*_missense_mappable.csv.",
    )
    p_ev.add_argument("--seed", type=int, default=42, help="Stratified split RNG seed.")
    p_ev.add_argument(
        "--min-label-rows",
        type=int,
        default=40,
        metavar="N",
        help="Minimum Path/LP + Ben/LB rows required to run evaluation.",
    )
    p_ev.add_argument(
        "--min-dated-rows",
        type=int,
        default=50,
        metavar="N",
        help="Minimum parseable germline_date_last_evaluated rows to use time-based split.",
    )
    p_ev.set_defaults(func=cmd_evaluate)

    p_db = sub.add_parser("db", help="SQLite index for proteins/runs (for UI/search).")
    dbsub = p_db.add_subparsers(dest="db_cmd", required=True)
    p_db_init = dbsub.add_parser("init", help="Create/initialize the SQLite database.")
    p_db_init.add_argument("--db", default=str(DEFAULT_DB_PATH), help="Path to avi.db (default: ./avi.db).")
    p_db_init.set_defaults(func=cmd_db_init)

    p_db_imp = dbsub.add_parser("import-presets", help="Import avi/presets.json into the DB.")
    p_db_imp.add_argument("--db", default=str(DEFAULT_DB_PATH), help="Path to avi.db (default: ./avi.db).")
    p_db_imp.set_defaults(func=cmd_db_import_presets)

    p_db_m = dbsub.add_parser("import-manifest", help="Import a proteins manifest CSV/TSV into the DB.")
    p_db_m.add_argument("--db", default=str(DEFAULT_DB_PATH), help="Path to avi.db (default: ./avi.db).")
    p_db_m.add_argument(
        "--manifest",
        required=True,
        help="CSV/TSV with header: uniprot_id (required), gene_symbol, preset_key, clinvar_esearch_term, output_basename, alphafold_fragment, synonyms (semicolon-separated).",
    )
    p_db_m.set_defaults(func=cmd_db_import_manifest)

    p_db_idx = dbsub.add_parser("index-runs", help="Index existing run directories into the DB.")
    p_db_idx.add_argument("--db", default=str(DEFAULT_DB_PATH), help="Path to avi.db (default: ./avi.db).")
    p_db_idx.add_argument("--runs-root", default="runs", help="Runs root directory to scan (default: runs).")
    p_db_idx.set_defaults(func=cmd_db_index_runs)

    p_ui = sub.add_parser("ui", help="Start local search UI (reads avi.db).")
    p_ui.add_argument("--db", default=str(DEFAULT_DB_PATH), help="Path to avi.db (default: ./avi.db).")
    p_ui.add_argument("--host", default="127.0.0.1", help="Bind host (default: 127.0.0.1).")
    p_ui.add_argument("--port", type=int, default=8000, help="Bind port (default: 8000).")
    p_ui.add_argument("--open", action="store_true", help="Open the UI in your browser.")
    p_ui.set_defaults(func=cmd_ui)

    p_ds = sub.add_parser(
        "dataset",
        help="Run pipeline across many proteins from manifest/DB (updates avi.db for UI).",
    )
    dssub = p_ds.add_subparsers(dest="ds_cmd", required=True)
    p_ds_run = dssub.add_parser("run", help="Run pipeline over a set of proteins.")
    p_ds_run.add_argument("--db", default=str(DEFAULT_DB_PATH), help="Path to avi.db (default: ./avi.db).")
    p_ds_run.add_argument(
        "--manifest",
        default=None,
        help="Optional manifest CSV/TSV to import before running (see data/proteins_manifest_example.csv).",
    )
    p_ds_run.add_argument(
        "--presets",
        default=None,
        help="Optional comma-separated preset_key filters (runs all proteins if omitted).",
    )
    p_ds_run.add_argument(
        "--protein-ids",
        default=None,
        help="Optional comma-separated proteins.id values from avi.db (overrides --presets).",
    )
    p_ds_run.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Run only first N selected proteins (useful for smoke runs).",
    )
    p_ds_run.add_argument(
        "--runs-parent",
        default="runs_dataset",
        help="Output parent directory under repo (default: runs_dataset).",
    )
    p_ds_run.add_argument("--force-download", action="store_true", help="Force refresh of raw downloads.")
    p_ds_run.add_argument("--continue-on-error", action="store_true", help="Continue after failures.")
    p_ds_run.add_argument("--report", action="store_true", help="Export executed HTML report per run.")
    p_ds_run.add_argument("--evaluate", action="store_true", help="Write evaluation_metrics.json per run.")
    p_ds_run.add_argument("--evaluate-seed", type=int, default=42)
    p_ds_run.add_argument("--evaluate-min-label-rows", type=int, default=40, metavar="N")
    p_ds_run.add_argument("--evaluate-min-dated-rows", type=int, default=50, metavar="N")
    p_ds_run.set_defaults(func=cmd_dataset_run)

    p_list = sub.add_parser(
        "list-runs",
        help="List recent run directories under runs/.",
    )
    p_list.add_argument(
        "--runs-root",
        default=None,
        help="Override runs root directory (default: ./runs).",
    )
    p_list.add_argument(
        "--limit",
        type=int,
        default=20,
        help="Maximum number of runs to list (default: 20).",
    )
    p_list.set_defaults(func=cmd_list_runs)

    return p


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    ns = parser.parse_args(argv)

    assert ns.func is not None
    return int(ns.func(ns))


if __name__ == "__main__":
    raise SystemExit(main())

