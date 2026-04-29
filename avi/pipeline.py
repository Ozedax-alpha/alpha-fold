"""Pipeline orchestration helpers shared by `avi` and `scripts/*`.

This module intentionally stays lightweight and offline-friendly: it defines
stage ordering, "done" checks via artifact existence, and optional resume state.
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Callable

from avi.paths import artifact_paths

STAGE_ORDER: list[str] = ["download", "variants", "subset", "report", "repro"]


def run_state_path(run_dir: Path) -> Path:
    return run_dir / "run_state.json"


def load_run_state(run_dir: Path) -> dict[str, Any]:
    p = run_state_path(run_dir)
    if not p.is_file():
        return {"stages": {}}
    try:
        doc: Any = json.loads(p.read_text(encoding="utf-8-sig"))
    except Exception:
        return {"stages": {}}
    if not isinstance(doc, dict):
        return {"stages": {}}
    if "stages" not in doc or not isinstance(doc.get("stages"), dict):
        doc["stages"] = {}
    return doc


def save_run_state(run_dir: Path, state: dict[str, Any], *, atomic_write_text: Callable[[Path, str], None]) -> None:
    atomic_write_text(run_state_path(run_dir), json.dumps(state, indent=2, ensure_ascii=False) + "\n")


def artifact_done_for_stage(stage: str, *, run_dir: Path, cfg: dict[str, Any], report_html: Path) -> bool:
    data_root = run_dir / "data"
    art = artifact_paths(cfg, data_root)
    if stage == "download":
        return art["alphafold_pdb"].is_file() and art["clinvar_json"].is_file()
    if stage == "variants":
        return art["variants_basic_csv"].is_file()
    if stage == "subset":
        return art["missense_mappable_csv"].is_file()
    if stage == "report":
        return report_html.is_file()
    if stage == "repro":
        return art["repro_manifest"].is_file()
    return False


@dataclass(frozen=True)
class StageContext:
    run_dir: Path
    env: dict[str, str] | None
    cfg: dict[str, Any] | None
    resume: bool
    report_html: Path


def run_stages(
    *,
    ctx: StageContext,
    stages: list[str],
    stage_runner: Callable[[str], int],
    atomic_write_text: Callable[[Path, str], None],
) -> int:
    """Run stages in canonical order with optional resume + run_state.json updates."""
    if any(s not in STAGE_ORDER for s in stages):
        unknown = [s for s in stages if s not in STAGE_ORDER]
        raise SystemExit(f"Unknown stage(s): {', '.join(unknown)}. Choose from: {', '.join(STAGE_ORDER)}")

    wanted = []
    for s in stages:
        if s not in wanted:
            wanted.append(s)

    cfg = ctx.cfg
    state = load_run_state(ctx.run_dir) if ctx.run_dir else None

    for s in STAGE_ORDER:
        if s not in wanted:
            continue

        if ctx.resume and cfg is not None:
            if artifact_done_for_stage(s, run_dir=ctx.run_dir, cfg=cfg, report_html=ctx.report_html):
                if state is not None:
                    stages_state = state.setdefault("stages", {})
                    if isinstance(stages_state, dict):
                        stages_state[str(s)] = {
                            "skipped": True,
                            "reason": "already complete",
                            "finished_utc": datetime.now(UTC).replace(microsecond=0).isoformat(),
                            "return_code": 0,
                        }
                        save_run_state(ctx.run_dir, state, atomic_write_text=atomic_write_text)
                print(f"= Skipping stage {s!r} (already complete)")
                continue

        if state is not None:
            stages_state = state.setdefault("stages", {})
            if isinstance(stages_state, dict):
                stages_state[str(s)] = {"started_utc": datetime.now(UTC).replace(microsecond=0).isoformat()}
                save_run_state(ctx.run_dir, state, atomic_write_text=atomic_write_text)

        rc = int(stage_runner(s))

        if state is not None:
            stages_state = state.setdefault("stages", {})
            if isinstance(stages_state, dict):
                prev = stages_state.get(str(s)) if isinstance(stages_state.get(str(s)), dict) else {}
                prev.update(
                    {
                        "finished_utc": datetime.now(UTC).replace(microsecond=0).isoformat(),
                        "return_code": rc,
                    }
                )
                stages_state[str(s)] = prev
                save_run_state(ctx.run_dir, state, atomic_write_text=atomic_write_text)

        if rc != 0:
            return rc

    return 0

