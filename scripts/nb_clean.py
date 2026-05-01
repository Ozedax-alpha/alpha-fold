"""
Normalize notebooks for clean diffs.

- Ensure each cell has an id
- Clear execution counts and outputs
- Drop transient metadata (where safe)

Designed to be used as a pre-commit hook.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from uuid import uuid4

import nbformat


def clean_notebook(path: Path) -> bool:
    nb = nbformat.read(str(path), as_version=4)
    changed = False

    # Ensure stable cell ids (required by newer nbformat).
    for cell in nb.cells:
        if not cell.get("id"):
            cell["id"] = str(uuid4())
            changed = True

        # Strip execution artifacts.
        if cell.get("cell_type") == "code":
            if cell.get("execution_count") is not None:
                cell["execution_count"] = None
                changed = True
            if cell.get("outputs"):
                cell["outputs"] = []
                changed = True

        # Notebook UIs often add per-cell metadata; keep only what we need.
        meta = cell.get("metadata") or {}
        if isinstance(meta, dict) and meta:
            keep = {}
            # Keep tags if present (used by some workflows).
            if "tags" in meta:
                keep["tags"] = meta["tags"]
            if keep != meta:
                cell["metadata"] = keep
                changed = True

    # Top-level metadata: keep kernelspec/language_info when present.
    nb_meta = nb.get("metadata") or {}
    if isinstance(nb_meta, dict) and nb_meta:
        keep_top = {}
        for k in ("kernelspec", "language_info"):
            if k in nb_meta:
                keep_top[k] = nb_meta[k]
        if keep_top != nb_meta:
            nb["metadata"] = keep_top
            changed = True

    if changed:
        nbformat.write(nb, str(path))
    return changed


def main() -> None:
    ap = argparse.ArgumentParser(description="Normalize .ipynb files for clean diffs.")
    ap.add_argument("paths", nargs="+", help="Notebook paths to clean.")
    args = ap.parse_args()

    any_changed = False
    for p in args.paths:
        path = Path(p)
        if path.suffix.lower() != ".ipynb":
            continue
        if not path.is_file():
            raise SystemExit(f"Missing notebook: {path}")
        any_changed |= clean_notebook(path)

    raise SystemExit(1 if any_changed else 0)


if __name__ == "__main__":
    main()

