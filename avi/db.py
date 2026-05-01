"""SQLite-backed index for proteins, runs, and metrics (UI/search)."""

from __future__ import annotations

import json
import sqlite3
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable


SCHEMA_VERSION = 1


@dataclass(frozen=True)
class DBPaths:
    db_path: Path


def default_db_path(base_dir: Path) -> Path:
    return base_dir / "avi.db"


def connect(db_path: Path) -> sqlite3.Connection:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(str(db_path))
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA foreign_keys = ON;")
    return con


def init_db(con: sqlite3.Connection) -> None:
    con.execute(
        """
        CREATE TABLE IF NOT EXISTS meta (
          key TEXT PRIMARY KEY,
          value TEXT NOT NULL
        );
        """
    )
    con.execute(
        """
        CREATE TABLE IF NOT EXISTS proteins (
          id INTEGER PRIMARY KEY,
          preset_key TEXT UNIQUE,
          uniprot_id TEXT NOT NULL,
          gene_symbol TEXT,
          clinvar_esearch_term TEXT,
          output_basename TEXT,
          alphafold_fragment TEXT,
          synonyms_json TEXT NOT NULL DEFAULT '[]',
          synonyms_text TEXT NOT NULL DEFAULT '',
          created_utc TEXT NOT NULL,
          updated_utc TEXT NOT NULL
        );
        """
    )
    # Ensure we can upsert by UniProt even when preset_key is absent.
    con.execute(
        "CREATE UNIQUE INDEX IF NOT EXISTS proteins_uniprot_id_uq ON proteins(uniprot_id);"
    )
    con.execute(
        """
        CREATE TABLE IF NOT EXISTS runs (
          id INTEGER PRIMARY KEY,
          protein_id INTEGER,
          run_dir TEXT NOT NULL UNIQUE,
          created_utc TEXT,
          status TEXT NOT NULL DEFAULT 'unknown',
          report_html_path TEXT,
          evaluation_json_path TEXT,
          error TEXT,
          updated_utc TEXT NOT NULL,
          FOREIGN KEY (protein_id) REFERENCES proteins(id) ON DELETE SET NULL
        );
        """
    )
    con.execute(
        """
        CREATE TABLE IF NOT EXISTS run_metrics (
          run_id INTEGER PRIMARY KEY,
          metrics_json TEXT NOT NULL,
          updated_utc TEXT NOT NULL,
          FOREIGN KEY (run_id) REFERENCES runs(id) ON DELETE CASCADE
        );
        """
    )
    con.execute(
        "INSERT OR REPLACE INTO meta(key, value) VALUES(?, ?);",
        ("schema_version", str(SCHEMA_VERSION)),
    )
    con.commit()

    # Lightweight "migration" for older DBs: ensure synonyms_text exists.
    cols = {r[1] for r in con.execute("PRAGMA table_info(proteins);").fetchall()}
    if "synonyms_text" not in cols:
        con.execute("ALTER TABLE proteins ADD COLUMN synonyms_text TEXT NOT NULL DEFAULT '';")
        con.commit()

    # Best-effort FTS5: works on most modern Python builds; if missing, UI will fall back.
    try:
        con.execute(
            """
            CREATE VIRTUAL TABLE IF NOT EXISTS proteins_fts USING fts5(
              preset_key,
              uniprot_id,
              gene_symbol,
              synonyms,
              content='proteins',
              content_rowid='id'
            );
            """
        )
        con.executescript(
            """
            CREATE TRIGGER IF NOT EXISTS proteins_ai AFTER INSERT ON proteins BEGIN
              INSERT INTO proteins_fts(rowid, preset_key, uniprot_id, gene_symbol, synonyms)
              VALUES (new.id, new.preset_key, new.uniprot_id, new.gene_symbol, new.synonyms_text);
            END;
            CREATE TRIGGER IF NOT EXISTS proteins_ad AFTER DELETE ON proteins BEGIN
              INSERT INTO proteins_fts(proteins_fts, rowid, preset_key, uniprot_id, gene_symbol, synonyms)
              VALUES('delete', old.id, old.preset_key, old.uniprot_id, old.gene_symbol, old.synonyms_text);
            END;
            CREATE TRIGGER IF NOT EXISTS proteins_au AFTER UPDATE ON proteins BEGIN
              INSERT INTO proteins_fts(proteins_fts, rowid, preset_key, uniprot_id, gene_symbol, synonyms)
              VALUES('delete', old.id, old.preset_key, old.uniprot_id, old.gene_symbol, old.synonyms_text);
              INSERT INTO proteins_fts(rowid, preset_key, uniprot_id, gene_symbol, synonyms)
              VALUES (new.id, new.preset_key, new.uniprot_id, new.gene_symbol, new.synonyms_text);
            END;
            """
        )
        con.commit()
    except sqlite3.OperationalError:
        # FTS5 or JSON1 may be unavailable; proceed without it.
        pass


def _now_utc() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat()


def upsert_protein(
    con: sqlite3.Connection,
    *,
    preset_key: str | None,
    uniprot_id: str,
    gene_symbol: str | None,
    clinvar_esearch_term: str | None,
    output_basename: str | None,
    alphafold_fragment: str | None,
    synonyms: list[str] | None = None,
) -> int:
    now = _now_utc()
    syn = synonyms or []
    syn_clean = [s.strip() for s in syn if isinstance(s, str) and s.strip()]
    syn_json = json.dumps(syn_clean, ensure_ascii=False)
    syn_text = " ".join(syn_clean)

    # If preset_key exists, use it as the stable unique key.
    if preset_key:
        con.execute(
            """
            INSERT INTO proteins(
              preset_key, uniprot_id, gene_symbol, clinvar_esearch_term, output_basename, alphafold_fragment,
              synonyms_json, synonyms_text, created_utc, updated_utc
            ) VALUES(?,?,?,?,?,?,?,?,?,?)
            ON CONFLICT(preset_key) DO UPDATE SET
              uniprot_id=excluded.uniprot_id,
              gene_symbol=excluded.gene_symbol,
              clinvar_esearch_term=excluded.clinvar_esearch_term,
              output_basename=excluded.output_basename,
              alphafold_fragment=excluded.alphafold_fragment,
              synonyms_json=excluded.synonyms_json,
              synonyms_text=excluded.synonyms_text,
              updated_utc=excluded.updated_utc;
            """,
            (
                preset_key,
                uniprot_id,
                gene_symbol,
                clinvar_esearch_term,
                output_basename,
                alphafold_fragment,
                syn_json,
                syn_text,
                now,
                now,
            ),
        )
        row = con.execute("SELECT id FROM proteins WHERE preset_key = ?;", (preset_key,)).fetchone()
        assert row is not None
        con.commit()
        return int(row["id"])

    # Otherwise, fall back to UniProt stable unique.
    con.execute(
        """
        INSERT INTO proteins(
          preset_key, uniprot_id, gene_symbol, clinvar_esearch_term, output_basename, alphafold_fragment,
          synonyms_json, synonyms_text, created_utc, updated_utc
        ) VALUES(NULL,?,?,?,?,?,?,?,?)
        ON CONFLICT(uniprot_id) DO UPDATE SET
          gene_symbol=excluded.gene_symbol,
          clinvar_esearch_term=excluded.clinvar_esearch_term,
          output_basename=excluded.output_basename,
          alphafold_fragment=excluded.alphafold_fragment,
          synonyms_json=excluded.synonyms_json,
          synonyms_text=excluded.synonyms_text,
          updated_utc=excluded.updated_utc;
        """,
        (
            uniprot_id,
            gene_symbol,
            clinvar_esearch_term,
            output_basename,
            alphafold_fragment,
            syn_json,
            syn_text,
            now,
            now,
        ),
    )
    row = con.execute("SELECT id FROM proteins WHERE uniprot_id = ?;", (uniprot_id,)).fetchone()
    assert row is not None
    con.commit()
    return int(row["id"])


def upsert_run(
    con: sqlite3.Connection,
    *,
    run_dir: Path,
    protein_id: int | None,
    created_utc: str | None,
    status: str,
    report_html_path: Path | None,
    evaluation_json_path: Path | None,
    error: str | None,
) -> int:
    now = _now_utc()
    con.execute(
        """
        INSERT INTO runs(
          protein_id, run_dir, created_utc, status, report_html_path, evaluation_json_path, error, updated_utc
        ) VALUES(?,?,?,?,?,?,?,?)
        ON CONFLICT(run_dir) DO UPDATE SET
          protein_id=excluded.protein_id,
          created_utc=excluded.created_utc,
          status=excluded.status,
          report_html_path=excluded.report_html_path,
          evaluation_json_path=excluded.evaluation_json_path,
          error=excluded.error,
          updated_utc=excluded.updated_utc;
        """,
        (
            protein_id,
            str(run_dir),
            created_utc,
            status,
            str(report_html_path) if report_html_path else None,
            str(evaluation_json_path) if evaluation_json_path else None,
            error,
            now,
        ),
    )
    row = con.execute("SELECT id FROM runs WHERE run_dir = ?;", (str(run_dir),)).fetchone()
    assert row is not None
    con.commit()
    return int(row["id"])


def upsert_run_metrics(con: sqlite3.Connection, *, run_id: int, metrics: dict[str, Any]) -> None:
    now = _now_utc()
    con.execute(
        """
        INSERT INTO run_metrics(run_id, metrics_json, updated_utc)
        VALUES(?,?,?)
        ON CONFLICT(run_id) DO UPDATE SET
          metrics_json=excluded.metrics_json,
          updated_utc=excluded.updated_utc;
        """,
        (run_id, json.dumps(metrics, ensure_ascii=False, indent=2), now),
    )
    con.commit()


def fts_available(con: sqlite3.Connection) -> bool:
    try:
        con.execute("SELECT 1 FROM proteins_fts LIMIT 1;").fetchone()
        return True
    except sqlite3.OperationalError:
        return False


def search_proteins(con: sqlite3.Connection, query: str, *, limit: int = 50) -> list[sqlite3.Row]:
    q = (query or "").strip()
    if not q:
        return con.execute(
            "SELECT * FROM proteins ORDER BY updated_utc DESC LIMIT ?;",
            (int(limit),),
        ).fetchall()
    if fts_available(con):
        # FTS5 query syntax; simplest: prefix match on tokens.
        fts_q = " ".join([f"{t}*" for t in q.split() if t])
        return con.execute(
            """
            SELECT p.* FROM proteins_fts f
            JOIN proteins p ON p.id = f.rowid
            WHERE proteins_fts MATCH ?
            ORDER BY p.updated_utc DESC
            LIMIT ?;
            """,
            (fts_q, int(limit)),
        ).fetchall()
    like = f"%{q.lower()}%"
    return con.execute(
        """
        SELECT * FROM proteins
        WHERE lower(uniprot_id) LIKE ?
           OR lower(coalesce(gene_symbol,'')) LIKE ?
           OR lower(coalesce(preset_key,'')) LIKE ?
           OR lower(synonyms_json) LIKE ?
        ORDER BY updated_utc DESC
        LIMIT ?;
        """,
        (like, like, like, like, int(limit)),
    ).fetchall()


def latest_run_for_protein(con: sqlite3.Connection, protein_id: int) -> sqlite3.Row | None:
    return con.execute(
        """
        SELECT * FROM runs
        WHERE protein_id = ?
        ORDER BY updated_utc DESC
        LIMIT 1;
        """,
        (int(protein_id),),
    ).fetchone()


def get_run_metrics(con: sqlite3.Connection, run_id: int) -> dict[str, Any] | None:
    row = con.execute("SELECT metrics_json FROM run_metrics WHERE run_id = ?;", (int(run_id),)).fetchone()
    if row is None:
        return None
    try:
        return json.loads(str(row["metrics_json"]))
    except Exception:
        return None


def iter_run_dirs(runs_root: Path) -> Iterable[Path]:
    if not runs_root.is_dir():
        return []
    for protein_dir in sorted([p for p in runs_root.iterdir() if p.is_dir()]):
        for run_dir in sorted([p for p in protein_dir.iterdir() if p.is_dir()]):
            yield run_dir

