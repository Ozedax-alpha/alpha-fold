"""Lightweight local web UI (search + protein detail) backed by SQLite."""

from __future__ import annotations

import html
import json
import os
import subprocess
import sys
import threading
import socketserver
from http.server import BaseHTTPRequestHandler
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlparse

from avi import db as avi_db


CSS = """
body{font-family:ui-sans-serif,system-ui,-apple-system,Segoe UI,Roboto,Arial;margin:0;background:#0b0f19;color:#e6e9ef}
a{color:#9cd4ff;text-decoration:none} a:hover{text-decoration:underline}
.wrap{max-width:980px;margin:0 auto;padding:20px}
.top{position:sticky;top:0;background:rgba(11,15,25,.92);backdrop-filter:blur(8px);border-bottom:1px solid #1d2640}
.row{display:flex;gap:12px;align-items:center}
input[type=text]{flex:1;padding:10px 12px;border-radius:10px;border:1px solid #2a355a;background:#0f1630;color:#e6e9ef}
button{padding:10px 12px;border-radius:10px;border:1px solid #2a355a;background:#15204a;color:#e6e9ef;cursor:pointer}
.card{border:1px solid #1d2640;background:#0f1630;border-radius:14px;padding:14px;margin:12px 0}
.muted{color:#aab3c5}
.pill{display:inline-block;padding:3px 8px;border-radius:999px;border:1px solid #2a355a;font-size:12px;color:#cbd5e1;margin-right:6px}
.h1{font-size:20px;font-weight:700}
.grid{display:grid;grid-template-columns:1fr 1fr;gap:10px}
.mono{font-family:ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,monospace;font-size:12px}
"""


def _page(title: str, body: str) -> bytes:
    t = html.escape(title)
    return (
        f"<!doctype html><html><head><meta charset='utf-8'/>"
        f"<meta name='viewport' content='width=device-width, initial-scale=1'/>"
        f"<title>{t}</title><style>{CSS}</style></head><body>{body}</body></html>"
    ).encode("utf-8")


def _topbar(q: str) -> str:
    qh = html.escape(q or "")
    return f"""
<div class="top">
  <div class="wrap">
    <div class="row">
      <div class="h1">AVI</div>
      <form action="/" method="get" style="flex:1" class="row">
        <input type="text" name="q" value="{qh}" placeholder="Search proteins (gene, UniProt, preset)…" />
        <button type="submit">Search</button>
      </form>
    </div>
    <div class="muted" style="margin-top:8px">Local search over indexed proteins and latest runs.</div>
  </div>
</div>
"""


def _protein_card(p: Any) -> str:
    gene = (p["gene_symbol"] or "").strip()
    uid = (p["uniprot_id"] or "").strip()
    preset = (p["preset_key"] or "").strip()
    title = gene or preset or uid
    sub = " · ".join([x for x in [preset and f"preset:{preset}", uid and f"UniProt:{uid}"] if x])
    return f"""
<div class="card">
  <div class="row" style="justify-content:space-between">
    <div>
      <div class="h1"><a href="/protein?id={int(p['id'])}">{html.escape(title)}</a></div>
      <div class="muted mono">{html.escape(sub)}</div>
    </div>
    <div class="pill">updated {html.escape(str(p['updated_utc']))}</div>
  </div>
</div>
"""


def _run_button(protein_id: int) -> str:
    return f"""
<form action="/run" method="post" style="margin-top:10px">
  <input type="hidden" name="id" value="{int(protein_id)}"/>
  <button type="submit">Run pipeline</button>
</form>
"""


def _metrics_summary(doc: dict[str, Any]) -> str:
    status = doc.get("status")
    split_mode = doc.get("split_mode")
    sm = doc.get("split_metadata") or {}
    pills = []
    if status:
        pills.append(f'<span class="pill">{html.escape(str(status))}</span>')
    if split_mode:
        pills.append(f'<span class="pill">{html.escape(str(split_mode))}</span>')
    if sm:
        np_ = sm.get("n_pos_all")
        nn_ = sm.get("n_neg_all")
        if np_ is not None and nn_ is not None:
            pills.append(f'<span class="pill">pos/neg: {int(np_)}/{int(nn_)}</span>')
    feats = doc.get("features") or {}
    feat_names = ", ".join(sorted([str(k) for k in feats.keys()])[:6])
    if feat_names:
        pills.append(f'<span class="pill">features: {html.escape(feat_names)}</span>')
    return " ".join(pills) if pills else '<span class="muted">No metrics.</span>'


class Handler(BaseHTTPRequestHandler):
    server_version = "avi-webui/0.1"

    @property
    def base_dir(self) -> Path:
        return Path(getattr(self.server, "base_dir"))

    @property
    def db_path(self) -> Path:
        return Path(getattr(self.server, "db_path"))

    def _send(self, status: int, data: bytes, *, content_type: str = "text/html; charset=utf-8") -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self) -> None:  # noqa: N802
        u = urlparse(self.path)
        path = u.path
        qs = parse_qs(u.query or "")

        if path == "/":
            q = (qs.get("q") or [""])[0]
            con = avi_db.connect(self.db_path)
            try:
                avi_db.init_db(con)
                rows = avi_db.search_proteins(con, q, limit=100)
            finally:
                con.close()
            cards = "".join([_protein_card(r) for r in rows]) or '<div class="card">No results.</div>'
            body = _topbar(q) + f'<div class="wrap">{cards}</div>'
            self._send(200, _page("AVI search", body))
            return

        if path == "/protein":
            pid = (qs.get("id") or [""])[0]
            try:
                protein_id = int(pid)
            except Exception:
                self._send(400, _page("Bad request", _topbar("") + '<div class="wrap"><div class="card">Missing protein id.</div></div>'))
                return

            con = avi_db.connect(self.db_path)
            try:
                p = con.execute("SELECT * FROM proteins WHERE id = ?;", (protein_id,)).fetchone()
                if p is None:
                    self._send(404, _page("Not found", _topbar("") + '<div class="wrap"><div class="card">Protein not found.</div></div>'))
                    return
                r = avi_db.latest_run_for_protein(con, protein_id)
                metrics = None
                if r is not None:
                    metrics = avi_db.get_run_metrics(con, int(r["id"]))
            finally:
                con.close()

            title = (p["gene_symbol"] or p["preset_key"] or p["uniprot_id"] or "protein").strip()
            info = f"""
<div class="card">
  <div class="h1">{html.escape(title)}</div>
  <div class="muted mono">UniProt: {html.escape(str(p['uniprot_id']))} · preset: {html.escape(str(p['preset_key'] or ''))}</div>
  {_run_button(int(p['id']))}
</div>
"""
            run_block = '<div class="card">No runs indexed for this protein.</div>'
            if r is not None:
                report = r["report_html_path"]
                report_link = ""
                if report and Path(str(report)).is_file():
                    report_link = f'<a class="pill" href="/file?path={html.escape(str(report))}">Open report HTML</a>'
                ev = r["evaluation_json_path"]
                ev_link = ""
                if ev and Path(str(ev)).is_file():
                    ev_link = f'<a class="pill" href="/file?path={html.escape(str(ev))}">Open evaluation JSON</a>'
                ms = _metrics_summary(metrics or {})
                run_block = f"""
<div class="card">
  <div class="row" style="justify-content:space-between">
    <div>
      <div class="h1">Latest run</div>
      <div class="muted mono">{html.escape(str(r['run_dir']))}</div>
      <div style="margin-top:8px">{ms}</div>
    </div>
    <div style="text-align:right">
      <div class="pill">status {html.escape(str(r['status']))}</div><br/>
      {report_link}<br/>{ev_link}
    </div>
  </div>
</div>
"""
            body = _topbar("") + f'<div class="wrap">{info}{run_block}<div class="card muted">Tip: run indexing from CLI after new runs.</div></div>'
            self._send(200, _page(f"Protein {title}", body))
            return

        if path == "/file":
            p = (qs.get("path") or [""])[0]
            if not p:
                self._send(400, _page("Bad request", _topbar("") + '<div class="wrap"><div class="card">Missing path.</div></div>'))
                return
            target = Path(p)
            try:
                target = target.resolve()
            except Exception:
                self._send(400, _page("Bad request", _topbar("") + '<div class="wrap"><div class="card">Bad path.</div></div>'))
                return
            # Prevent arbitrary file reads: only allow under repo base_dir or under runs roots.
            base = self.base_dir.resolve()
            allowed = False
            for root in (
                base / "runs",
                base / "runs_smoke",
                base / "runs_batch_smoke",
                base / "runs_batch_eval_smoke",
                base / "runs_dataset",
                base / "runs_dataset_smoke",
            ):
                try:
                    if target.is_relative_to(root.resolve()):  # py3.9+? (we're 3.11+)
                        allowed = True
                        break
                except Exception:
                    pass
            if not allowed:
                self._send(403, _page("Forbidden", _topbar("") + '<div class="wrap"><div class="card">File not allowed.</div></div>'))
                return
            if not target.is_file():
                self._send(404, _page("Not found", _topbar("") + '<div class="wrap"><div class="card">File not found.</div></div>'))
                return
            data = target.read_bytes()
            ct = "text/plain; charset=utf-8"
            if target.suffix.lower() == ".html":
                ct = "text/html; charset=utf-8"
            self._send(200, data, content_type=ct)
            return

        self._send(404, _page("Not found", _topbar("") + '<div class="wrap"><div class="card">Not found.</div></div>'))

    def do_POST(self) -> None:  # noqa: N802
        u = urlparse(self.path)
        if u.path != "/run":
            self._send(404, _page("Not found", _topbar("") + '<div class="wrap"><div class="card">Not found.</div></div>'))
            return

        length = int(self.headers.get("Content-Length") or "0")
        raw = self.rfile.read(length) if length > 0 else b""
        form = parse_qs(raw.decode("utf-8", errors="ignore"))
        pid_s = (form.get("id") or [""])[0]
        try:
            protein_id = int(pid_s)
        except Exception:
            self._send(400, _page("Bad request", _topbar("") + '<div class="wrap"><div class="card">Bad protein id.</div></div>'))
            return

        # Kick off a background run (no blocking the UI thread).
        base = self.base_dir.resolve()
        db_path = self.db_path.resolve()

        def _bg() -> None:
            cmd = [
                sys.executable,
                "-m",
                "avi",
                "dataset",
                "run",
                "--db",
                str(db_path),
                "--protein-ids",
                str(protein_id),
                "--runs-parent",
                "runs_dataset",
                "--report",
                "--evaluate",
                "--continue-on-error",
            ]
            subprocess.call(cmd, cwd=str(base))

        threading.Thread(target=_bg, daemon=True).start()

        # Redirect back to protein page.
        self.send_response(303)
        self.send_header("Location", f"/protein?id={protein_id}")
        self.end_headers()

    def log_message(self, fmt: str, *args: Any) -> None:  # noqa: D401
        if os.environ.get("AVI_WEBUI_QUIET", "").strip() in ("1", "true", "yes"):
            return
        super().log_message(fmt, *args)


def serve(*, base_dir: Path, db_path: Path, host: str, port: int) -> None:
    class _Server(socketserver.ThreadingMixIn, socketserver.TCPServer):
        daemon_threads = True
        allow_reuse_address = True

    with _Server((host, int(port)), Handler) as httpd:
        httpd.base_dir = str(base_dir)
        httpd.db_path = str(db_path)
        print(f"AVI UI: http://{host}:{port}/")
        httpd.serve_forever()

