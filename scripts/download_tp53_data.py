"""
Download human TP53 AlphaFold structure and ClinVar-associated records.

The Variation API gene endpoint used in early drafts is often unavailable (404).
This script uses NCBI E-utilities (esearch + batched esummary) instead, which is
stable for programmatic ClinVar access.
"""

from __future__ import annotations

import argparse
import json
import os
import time

import requests

from pipeline_runtime import RAW_DIR, load_config

RAW_DIR.mkdir(parents=True, exist_ok=True)

_CFG = load_config()
UNIPROT_ID = str(_CFG["uniprot_id"])
_FRAG = str(_CFG["alphafold_fragment"])
_OUT_BASE = str(_CFG["output_basename"])
# AlphaFold DB model versions change over time; resolve the PDB URL from the API.
AF_PREDICTION_API = f"https://alphafold.ebi.ac.uk/api/prediction/{UNIPROT_ID}"
AF_ENTRY_ID = f"AF-{UNIPROT_ID}-{_FRAG}"
AFDB_OUT = RAW_DIR / f"AF-{UNIPROT_ID}-{_FRAG}-alphafold.pdb"

CLINVAR_OUT = RAW_DIR / f"clinvar_{_OUT_BASE}_variants.json"

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "alphafold-variant-interpretation"
ESEARCH_TERM = str(_CFG["clinvar_esearch_term"])
ESUMMARY_BATCH = 200
REQUEST_PAUSE_SEC = 0.35  # stay under ~3 req/s without an API key


def _eutils_params(extra: dict) -> dict:
    params = {"tool": TOOL, **extra}
    email = os.environ.get("ENTREZ_EMAIL")
    if email:
        params["email"] = email
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        params["api_key"] = api_key
    return params


def download_alphafold_structure() -> None:
    if AFDB_OUT.exists():
        print(f"AlphaFold structure already exists at {AFDB_OUT}")
        return
    print(f"Resolving AlphaFold PDB URL for {UNIPROT_ID} ({AF_ENTRY_ID})...")
    meta = requests.get(AF_PREDICTION_API, timeout=60)
    meta.raise_for_status()
    predictions = meta.json()
    if not isinstance(predictions, list):
        raise RuntimeError("Unexpected AlphaFold API response (expected a JSON list).")
    chosen = None
    for p in predictions:
        if not isinstance(p, dict):
            continue
        if p.get("entryId") == AF_ENTRY_ID and p.get("uniprotAccession") == UNIPROT_ID:
            chosen = p
            break
    if chosen is None:
        raise RuntimeError(
            f"Could not find canonical entry {AF_ENTRY_ID} in AlphaFold API response."
        )
    pdb_url = chosen.get("pdbUrl")
    if not pdb_url:
        raise RuntimeError("AlphaFold API record has no pdbUrl.")
    ver = chosen.get("latestVersion", "?")
    print(f"Downloading AlphaFold PDB (model version {ver})...")
    resp = requests.get(pdb_url, timeout=120)
    resp.raise_for_status()
    AFDB_OUT.write_bytes(resp.content)
    print(f"Saved AlphaFold PDB to {AFDB_OUT}")


def _esearch_clinvar_tp53_ids(session: requests.Session) -> list[str]:
    params = _eutils_params(
        {
            "db": "clinvar",
            "term": ESEARCH_TERM,
            "retmax": 20000,
            "retmode": "json",
        }
    )
    r = session.get(f"{EUTILS_BASE}/esearch.fcgi", params=params, timeout=60)
    r.raise_for_status()
    data = r.json()
    idlist = data.get("esearchresult", {}).get("idlist", [])
    if not idlist:
        raise RuntimeError(f"ESearch returned no ClinVar IDs for {ESEARCH_TERM!r}.")
    return idlist


def _esummary_batch(
    session: requests.Session, ids: list[str]
) -> dict[str, dict]:
    params = _eutils_params(
        {
            "db": "clinvar",
            "id": ",".join(ids),
            "retmode": "json",
        }
    )
    r = session.get(f"{EUTILS_BASE}/esummary.fcgi", params=params, timeout=120)
    r.raise_for_status()
    payload = r.json()
    result = payload.get("result", {})
    uids = result.get("uids", [])
    out: dict[str, dict] = {}
    for uid in uids:
        rec = result.get(uid)
        if isinstance(rec, dict):
            out[str(uid)] = rec
    return out


def download_clinvar_tp53() -> None:
    if CLINVAR_OUT.exists():
        print(f"ClinVar TP53 data already exists at {CLINVAR_OUT}")
        return
    print(f"Downloading ClinVar records ({ESEARCH_TERM}) via E-utilities...")
    with requests.Session() as session:
        ids = _esearch_clinvar_tp53_ids(session)
        summaries: dict[str, dict] = {}
        for i in range(0, len(ids), ESUMMARY_BATCH):
            batch = ids[i : i + ESUMMARY_BATCH]
            summaries.update(_esummary_batch(session, batch))
            time.sleep(REQUEST_PAUSE_SEC)
    doc = {
        "source": "ncbi_eutils_esummary",
        "db": "clinvar",
        "gene_query": ESEARCH_TERM,
        "id_count": len(ids),
        "summaries": summaries,
    }
    CLINVAR_OUT.write_text(json.dumps(doc, indent=2), encoding="utf-8")
    print(f"Saved ClinVar summaries for {len(ids)} records to {CLINVAR_OUT}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--force",
        action="store_true",
        help="Redownload even if raw files already exist (overwrites PDB/ClinVar JSON).",
    )
    args = ap.parse_args()

    if args.force:
        if AFDB_OUT.exists():
            AFDB_OUT.unlink()
            print(f"Removed existing {AFDB_OUT} (--force)")
        if CLINVAR_OUT.exists():
            CLINVAR_OUT.unlink()
            print(f"Removed existing {CLINVAR_OUT} (--force)")

    download_alphafold_structure()
    download_clinvar_tp53()
    print("Done.")


if __name__ == "__main__":
    main()
