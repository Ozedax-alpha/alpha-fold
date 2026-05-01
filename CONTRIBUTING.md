# Contributing

## Environment

- Python **3.11+** (CI uses 3.13).
- From the repo root:

```bash
pip install -e ".[dev]"
```

Optional HTTP cache for repeated downloads: `pip install -e ".[cache]"` and set `AVI_USE_HTTP_CACHE=1` (see `.env.example`).

## Pre-commit

Format and lightweight checks before commit:

```bash
pre-commit install
pre-commit run --all-files
```

## Tests

```bash
pytest -q
```

## NCBI E-utilities

For any command that downloads ClinVar data, set **`ENTREZ_EMAIL`** to a real address (see [NCBI usage guidelines](https://www.ncbi.nlm.nih.gov/books/NBK25497/)). **`NCBI_API_KEY`** is optional and increases allowed request rates.

## Minimal repro (one gene)

1. `py -3 -m avi doctor`
2. `py -3 -m avi init --preset tp53 --run` (or another preset from `avi/presets.json`)
3. Expected artifacts under `<run-dir>/data/raw/` and `data/processed/`:
   - `AF-<UniProt>-F1-alphafold.pdb`, confidence JSON when available
   - `clinvar_<basename>_variants.json`
   - `<basename>_variants_basic.csv`, `<basename>_missense_mappable.csv`
4. Optional: `py -3 -m avi report --run-dir <run-dir>` and `py -3 -m avi evaluate --run-dir <run-dir>`

## CI

- **tests.yml** — `pytest` on push/PR (`pip install -e .`).
- **tests-network.yml** — optional full download for the `insulin` preset (`workflow_dispatch` or weekly cron).
