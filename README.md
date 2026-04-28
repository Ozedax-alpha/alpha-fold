# AlphaFold variant interpretation (TP53)

Small research scaffold: download human **TP53** (UniProt `P04637`) from AlphaFold DB, pull **ClinVar** records for the gene via NCBI E-utilities, and produce a tidy table of variant fields plus protein-position checks against the PDB.

## Recommended UX: one front door (`avi`)

Run everything via the CLI:

```bash
py -3 -m avi doctor
py -3 -m avi init --preset tp53
py -3 -m avi init --from-uniprot P38398 --run --report --open-report
py -3 -m avi add-preset brca1 --uniprot-id P38398 --gene-symbol BRCA1 --clinvar-esearch-term "BRCA1[gene]" --output-basename brca1 --alphafold-fragment F1
py -3 -m avi run --run-dir runs/tp53/<timestamp>
py -3 -m avi explain --run-dir runs/tp53/<timestamp>
py -3 -m avi notebook --run-dir runs/tp53/<timestamp>
py -3 -m avi report --run-dir runs/tp53/<timestamp>
py -3 -m avi open-report --run-dir runs/tp53/<timestamp>
py -3 -m avi list-runs
py -3 -m avi repro --run-dir runs/tp53/<timestamp>
py -3 -m avi clean --keep 20
py -3 -m avi batch --presets tp53,insulin --continue-on-error
py -3 -m avi add-preset brca1 --from-uniprot P38398 --auto --overwrite
py -3 -m avi batch --from-uniprot P04637,P38398 --continue-on-error
```

Notes:
- `avi init` creates a timestamped run directory under `runs/` containing `pipeline_config.json`.
- `avi run --run-dir ...` writes outputs to `<run-dir>/data/{raw,processed}/` via `PIPELINE_CONFIG_PATH` / `PIPELINE_OUTPUT_DIR`.
- Presets live in `avi/presets.json` (edit to add targets).

## Layout

- `data/raw/` — AlphaFold PDB and ClinVar JSON from the download script
- `data/processed/` — CSV produced by the processing script
- `scripts/` — download and processing entry points
- `notebooks/` — exploratory notebook starter

## Setup

Create a virtual environment (recommended), then:

```bash
pip install -r requirements.txt
```

For a **frozen** environment matching a prior run: `pip install -r requirements.lock.txt`.

Optional: set `ENTREZ_EMAIL` (and `NCBI_API_KEY` if you have one) for [NCBI E-utilities etiquette](https://www.ncbi.nlm.nih.gov/books/NBK25497/).

Gene / UniProt / output names are set in **`pipeline_config.json`** (defaults target TP53).

## Run

From the project root, **all stages** (download → process → missense subset):

```bash
python scripts/run_pipeline.py
```

Or step-by-step:

```bash
python scripts/download_tp53_data.py
python scripts/process_tp53_variants.py
python scripts/build_missense_subset.py
```

Outputs (with default `output_basename` `tp53`):

- `data/raw/AF-P04637-F1-alphafold.pdb` (PDB URL from AlphaFold API)
- `data/raw/clinvar_tp53_variants.json`
- `data/processed/tp53_variants_basic.csv`
- `data/processed/tp53_missense_mappable.csv` (adds pLDDT columns if `data/raw/AF-P04637-F1-confidence*.json` exists)

**Note:** An older ClinVar Variation API path (`.../beta/clinvar/gene/TP53`) often returns 404. The downloader uses **esearch + esummary** instead, which matches the JSON shape expected by `process_tp53_variants.py`.

## Notebook

Open `notebooks/01_tp53_exploration.ipynb` after the two scripts have run to explore the CSV and structure locally.
