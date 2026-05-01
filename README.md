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
py -3 -m avi report --run-dir runs/tp53/<timestamp> --notebook notebooks/01_target_exploration.ipynb
py -3 -m avi open-report --run-dir runs/tp53/<timestamp>
py -3 -m avi list-runs
py -3 -m avi repro --run-dir runs/tp53/<timestamp>
py -3 -m avi clean --keep 20
py -3 -m avi gc --delete-partial --older-than-days 7 --dry-run
py -3 -m avi gc --delete-partial --older-than-days 7 --yes
py -3 -m avi batch --presets tp53,insulin --continue-on-error
py -3 -m avi add-preset brca1 --from-uniprot P38398 --auto --overwrite
py -3 -m avi batch --from-uniprot P04637,P38398 --continue-on-error
```

Notes:
- `avi init` creates a timestamped run directory under `runs/` containing `pipeline_config.json`.
- `avi run --run-dir ...` writes outputs to `<run-dir>/data/{raw,processed}/` via `PIPELINE_CONFIG_PATH` / `PIPELINE_OUTPUT_DIR`.
- `avi report` defaults to `notebooks/01_target_exploration.ipynb`; use `--notebook` to override.
- `avi gc` safely garbage-collects broken/partial runs (`--dry-run` by default unless `--yes`).
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

Install the `avi` command on your PATH (editable install):

```bash
pip install -e .
# optional HTTP cache for iterative downloads:
pip install -e ".[cache]"
```

### Environment variables

Copy `.env.example` to `.env` (or export in your shell). Important values:

| Variable | Purpose |
|----------|---------|
| `ENTREZ_EMAIL` | Strongly recommended for [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25497/). |
| `NCBI_API_KEY` | Optional; raises allowed request rate. |
| `AVI_USE_HTTP_CACHE` | Set to `1` to cache AlphaFold / NCBI GETs under `.cache/avi-http/` (requires `pip install -e ".[cache]"`). |

### Refresh older run directories

If you created runs **before** confidence JSON / SASA / evaluation columns existed, re-download and reprocess in place:

```bash
py -3 -m avi run --run-dir runs/brca1/<timestamp> --force-download
py -3 -m avi report --run-dir runs/brca1/<timestamp>
```

Gene / UniProt / output names are set in **`pipeline_config.json`** (defaults target TP53). Optional keys:

- `preferred_transcript_prefix` — e.g. `NM_000546` so multi-transcript `variation_name` strings are sliced before HGVS parsing.
- `uniprot_residue_offset` — integer added to ClinVar protein index before matching the AlphaFold PDB (rare isoform cases).
- `analysis_region_start` / `analysis_region_end` / `analysis_region_label` — notebook “paper slice” filters (TP53 defaults approximate the DNA-binding domain, residues 102–292).

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
- `data/processed/tp53_missense_mappable.csv` (per-residue **pLDDT**, **±5-residue mean pLDDT**, **Shrake–Rupley SASA (Å²)**, and **CA distance to centroid** when PDB/confidence are present)

**Note:** An older ClinVar Variation API path (`.../beta/clinvar/gene/TP53`) often returns 404. The downloader uses **esearch + esummary** instead, which matches the JSON shape expected by `process_tp53_variants.py`.

## Held-out evaluation (not clinical)

After `*_missense_mappable.csv` exists:

```bash
py -3 -m avi evaluate --run-dir runs/tp53/<timestamp>
```

Writes `data/processed/evaluation_metrics.json` with Cohen’s *d* (train) and a rank AUC (test) for each numeric feature. Uses **time-based** holdout when `germline_date_last_evaluated` is populated for enough rows; otherwise a **stratified random** split (`--seed`).

## Notebook

Open `notebooks/01_target_exploration.ipynb` after the two scripts have run to explore the CSV and structure locally. It includes optional **analysis region** figures when `analysis_region_*` is set in config (TP53 defaults target the DNA-binding domain slice).

## Quality checks

- Local end-to-end smoke (isolated under `runs_smoke/`):

```powershell
powershell -ExecutionPolicy Bypass -File scripts/smoke_local.ps1
```

To skip the report/notebook execution (useful if nbconvert/Jupyter isn’t installed):

```powershell
powershell -ExecutionPolicy Bypass -File scripts/smoke_local.ps1 -SkipReport
```

- CI runs `pytest` automatically on pushes to `main` and pull requests via `.github/workflows/tests.yml`.
- Optional **network** smoke (AlphaFold + ClinVar download for the `insulin` preset): `.github/workflows/tests-network.yml` (`workflow_dispatch` or weekly cron). Uses a placeholder `ENTREZ_EMAIL`; set a real address in the workflow or fork settings if NCBI tightens checks.

## Roadmap checkpoints (recommended order)

1. **Golden TP53 run** — one reproducible run directory plus HTML report:
   - `py -3 -m avi doctor`
   - `py -3 -m avi init --preset tp53 --run --report` (add `--open-report` to open the browser)
2. **Other genes** — confirm presets / `--from-uniprot` paths, e.g. `brca1` and `py -3 -m avi init --from-uniprot P38398 --run --report`.
3. **Structure-backed columns** — downloads now pull AlphaFold **per-residue confidence** (`plddtDocUrl`) when available; the missense subset adds **pLDDT**, **±5-residue mean pLDDT**, and **CA distance to the global centroid** (crude burial proxy).
4. **Hardening** — E-utilities pacing speeds up when `NCBI_API_KEY` is set; download errors include short “what to do next” hints; `avi batch` prints a short failure summary when any run fails.
5. **Evaluation** — see `notebooks/01_target_exploration.ipynb` (cells after the missense preview): exploratory **pathogenic + likely pathogenic** vs **benign + likely benign** boxplots on pLDDT and centroid distance. Re-run `avi report --run-dir ...` after refreshing data so the executed HTML includes those figures.

**Hypothesis (exploratory, not clinical):** pathogenic missense variants may skew toward different AlphaFold confidence or geometry than benign variants for some genes; treat any pattern as hypothesis-generating unless validated on held-out data and external evidence.
