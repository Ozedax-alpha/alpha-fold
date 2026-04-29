param(
    [string]$RunsRoot = "runs_smoke",
    [switch]$SkipReport
)

$ErrorActionPreference = "Stop"

if ($SkipReport) {
    Write-Host "== Local smoke: init -> run(resume) -> (skip report) -> gc(dry-run) =="
} else {
    Write-Host "== Local smoke: init -> run(resume) -> report -> gc(dry-run) =="
}

# 1) Create an isolated run dir (no network required).
py -3 -m avi init --preset tp53 --runs-parent $RunsRoot | Out-Host

$targetRoot = Join-Path $RunsRoot "tp53"
if (-not (Test-Path $targetRoot)) {
    throw "Missing expected target root: $targetRoot"
}

$runDir = Get-ChildItem $targetRoot -Directory | Sort-Object LastWriteTime -Descending | Select-Object -First 1
if ($null -eq $runDir) {
    throw "Could not resolve newly created run directory under $targetRoot"
}
$run = $runDir.FullName
Write-Host "Using run dir: $run"

# 2) Seed minimal processed artifact so report preflight passes.
$processed = Join-Path $run "data\processed"
New-Item -ItemType Directory -Path $processed -Force | Out-Null
$csvPath = Join-Path $processed "tp53_variants_basic.csv"
"id,molecular_consequence,protein_position,in_alphafold_structure" | Set-Content -Path $csvPath -Encoding UTF8
"1,missense variant,10,True" | Add-Content -Path $csvPath -Encoding UTF8

# 3) Run resume on repro stage only (offline-safe).
py -3 -m avi run --run-dir $run --resume --stage repro | Out-Host

# 4) Optional: create a tiny notebook for fast report execution.
if (-not $SkipReport) {
    $nbPath = Join-Path $run "smoke_report.ipynb"
    $nbJson = @"
{
  "cells": [
    {
      "cell_type": "code",
      "execution_count": null,
      "metadata": {},
      "outputs": [],
      "source": [
        "print('smoke report ok')\n"
      ]
    }
  ],
  "metadata": {
    "kernelspec": {
      "display_name": "Python 3",
      "language": "python",
      "name": "python3"
    },
    "language_info": {
      "name": "python"
    }
  },
  "nbformat": 4,
  "nbformat_minor": 5
}
"@
    $nbJson | Set-Content -Path $nbPath -Encoding UTF8

    py -3 -m avi report --run-dir $run --notebook $nbPath | Out-Host
} else {
    Write-Host "Skipping report step (-SkipReport provided)."
}

# 5) GC preview only.
py -3 -m avi gc --runs-root $RunsRoot --delete-partial --older-than-days 0 --dry-run | Out-Host

Write-Host "Smoke completed successfully."
