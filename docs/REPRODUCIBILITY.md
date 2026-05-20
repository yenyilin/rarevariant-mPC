# Reproducibility

Two reproducibility paths are documented below.

1. **Verify the pipeline runs** — uses the synthetic cohort under
   `data/synthetic/`. Quick smoke test; no real-data dependencies. Numerical
   outputs will **not** match the published values, by design.

2. **Reproduce the published numerical results** — uses the real cohort
   tables that ship with this repository (`data/clinical/` and `data/raw/`).
   No controlled-access download is needed for any of the listed commands.

The only stages that genuinely require controlled-access data live
*upstream* of this repository (primary variant calling from FASTQ; raw
wet-lab image stacks). See the "Out of scope" section below.

## Verify the pipeline runs (synthetic data)

```bash
# Generate the schema-matched synthetic cohort
python data/synthetic/generate_synthetic_clinical.py

# Build the analysis-ready table from the synthetic inputs
python code/01_cohort_prep/prepare_epc_cohort.py --synthetic

# Run every analysis (numerical outputs will NOT match published values)
python code/03_ddr_enrichment/ddr_burden_bootstrap.py --arm both
python code/03_ddr_enrichment/ddr_arm_compare.py
python code/03_ddr_enrichment/ddr_per_patient_test.py
python code/04_firth_regression/ddr_logistic_regression.py
python code/05_time_bias/ddr_followup_plot.py
python code/05_time_bias/ddr_burden_vs_time.py
python code/06_replication/australian_binomial.py
```

What you get:
- A working environment confirmation
- Plots in `figures/` (synthetic-data versions)

## Reproduce the published numerical results

The default invocation of every script reads the real cohort tables
committed under `data/clinical/` and `data/raw/`:

```bash
python code/01_cohort_prep/prepare_epc_cohort.py
python code/03_ddr_enrichment/ddr_burden_bootstrap.py --arm both
python code/03_ddr_enrichment/ddr_arm_compare.py
python code/03_ddr_enrichment/ddr_per_patient_test.py
python code/04_firth_regression/ddr_logistic_regression.py
python code/04_firth_regression/ddr_logistic_regression.py --extended
python code/05_time_bias/ddr_followup_plot.py
python code/05_time_bias/ddr_burden_vs_time.py
python code/06_replication/australian_binomial.py
```

Expected output values are listed in `docs/METHODS_CROSSREF.md`.

## Out of scope for this repository

| Stage | Why it is not here |
|---|---|
| Primary variant calling from FASTQ | Uses standard tools (BWA / Strelka2 / Picard) with the parameters documented in `code/02_variant_calling/README.md`. Re-running it requires controlled-access raw data, deposited to EGA. |
| Raw wet-lab image stacks | Protocols are in Supplementary Methods §5. The *quantified* assay readings (CCK-8, scratch wound, Boyden chamber, olaparib dose–response) are committed under `data/raw/functional_assays/` and consumed by `code/07_functional/`. The raw image stacks behind these quantifications are deposited with the controlled-access archive. |
| Manuscript text | Not part of the analytical reproducibility surface. |

## Software environment

- Python 3.11 (versions pinned in `environment.yml` / `requirements.txt`)
- Tested on macOS 14.x–15.x; should work on any POSIX-compatible system with Python 3.11. Linux/WSL users have run earlier versions of the
  pipeline.
- A Docker image is **not** currently provided. If you encounter
  environment-specific issues, please open a GitHub issue.

## Why expected values may differ slightly

- Bootstrap p-values vary between runs unless the seed is fixed; we use
  `--seed 42` throughout. Set the same seed to reproduce exactly.
- Profile-likelihood CI bisection terminates at 1e-4 tolerance; the
  reported CI bounds (3.07 / 3528.42) are stable to two decimals.
- BLAS implementation differences can change the last digit of OR.
  Manuscript values were computed on macOS with `numpy=2.3.3` (installed via
  pip into python.org Python 3.11; Accelerate BLAS).
