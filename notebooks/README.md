# Notebooks (marimo)

Reviewer-facing interactive walkthroughs that complement the production
scripts in `code/`. Each notebook tells a self-contained story for one
analysis and lets reviewers perturb a key parameter to see the result
move.

We use [marimo](https://marimo.io/) instead of Jupyter for two reasons:

1. **Reviewable diffs.** Marimo notebooks are stored as plain `.py`
   files, so `git diff` and code review actually work. Jupyter
   `.ipynb` JSON diffs are notoriously unreadable.
2. **Reactive execution.** Cells re-run automatically when their
   inputs change; there is no hidden out-of-order state.

## Notebooks

| File | What it shows | Manuscript reference |
|---|---|---|
| `01_ddr_convergent_walkthrough.py` | The central finding, end-to-end on the real committed cohort. Walks through three statistically independent tests (per-patient Mann-Whitney burden, gene-set bootstrap competitive null, Firth penalized regression) and shows how they converge on the same pathway-level conclusion. Reproduces the published p = 4.57 × 10⁻⁶, p < 1 × 10⁻⁴, and OR = 26.80 by calling the production scripts in `code/` and reading their outputs. | Methods Methods §§2.3.1–2.3.2; Supplementary Methods §§4.2 and 4.7; Results ¶2–3; Abstract |

## Running

```bash
# Edit interactively (auto-reactive, with the cell graph)
marimo edit notebooks/01_ddr_convergent_walkthrough.py

# Read-only app for reviewers
marimo run  notebooks/01_ddr_convergent_walkthrough.py

# Or open as a regular Python script — this works because marimo files
# are valid Python:
python  notebooks/01_ddr_convergent_walkthrough.py
```

`marimo run` produces a web app at `http://localhost:2718`. To share
with a reviewer over the network, add `--host 0.0.0.0 --port 2718` and
expose the port — or deploy to marimo's hosted runner / Hugging Face
Spaces / a static export (`marimo export html-wasm`).

## Reproducing the manuscript values

The notebook reproduces the manuscript's published values out of the
box — the real, de-identified cohort tables it reads (`data/clinical/`,
`data/raw/`) are already committed to this repository. The only
controlled-access pieces are the raw FASTQ/BAM/VCF *upstream* of these
aggregate matrices (deposited to a public controlled-access repository 
before manuscript publication; see `docs/DATA_AVAILABILITY.md`), which 
the notebook does not need.

A `--synthetic` path is also available end-to-end via the same
production scripts; numerical results will differ from the published
values by design (the synthetic cohort encodes a deliberately weaker
effect). See `docs/REPRODUCIBILITY.md`.
