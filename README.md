# Pathway-Level Germline DDR Rare Variants in Metastatic Prostate Cancer — Analysis Code

[![DOI](https://zenodo.org/badge/1235749000.svg)](https://doi.org/10.5281/zenodo.20317985)

This repository contains the analysis code and de-identified derived data
supporting the manuscript:

> **Pathway-Level Germline Rare Variant Enrichment Across DNA Damage
> Response Genes in Metastatic versus Non-Metastatic Prostate Cancer:
> An Extreme Phenotype Study.**
> Lin Y-Y et al., submitted to *European Urology* (2026).
> Preprint: https://doi.org/10.1101/2025.04.28.25326584
> *(v1 — an updated preprint matching this submission will be posted to medRxiv on acceptance.)*

Most genetic studies of prostate cancer have addressed *cancer initiation*;
this study targets a distinct phenotype of **metastatic progression after
diagnosis** using an extreme-phenotype cohort design that sequences
patients at opposite ends of the post-diagnosis course. The finding is a
**pathway-level** germline rare-variant burden across DNA Damage Response
(DDR) genes that distinguishes the two arms, with convergent
contributions from *multiple DDR genes rather than a single dominant
variant*. The pathway-level signal replicates in an independent Australian
extreme-phenotype cohort and at the pathway level in the Pan Prostate
Cancer Group (n = 976). In vitro functional studies of variants currently
classified as uncertain or benign in clinical databases reveal measurable
phenotypic effects, pointing to patients who may benefit from DDR-targeted
therapy but are not flagged by current variant-classification criteria.

A persistent DOI for the code and derived data is minted via Zenodo and
listed in the manuscript's Data Availability statement.

## Repository contents

### In this repo

- All scripts used to generate the statistical results, tables, and figures
  reported in the manuscript.
- The **de-identified clinical cohort table** (`data/clinical/clinical.tsv`,
  52 patients) — the real Table 1 cohort, cleared for public release.
- **Primary pipeline input tables** (`data/raw/`) — per-patient and per-gene
  aggregate matrices (gnsRV burden, oncoplot matrices, cBioPortal exports,
  quantified wet-lab assay readings).
- A **synthetic cohort** (`data/synthetic/`) matching the schema of the
  real data, for end-to-end smoke testing without controlled-access
  credentials.
- **De-identified summary tables** (`data/derived/`): the prespecified DDR
  gene panel, the 56-variant genome-wide screen list (with gnomAD MAFs,
  FATHMM-MKL, EVEE, and ClinVar annotations), Australian EPC carrier
  counts, and PPCG gene-level association summary.
- A **claim-to-script crossreference** (`docs/METHODS_CROSSREF.md`) that
  maps each numerical result in the manuscript to the file that produces it.

### Not in this repo

Raw FASTQ, BAM, VCF, and raw wet-lab image stacks are **not** in this
repository. They remained controlled-access and are being deposited to 
a public controlled-access repository before manuscript publication (see 
  `docs/DATA_AVAILABILITY.md`). The aggregate matrices and quantified 
readings that the analysis pipeline actually *consumes* are committed 
under `data/clinical/` and `data/raw/`

### A note on tools

Every stage uses standard, widely validated tools  (e.g. BWA-MEM, 
Strelka2, ANNOVAR, Firth penalized regression, SKAT) so the findings do 
not depend on bespoke, unvalidated code and can be reproduced with any
equivalent workflow. The methodological contribution is the extreme-phenotype 
cohort design and the rare-variant prioritization strategy, not new software.


## Quickstart

Python 3.11 with `pandas`, `numpy`, `scipy`, `matplotlib`, `openpyxl`,
`marimo` (for `notebooks/`), and `pytest` (for the smoke tests). Exact
versions are pinned in `environment.yml` / `requirements.txt`. Nothing
in `code/` requires R: the PPCG gene-level SKAT tests are run by the
PPCG consortium, not from this repo.

```bash
git clone https://github.com/yenyilin/rarevariant-mPC.git
cd rarevariant-mPC

# Option A — conda (canonical, matches environment.yml exactly)
conda env create -f environment.yml
conda activate rarevariant-mpc

# Option B — pip on python.org Python 3.11 (the setup the maintainer actually used)
python3.11 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# Option C — uv (much faster than pip; same requirements.txt)
uv venv --python 3.11
source .venv/bin/activate
uv pip install -r requirements.txt

# 1. Generate the synthetic cohort
python data/synthetic/generate_synthetic_clinical.py

# 2. Build the analysis-ready table
python code/01_cohort_prep/prepare_epc_cohort.py --synthetic

# 3. Reproduce Firth regression (synthetic OR will not match published OR=26.80)
python code/04_firth_regression/ddr_logistic_regression.py

# 4. Run the test suite
pytest

# 5. Open an interactive walkthrough (marimo)
marimo edit notebooks/01_ddr_convergent_walkthrough.py
```

To reproduce the published numerical results, just re-run steps 2 onward
without the `--synthetic` flag — the real cohort tables are already
committed under `data/clinical/` and `data/raw/`. See
`docs/REPRODUCIBILITY.md` for the full command list.

## Reproducing specific manuscript claims

See `docs/METHODS_CROSSREF.md` for the full table. Examples:

| Manuscript claim | Script |
|---|---|
| DDR gnsRV enrichment, p = 4.57 × 10⁻⁶ | `code/03_ddr_enrichment/ddr_arm_compare.py` |
| Firth OR = 26.80, 95% CI 3.07–3528.42 | `code/04_firth_regression/ddr_logistic_regression.py` |
| Mann-Whitney follow-up comparison, p < 0.001 | `code/05_time_bias/ddr_followup_plot.py` |
| Australian EPC binomial, p = 0.027–0.050 | `code/06_replication/australian_binomial.py` |

## Repository layout

```
.
├── code/
│   ├── 01_cohort_prep/        Build analysis table from clinical + genotype inputs
│   ├── 02_variant_calling/    BWA / Strelka2 / annotation pipeline (documented)
│   ├── 03_ddr_enrichment/     DDR pathway-burden tests (Mann-Whitney, gene-set bootstrap, per-patient p)
│   ├── 04_firth_regression/   Firth penalized logistic regression (OR=26.80)
│   ├── 05_time_bias/          Mann-Whitney + Spearman time-bias analyses
│   ├── 06_replication/        Australian EPC binomial (PPCG SKAT done by consortium)
│   ├── 07_functional/         Quantification of in vitro assays
│   └── 08_figures/            One script per main and supplementary figure
├── data/
│   ├── clinical/              De-identified clinical cohort table (52 patients)
│   ├── raw/                   Primary pipeline inputs (per-patient aggregate matrices)
│   ├── derived/               Aggregate / de-identified summary tables
│   ├── synthetic/             Simulated cohort matching the real schema
│   └── processed/             Pipeline outputs (gitignored; regenerated at runtime)
├── docs/
│   ├── DATA_AVAILABILITY.md   Accession numbers and access procedure
│   ├── PRECOMMIT_CHECKLIST.md Pre-release checks (PHI / controlled-data scan)
│   └── METHODS_CROSSREF.md    Manuscript claim → script crossreference
├── notebooks/                 Reviewer-facing marimo walkthroughs (.py)
├── tests/                     Sanity checks runnable on synthetic data
└── .github/workflows/         CI (pytest on synthetic data)
```

## Citation

```
@article{rarevariant_mpca_2026,
  title  = {Pathway-Level Germline Rare Variant Enrichment Across DNA
            Damage Response Genes in Metastatic versus Non-Metastatic
            Prostate Cancer: An Extreme Phenotype Study},
  author = {Lin, Y.-Y. and others},
  year   = {2026},
  note   = {Preprint: doi:10.1101/2025.04.28.25326584;
            Code and derived data: doi:10.5281/zenodo.20317985}
}
```

## License

Code is released under the MIT License (`LICENSE`). Documentation and the
public data tables (`data/derived/`, `data/clinical/`, `data/synthetic/`)
are released under CC-BY-4.0. Raw sequencing data deposited externally
remains under controlled access governed by the original consent.

## Contact

Issues and questions: please open a GitHub issue.
Data access: see `docs/DATA_AVAILABILITY.md`.
