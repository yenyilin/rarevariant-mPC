# Data

Four committed subdirectories (plus `processed/`, generated at runtime),
with strict rules about what may be committed to this repo.

## `derived/` — committed

Aggregate tables that contain no individual-level genotype or clinical data.
Safe to redistribute under CC-BY-4.0.

| File | Description | Source in manuscript |
|---|---|---|
| `ddr_panel_25genes.tsv` | The 25-gene prespecified DDR panel | Methods §2.3, Fig. 1, Supp Table S4 |
| `gnsrv_56_screen.tsv` | The 56 metastatic-exclusive variants from the genome-wide screen, with gnomAD MAFs, FATHMM-MKL scores, and ClinVar annotations | Results, Fig. 2, Supp Table S4 |
| `australian_epc_counts.tsv` | Per-subgroup carrier counts in the Australian EPC | Results, Fig. 3A  |
| `ppcg_gene_associations.tsv` | Gene-level Wald + Burden + SKAT statistics in PPCG | Results, Fig. 3B |

## `clinical/` — committed

The de-identified clinical cohort table for the 52 EPC patients. All direct
identifiers have been removed; the diagnosis year was dropped, and the
remaining quasi-identifier (metastasis location, anatomical-site level)
was reviewed and cleared for release. Published as an explicit exception
to the `clinical*.tsv` ignore rule. See `data/clinical/README.md`.

- `clinical.tsv` — 52 patients, the real clinical table used by the
  manuscript (Table 1, `prepare_epc_cohort.py` real mode)

## `synthetic/` — committed

A simulated cohort with the same column schema as the real clinical data.
Lets reviewers run every script end-to-end without applying for access.
**Numerical results from synthetic data will not match the published
results**: that is by design.

`generate_synthetic_clinical.py` produces:

- `clinical_synthetic.tsv` — 52 simulated patients, schema-matched
- `gnsrv_burden_per_patient_synthetic.tsv` — per-patient DDR gnsRV burden, schema-matched

## `raw/` — primary pipeline inputs (committed)

The source-of-truth tables fed into the analysis pipeline. We commit these
to make the analysis reproducible end-to-end as far as the underlying data
permits — the patient-level germline VCFs and BAMs that *produced* these
aggregates are not in this repo, and are being deposited to EGA in parallel
with manuscript submission (see `docs/DATA_AVAILABILITY.md`).

> The name "raw" is relative to this pipeline, not relative to sequencing —
> these files are the *primary inputs* the analysis scripts read; they are
> not raw FASTQ.

| File | Description | Used by |
|---|---|---|
| `gnsrv_burden_per_patient.tsv` | Per-patient DDR gnsRV burden | `code/01_cohort_prep` |
| `gnsrv_per_gene_per_patient.tsv` | Per-gene × per-patient gnsRV indicator | `code/04_firth_regression`, `code/08_figures/fig1*` |
| `synonymous_per_gene_per_patient.tsv` | Per-gene × per-patient synonymous indicator (negative control) | `code/03_ddr_enrichment/ddr_burden_bootstrap.py` |
| `oncoplot/*.tsv` | Pre-shaped oncoplot matrices (DDR + genome-wide; gnsRV, synonymous, CNV) | `code/08_figures/fig1*`, `fig2*` |
| `cbioportal/*.tsv` | cBioPortal exports for Fig 3C/D/E (TCGA-PRAD, MSKCC, SU2C-PCF) | `code/08_figures/fig3c*`, `fig3d*`, `fig3e*` |
| `functional_assays/*.tsv` | Quantified wet-lab assay readings (CCK-8, scratch wound, Boyden, olaparib dose–response) | `code/07_functional/`, `code/08_figures/fig4*` |

**Not in the repo, controlled access:**

- `vcf/*.vcf.gz` — per-patient germline variant calls
- `cnv/*.tsv` — Bionano Nexus per-patient CNV output
- Per-patient BAM/FASTQ

The analysis scripts read from `data/raw/` by default; pass `--synthetic` to
read the schema-matched simulated inputs from `data/synthetic/` instead.

## `processed/` — gitignored, generated at runtime

Per-patient derived tables written by the analysis scripts — most notably
`epc_cohort.tsv`, produced by `code/01_cohort_prep/prepare_epc_cohort.py`
from the clinical and gnsRV-burden inputs. **Never committed**; regenerated
from the inputs above whenever the pipeline is rerun.
