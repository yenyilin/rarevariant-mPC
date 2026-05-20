# Data Availability

Our goal with this release is to make the analysis reproducible from public
artifacts as far as the underlying data permits. Code, derived tables, and
the de-identified clinical cohort are committed here; the controlled-access
raw sequencing data is being deposited to public repository in parallel with 
manuscript submission. Three data access classes are listed below.

## Public (in this repository)

Aggregate / de-identified summary tables under `data/derived/` and a
synthetic cohort under `data/synthetic/`. Released under CC-BY-4.0.

| File | Description |
|---|---|
| `data/derived/ddr_panel_25genes.tsv` | The 25 prespecified DDR genes |
| `data/derived/gnsrv_56_screen.tsv` | Genome-wide screen variant list (gnomAD MAFs + FATHMM-MKL + ClinVar annotations) |
| `data/derived/australian_epc_counts.tsv` | Australian EPC carrier counts |
| `data/derived/ppcg_gene_associations.tsv` | PPCG gene-level associations |
| `data/clinical/clinical.tsv` | 52-patient de-identified clinical cohort table (real) |
| `data/synthetic/clinical_synthetic.tsv` | 52 simulated patients (schema match) |
| `data/synthetic/gnsrv_burden_per_patient_synthetic.tsv` | Per-patient synthetic gnsRV burden |

The de-identified clinical table is released publicly; see
`data/clinical/README.md` for the de-identification basis.

## Controlled access (request required)

Raw FASTQ, BAM, VCF, and per-patient genotype / CNV matrices for the
Vancouver Prostate Centre (VPC) discovery EPC (n = 52). The de-identified
clinical table for the same cohort is **not** controlled access — it is in
the public section above (`data/clinical/clinical.tsv`).

- **Deposition status:** Deposition to public repository was in progress at 
the time of manuscript submission. The accession below is provisional;  
it will be finalized, and the data released to controlled access, before publication.

## Third-party cohorts (separate access channels)

| Cohort | Source | Access |
|---|---|---|
| Australian EPC (n = 53) | Whole-genome sequencing[^auepc] | Original cohort PI |
| PPCG (n = 976) | Pan Prostate Cancer Group consortium | https://pancancer.org/ |
| TCGA-PRAD, MSKCC, SU2C-PCF, etc. | Public somatic CNV / outcome data | See `data/raw/cbioportal/README.md` for study IDs, gene set, and the exact cBioPortal query URL. |



## Software

A static, citable archive of the code in this repository is minted on
Zenodo for each release:

- **Zenodo concept DOI:** 10.5281/zenodo.<TBD>
- **Release v1.0 (resubmission snapshot):** 10.5281/zenodo.<TBD>

## Contact

- Code and reproducibility issues: GitHub Issues on this repository.
- Controlled-access data (VPC EPC): contact the corresponding author or the
  Vancouver Prostate Centre Data Access Committee directly.

## References
[^auepc]: Jaratlerdsiri W, Jiang J, Gong T, Patrick SM, Willet C, Chew T, et al.
  African-specific molecular taxonomy of prostate cancer. *Nature*
  2022;609:552–9. doi:10.1038/s41586-022-05154-6
