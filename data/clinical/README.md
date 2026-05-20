# `data/clinical/` — de-identified clinical cohort

`clinical.tsv` is the real clinical table for the 52 EPC patients used
throughout the manuscript (Table 1; the real-mode input to
`code/01_cohort_prep/prepare_epc_cohort.py`).

## Why this is committed

The `.gitignore` blocks every `clinical*.tsv` to keep PHI out of the repo.
This file is the one allowed exception (`!data/clinical/clinical.tsv`)
because it's been de-identified.

- All direct identifiers removed (no names, medical record numbers, dates of
  birth, or institution-specific accession IDs).
- Quasi-identifiers reviewed before release: `Diagnosis year` was dropped
  (it touched no analytical result, so retaining it would only have added
  re-identification risk). `Metastasis location` is recorded at the
  anatomical-site level and was cleared as non-identifying for this cohort.
- Patient IDs are study-internal sequential labels, not institutional IDs.

## Schema

One row per patient, tab-separated, sentence-case headers. The columns
consumed by the analysis pipeline (`patient_id`, `metastatic`, `ancestry`,
`age_at_diagnosis`, `biopsy_gleason`, `surgical_margins`,
`time_to_metastasis`, `followup_duration_months`) are mapped from these
headers by the `COLUMN_RENAME` table in `prepare_epc_cohort.py`. Columns not
in that table are carried for reference and ignored by the pipeline.

## Relationship to the other tiers

- `data/synthetic/clinical_synthetic.tsv` mirrors this file's schema with
  simulated values, for users who want to dry-run the pipeline.
- This file is the genuine cohort; results computed from it reproduce the
  published numbers.
