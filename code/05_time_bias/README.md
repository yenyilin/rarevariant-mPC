# `code/05_time_bias/` — Time-dependent selection bias analyses

These scripts test for **time-dependent selection bias** — the concern
that patients who metastasize earlier inherently differ from those who
remain non-metastatic for long periods, which could confound the
appearance of metastatic-exclusive variants.

Two complementary tests rule out the most concerning time-bias scenarios:

| Script | Test | Manuscript result |
|---|---|---|
| `ddr_followup_plot.py` | Mann-Whitney U: non-metastatic follow-up duration vs. metastatic time-to-metastasis | **p < 0.001** — non-metastatic patients have substantially longer follow-up than the median time-to-metastasis in the metastatic arm, arguing against censoring-driven bias |
| `ddr_burden_vs_time.py` | Spearman correlation: DDR gnsRV burden vs. time-to-metastasis (within the metastatic arm) | **ρ = −0.152, p = 0.458, n = 26** — no evidence that high-burden patients metastasize earlier; null result rules out timing-bias explanation |

## Inputs

Both scripts read:
- `data/processed/epc_cohort.tsv` (52 patients × clinical fields)

Required columns:
- `metastatic` (0/1)
- `time_to_metastasis` (months; NA for non-metastatic)
- `followup_duration_months` (months; for non-metastatic patients)
- `ddr_burden` (count of DDR gnsRVs per patient)

The `epc_cohort.tsv` file is produced by `code/01_cohort_prep/prepare_epc_cohort.py` from controlled-access clinical data. For end-to-end testing with synthetic data, see `data/synthetic/` and `tests/test_pipeline_smoke.py`.

## Outputs

| Script | Output | Manuscript figure |
|---|---|---|
| `ddr_followup_plot.py` | `figures/followup_comparison.{pdf,png}` | Supplementary Methods Fig. 5 |
| `ddr_burden_vs_time.py` | `figures/ddr_burden_vs_time.{pdf,png}` | Supplementary Methods Fig. 6 |

## Usage

```bash
python3 code/05_time_bias/ddr_followup_plot.py
python3 code/05_time_bias/ddr_burden_vs_time.py
```

Both scripts print test statistics to stdout and save figures to
`figures/`.

## Note on file history

These scripts were previously located in `code/08_figures/` under names
`figS_MX_ddr_followup_plot.py` and `figS_M6_ddr_burden_vs_time.py`.
They were moved here to match the documented analytic stage 
 (`code/05_time_bias/`), and it groups the two time-bias analyses by 
 analytical purpose (time-dependent selection bias) rather than by 
 output type (figure generation).
