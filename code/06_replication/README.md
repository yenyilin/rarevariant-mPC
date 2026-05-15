# `code/06_replication/` — Independent replication analyses

Two independent cohorts test whether the metastatic-exclusive candidate
gnsRVs identified in the VPC EPC discovery cohort replicate beyond the
discovery data.

## Australian EPC — binomial carrier-frequency test

`australian_binomial.py` — runnable from this repository.

A one-sided binomial test of carrier frequency in the Australian EPC
against a gnomAD-derived population expectation. For each subgroup of
size *n* with *k* observed carriers, the test asks whether *k* exceeds
the count expected under `p_carrier`, the probability that a random
individual carries at least one candidate variant. `p_carrier` is
computed on the fly from the 56-variant screen table (see the script
header for the filter rule and gnomAD column used).

```bash
python code/06_replication/australian_binomial.py
```

- **Inputs:** `data/derived/gnsrv_56_screen.tsv`, `data/derived/australian_epc_counts.tsv`
- **Output:** `data/processed/australian_binomial_results.csv`

## PPCG — gene-level association testing

**The PPCG analyses were performed by the Pan Prostate Cancer Group
consortium on consortium-internal data. The authors received summary
results only; no PPCG individual-level genotype data was held by the
authors, and these analyses are not runnable from this repository.**

Gene-level association testing in the PPCG cohort (*n* = 976) comprised:

- **Wald test** — per-variant association (logistic or linear, by endpoint).
- **Burden test** — all gnsRVs within a gene collapsed to a single
  per-gene carrier indicator.
- **SKAT** — a single pre-planned aggregate test across all 476 candidate
  gnsRVs (one test; no multiple-testing correction required).

Covariates were age at radical prostatectomy and genetic ancestry. The
reported aggregate SKAT result is *p* = 0.03 (Methods §2.3; Results;
Fig. 3B; Supplementary Methods §4.5).

The summary results are committed as
`data/derived/ppcg_gene_associations.tsv` (per-gene endpoint, Wald,
burden, and SKAT statistics). PPCG cohort access is governed by the
consortium — see `docs/DATA_AVAILABILITY.md` (Tier 3) and
https://panprostate.org/.
