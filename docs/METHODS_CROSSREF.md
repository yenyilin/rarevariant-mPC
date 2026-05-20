# Manuscript claim → script crossreference

Every numerical result in the manuscript is produced by a specific
script in this repository. This table maps the two together so reviewers
can verify any reported value.

## Discovery cohort (EPC, n = 52)

| Manuscript location | Claim | Script | Output |
|---|---|---|---|
| Results, ¶1 | 10,854 metastatic-exclusive gnsRVs; 6,080 fathmm-MKL deleterious | upstream variant filtering — see `code/02_variant_calling/README.md` | input matrix |
| Results, ¶2; Fig. 1A | 71 unique gnsRVs (81 carrier-variant observations) exclusive to the metastatic arm vs. 21 exclusive to the non-metastatic arm | upstream variant filtering — see `code/02_variant_calling/README.md` | input to bootstrap |
| Results, ¶2 | DDR per-patient burden: met > nm; **p = 4.57 × 10⁻⁶** (Mann-Whitney two-sided) | `code/03_ddr_enrichment/ddr_arm_compare.py` | stdout |
| Results, ¶2; Fig. 1D, E | Synonymous control: 15 vs. 17, p = n.s. | upstream variant filtering — see `code/02_variant_calling/README.md` | input matrix |
| Results, ¶2 | Bootstrap vs. random 25-gene panels: p < 1 × 10⁻⁴ (within-arm competitive null) | `code/03_ddr_enrichment/ddr_burden_bootstrap.py --arm met` | stdout |
| Results, ¶2 | Per-patient empirical p + Fisher's combined: p = 3.9 × 10⁻¹⁹ (met) vs p = 0.89 (nm) | `code/03_ddr_enrichment/ddr_per_patient_test.py` | stdout |
| Results, ¶2; Fig. 1B, C | Mean HRD: 21.58 (met) vs. 14.62 (nm); p = 0.04 | `code/08_figures/fig1_ddr_oncoplot.py` (renders the HRD column of the oncoplot) | manuscript table |
| Methods §2.3.2; Results, ¶3; Abstract | Firth primary: OR = 26.80, 95% CI 3.07–3528.42, LRT p < 0.001 | `code/04_firth_regression/ddr_logistic_regression.py` | `data/processed/firth_results.csv` |
| Methods §2.3.2; Results, ¶3 | Firth extended (+Gleason+margins): OR = 24.93, 95% CI 2.91–3261.56, p < 0.001 | `code/04_firth_regression/ddr_logistic_regression.py --extended` | same |
| Results, ¶3; Methods §2.3.2; Supplementary Methods Fig. 5 | Mann-Whitney follow-up vs. time-to-mets: p < 0.001 | `code/05_time_bias/ddr_followup_plot.py` | `figures/followup_comparison.pdf` |
| Results, ¶3; Methods §2.3.2; Supplementary Methods Fig. 6 | Spearman burden vs. time-to-mets: ρ = −0.152, p = 0.458, n = 26 | `code/05_time_bias/ddr_burden_vs_time.py` | `figures/ddr_burden_vs_time.pdf` |
| Results, ¶4; Fig. 2 | Genome-wide screen: 56 gnsRVs in 53 genes, ≥ 3 met carriers, 0 nm | upstream variant filtering | `data/derived/gnsrv_56_screen.tsv` |
| Results, ¶4 | Total gnsRV burden across 53 screen genes: p = 2.98 × 10⁻⁸ | extension of `ddr_burden_bootstrap.py` (panel = 53 screen genes) | stdout |

## Replication

| Manuscript location | Claim | Script | Output |
|---|---|---|---|
| Results, ¶5; Fig. 3A; Table 2 | Australian EPC bone/visceral 9/17, p = 0.0506 | `code/06_replication/australian_binomial.py` | stdout |
| Results, ¶5; Fig. 3A | Australian EPC met+BCR 14/29, p = 0.0400 | same | stdout |
| Results, ¶5; Fig. 3A | Australian EPC met+BCR+nodes 16/33, p = 0.0278 | same | stdout |
| Results, ¶5; Fig. 3A | Australian EPC non-met 5/14, p = 0.4534 | same | stdout |
| Methods §2.4 | p_carrier = 0.311 (computed on-the-fly) | same (`compute_p_carrier`) | stdout |
| Fig. 3A (panel) | Australian EPC carrier-frequency table — rendered figure panel | `code/08_figures/fig3a_australian_replication.py` | `figures/fig3a_australian_replication.pdf` |
| Results, ¶6; Fig. 3B; Table 3 | PPCG aggregate SKAT: p = 0.03 (n = 976; 476 candidate gnsRVs) | PPCG consortium — analysis performed by the consortium; summary results provided | `data/derived/ppcg_gene_associations.tsv` |
| Results, ¶6; Fig. 3B | Per-gene burden hits: FOCAD, CHEK2, ACSF3, SIVA1, SPATA9, ZSWIM4, DNAJC10, HPGDS | same | `data/derived/ppcg_gene_associations.tsv` |
| Discussion | CHEK2 I157T age-at-diagnosis: p = 1.0 × 10⁻⁶ | same | `data/derived/ppcg_gene_associations.tsv` |
| Fig. 3B (panel) | PPCG gene-level association table rendered figure panel | `code/08_figures/fig3b_ppcg_associations.py` | `figures/fig3b_ppcg_associations.pdf` |

## Functional

| Manuscript location | Claim | Script | Output |
|---|---|---|---|
| Results, ¶7; Fig. 4A | KDM6B K973Q proliferation: parental > C6 > C3 | `code/07_functional/fig4a_proliferation.py` | `figures/fig4a_proliferation.pdf` |
| Results, ¶7; Fig. 4B | Scratch wound: enhanced migration in both knock-ins, C3 > C6 | `code/07_functional/fig4b_scratch_wound.py` | `figures/fig4b_scratch_wound.pdf` |
| Results, ¶7; Fig. 4C | Boyden chamber: ~4-fold migration/invasion, C3 > C6 | `code/07_functional/fig4c_boyden_chamber.py` | `figures/fig4c_boyden_chamber.pdf` |
| Results, ¶7; Fig. 4D | BRCA2 I1962T Olaparib: reduced viability vs. parental, less than BRCA1 Q1200X | `code/07_functional/fig4d_olaparib.py` | `figures/fig4d_olaparib.pdf` |

## Cohort and supplementary

| Manuscript location | Claim | Script | Output |
|---|---|---|---|
| Methods §2.1; Table 1 | SMD ≤ 0.1 for most variables; Gleason and margins SMD ~0.20 | `code/08_figures/generate_table1.py` (`plot_love`) | `figures/table1_love.pdf` |
| Methods §2.1; Supplementary Methods Fig. 5 | Swimmer plot of follow-up | `code/08_figures/generate_table1.py` (`plot_swimmer`) | `figures/table1_swimmer.pdf` |
| Supplementary Methods Fig. 1 | Two-stage discovery-replication study-design schematic | `code/08_figures/figS_study_design.py` | `figures/figS_study_design.pdf` |
| Supplementary Methods Fig. 4; Supplementary Methods §4.2.4 | Per-patient empirical p-value forest plot | `code/08_figures/figS_per_patient_pvalues.py` | `figures/figS_per_patient_pvalues.pdf` |
| Supplementary Methods §1.2 | Continuous vs. categorical variable comparisons (age p = 0.359; biopsy GS p = 0.367; etc.) | `code/01_cohort_prep/prepare_epc_cohort.py` (descriptive stage) | stdout |


## Notes for reviewers

- Path conventions: `data/clinical/` holds the de-identified clinical table (committed); `data/raw/` holds primary pipeline inputs as aggregate per-patient matrices (committed; the raw VCFs/BAMs that produced them are controlled-access and will be released to a public controlled-access repository before manuscript publication — see `docs/DATA_AVAILABILITY.md`); `data/derived/` holds aggregate summary tables; `data/synthetic/` mirrors the schema for end-to-end smoke testing; `data/processed/` is generated at runtime and gitignored.
- Reseed: every stochastic script accepts `--seed`; default is 42. Use the same seed to reproduce exactly.
