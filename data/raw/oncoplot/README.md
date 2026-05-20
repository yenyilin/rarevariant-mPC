# Oncoplot input matrices (committed)

Per-gene × per-patient matrices used by the DDR and genome-wide oncoplots
(Fig. 1 and Fig. 2). These are aggregate indicator and CNV matrices, and they
do not contain raw genotypes, VCFs, or BAMs. They are committed to keep
the oncoplot figures fully reproducible from this repository.

The patient-level germline VCFs and CNV calls that *produced* these
matrices remain controlled access and will be released to a public 
controlled-access repository before manuscript publication (see `docs/DATA_AVAILABILITY.md`).

| File | Used by |
|---|---|
| `met_ddr_gnsrv.tsv` (nonsynonymous, metastatic arm) | `code/08_figures/fig1_ddr_oncoplot.py` |
| `nonmet_ddr_gnsrv.tsv` (nonsynonymous, non-metastatic arm) | `code/08_figures/fig1_ddr_oncoplot.py` |
| `met_ddr_synonymous_rv.tsv` (synonymous negative control, met arm) | `code/08_figures/fig1_ddr_oncoplot.py` |
| `nonmet_ddr_synonymous_rv.tsv` (synonymous negative control, non-met arm) | `code/08_figures/fig1_ddr_oncoplot.py` |
| `cnv_ddr.tsv` | `code/08_figures/fig1_ddr_oncoplot.py` |
| `met_genomewide.tsv` | `code/08_figures/fig2_genomewide_oncoplot.py` |
| `cnv_genomewide.tsv` | `code/08_figures/fig2_genomewide_oncoplot.py` |

Each TSV is a wide gene × patient matrix; the column layout is documented
inline in `code/08_figures/_oncoplot_helpers.py` (`parse_tsv`).
