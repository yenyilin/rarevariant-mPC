# Figure-composition scripts

These scripts assemble the manuscript's main and supplementary figures
from results produced by the analysis stages in `code/`. All write their
output to `figures/` at the repository root.

The Figure 4 functional panels are **not** here — they live in
`code/07_functional/`, because those scripts perform assay analysis, not
just figure composition.

Every script below now runs end-to-end from data committed to this repo
(`data/raw/`, `data/clinical/`, `data/synthetic/`, `data/derived/`). The
controlled-access archive on EGA holds the patient-level VCFs and BAMs that
*produced* the aggregate matrices in `data/raw/` — none of the scripts here
need that level.

| Script | Manuscript | Inputs |
|---|---|---|
| `generate_table1.py` | Table 1 + love plot + swimmer plot (Methods §2.1) | Clinical TSV (CLI argument) |
| `fig1_ddr_oncoplot.py` | Fig. 1 — 5-panel DDR oncoplot (gnsRV met/nm; synonymous-variant control met/nm; HRD track; CNV side bar) | `data/raw/oncoplot/*` |
| `fig2_genomewide_oncoplot.py` | Fig. 2 — genome-wide screen oncoplot | `data/raw/oncoplot/*` |
| `_oncoplot_helpers.py` | Shared helper imported by the oncoplot scripts — not run directly | — |
| `fig3a_australian_replication.py` | Fig. 3A — Australian EPC carrier-frequency table | None (values hard-coded; author-confirmed) |
| `fig3b_ppcg_associations.py` | Fig. 3B — PPCG gene-level association table | None (values hard-coded; author-confirmed) |
| `fig3c_cnv_oncoprint.py` | Fig. 3C — CNV oncoprint, 12 candidate genes across 4 public PCa cohorts | `data/raw/cbioportal/fig3c_cnv_oncoprint.tsv` |
| `fig3d_cnv_barplot.py` | Fig. 3D — CNV alteration-frequency stacked barplot | `data/raw/cbioportal/fig3d_cnv_frequency.tsv` |
| `fig3e_overall_survival.py` | Fig. 3E — Kaplan-Meier overall survival (Altered vs Unaltered) | `data/raw/cbioportal/fig3e_overall_survival.tsv` |
| `figS_per_patient_pvalues.py` | Supplementary Methods Fig.4 — per-patient empirical p-value forest plot (Supp. Methods §4.2.4) | `data/raw/gnsrv_per_gene_per_patient.tsv` (+ DDR panel & FLAGS list from `data/derived/`) |
| `figS_study_design.py` | Fig. S1 — two-stage discovery-replication schematic | None (pure diagram) |

## Usage

```bash
# No external data — pure diagrams or hard-coded values
python code/08_figures/figS_study_design.py
python code/08_figures/fig3a_australian_replication.py
python code/08_figures/fig3b_ppcg_associations.py

# Per-patient aggregate matrices (committed under data/raw/)
python code/08_figures/fig1_ddr_oncoplot.py
python code/08_figures/fig2_genomewide_oncoplot.py
python code/08_figures/figS_per_patient_pvalues.py

# Public cBioPortal exports (committed under data/raw/cbioportal/)
python code/08_figures/fig3c_cnv_oncoprint.py
python code/08_figures/fig3d_cnv_barplot.py
python code/08_figures/fig3e_overall_survival.py

# Table 1 — pass a clinical TSV (the real one in data/clinical/, or synthetic)
python code/08_figures/generate_table1.py data/clinical/clinical.tsv -o figures/
```

## About these scripts

Paths resolve through the repo's `data/` and `figures/` directories
(`REPO = Path(__file__).resolve().parents[2]`); each script reads committed
or controlled-access data and writes figures to `figures/`. The scientific
code is unchanged from the manuscript working tree.
