# `code/07_functional/` — Functional characterization (Figure 4)

Analysis-and-plotting scripts for the in vitro functional assays. Each
script reads its raw assay data, computes the reported statistics, and
renders one Figure 4 panel. Wet-lab assay protocols (CRISPR/Cas9 prime
editing, scratch-wound migration, Boyden-chamber invasion, CCK-8
proliferation, Olaparib CCK-8 viability) are documented in the
manuscript's Supplementary Methods.

| Script | Panel | Assay |
|---|---|---|
| `fig4a_proliferation.py` | Fig. 4A | CCK-8 proliferation |
| `fig4b_scratch_wound.py` | Fig. 4B | Scratch-wound migration |
| `fig4c_boyden_chamber.py` | Fig. 4C | Boyden-chamber invasion + migration |
| `fig4d_olaparib.py` | Fig. 4D | Olaparib CCK-8 viability |

## Inputs — raw assay data

Located in `data/raw/functional_assays/`:

| Assay | File |
|---|---|
| CCK-8 proliferation | `fig4a_cck8_proliferation.tsv` |
| Scratch-wound migration | `fig4b_scratch_wound.tsv` |
| Boyden-chamber invasion + migration | `fig4c_boyden_chamber.tsv` |
| Olaparib dose-response | `fig4d_olaparib_dose_response.tsv` |

## Outputs

Figure panels are written to `figures/` at the repository root.

The scratch-wound and Boyden scripts also write step-by-step
quantification — technical-replicate means, log-transform, paired log
differences — and a statistics summary to `data/processed/functional/`:

- `scratch_step1_br_means.tsv`, `scratch_step2_log_br_means.tsv`,
  `scratch_step3_log_differences.tsv`, `scratch_stats_summary.txt`
- `boyden_step1_br_means.tsv`, `boyden_step2_log_br_means.tsv`,
  `boyden_step3_log_differences.tsv`, `boyden_stats_summary.txt`

The proliferation (4A) and Olaparib (4D) scripts compute their statistics
in-script and do not write intermediate tables.

## Usage

```bash
python code/07_functional/fig4a_proliferation.py
python code/07_functional/fig4b_scratch_wound.py
python code/07_functional/fig4c_boyden_chamber.py
python code/07_functional/fig4d_olaparib.py
```
