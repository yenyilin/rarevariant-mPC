# cBioPortal-downloaded data (Figure 3 C, D, E)

External public-data downloads from [cBioPortal](https://www.cbioportal.org/)
used to generate panels C, D, and E of Figure 3.

## Files

| File | Figure | Description | Used by script |
|------|--------|-------------|----------------|
| `fig3c_cnv_oncoprint.tsv` | Fig 3C | Patient-level CNV alterations (gene × sample matrix of amplifications / deletions) for the candidate gnsRV gene set. | `code/08_figures/fig3c_cnv_oncoprint.py` |
| `fig3d_cnv_frequency.tsv` | Fig 3D | Per-cancer-type CNV alteration frequency summary (study-level aggregate counts). | `code/08_figures/fig3d_cnv_barplot.py` |
| `fig3e_overall_survival.tsv` | Fig 3E | Clinical attributes including overall survival status and time-to-event (months), used for Kaplan–Meier analysis. | `code/08_figures/fig3e_overall_survival.py` |

## Provenance

Files were downloaded from cBioPortal on **2026-04-26**.

### Cohorts queried (four public PCa studies, 1,282 samples total)

| Cohort (display name) | cBioPortal study ID | Reference |
|------------------------|---------------------|-----------|
| Metastatic Prostate Cancer (SU2C/PCF Dream Team, Cell 2015) | `prad_su2c_2015` | Robinson et al., *Cell* 2015 |
| Metastatic Prostate Adenocarcinoma (SU2C/PCF Dream Team, PNAS 2019) | `prad_su2c_2019` | Abida et al., *PNAS* 2019 |
| Prostate Adenocarcinoma (MSK, Cancer Cell 2010) | `prad_mskcc` | Taylor et al., *Cancer Cell* 2010 |
| Prostate Adenocarcinoma (TCGA, PanCancer Atlas) | `prad_tcga_pan_can_atlas_2018` | TCGA Research Network, *Cell* 2015 / PanCancer Atlas 2018 |

### Gene set queried (12 genes from the discovery EPC genome-wide screen; Fig. 2)

```
CRISPLD1, ERBIN, KDM6B, LAMA5, MASP1, MCM2,
PKD1L2,  POLB,  RP1L1, SPRY3, SUSD1, ZC3H7A
```

### Canonical cBioPortal query URL

The exact URL that produced these downloads (encodes all four study IDs,
the twelve-gene panel, and the alteration-type filters used):

```
https://www.cbioportal.org/results/oncoprint?cancer_study_list=prad_su2c_2015%2Cprad_su2c_2019%2Cprad_mskcc%2Cprad_tcga_pan_can_atlas_2018&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&profileFilter=mutations%2Cstructural_variants%2Cgistic%2Ccna&case_set_id=all&gene_list=CRISPLD1%250AERBIN%250AKDM6B%250ALAMA5%250AMASP1%250AMCM2%250APKD1L2%250APOLB%250ARP1L1%250ASPRY3%250ASUSD1%250AZC3H7A&geneset_list=%20&tab_index=tab_visualize&Action=Submit
```

Query parameters (URL-decoded):
- `cancer_study_list` — the four study IDs above, comma-separated.
- `gene_list` — the twelve genes, newline-separated.
- `profileFilter=mutations,structural_variants,gistic,cna` — selects four
  DNA-level molecular profile types: somatic mutations, structural variants,
  GISTIC2 discretized CNAs, and raw / continuous CNA scores. The CNV
  analyses in Fig 3 C/D/E use only the GISTIC / CNA portion of this filter;
  the mutations and structural-variants portions are present in the URL but
  do not contribute to the panels shown.
- `case_set_id=all` — uses cBioPortal's "All samples" case list for each
  selected study (no case-level filtering beyond what each study itself
  defines as "all").
- `Z_SCORE_THRESHOLD=2.0`, `RPPA_SCORE_THRESHOLD=2.0` — mRNA z-score and
  RPPA protein-score thresholds. These thresholds govern expression-based 
  alteration calls and are inert when `profileFilter` does not include 
  mRNA / RPPA profiles (as is the case here); they are preserved in the 
  URL for byte-identical reproducibility.

The original cBioPortal export filenames (`PATIENT_DATA_CNV_oncoprint.tsv`,
`cancer_types_summary.txt`, `Overall.txt`) were renamed to follow the
repository's `fig{panel}_{description}.tsv` convention for reproducibility
and grep-friendliness.

## To re-download

For an independent re-download from cBioPortal:

1. **Open the canonical query URL** above in a web browser. This loads the
   pre-configured OncoPrint view with all four studies and twelve genes.
2. **OncoPrint export** (Fig 3C):
   - The OncoPrint tab opens by default → click the download icon
     (top-right of the panel) → select the patient-level CNV matrix export.
   - Save as `fig3c_cnv_oncoprint.tsv` in this directory.
3. **Cancer Types Summary export** (Fig 3D):
   - Switch to the *Cancer Types Summary* tab → download the per-cancer-type
     CNV frequency table.
   - Save as `fig3d_cnv_frequency.tsv` in this directory.
4. **Clinical / Overall Survival export** (Fig 3E):
   - Switch to the *Comparison/Survival* tab → download the survival table
     containing patient ID, study ID, time-to-event (months), event status,
     and KM survival probability.
   - Save as `fig3e_overall_survival.tsv` in this directory.

Numerical results in the published Fig 3 reflect cBioPortal data as of the
download date (**2026-04-26**). cBioPortal updates its underlying data
periodically, so values obtained after re-download may differ from the
published figures if any of the four underlying studies have since been
updated.

## Format notes

- All three files are tab-separated. Original `.txt` extensions from
  cBioPortal were renamed to `.tsv` for consistency with the rest of the
  repository.
- File headers vary by export type; see the docstring of each consuming
  script (`fig3{c,d,e}_*.py`) for the expected column structure.
- cBioPortal occasionally emits BOM characters or non-UTF-8 encoding;
  if a script fails on read, open the file in a text editor and confirm
  UTF-8 encoding without BOM.

## License

Data is publicly available under cBioPortal's terms of use (typically
CC-BY 4.0 for derived analyses; the underlying primary studies have their
own data-use agreements — verify per study before redistribution).
This README does not relicense the data; the files themselves remain
subject to their original publication's data-use terms.
