# Variant calling pipeline

The pipeline that produces germline VCFs from raw FASTQ is a thin wrapper
around standard tools and does not redistribute the underlying sequencing
data. It is documented here so reviewers can verify the parameters reported
in the manuscript Methods §2.2. The use of standard tools throughout this
stage is deliberate: see the top-level README ("Scope") for why the study
relies on established software rather than bespoke methods.

## Inputs

- Paired-end Illumina NextSeq500 reads (FASTQ.gz) from FFPE tumor and
  matched distant benign prostate tissue. Deposition to a controlled-access
  archive was in progress at the time of manuscript submission; see
  `docs/DATA_AVAILABILITY.md` for the accession number and current status.

## Stages and parameters

| Stage | Tool | Version | Key parameters |
|---|---|---|---|
| Alignment | BWA-MEM | 0.7.17 | reference: GRCh38, Ensembl Release 92 |
| Duplicate marking | Picard | 2.18.11 | `MarkDuplicates` |
| Filter mismatched reads | custom | — | drop a read pair if either 150 bp mate has ≥5 mismatched bases (≈3%) |
| Germline calling | Strelka2 | 2.9.10 | germline mode |
| Post-call filter | bcftools | 1.18 | DP ≥ 15; AAF ≥ 0.20 in tumor and matched benign |
| Annotation | ANNOVAR | 2017-07-17 | GENCODE v27, gnomAD v2.1.1 |
| Function prediction | FATHMM-MKL | dbNSFP 3.0a | deleterious threshold ≥ 0.5 |
| CNV calling | Bionano Nexus | 11.0 | min depth 20, MQ 30, BQ 20, sig 1e-6, min calls 3, max distance 1Mb |
| HRD score | scarHRD | 0.1.1 |  |

## Reference filters

- **Rarity:** MAF ≤ 2% in gnomAD v2.1.1 (both WGS and WES sub-cohorts).
  gnomAD-absent variants are assigned MAF = 0 and retained.
- **Consequence:** missense, nonsense, splice-site (±2 bp from exon boundary).
- **AAF threshold rationale:** the ≥ 0.20 threshold also excludes FFPE
  deamination artifacts, which manifest at AAF < 0.10. FF/FFPE concordance
  was 97.72% in seven matched-pair patients (Methods §2.2;
  Supplementary Methods §3.4).

## What is **not** here

- The raw FASTQ files and the intermediate BAM/VCF files — controlled-access
  data; see `docs/DATA_AVAILABILITY.md` (Controlled access) for deposition
  status and the data-access procedure.
- Wrapper shell scripts for the sequencing core's specific cluster. These
  are environment-specific and are not part of the analytical method.

If you need to rerun the variant-calling pipeline on your own data, the
parameters above are sufficient to reproduce it with any standard germline
WES workflow.
