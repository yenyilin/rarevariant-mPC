"""
Build data/processed/epc_cohort.tsv (the analysis-ready cohort table) from a
clinical TSV and a per-patient gnsRV burden TSV.

Default inputs:
  1. Clinical TSV
       real:      data/clinical/clinical.tsv
       synthetic: data/synthetic/clinical_synthetic.tsv
  2. Per-patient gnsRV burden TSV
       real:      data/raw/gnsrv_burden_per_patient.tsv
       synthetic: data/synthetic/gnsrv_burden_per_patient_synthetic.tsv

────────────────────────────────────────────────────────────────────────────
INPUT CONTRACT
────────────────────────────────────────────────────────────────────────────
Clinical TSV — sentence-case headers (the "new" schema; the previous ALL-CAPS
schema with explicit DATE columns is no longer used — date arithmetic was
dropped). Time fields must be supplied pre-derived in "(months)".

  Required — the script aborts with an error if any are missing:
      Patient ID
      Metastasis (1=yes, 0=no)
      Ancestry (inferred)
      Age at diagnosis (years)

  Optional — used by downstream scripts; if a column is absent the matching
  output column is emitted empty (NA) instead of raising an error:
      Biopsy total Gleason score
      Surgical margins                ("POSITIVE" / "NEGATIVE")
      Time to metastasis (months)     (blank for non-metastatic patients)
      Follow-up duration (months)

  Any other sentence-case column is recognised only if it appears in the
  COLUMN_RENAME map below; unmapped columns are ignored.

Per-patient gnsRV burden TSV:
      column 1          per-patient identifier (any header; matched by
                        position and coerced to align with Patient ID)
      DDR_total_gnsrv          non-synonymous DDR rare-variant burden
      DDR_total_synonymous_rv  synonymous DDR rare-variant burden
  If this file is absent, the ddr_burden columns are left empty with a warning.

────────────────────────────────────────────────────────────────────────────
OUTPUT
────────────────────────────────────────────────────────────────────────────
Columns consumed by the downstream analysis scripts:
    patient_id, metastatic, ancestry, age_at_diagnosis,
    biopsy_gleason, surgical_margins_binary,
    ddr_burden, ddr_synonymous_burden,
    time_to_metastasis, followup_duration_months

Usage:
    python code/01_cohort_prep/prepare_epc_cohort.py            # real defaults
    python code/01_cohort_prep/prepare_epc_cohort.py --synthetic
    python code/01_cohort_prep/prepare_epc_cohort.py \\
        --clinical /path/to/clinical.tsv \\
        --extra    /path/to/gnsrv_burden_per_patient.tsv
"""

import argparse
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[2]


# Map new sentence-case headers (left) to the internal column names used
# by every downstream analysis script (right). Add entries here if the
# upstream TSV grows new columns; do not rename the values without also
# updating the consumer scripts (Firth regression, time-bias plots, etc.).
COLUMN_RENAME = {
    "Patient ID":                      "patient_id",
    "Metastasis (1=yes, 0=no)":        "metastatic",
    "Ancestry (inferred)":             "ancestry",
    "Age at diagnosis (years)":        "age_at_diagnosis",
    "PSA at diagnosis (ng/mL)":        "psa_at_diagnosis",
    "Biopsy total Gleason score":      "biopsy_gleason",
    "Biopsy positive cores":           "biopsy_positive_cores",
    "Biopsy total cores":              "biopsy_total_cores",
    "Clinical stage":                  "clinical_stage",
    "Pathological stage":              "pathological_stage",
    "Surgical margins":                "surgical_margins",
    "Lymph node status":               "lymph_node_status",
    "Status at follow-up":             "status_at_followup",
    "Follow-up duration (months)":     "followup_duration_months",
    "Time to PSA recurrence (months)": "time_to_psa_recurrence",
    "Time to metastasis (months)":     "time_to_metastasis",
    "Time to castration resistance (months)": "time_to_crpc",
    "Metastasis location":             "metastasis_location",
    "Metastasis treatment":            "metastasis_treatment",
    "Castration resistance treatment": "castration_resistance_treatment",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument(
        "--synthetic",
        action="store_true",
        help="Read from data/synthetic/ instead of data/raw/.",
    )
    p.add_argument(
        "--clinical", type=Path, default=None,
        help="Override path to the clinical TSV.",
    )
    p.add_argument(
        "--extra", type=Path, default=None,
        help="Override path to the per-patient gnsRV burden TSV "
             "(if absent, ddr_burden columns are left empty with a warning).",
    )
    p.add_argument(
        "--out", type=Path,
        default=REPO / "data" / "processed" / "epc_cohort.tsv",
        help="Output CSV path.",
    )
    return p.parse_args()


def standardise_ancestry(val) -> str:
    if pd.isna(val):
        return "Other"
    v = str(val).upper()
    if "EUROPEAN" in v:
        return "European"
    if "EAST ASIAN" in v:
        return "East Asian"
    return "Other"


def main() -> None:
    args = parse_args()

    if args.clinical is not None:
        clinical_tsv = args.clinical
    elif args.synthetic:
        clinical_tsv = REPO / "data" / "synthetic" / "clinical_synthetic.tsv"
    else:
        clinical_tsv = REPO / "data" / "clinical" / "clinical.tsv"

    if args.extra is not None:
        extra_tsv = args.extra
    elif args.synthetic:
        extra_tsv = REPO / "data" / "synthetic" / "gnsrv_burden_per_patient_synthetic.tsv"
    else:
        extra_tsv = REPO / "data" / "raw" / "gnsrv_burden_per_patient.tsv"

    out_csv = args.out
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    # ── Load clinical and rename columns ─────────────────────────────────────
    clin = pd.read_csv(clinical_tsv, sep="\t")
    clin.columns = [c.strip().strip('"') for c in clin.columns]
    clin = clin.rename(columns=COLUMN_RENAME)

    required = {"patient_id", "metastatic", "ancestry", "age_at_diagnosis"}
    missing = required - set(clin.columns)
    if missing:
        raise SystemExit(
            f"Missing required columns after rename: {sorted(missing)}.\n"
            f"Found columns: {sorted(clin.columns)}\n"
            f"Update COLUMN_RENAME in this script to map your headers."
        )

    clin["patient_id"] = pd.to_numeric(clin["patient_id"], errors="coerce").astype("Int64")
    clin["metastatic"] = pd.to_numeric(clin["metastatic"], errors="coerce").astype(int)
    print(
        f"Clinical: {len(clin)} patients "
        f"({int((clin['metastatic']==1).sum())} metastatic, "
        f"{int((clin['metastatic']==0).sum())} non-metastatic) "
        f"from {clinical_tsv}"
    )

    # Pre-derived months columns are read directly; coerce to float and
    # leave NaN where the source is blank (e.g. time_to_metastasis is
    # absent for non-metastatic patients).
    for col in ("time_to_metastasis", "followup_duration_months",
                "age_at_diagnosis", "biopsy_gleason"):
        if col in clin.columns:
            clin[col] = pd.to_numeric(clin[col], errors="coerce")

    # ── Merge per-patient DDR burden (optional) ──────────────────────────────
    if extra_tsv.exists():
        extra = pd.read_csv(extra_tsv, sep="\t")
        extra = extra.rename(
            columns={
                extra.columns[0]: "patient_id",
                "DDR_total_gnsrv": "ddr_burden",
                "DDR_total_synonymous_rv": "ddr_synonymous_burden",
            }
        )
        extra["patient_id"] = pd.to_numeric(extra["patient_id"], errors="coerce").astype("Int64")
        clin = clin.merge(
            extra[["patient_id", "ddr_burden", "ddr_synonymous_burden"]],
            on="patient_id", how="left",
        )
        n_missing = int(clin["ddr_burden"].isna().sum())
        print(f"DDR burden merged from {extra_tsv} "
              f"({len(clin) - n_missing}/{len(clin)} patients)")
    else:
        print(f"WARNING: {extra_tsv} not found — ddr_burden columns left empty.")
        clin["ddr_burden"] = pd.NA
        clin["ddr_synonymous_burden"] = pd.NA

    # ── Derived columns ──────────────────────────────────────────────────────
    clin["ancestry"] = clin["ancestry"].apply(standardise_ancestry)
    if "surgical_margins" in clin.columns:
        clin["surgical_margins_binary"] = clin["surgical_margins"].map(
            {"POSITIVE": 1, "NEGATIVE": 0}
        )
    else:
        clin["surgical_margins_binary"] = pd.NA

    # ── Select and save ──────────────────────────────────────────────────────
    keep = [
        "patient_id",
        "metastatic",
        "ancestry",
        "age_at_diagnosis",
        "biopsy_gleason",
        "surgical_margins_binary",
        "ddr_burden",
        "ddr_synonymous_burden",
        "time_to_metastasis",
        "followup_duration_months",
    ]
    for col in keep:
        if col not in clin.columns:
            clin[col] = pd.NA

    out = clin[keep].copy()
    out.to_csv(out_csv, sep="\t", index=False)
    print(f"\nSaved: {out_csv}")
    print(out.describe(include="all").round(2))


if __name__ == "__main__":
    main()
