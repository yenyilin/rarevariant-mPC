"""
Generate a synthetic cohort that matches the schema of data/clinical/clinical.tsv
and data/raw/gnsrv_burden_per_patient.tsv. Lets reviewers run every analysis
script end-to-end without controlled-access credentials.

Numerical results from this synthetic cohort will NOT match the published
results — the simulated effect size is intentionally weaker than the real
OR = 26.80, so that the test suite can verify the pipeline runs without
implying that the synthetic data reproduces the biology.

Schema: 27-column sentence-case headers with pre-derived "(months)" time
fields (no date arithmetic) — matches data/clinical/clinical.tsv as of the
schema migration.

Usage:
    python data/synthetic/generate_synthetic_clinical.py [--seed 42]

Outputs:
    data/synthetic/clinical_synthetic.tsv
    data/synthetic/gnsrv_burden_per_patient_synthetic.tsv
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd

OUT_DIR = Path(__file__).resolve().parent
N_TOTAL = 52
N_MET = 26


def _fmt_months(arr: np.ndarray) -> list[str]:
    """Format a months array, emitting empty string for NaN (matches real file)."""
    return [f"{x:.1f}" if not np.isnan(x) else "" for x in arr]


def simulate(seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    pid = np.arange(1, N_TOTAL + 1)
    metastatic = np.array([1] * N_MET + [0] * (N_TOTAL - N_MET))

    # ── Demographics ─────────────────────────────────────────────────────────
    ancestry = rng.choice(
        ["EUROPEAN", "EAST ASIAN", "SOUTH ASIAN", "HISPANIC", "WESTERN ASIAN"],
        size=N_TOTAL,
        p=[0.62, 0.25, 0.04, 0.04, 0.05],
    )
    age = np.round(rng.normal(loc=65, scale=6, size=N_TOTAL)).clip(45, 85).astype(int)
    psa = np.round(rng.lognormal(mean=2.0, sigma=0.6, size=N_TOTAL), 2)

    # ── Biopsy Gleason ───────────────────────────────────────────────────────
    biopsy_total = rng.choice([7, 8, 9], size=N_TOTAL, p=[0.15, 0.55, 0.30])
    biopsy_g1 = np.where(biopsy_total >= 8, 4, rng.choice([3, 4], size=N_TOTAL))
    biopsy_g2 = biopsy_total - biopsy_g1
    biopsy_pos = rng.integers(1, 8, size=N_TOTAL)
    biopsy_tot = biopsy_pos + rng.integers(0, 7, size=N_TOTAL)  # always ≥ positive

    # ── Clinical / treatment ─────────────────────────────────────────────────
    clinical_stage = rng.choice(
        ["T1c", "T2a", "T2b", "T2c", "T3a", "T3b"],
        size=N_TOTAL,
        p=[0.15, 0.20, 0.15, 0.20, 0.20, 0.10],
    )
    main_treatment = rng.choice(
        ["RADICAL PROSTATECTOMY", "EBRT+ADT", "EBRT"],
        size=N_TOTAL,
        p=[0.80, 0.15, 0.05],
    )
    adjuvant = rng.choice(["", "HORMONAL", "RADIATION"], size=N_TOTAL, p=[0.65, 0.25, 0.10])

    # ── Pathological Gleason / stage ─────────────────────────────────────────
    patho_total = rng.choice([7, 8, 9, 10], size=N_TOTAL, p=[0.30, 0.40, 0.25, 0.05])
    patho_g1 = np.where(patho_total >= 8, 4, rng.choice([3, 4], size=N_TOTAL))
    patho_g2 = patho_total - patho_g1
    margins = rng.choice(["POSITIVE", "NEGATIVE"], size=N_TOTAL, p=[0.42, 0.58])
    lymph_node = rng.choice(["NEGATIVE", "POSITIVE"], size=N_TOTAL, p=[0.92, 0.08])
    patho_stage = rng.choice(
        ["pT2", "pT2a", "pT2b", "pT2c", "pT3a", "pT3b", "pT4"],
        size=N_TOTAL,
        p=[0.05, 0.10, 0.10, 0.20, 0.30, 0.20, 0.05],
    )

    # ── Pre-derived (months) time fields — NO date arithmetic ───────────────
    time_to_met = np.where(
        metastatic == 1,
        np.round(rng.uniform(15, 110, size=N_TOTAL), 1),
        np.nan,
    )
    followup = np.where(
        metastatic == 1,
        np.round(time_to_met + rng.uniform(20, 100, size=N_TOTAL), 1),
        np.round(rng.uniform(95, 165, size=N_TOTAL), 1),
    )
    time_to_psa = np.where(
        rng.random(N_TOTAL) < 0.6,
        np.round(rng.uniform(6, 80, size=N_TOTAL), 1),
        np.nan,
    )
    time_to_crpc = np.where(
        metastatic == 1,
        np.round(time_to_met + rng.uniform(10, 40, size=N_TOTAL), 1),
        np.nan,
    )

    status = np.where(
        metastatic == 1,
        rng.choice(
            ["ALIVE METASTATIC DISEASE", "DEAD METASTATIC DISEASE"],
            size=N_TOTAL, p=[0.6, 0.4],
        ),
        rng.choice(
            ["ALIVE NO EVIDENCE OF DISEASE", "ALIVE PSA RECURRENCE"],
            size=N_TOTAL, p=[0.75, 0.25],
        ),
    )
    met_location = np.where(
        metastatic == 1,
        rng.choice(
            ["BONE", "LYMPH NODE", "LUNG", "LIVER", "BONE AND LYMPH NODE"],
            size=N_TOTAL,
        ),
        "",
    )
    met_treatment = np.where(
        metastatic == 1,
        rng.choice(
            ["HORMONE THERAPY", "RADIATION HORMONAL", "ADT + DOCETAXEL"],
            size=N_TOTAL,
        ),
        "",
    )
    crpc_treatment = np.where(
        metastatic == 1,
        rng.choice(
            ["", "DOCETAXEL", "ENZALUTAMIDE", "ABIRATERONE"],
            size=N_TOTAL,
        ),
        "",
    )

    clinical = pd.DataFrame({
        "Patient ID":                              pid,
        "Metastasis (1=yes, 0=no)":                metastatic,
        "Ancestry (inferred)":                     ancestry,
        "Age at diagnosis (years)":                age,
        "PSA at diagnosis (ng/mL)":                psa,
        "Biopsy Gleason score 1":                  biopsy_g1,
        "Biopsy Gleason score 2":                  biopsy_g2,
        "Biopsy total Gleason score":              biopsy_total,
        "Biopsy positive cores":                   biopsy_pos,
        "Biopsy total cores":                      biopsy_tot,
        "Clinical stage":                          clinical_stage,
        "Main treatment":                          main_treatment,
        "Adjuvant therapy":                        adjuvant,
        "Pathological Gleason score 1":            patho_g1,
        "Pathological Gleason score 2":            patho_g2,
        "Pathological total Gleason score":        patho_total,
        "Surgical margins":                        margins,
        "Lymph node status":                       lymph_node,
        "Pathological stage":                      patho_stage,
        "Follow-up duration (months)":             _fmt_months(followup),
        "Status at follow-up":                     status,
        "Time to PSA recurrence (months)":         _fmt_months(time_to_psa),
        "Time to metastasis (months)":             _fmt_months(time_to_met),
        "Metastasis location":                     met_location,
        "Metastasis treatment":                    met_treatment,
        "Time to castration resistance (months)":  _fmt_months(time_to_crpc),
        "Castration resistance treatment":         crpc_treatment,
    })

    # ── gnsRV burden per patient: weak metastatic association ───────────────
    # (effect size intentionally smaller than the real OR = 26.80)
    base_rate = 0.6 + 0.7 * metastatic
    ddr_total = rng.poisson(lam=base_rate)
    ddr_syn = rng.poisson(lam=0.5, size=N_TOTAL)
    extra = pd.DataFrame(
        {
            "metastatic":                metastatic,
            "DDR_total_gnsrv":           ddr_total,
            "DDR_total_synonymous_rv":   ddr_syn,
        },
        index=pid,  # patient_id as index → emitted as unnamed first column
    )
    return clinical, extra


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    clinical, extra = simulate(args.seed)
    clinical_out = OUT_DIR / "clinical_synthetic.tsv"
    extra_out = OUT_DIR / "gnsrv_burden_per_patient_synthetic.tsv"
    clinical.to_csv(clinical_out, sep="\t", index=False)
    extra.to_csv(extra_out, sep="\t")  # index=True so patient_id is the leading unnamed col
    print(f"Wrote {clinical_out.name}                ({len(clinical)} rows, {len(clinical.columns)} cols)")
    print(f"Wrote {extra_out.name}    ({len(extra)} rows, {len(extra.columns) + 1} cols incl. patient_id)")


if __name__ == "__main__":
    main()
