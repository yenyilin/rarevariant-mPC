"""
Smoke tests that exercise the full pipeline on the synthetic cohort.
These do not assert manuscript values (synthetic effect is intentionally
weaker than the real OR=26.80) — they assert the pipeline runs and
produces sane outputs.
"""

from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]

# Published one-sided binomial p-values from Table 2 / Fig. 3A (Results, ¶5).
# Stored here, in the test, rather than in the input counts TSV — see the
# rationale in data/derived/australian_epc_counts.tsv ("Expected carrier
# counts and one-sided binomial p-values are computed … they are not
# stored here to avoid drift between source data and derived values").
PUBLISHED_AUSTRALIAN_P = {
    "metastatic_bone_visceral":        0.049,
    "metastatic_plus_BCR":             0.038,
    "metastatic_plus_BCR_plus_nodes":  0.027,
    "non_metastatic_control":          0.449,
}


def run(cmd: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, check=True, cwd=REPO)


@pytest.fixture(scope="session", autouse=True)
def synthetic_cohort():
    run([sys.executable, "data/synthetic/generate_synthetic_clinical.py", "--seed", "42"])
    run([sys.executable, "code/01_cohort_prep/prepare_epc_cohort.py", "--synthetic"])
    return REPO / "data" / "processed" / "epc_cohort.tsv"


def test_synthetic_files_exist():
    for name in ("clinical_synthetic.tsv", "gnsrv_burden_per_patient_synthetic.tsv"):
        assert (REPO / "data" / "synthetic" / name).exists()


def test_processed_table_shape(synthetic_cohort):
    df = pd.read_csv(synthetic_cohort, sep="\t")
    assert len(df) == 52
    assert set(df["metastatic"].unique()) == {0, 1}
    assert df["metastatic"].sum() == 26
    expected_cols = {
        "patient_id", "metastatic", "ancestry", "age_at_diagnosis",
        "ddr_burden", "time_to_metastasis", "followup_duration_months",
    }
    assert expected_cols.issubset(df.columns)


def test_firth_regression_runs(synthetic_cohort):
    out = REPO / "data" / "processed" / "firth_results.csv"
    if out.exists():
        out.unlink()
    run([sys.executable, "code/04_firth_regression/ddr_logistic_regression.py"])
    assert out.exists()
    res = pd.read_csv(out).iloc[0]
    assert res["OR"] > 0
    assert res["CI_lo"] < res["OR"] < res["CI_hi"]
    assert 0 <= res["p_LRT"] <= 1


def test_australian_binomial_recovers_published_pvalues():
    """Running the binomial test against the auto-computed p_carrier (≈0.311)
    must recover the published Table 2 / Fig. 3A one-sided p-values."""
    out = REPO / "data" / "processed" / "australian_binomial_results.csv"
    if out.exists():
        out.unlink()
    run([sys.executable, "code/06_replication/australian_binomial.py"])
    res = pd.read_csv(out)
    for _, row in res.iterrows():
        subgroup = row["subgroup"]
        published = PUBLISHED_AUSTRALIAN_P[subgroup]
        recomputed = float(row["p_recomputed"])
        assert abs(published - recomputed) < 0.01, (
            f"{subgroup}: published {published} vs recomputed {recomputed}"
        )


def test_followup_plot_generates_figure(synthetic_cohort):
    pdf = REPO / "figures" / "followup_comparison.pdf"
    if pdf.exists():
        pdf.unlink()
    run([sys.executable, "code/05_time_bias/ddr_followup_plot.py"])
    assert pdf.exists() and pdf.stat().st_size > 0
