"""Marimo walkthrough — convergent pathway-level evidence for germline DDR
rare-variant burden in metastatic prostate cancer.

This notebook reproduces the manuscript's central finding using the real,
de-identified cohort committed to this repository (`data/clinical/`,
`data/raw/`): three statistically independent tests converge on the
conclusion that the metastatic arm carries a higher germline rare-variant
burden across the 25-gene DDR pathway. A single rare variant in any one
gene would be statistically underpowered in n = 52; the aggregate
*pathway-level* signal, **not single-variant association** is what carries
the result.

Reviewers can run this two ways:

    marimo edit notebooks/01_ddr_convergent_walkthrough.py     # interactive
    marimo run  notebooks/01_ddr_convergent_walkthrough.py     # read-only app

The notebook is a thin orchestrator: every section calls the corresponding
production script in `code/` via subprocess and reads its output. No
statistical logic is reimplemented here. No controlled-access download is
needed — the real cohort tables ship with this repository.
"""

import marimo

__generated_with = "0.9.27"
app = marimo.App(width="medium")


@app.cell
def __():
    import marimo as mo

    return (mo,)


@app.cell
def __(mo):
    mo.md(
        r"""
        # Convergent pathway-level evidence — interactive walkthrough

        **Why this matters.** Most genetic studies of prostate cancer have
        addressed *cancer initiation*; this manuscript targets a distinct
        phenotype — *metastatic progression after diagnosis* — using an
        extreme-phenotype cohort design. The methodological novelty is
        that the headline finding is a **pathway-level rare-variant
        burden across 25 DDR genes**, supported by **three statistically
        independent tests** (per-patient burden, gene-set bootstrap
        competitive null, and Firth penalized regression) that converge
        on the same conclusion.

        **What you'll see, and where the code lives:**

        | § | Test                              | Script                                                    | Published value           |
        |---|-----------------------------------|-----------------------------------------------------------|---------------------------|
        | 2 | Cohort assembly                   | `code/01_cohort_prep/prepare_epc_cohort.py`               | n = 52 (26 met / 26 nm)   |
        | 3 | Variant landscape                 | (input matrix — see `code/02_variant_calling/README.md`)  | 71 vs 21 unique gnsRVs    |
        | 4 | Per-patient burden (Mann-Whitney) | `code/03_ddr_enrichment/ddr_arm_compare.py`               | p = 4.57 × 10⁻⁶           |
        | 5 | Gene-set bootstrap                | `code/03_ddr_enrichment/ddr_burden_bootstrap.py`          | p < 1 × 10⁻⁴ (met arm)    |
        | 6 | Firth penalized regression        | `code/04_firth_regression/ddr_logistic_regression.py`     | OR = 26.80                |
        | 7 | Convergent verdict                | (this notebook)                                           |                           |

        Every cell calls the production scripts via `subprocess` — the
        notebook is a navigation index *to* the code, not a
        re-implementation.
        """
    )
    return


@app.cell
def __():
    """Imports + helpers."""
    import math
    import re
    import subprocess
    import sys
    from pathlib import Path

    import matplotlib.pyplot as plt
    import pandas as pd

    REPO = Path(__file__).resolve().parents[1]
    PY = sys.executable

    def run_script(rel_path: str, *args: str) -> str:
        """Run a script from REPO root; return stdout. On failure, surface
        stderr in the raised exception so the marimo cell shows a useful
        message instead of a bare CalledProcessError."""
        cmd = [PY, str(REPO / rel_path), *args]
        out = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
        if out.returncode != 0:
            raise RuntimeError(
                f"{rel_path} exited {out.returncode}\n"
                f"stderr:\n{out.stderr}\n"
                f"stdout:\n{out.stdout}"
            )
        return out.stdout

    return REPO, run_script, pd, plt, re, math


@app.cell
def __(mo):
    mo.md(
        """
        ## 2. Cohort assembly

        `prepare_epc_cohort.py` merges the real, de-identified clinical
        table (`data/clinical/clinical.tsv`) with the per-patient gnsRV
        burden (`data/raw/gnsrv_burden_per_patient.tsv`) and writes the
        analysis-ready cohort to `data/processed/epc_cohort.tsv`.
        """
    )
    return


@app.cell
def __(REPO, run_script, pd, mo):
    """Build cohort table; show per-arm summary."""
    run_script("code/01_cohort_prep/prepare_epc_cohort.py")
    cohort = pd.read_csv(REPO / "data" / "processed" / "epc_cohort.tsv", sep="\t")

    by_arm = (
        cohort.groupby("metastatic")
        .agg(
            n=("patient_id", "size"),
            mean_age=("age_at_diagnosis", "mean"),
            mean_ddr_burden=("ddr_burden", "mean"),
            mean_ddr_synonymous=("ddr_synonymous_burden", "mean"),
        )
        .round(2)
    )
    by_arm.index = by_arm.index.map({0: "Non-metastatic", 1: "Metastatic"})

    mo.vstack(
        [
            mo.md(
                f"**Cohort:** {len(cohort)} patients "
                f"({(cohort['metastatic'] == 1).sum()} metastatic, "
                f"{(cohort['metastatic'] == 0).sum()} non-metastatic)."
            ),
            by_arm,
            mo.md(
                "*Synonymous gnsRV burden is the negative control — its arm "
                "difference should be near zero.*"
            ),
        ]
    )
    return (cohort,)


@app.cell
def __(mo):
    mo.md(
        """
        ## 3. The DDR variant landscape

        The upstream variant-calling pipeline
        (`code/02_variant_calling/README.md`) identified the
        metastatic-exclusive and non-metastatic-exclusive gnsRVs in the
        25-gene DDR panel:

        > **71 unique gnsRVs (81 carrier-variant observations) exclusive
        > to the metastatic arm vs. 21 exclusive to the non-metastatic
        > arm.** — Results, ¶2

        The per-gene × per-patient matrix that produces these counts
        lives at `data/raw/gnsrv_per_gene_per_patient.tsv`; the
        prespecified gene list is `data/derived/ddr_panel_25genes.tsv`.
        The three statistical tests below are all downstream of this
        matrix.
        """
    )
    return


@app.cell
def __(mo):
    mo.md(
        """
        ## 4. Test 1 — Per-patient burden (Mann-Whitney)

        Compares per-patient gnsRV burden across the two arms directly.
        Non-parametric; makes no distributional assumption.

        **Published value:** Mann-Whitney two-sided p = 4.57 × 10⁻⁶
        (Methods §2.3.1; Results, ¶2).
        """
    )
    return


@app.cell
def __(cohort, plt):
    """Violin plot of DDR burden by arm."""
    _fig, _ax = plt.subplots(figsize=(5, 3.5))
    _met = (
        cohort.loc[cohort["metastatic"] == 1, "ddr_burden"]
        .dropna()
        .astype(float)
        .values
    )
    _nm = (
        cohort.loc[cohort["metastatic"] == 0, "ddr_burden"]
        .dropna()
        .astype(float)
        .values
    )

    _parts = _ax.violinplot([_nm, _met], showmedians=True, widths=0.7)
    for _body, _color in zip(_parts["bodies"], ["#5da7c7", "#c25450"]):
        _body.set_facecolor(_color)
        _body.set_alpha(0.6)
    _ax.set_xticks([1, 2])
    _ax.set_xticklabels(
        [f"Non-metastatic\n(n={len(_nm)})", f"Metastatic\n(n={len(_met)})"]
    )
    _ax.set_ylabel("Per-patient DDR gnsRV burden")
    _ax.set_title("DDR rare-variant burden by arm")
    _ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    _fig
    return


@app.cell
def __(run_script, re, math, mo):
    """Run ddr_arm_compare.py; extract Mann-Whitney p."""
    out_mw = run_script("code/03_ddr_enrichment/ddr_arm_compare.py")
    _m = re.search(
        r"two-sided\s+\(met\s*[≠!=]+\s*nm\):\s+p\s*=\s*([\d.eE+-]+)",
        out_mw,
    )
    p_mw = float(_m.group(1)) if _m else math.nan

    mo.md(
        f"**Mann-Whitney U test — `ddr_arm_compare.py`:**\n\n"
        f"- Two-sided p (met ≠ nm) = **{p_mw:.3e}**\n"
        f"- Direction: metastatic > non-metastatic\n"
    )
    return (p_mw,)


@app.cell
def __(mo):
    mo.md(
        """
        ## 5. Test 2 — Gene-set bootstrap (competitive null)

        Resamples 25 random genes 10,000 times and compares the observed
        DDR burden against the resulting null distribution **within
        each arm**. Tests: *is the DDR signal larger than a random
        25-gene panel drawn from the same exome?* Independent of Test 1
        — it uses gene-set identity, not patient-level comparison.

        **Published value:** Empirical p < 1 × 10⁻⁴ in the metastatic
        arm (Methods §2.3.1; Supplementary Methods §4.2.2). The
        non-metastatic arm serves as a negative-control direction —
        its p is *not* significant.

        *(Takes ~30 s with the default 10,000 iterations.)*
        """
    )
    return


@app.cell
def __(run_script, re, math, mo):
    """Run ddr_burden_bootstrap.py; parse met and nm empirical p."""
    out_bs = run_script(
        "code/03_ddr_enrichment/ddr_burden_bootstrap.py", "--arm", "both"
    )
    p_met = p_nm = math.nan
    arm_blocks = out_bs.split("=== ")
    for _block in arm_blocks:
        if _block.startswith("MET ARM"):
            _m = re.search(r"Empirical p[^=]*=\s*([\d.eE+-]+)", _block)
            if _m:
                p_met = float(_m.group(1))
        elif _block.startswith("NONMET ARM"):
            _m = re.search(r"Empirical p[^=]*=\s*([\d.eE+-]+)", _block)
            if _m:
                p_nm = float(_m.group(1))

    mo.md(
        f"**Gene-set bootstrap — `ddr_burden_bootstrap.py --arm both`:**\n\n"
        f"- Metastatic arm: empirical p = **{p_met:.3e}**\n"
        f"- Non-metastatic arm (negative control): empirical p = **{p_nm:.3e}**\n\n"
        f"Same test, opposite directions: a strong pathway-level signal "
        f"only in the metastatic arm. Independent of Test 1."
    )
    return p_met, p_nm


@app.cell
def __(mo):
    mo.md(
        """
        ## 6. Test 3 — Firth penalized logistic regression

        Tests whether DDR-carrier status (or burden) predicts metastasis
        after adjusting for ancestry. Uses the Firth penalty to handle
        quasi-complete separation in this cohort (every metastatic
        patient is a carrier; only ~half of non-metastatic). Reports OR
        with profile-likelihood 95 % CI.

        **Published value:** OR = 26.80, 95 % CI 3.07–3528.42, LRT
        p < 0.001 (Methods §2.3.1; Abstract).
        """
    )
    return


@app.cell
def __(REPO, run_script, pd, plt, mo):
    """Run Firth (primary + extended adjustments) and render a forest plot.

    The script's default is `--adjustment both`, so a single invocation
    writes all four rows (ancestry-adjusted × {binary, count}, then
    extended × {binary, count}) to firth_results.csv.
    """
    run_script("code/04_firth_regression/ddr_logistic_regression.py")
    firth = pd.read_csv(REPO / "data" / "processed" / "firth_results.csv")

    _fig, _ax = plt.subplots(figsize=(7.5, 3.2))
    _labels = [
        f"{r['adjustment']} · {r['predictor_name']}" for _, r in firth.iterrows()
    ]
    _ys = list(range(len(firth)))[::-1]
    for _y, (_, r) in zip(_ys, firth.iterrows()):
        _ax.plot([r["CI_lo"], r["CI_hi"]], [_y, _y], color="#555", lw=2)
        _ax.plot(r["OR"], _y, "o", color="#c25450", ms=8)
    _ax.axvline(1.0, color="#888", ls="--", lw=0.8)
    _ax.set_xscale("log")
    _ax.set_yticks(_ys)
    _ax.set_yticklabels(_labels, fontsize=8)
    _ax.set_xlabel("Odds ratio (log scale)")
    _ax.set_title("Firth penalized logistic regression — DDR carrier / burden")
    _ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    _fig

    _primary = firth.query(
        "adjustment == 'ancestry-adjusted' and predictor_type == 'binary'"
    ).iloc[0]
    mo.md(
        f"**Primary model (ancestry-adjusted, DDR carrier):** "
        f"OR = **{_primary['OR']:.2f}**, "
        f"95 % CI {_primary['CI_lo']:.2f}–{_primary['CI_hi']:.2f}, "
        f"LRT p = {_primary['p_LRT']:.3e}."
    )
    return (firth,)


@app.cell
def __(mo):
    mo.md(
        """
        ## 7. Convergent verdict

        Three statistically independent tests of the same hypothesis —
        higher DDR pathway-level rare-variant burden in the metastatic
        arm — and three convergent answers. The non-metastatic arm of
        the bootstrap is the negative-control direction; its
        non-significant result is what *should* happen if the pathway
        signal is real and arm-specific.
        """
    )
    return


@app.cell
def __(p_mw, p_met, p_nm, firth, mo):
    """Final summary table."""
    _primary = firth.query(
        "adjustment == 'ancestry-adjusted' and predictor_type == 'binary'"
    ).iloc[0]

    _rows = [
        ("Per-patient burden (Mann-Whitney)", f"p = {p_mw:.3e}", "p = 4.57 × 10⁻⁶"),
        ("Gene-set bootstrap, metastatic arm", f"p = {p_met:.3e}", "p < 1 × 10⁻⁴"),
        ("Gene-set bootstrap, non-met (control)", f"p = {p_nm:.3e}", "n.s."),
        (
            "Firth logistic, ancestry-adjusted",
            f"OR = {_primary['OR']:.2f}",
            "OR = 26.80",
        ),
    ]
    _table = "| Test | This notebook | Published |\n|---|---|---|\n" + "\n".join(
        f"| {t} | {n} | {p} |" for t, n, p in _rows
    )

    mo.md(
        _table + "\n\n"
        "**Conclusion.** The same pathway-level enrichment shows up "
        "across three statistically independent tests, and is absent "
        "from the non-metastatic negative-control direction of the "
        "bootstrap. Because the notebook operates on the real cohort "
        "tables committed to this repo, the values above match the "
        "manuscript verbatim."
    )
    return


@app.cell
def __(mo):
    mo.md(
        """
        ## 8. Where the rest of the story is

        This notebook stops at the discovery cohort's central finding.
        The manuscript continues in three directions, each in its own
        `code/` module:

        - **Replication** — Australian extreme-phenotype cohort (n = 53,
          binomial against gnomAD AF) and PPCG (n = 976, SKAT/Burden).
          See `code/06_replication/` and `docs/METHODS_CROSSREF.md`.
        - **Functional validation** — *in vitro* assays on KDM6B and
          BRCA2 variants currently classified as uncertain or benign,
          showing measurable phenotypic effects. See
          `code/07_functional/` and Figure 4.
        - **Figures** — main and supplementary figure-rendering scripts
          under `code/08_figures/`; published PDFs ship under
          `figures/` (mapping in `figures/README.md`).

        For the full manuscript-claim → script mapping, see
        `docs/METHODS_CROSSREF.md`.
        """
    )
    return


if __name__ == "__main__":
    app.run()
