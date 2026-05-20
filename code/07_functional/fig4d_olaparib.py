"""
Figure 4D — BRCA2 I1962T Olaparib sensitivity
==============================================

Bar chart of olaparib IC50 (µM) across three HEK293T-derived cell lines,
with per-biological-replicate scatter dots, SEM error bars, one-way ANOVA
omnibus annotation, and pairwise Welch t-test significance brackets
(Bonferroni × 3) (Methods §2.5; Results, functional paragraph; Fig. 4D
caption).

Cell lines:
  - HEK293T parental (wild-type BRCA2)
  - HEK293T BRCA2 I1962T knock-in (the candidate variant being tested)
  - HEK293T BRCA1 Q1200X knock-in (positive control for PARP-inhibitor sensitivity)

Reported IC50 values (mean ± SEM, n = 3 biological replicates):
  HEK293T parental:        4.10 ± 0.05 µM
  HEK293T BRCA2 I1962T:    2.50 ± 0.12 µM
  HEK293T BRCA1 Q1200X:    2.10 ± 0.27 µM

IC50 provenance — IMPORTANT NOTE
--------------------------------
IC50 values shown in Figure 4D are computed in GraphPad Prism (4-parameter
logistic, variable slope) from the raw 8-point dose-response data and are
stored in the input TSV as precomputed per-biological-replicate values
(`measurement_type == "IC50_uM_per_BR"`). This script reads those Prism-fit
values DIRECTLY rather than re-fitting in Python.

Rationale:
  1. The published figure uses the Prism-fit IC50s; using them here
     ensures the public-repo figure matches the publication exactly.
  2. Prism's 4PL implementation (constrained nonlinear regression with
     proprietary optimisation) yields numerically different IC50 estimates
     from scipy.optimize.curve_fit for low-curvature dose-response data,
     and we do not want apparent reproducibility discrepancies that reflect
     tool choice rather than the underlying measurements.
  3. The raw per-concentration absorbance measurements (`measurement_type
     == "absorbance"`) are still included in the input TSV for transparency
     and so readers can independently re-fit if they wish, but they are
     not consumed by this script's IC50 plotting code.

To independently reproduce the IC50 fits from raw dose-response data:
  - Open the dose-response file in GraphPad Prism (v10+).
  - XY analysis → Nonlinear regression → Dose-response (variable slope,
    four parameters), constrain top = 100, bottom = 0.
  - Read off IC50, 95% CI, and LogIC50 SEM per cell line.
  - Compare against the values in the TSV's IC50_uM_per_BR rows.

Statistics shown in the figure:
  - One-way ANOVA (omnibus across the three cell lines)
  - Pairwise Welch's t-test, Bonferroni-corrected × 3 comparisons
  - Significance stars: *** p<0.001, ** p<0.01, * p<0.05, ns otherwise

Output:
  figures/fig4d_olaparib.pdf
  figures/fig4d_olaparib.png

Usage:
  python code/07_functional/fig4d_olaparib.py
"""

from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
DATA_FILE = (
    REPO / "data" / "raw" / "functional_assays" / "fig4d_olaparib_dose_response.tsv"
)
FIGURES_DIR = REPO / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

# Cell-line order and display labels for Fig 4D
CELL_LINES = ["HEK293T_Parental", "BRCA2_I1962T", "BRCA1_Q1200X"]
LABELS = ["Parental\n(HEK293T)", "BRCA2\nI1962T", "BRCA1\nQ1200X\n(+ctrl)"]

# Editorial palette: neutral grey (parental), navy (candidate variant),
# dark cyan (positive control)
COLORS = ["#9E9E9E", "#1565C0", "#00838F"]
COLORS_LITE = ["#E0E0E0", "#90CAF9", "#80DEEA"]

FONT = "Arial"


def load_ic50_per_BR(path: Path) -> list[np.ndarray]:
    """Read the long-format dose-response file and return per-BR IC50 values
    for each cell line (in CELL_LINES order). IC50 values stored under
    `measurement_type == "IC50_uM_per_BR"` are the Prism 4PL fits — see the
    module docstring for provenance.

    IC50 values are pre-fit in Prism. This script consumes them directly;
    do NOT re-fit in Python for the published figure (the constrained 4PL
    optimisation differs numerically across tools).
    """
    df = pd.read_csv(path, sep="\t", comment="#")
    per_br = df[df["measurement_type"] == "IC50_uM_per_BR"]
    values = []
    for cell in CELL_LINES:
        vals = per_br.loc[per_br["cell_line"] == cell, "value"].to_numpy(dtype=float)
        if len(vals) == 0:
            raise SystemExit(
                f"No IC50_uM_per_BR rows found for cell line '{cell}' in {path}"
            )
        values.append(vals)
    return values


def fmt_star(p: float) -> str:
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def fmt_p_num(p: float) -> str:
    if p >= 0.001:
        return f"p = {p:.3f}"
    exp = int(np.floor(np.log10(p)))
    coef = p / 10**exp
    return f"p = {coef:.2f}×10$^{{{exp}}}$"


def print_summary(
    VALUES, means, sems, f_stat, p_anova, key_pairs, raw_p, bonf_p
) -> None:
    SEP = "─" * 65
    print(SEP)
    print("Olaparib IC50 — Statistical Summary (Fig. 4D)")
    print(SEP)
    print("Design : 3 independent biological replicates per cell line")
    print("IC50s  : Prism 4PL fits (read from input TSV; not re-fit in Python)")
    print("Test   : One-way ANOVA (omnibus) + pairwise Welch t-test, Bonferroni ×3")
    print()
    print("Group means ± SEM (µM Olaparib):")
    for lab, m, s, v in zip(LABELS, means, sems, VALUES):
        lab_flat = lab.replace("\n", " ")
        print(
            f"  {lab_flat:<28}: {m:.3f} ± {s:.3f}  per-BR: {np.array2string(v, precision=2)}"
        )
    print()
    print(f"One-way ANOVA: F={f_stat:.3f}, {fmt_p_num(p_anova)}")
    print()
    print(f"{'Comparison':<40} {'raw p':>8}  {'Bonf p':>8}  sig")
    print("─" * 62)
    for i, j in key_pairs:
        a = LABELS[i].replace("\n", " ")
        b = LABELS[j].replace("\n", " ")
        rp = raw_p[(i, j)]
        bp = bonf_p[(i, j)]
        print(f"  {a} vs {b:<16}  {rp:>8.4f}    {bp:>8.4f}  {fmt_star(bp)}")
    print(SEP)


def main() -> None:
    VALUES = load_ic50_per_BR(DATA_FILE)

    # ── Statistics ──────────────────────────────────────────────────────────
    f_stat, p_anova = stats.f_oneway(*VALUES)
    key_pairs = list(combinations(range(len(LABELS)), 2))
    raw_p, bonf_p = {}, {}
    for i, j in key_pairs:
        _, p = stats.ttest_ind(VALUES[i], VALUES[j], equal_var=False)
        raw_p[(i, j)] = float(p)
        bonf_p[(i, j)] = min(float(p) * len(key_pairs), 1.0)

    means = [v.mean() for v in VALUES]
    sems = [stats.sem(v) for v in VALUES]

    print_summary(VALUES, means, sems, f_stat, p_anova, key_pairs, raw_p, bonf_p)

    # ── Plot ────────────────────────────────────────────────────────────────
    plt.rcParams.update({"font.family": FONT, "font.size": 11})
    fig, ax = plt.subplots(figsize=(5.5 * 0.8, 5.5))
    x = np.arange(len(LABELS))
    bar_width = 0.52

    # Bars (mean IC50)
    ax.bar(
        x,
        means,
        width=bar_width,
        color=COLORS,
        edgecolor="white",
        linewidth=0.8,
        alpha=0.88,
        zorder=2,
    )

    # Per-BR scatter dots
    rng = np.random.default_rng(42)
    for i, v in enumerate(VALUES):
        jitter = rng.uniform(-0.09, 0.09, len(v))
        ax.scatter(
            x[i] + jitter,
            v,
            color=COLORS_LITE[i],
            s=48,
            zorder=3,
            edgecolors=COLORS[i],
            linewidths=0.8,
        )

    # SEM error bars (on top of scatter)
    ax.errorbar(
        x,
        means,
        yerr=sems,
        fmt="none",
        color="#222222",
        capsize=4,
        capthick=1.1,
        linewidth=1.1,
        zorder=4,
    )

    # ── Significance brackets ───────────────────────────────────────────────
    y_top = max(m + s for m, s in zip(means, sems))
    h = y_top * 0.020
    bracket_specs = [
        (0, 1, y_top * 1.10),
        (0, 2, y_top * 1.22),
        (1, 2, y_top * 1.34),
    ]
    for i, j, yb in bracket_specs:
        xa, xb = x[i], x[j]
        bp = bonf_p[(i, j)]
        ax.plot(
            [xa, xa, xb, xb],
            [yb, yb + h, yb + h, yb],
            color="#444444",
            linewidth=0.9,
        )
        ax.text(
            (xa + xb) / 2,
            yb + h * 1.4,
            f"{fmt_star(bp)}  ({fmt_p_num(bp)})",
            ha="center",
            va="bottom",
            fontsize=14,  # fontsize=9,
            fontweight="bold",
            color="#222222",
        )

    # ── Axes ────────────────────────────────────────────────────────────────
    ax.set_xticks(x)
    ax.set_xticklabels(LABELS, fontsize=10)
    ax.set_ylabel("IC$_{50}$ (µM Olaparib)", fontsize=11)
    ax.set_ylim(0, y_top * 1.58)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(1))
    ax.yaxis.grid(True, linestyle="--", linewidth=0.45, alpha=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", labelsize=10)

    ax.set_title(
        "BRCA2 I1962T: Olaparib sensitivity\n(HEK293T prime-edited clones)",
        fontsize=12,
        fontweight="bold",
        pad=8,
    )

    ax.text(
        0.97,
        0.97,
        f"ANOVA {fmt_p_num(p_anova)}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=14,
        color="#666666",
        style="italic",  # fontsize=9
    )

    plt.tight_layout()
    pdf = FIGURES_DIR / "fig4d_olaparib.pdf"
    png = FIGURES_DIR / "fig4d_olaparib.png"
    plt.savefig(pdf, dpi=300, bbox_inches="tight")
    plt.savefig(png, dpi=300, bbox_inches="tight")
    print(f"\nSaved: {pdf}")
    print(f"Saved: {png}")


if __name__ == "__main__":
    main()
