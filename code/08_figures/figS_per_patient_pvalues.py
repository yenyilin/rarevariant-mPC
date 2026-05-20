"""
Supplementary Figure S4 — per-patient empirical p-values for DDR burden vs.
within-patient bootstrap null distributions, in the metastatic and
non-metastatic arms side by side.

The bootstrap null is computed on-the-fly from the gene × patient matrix
(the same input used by ddr_burden_bootstrap.py and ddr_per_patient_test.py);
with a fixed seed (default 42) the per-patient null distributions match
those used in Supplementary Methods §4.2.4 and the per-patient test script.

Default behaviour excludes both the DDR panel and the Shyr 2014 FLAGS
list (~150 hyper-variable genes) from the random-panel pool, matching
the canonical configuration reported in §4.2.3–§4.2.4. Pass --no-exclude
to disable FLAGS exclusion for sensitivity analysis (matches §4.2.5).

Computes per-patient empirical p-values with Laplace smoothing, and
renders a two-panel forest plot of -log10(p_emp) ranked by significance
within each arm. Vertical reference lines mark p = 0.05 and p = 0.001.

This is the figure complement to code/03_ddr_enrichment/ddr_per_patient_test.py.

Usage
-----
    python code/08_figures/figS_per_patient_pvalues.py
    python code/08_figures/figS_per_patient_pvalues.py --no-exclude
    python code/08_figures/figS_per_patient_pvalues.py --n-iterations 50000
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
FIGURES_DIR = REPO / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

DEFAULT_MET = (
    2,
    4,
    3,
    5,
    3,
    1,
    5,
    1,
    7,
    1,
    2,
    1,
    3,
    1,
    5,
    1,
    4,
    4,
    3,
    7,
    5,
    3,
    6,
    2,
    1,
    1,
)
DEFAULT_NONMET = (
    1,
    1,
    0,
    0,
    0,
    0,
    0,
    1,
    2,
    0,
    1,
    0,
    3,
    1,
    2,
    0,
    0,
    2,
    1,
    2,
    0,
    0,
    1,
    1,
    2,
    0,
)

PANEL_SIZE = 25  # number of genes per random panel (matches DDR panel size)
N_MET = 26  # first N_MET cols of the matrix  = metastatic arm
N_NONMET = 26  # next  N_NONMET cols              = non-metastatic arm

C_SIG_STRONG = "#a51e1e"  # p < 0.001 — dark red
C_SIG = "#d97766"  # 0.001 ≤ p < 0.05 — light red
C_NS = "#9aa3ad"  # p ≥ 0.05 — grey
C_REF = "#444444"  # reference lines


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument(
        "--matrix",
        type=Path,
        default=REPO / "data" / "raw" / "gnsrv_per_gene_per_patient.tsv",
        help=f"Wide gene × patient TSV (col 1 = gene name → row index). "
        f"First {N_MET} patient columns = met arm; next {N_NONMET} = non-met.",
    )
    p.add_argument(
        "--ddr-panel",
        type=Path,
        default=REPO / "data" / "derived" / "ddr_panel_25genes.tsv",
        help="TSV with a 'gene' column listing the 25 DDR genes. "
        "These are excluded from the random pool.",
    )
    p.add_argument(
        "--include-ddr-in-pool",
        action="store_true",
        help="Allow DDR genes to be drawn into random panels (default: exclude them).",
    )
    p.add_argument(
        "--exclude-genes",
        type=Path,
        default=REPO / "data" / "derived" / "flags_shyr2014.tsv",
        help="TSV (header `gene`) of additional genes to exclude from the "
        "random pool — typically a FLAGS list (Shyr 2014). "
        "Lines starting with `#` are ignored. Use --no-exclude to disable.",
    )
    p.add_argument(
        "--no-exclude",
        action="store_true",
        help="Disable the FLAGS-based exclusion (use the full gene pool minus DDR).",
    )
    p.add_argument("--n-iterations", type=int, default=10_000)
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def load_exclude_list(path: Path) -> set[str]:
    """Load a one-column TSV of gene symbols to exclude. Tolerates `#` comments
    and either a header line `gene` or no header."""
    if not path.exists():
        print(
            f"WARNING: exclusion list {path} not found — proceeding with no exclusions."
        )
        return set()
    genes: set[str] = set()
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#") or s.lower() == "gene":
            continue
        genes.add(s.split("\t")[0])
    return genes


def build_null_matrix(
    matrix: pd.DataFrame,
    pool: list[str],
    n_iterations: int,
    seed: int,
) -> np.ndarray:
    """Sample n_iterations random PANEL_SIZE-gene panels and compute per-patient
    burdens. Returns an (n_iterations, n_patients) integer array. Same RNG
    construction as ddr_burden_bootstrap.py and ddr_per_patient_test.py —
    seed=42 produces identical gene draws across all three scripts."""
    rng = np.random.default_rng(seed)
    n_patients = matrix.shape[1]
    null_matrix = np.empty((n_iterations, n_patients), dtype=np.int64)
    for i in range(n_iterations):
        sampled = rng.choice(pool, size=PANEL_SIZE, replace=False)
        null_matrix[i] = matrix.loc[sampled].to_numpy().sum(axis=0)
    return null_matrix


def per_patient_p(observed: np.ndarray, null_matrix: np.ndarray) -> np.ndarray:
    """Empirical one-sided p_i = (#null_rows where null_i >= observed_i + 1) / (N+1)."""
    n_iter = null_matrix.shape[0]
    n_ge = (null_matrix >= observed[None, :]).sum(axis=0).astype(int)
    return (n_ge + 1) / (n_iter + 1)


def color_for(p: float) -> str:
    if p < 0.001:
        return C_SIG_STRONG
    if p < 0.05:
        return C_SIG
    return C_NS


def panel(
    ax,
    label: str,
    observed: np.ndarray,
    p_emp: np.ndarray,
    n_iter: int,
    pids_in: list[str],
) -> None:
    # Sort patients by p-value (most significant at the top of the plot)
    order = np.argsort(p_emp)
    pids = np.array(pids_in)[order]
    obs = observed[order]
    p = p_emp[order]
    nlp = -np.log10(p)
    floor = -np.log10(1.0 / (n_iter + 1))

    y = np.arange(len(observed))
    colors = [color_for(pp) for pp in p]
    ax.barh(y, nlp, color=colors, edgecolor="white", linewidth=0.6, height=0.7)

    # Annotate with observed burden + p
    for yi, (oi, pi, npi) in enumerate(zip(obs, p, nlp)):
        # marker if at empirical floor
        suffix = "  (floor)" if pi <= 1.0 / (n_iter + 1) + 1e-12 else ""
        ax.text(
            npi + 0.05,
            yi,
            f"obs={int(oi)}, p={pi:.1e}{suffix}",
            va="center",
            fontsize=7.5,
            color="#222",
        )

    # Reference lines
    for thr, txt in [(0.05, "p = 0.05"), (0.001, "p = 0.001")]:
        x = -np.log10(thr)
        ax.axvline(x, color=C_REF, linestyle="--", linewidth=0.8, alpha=0.6)
        ax.text(
            x,
            len(observed) - 0.4,
            txt,
            fontsize=7,
            color=C_REF,
            ha="center",
            va="bottom",
        )
    # Empirical floor line
    ax.axvline(floor, color="#666", linestyle=":", linewidth=0.8, alpha=0.5)
    ax.text(
        floor,
        -1.0,
        f"floor (1/{n_iter:,}+1)",
        fontsize=6.5,
        color="#666",
        ha="center",
        va="top",
    )

    n_lt_05 = int((p_emp < 0.05).sum())
    n_lt_001 = int((p_emp < 0.001).sum())
    ax.set_yticks(y)
    ax.set_yticklabels(pids, fontsize=7.5)
    ax.set_xlabel("$-\\log_{10}(p_{\\mathrm{emp}})$", fontsize=10)
    ax.set_title(
        f"{label}   ({n_lt_05}/{len(observed)} at p<0.05;  "
        f"{n_lt_001}/{len(observed)} at p<0.001)",
        fontsize=10,
        loc="left",
        weight="bold",
        pad=8,
    )
    ax.set_xlim(0, max(floor * 1.15, 5))
    ax.invert_yaxis()  # most-significant patient at top
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def main() -> None:
    args = parse_args()

    # ── Load matrix and build the random pool ──────────────────────────────
    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)
    matrix = matrix.groupby(level=0).sum()  # collapse duplicate gene rows
    if matrix.shape[1] < N_MET + N_NONMET:
        raise SystemExit(
            f"Matrix has {matrix.shape[1]} patient columns; "
            f"need at least {N_MET + N_NONMET}."
        )

    ddr_genes = set(pd.read_csv(args.ddr_panel, sep="\t", comment="#")["gene"].tolist())
    flags = set() if args.no_exclude else load_exclude_list(args.exclude_genes)
    excluded = (set() if args.include_ddr_in_pool else ddr_genes) | flags
    pool = [g for g in matrix.index if g not in excluded]

    if len(pool) < PANEL_SIZE:
        raise SystemExit(
            f"Random pool has only {len(pool)} genes; need at least {PANEL_SIZE}."
        )

    config_label = "FLAGS-excluded" if flags else "no FLAGS exclusion"
    print(
        f"Bootstrap: {args.n_iterations:,} random {PANEL_SIZE}-gene panels, seed = {args.seed}"
    )
    print(f"Pool:      {len(pool)} genes ({config_label})")
    print(f"Matrix:    {matrix.shape[0]} genes × {matrix.shape[1]} patients")

    # ── Single bootstrap covering both arms (52 patients) ───────────────────
    null_matrix_full = build_null_matrix(matrix, pool, args.n_iterations, args.seed)

    # ── Per-arm observed vectors and per-patient empirical p-values ─────────
    met_obs = np.asarray(DEFAULT_MET, dtype=int)
    nm_obs = np.asarray(DEFAULT_NONMET, dtype=int)

    met_null = null_matrix_full[:, :N_MET]
    nm_null = null_matrix_full[:, N_MET : N_MET + N_NONMET]
    met_ids = [str(c) for c in matrix.columns[:N_MET]]
    nm_ids = [str(c) for c in matrix.columns[N_MET : N_MET + N_NONMET]]

    p_met = per_patient_p(met_obs, met_null)
    p_nm = per_patient_p(nm_obs, nm_null)

    print(f"\nMet patient IDs (from matrix header):     {met_ids}")
    print(f"Non-met patient IDs (from matrix header): {nm_ids}")
    print(
        f"\nMet:     {int((p_met < 0.05).sum())}/{N_MET} at p<0.05, "
        f"{int((p_met < 0.001).sum())}/{N_MET} at p<0.001"
    )
    print(
        f"Non-met: {int((p_nm < 0.05).sum())}/{N_NONMET} at p<0.05, "
        f"{int((p_nm < 0.001).sum())}/{N_NONMET} at p<0.001"
    )

    # ── Render figure ───────────────────────────────────────────────────────
    fig, (ax_m, ax_n) = plt.subplots(1, 2, figsize=(13, 8), sharex=False)
    panel(ax_m, "A   Metastatic (n = 26)", met_obs, p_met, args.n_iterations, met_ids)
    panel(ax_n, "B   Non-metastatic (n = 26)", nm_obs, p_nm, args.n_iterations, nm_ids)

    # Shared legend
    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor=C_SIG_STRONG, label="p < 0.001"),
        Patch(facecolor=C_SIG, label="0.001 ≤ p < 0.05"),
        Patch(facecolor=C_NS, label="p ≥ 0.05"),
    ]
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        ncols=3,
        frameon=False,
        fontsize=9,
        bbox_to_anchor=(0.5, -0.02),
    )

    fig.suptitle(
        f"Per-patient empirical p-values for DDR gnsRV burden vs.\n"
        f"{args.n_iterations:,} random 25-gene panels "
        f"({config_label}; Laplace-smoothed, one-sided)",
        fontsize=11,
        weight="bold",
        y=0.99,
    )
    plt.tight_layout(rect=(0, 0.04, 1, 0.96))
    pdf = FIGURES_DIR / "figS_per_patient_pvalues.pdf"
    png = FIGURES_DIR / "figS_per_patient_pvalues.png"
    plt.savefig(pdf, dpi=300, bbox_inches="tight")
    plt.savefig(png, dpi=300, bbox_inches="tight")
    print(f"\nSaved: {pdf}")
    print(f"Saved: {png}")


if __name__ == "__main__":
    main()
