#!/usr/bin/env python3
"""
fig2_genomewide_oncoplot.py
───────────────────────────
Single-panel genome-wide gnsRV oncoplot (Fig. 2) for the metastatic arm.
Designed for the 53-gene exploratory screen result, but works for any
gene × patient TSV in the standard flag format.

Input files (in data/raw/oncoplot/):
  met_genomewide.tsv       gnsRV data  (same flag format as met_ddr_gnsrv.tsv)
  cnv_genomewide.tsv       optional somatic CNV per gene
                           columns: GENE | MET | NONMET

Configuration — edit the block below before running:
  SHOW_PATIENT_IDS   True / False
  SHOW_HRD           True / False  (requires HRD row in TSV)
  SHOW_CNV           True / False
  GENE_SORT          "original" | "frequency" | "alphabetical"

Usage:
    python code/08_figures/fig2_genomewide_oncoplot.py

Output:
    figures/fig2_genomewide_oncoplot.{pdf,png}
"""

import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parents[1]
ONCOPLOT_DATA = PROJECT_DIR / "data" / "raw" / "oncoplot"
FIGURES_DIR = PROJECT_DIR / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

sys.path.insert(0, str(SCRIPT_DIR))
from _oncoplot_helpers import (  # shared parsers + palette
    EXTRA_MULTI,
    EXTRA_PATHOGENIC,
    HRD_THRESHOLD,
    PALETTE,
    TILE_H,
    ZYGOSITY_BOOST,
    parse_ddr_tex,
    parse_tsv,
)

# ─────────────────────────────────────────────────────────────────────────────
#  USER CONFIGURATION  ← edit here
# ─────────────────────────────────────────────────────────────────────────────

SHOW_PATIENT_IDS = True  # True  → patient IDs on x-axis
# False → plain axis with descriptive label

SHOW_HRD = True  # True  → HRD colour track above the matrix
#         (needs an HRD row in the TSV)

SHOW_CNV = True  # True  → somatic CNV side bar on the right

GENE_SORT = "original"  # "original"    → input row order (default)
# "frequency"   → most-altered gene at top
# "alphabetical"→ A–Z

# ─────────────────────────────────────────────────────────────────────────────
#  Drawing function
# ─────────────────────────────────────────────────────────────────────────────


def draw_genomewide(
    fig,
    data_path,
    cnv_path=None,
    arm_label="metastatic arm",
    bar_color=None,
    zygosity_dict=None,
    show_patient_ids=SHOW_PATIENT_IDS,
    show_hrd=SHOW_HRD,
    show_cnv=SHOW_CNV,
):
    """Draw the full oncoplot directly onto *fig* using its own GridSpec."""

    if bar_color is None:
        bar_color = PALETTE["bar_met"]
    if zygosity_dict is None:
        zygosity_dict = ZYGOSITY_BOOST

    # ── Parse ──────────────────────────────────────────────────────────────
    if str(data_path).endswith(".tex"):
        patients, hrd, genes_raw, matrix, cnv, totals = parse_ddr_tex(data_path)
    else:
        patients, hrd, genes_raw, matrix, cnv, totals = parse_tsv(data_path, cnv_path)
    n_pat = len(patients)

    for g, p in EXTRA_MULTI:
        if g in matrix and p in matrix[g]:
            matrix[g][p]["tags"].add("multi")

    if GENE_SORT == "frequency":
        genes = sorted(genes_raw, key=lambda g: len(matrix[g]), reverse=True)
    elif GENE_SORT == "alphabetical":
        genes = sorted(genes_raw)
    else:
        genes = list(genes_raw)
    n_genes = len(genes)
    gene_freq = {g: len(matrix[g]) for g in genes}

    # ── GridSpec ───────────────────────────────────────────────────────────
    n_cols = 3 if show_cnv else 2
    n_rows = 3 if show_hrd else 2  # top_bar [+ hrd] + main
    h_ratios = [0.8] + ([0.2] if show_hrd else []) + [n_genes * 0.30]
    w_ratios = [n_pat * 0.36, 1.1] + ([1.1] if show_cnv else [])

    gs = GridSpec(
        nrows=n_rows,
        ncols=n_cols,
        figure=fig,
        height_ratios=h_ratios,
        width_ratios=w_ratios,
        hspace=0.06,
        wspace=0.10,
        left=0.13,
        right=0.95,
        top=0.91,
        bottom=0.03,
    )

    r_top = 0
    r_hrd = 1 if show_hrd else None
    r_main = 2 if show_hrd else 1

    ax_top = fig.add_subplot(gs[r_top, 0])
    ax_hrd = fig.add_subplot(gs[r_hrd, 0]) if show_hrd else None
    ax_main = fig.add_subplot(gs[r_main, 0])
    ax_right = fig.add_subplot(gs[r_main, 1])
    ax_cnv = fig.add_subplot(gs[r_main, 2]) if show_cnv else None

    # ── Main matrix ────────────────────────────────────────────────────────
    cell_w = 0.90
    for j in range(n_pat):
        pat = patients[j]
        for i, gene in enumerate(genes):
            row_y = n_genes - 1 - i
            ax_main.add_patch(
                mpatches.FancyBboxPatch(
                    (j - cell_w / 2, row_y - 0.45),
                    cell_w,
                    0.90,
                    boxstyle="round,pad=0.02",
                    facecolor=PALETTE["bg_empty"],
                    edgecolor="white",
                    linewidth=0.6,
                )
            )
            if pat not in matrix[gene]:
                continue
            entry = matrix[gene][pat]
            tags = entry["tags"]
            color = PALETTE[entry["type"]]
            h = TILE_H[entry["type"]]
            ax_main.add_patch(
                mpatches.FancyBboxPatch(
                    (j - cell_w / 2 + 0.03, row_y - h / 2),
                    cell_w - 0.06,
                    h,
                    boxstyle="round,pad=0.01",
                    facecolor=color,
                    edgecolor="none",
                )
            )
            if "multi" in tags:
                ax_main.add_patch(
                    mpatches.FancyBboxPatch(
                        (j - cell_w / 2 + 0.03, row_y - h / 2),
                        cell_w - 0.06,
                        h,
                        boxstyle="round,pad=0.01",
                        facecolor="none",
                        edgecolor=PALETTE["multi_border"],
                        linewidth=2.0,
                    )
                )
            if (gene, pat) in zygosity_dict or "zygosity" in tags:
                ax_main.plot(
                    j,
                    row_y + h / 2 + 0.06,
                    marker="^",
                    color=PALETTE["zygosity_marker"],
                    markersize=5.5,
                    markeredgecolor="black",
                    markeredgewidth=0.6,
                    zorder=10,
                    clip_on=False,
                )

    ax_main.set_xlim(-0.5, n_pat - 0.5)
    ax_main.set_ylim(-0.5, n_genes - 0.5)
    ax_main.set_yticks(range(n_genes))
    ax_main.set_yticklabels(
        [genes[n_genes - 1 - i] for i in range(n_genes)],
        fontsize=9.5,
        fontstyle="italic",
    )
    ax_main.tick_params(axis="both", length=0)

    if show_patient_ids:
        ax_main.set_xticks(range(n_pat))
        ax_main.set_xticklabels(
            [str(p) for p in patients],
            fontsize=7.5,
            rotation=90,
        )
        ax_main.set_xlabel(f"Patient ID ({arm_label})", fontsize=10.5, labelpad=6)
    else:
        ax_main.set_xticks([])
        ax_main.set_xlabel(
            f"Patients ({arm_label}, n={n_pat})",
            fontsize=10.5,
            labelpad=6,
        )

    for sp in ax_main.spines.values():
        sp.set_visible(False)

    # ── Top bar — gnsRV count per patient ──────────────────────────────────
    top_vals = [totals.get(p, 0) for p in patients]
    ax_top.bar(
        range(n_pat),
        top_vals,
        color=PALETTE["bar_top"],
        width=0.72,
        edgecolor="white",
        linewidth=0.4,
    )
    ax_top.set_xlim(-0.5, n_pat - 0.5)
    ax_top.set_ylabel("gnsRV\ncount", fontsize=9.5)
    ax_top.set_xticks([])
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_top.spines["bottom"].set_visible(False)
    ax_top.tick_params(axis="y", labelsize=8.5)

    # ── HRD track (optional) ───────────────────────────────────────────────
    if ax_hrd is not None:
        hrd_vals = [hrd.get(p, 0) for p in patients]
        hrd_max = max(hrd_vals) if hrd_vals else 50
        hrd_cmap = mcolors.LinearSegmentedColormap.from_list(
            "hrd", ["#95a5a6", "#f2f3f4", "#1e8449"]
        )  # gray → white → green
        hrd_norm = mcolors.TwoSlopeNorm(
            vmin=0, vcenter=HRD_THRESHOLD, vmax=max(hrd_max, HRD_THRESHOLD + 1)
        )
        for j, (pat, val) in enumerate(zip(patients, hrd_vals)):
            fc = hrd_cmap(hrd_norm(val)) if hrd_max > 0 else "#f2f3f4"
            ax_hrd.add_patch(
                mpatches.Rectangle(
                    (j - 0.45, 0),
                    0.9,
                    1,
                    facecolor=fc,
                    edgecolor="white",
                    linewidth=0.5,
                )
            )
            ax_hrd.text(
                j,
                0.5,
                str(val),
                ha="center",
                va="center",
                fontsize=6.5,
                fontweight="bold",
                color="white" if val >= HRD_THRESHOLD else "#333333",
            )
        ax_hrd.set_xlim(-0.5, n_pat - 0.5)
        ax_hrd.set_ylim(0, 1)
        ax_hrd.set_yticks([0.5])
        ax_hrd.set_yticklabels(["HRD"], fontsize=9.5, fontweight="bold")
        ax_hrd.set_xticks([])
        ax_hrd.tick_params(axis="both", length=0)
        for sp in ax_hrd.spines.values():
            sp.set_visible(False)

    # ── Right bar — patients per gene ──────────────────────────────────────
    right_vals = [gene_freq[genes[n_genes - 1 - i]] for i in range(n_genes)]
    bars_r = ax_right.barh(
        range(n_genes),
        right_vals,
        color=bar_color,
        height=0.60,
        edgecolor="white",
        linewidth=0.4,
    )
    for bar, v in zip(bars_r, right_vals):
        if v > 0:
            ax_right.text(
                bar.get_width() + 0.12,
                bar.get_y() + bar.get_height() / 2,
                str(v),
                ha="left",
                va="center",
                fontsize=8,
            )
    ax_right.set_ylim(-0.5, n_genes - 0.5)
    ax_right.set_xlabel("# patients", fontsize=9.5)
    ax_right.set_yticks([])
    ax_right.spines["top"].set_visible(False)
    ax_right.spines["right"].set_visible(False)
    ax_right.tick_params(axis="x", labelsize=8.5)
    ax_right.set_title("EPC gnsRV\n variants", fontsize=9.5, pad=3)

    # ── CNV bar (optional) ─────────────────────────────────────────────────
    if ax_cnv is not None:
        y_pos = np.arange(n_genes)
        cnv_m = [cnv.get(genes[n_genes - 1 - i], (0, 0))[0] for i in range(n_genes)]
        cnv_n = [cnv.get(genes[n_genes - 1 - i], (0, 0))[1] for i in range(n_genes)]
        ax_cnv.barh(
            y_pos + 0.15,
            cnv_m,
            height=0.28,
            color=PALETTE["bar_met"],
            edgecolor="white",
            linewidth=0.3,
            label="Metastatic",
        )
        ax_cnv.barh(
            y_pos - 0.15,
            cnv_n,
            height=0.28,
            color=PALETTE["bar_nonmet"],
            edgecolor="white",
            linewidth=0.3,
            label="Non-metastatic",
        )
        ax_cnv.set_ylim(-0.5, n_genes - 0.5)
        ax_cnv.set_xlabel("# patients", fontsize=9.5)
        ax_cnv.set_yticks([])
        ax_cnv.spines["top"].set_visible(False)
        ax_cnv.spines["right"].set_visible(False)
        ax_cnv.tick_params(axis="x", labelsize=8.5)
        ax_cnv.legend(
            fontsize=10,  # 8.5,
            frameon=False,
            loc="lower right",
            bbox_to_anchor=(1.0, 1.05),
            bbox_transform=ax_cnv.transAxes,
        )
        ax_cnv.set_title("Somatic CNV\n (EPC)", fontsize=9.5, pad=3)


# ─────────────────────────────────────────────────────────────────────────────
#  Build figure
# ─────────────────────────────────────────────────────────────────────────────


def build_figure(
    data_path=None,
    cnv_path=None,
    out_dir=FIGURES_DIR,
):
    # ── Resolve paths ──────────────────────────────────────────────────────
    if data_path is None:
        data_path = ONCOPLOT_DATA / "met_genomewide.tsv"
    data_path = Path(data_path)
    if not data_path.exists():
        raise FileNotFoundError(f"Data file not found: {data_path}")

    if cnv_path is None:
        _cnv = ONCOPLOT_DATA / "cnv_genomewide.tsv"
        cnv_path = _cnv if _cnv.exists() else None

    # Pre-parse to get dimensions for figure sizing
    if str(data_path).endswith(".tex"):
        patients, _, genes_raw, _, _, _ = parse_ddr_tex(data_path)
    else:
        patients, _, genes_raw, _, _, _ = parse_tsv(data_path)
    n_genes = len(genes_raw)
    n_pat = len(patients)

    # ── Figure size ────────────────────────────────────────────────────────
    hrd_h = 0.2 if SHOW_HRD else 0
    fig_h = max(10, (n_genes * 0.30 + 0.8 + hrd_h) * 1.18 + 2.0)
    fig_w = 14
    fig = plt.figure(figsize=(fig_w, fig_h))

    fig.suptitle(
        f"Genome-wide screen gnsRVs — metastatic arm  "
        f"({n_pat} patients, {n_genes} genes)",
        fontsize=13,
        fontweight="bold",
        y=0.975,
    )

    # ── Draw ───────────────────────────────────────────────────────────────
    draw_genomewide(
        fig,
        data_path,
        cnv_path=cnv_path,
        arm_label="metastatic arm",
        bar_color=PALETTE["bar_met"],
        zygosity_dict=ZYGOSITY_BOOST,
        show_patient_ids=SHOW_PATIENT_IDS,
        show_hrd=SHOW_HRD,
        show_cnv=SHOW_CNV,
    )

    # ── Legend (top-right of figure) ───────────────────────────────────────
    legend_handles = [
        mpatches.Patch(
            facecolor=PALETTE["deleterious"], label="Deleterious (fathmm-MKL)"
        ),
        mpatches.Patch(facecolor=PALETTE["tolerated"], label="Tolerated (fathmm-MKL)"),
        mpatches.Patch(
            facecolor=PALETTE["bg_empty"],
            edgecolor=PALETTE["multi_border"],
            linewidth=2.0,
            label="Multiple gnsRVs in gene",
        ),
        plt.Line2D(
            [0],
            [0],
            marker="^",
            color="w",
            markerfacecolor=PALETTE["zygosity_marker"],
            markeredgecolor="black",
            markersize=9,
            label="Zygosity / dosage boost",
        ),
        mpatches.Patch(
            facecolor="#1e8449", label=f"HRD ≥ {HRD_THRESHOLD} (DDR deficient)"
        ),
        mpatches.Patch(facecolor="#95a5a6", label=f"HRD < {HRD_THRESHOLD}"),
        mpatches.Patch(
            facecolor=PALETTE["bar_met"], alpha=0.9, label="Metastatic (n = 26)"
        ),
        mpatches.Patch(
            facecolor=PALETTE["bar_nonmet"],
            alpha=0.9,
            label="Non-metastatic (n =26)",
        ),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper right",
        bbox_to_anchor=(0.975, 0.965),
        ncol=4,  # 2,
        fontsize=14,  # 10,
        frameon=True,
        fancybox=True,
        edgecolor="#cccccc",
        title="Legend",
        title_fontsize=11,
        handlelength=1.6,
        columnspacing=1.2,
    )

    # ── Save ───────────────────────────────────────────────────────────────
    out = Path(out_dir)
    for ext in ("pdf", "png"):
        p = out / f"fig2_genomewide_oncoplot.{ext}"
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(f"Saved: {p}")
    plt.close(fig)
    print(
        f"\nConfiguration used:  SHOW_PATIENT_IDS={SHOW_PATIENT_IDS}  "
        f"SHOW_HRD={SHOW_HRD}  SHOW_CNV={SHOW_CNV}  GENE_SORT={GENE_SORT!r}"
    )


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    build_figure()
