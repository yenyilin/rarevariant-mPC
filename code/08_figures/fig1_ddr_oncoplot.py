#!/usr/bin/env python3
"""
fig1_ddr_oncoplot.py
────────────────────
Five-panel DDR figure (Fig. 1) comparing germline nonsynonymous (gnsRV) and
synonymous rare-variant landscapes side by side as a negative control.

  Panel (A)  Bar chart : gnsRV vs synonymous-variant enrichment
  Panel (B)  Oncoplot  : gnsRV       — metastatic arm     (top-left)
  Panel (C)  Oncoplot  : gnsRV       — non-metastatic arm (bottom-left)
  Panel (D)  Oncoplot  : synonymous  — metastatic arm     (top-right, negative control of B)
  Panel (E)  Oncoplot  : synonymous  — non-metastatic arm (bottom-right, negative control of C)

Input files (in data/raw/oncoplot/):
  met_ddr_gnsrv.tsv             — nonsynonymous, metastatic arm
  nonmet_ddr_gnsrv.tsv          — nonsynonymous, non-metastatic arm
  met_ddr_synonymous_rv.tsv     — synonymous, metastatic arm (negative control)
  nonmet_ddr_synonymous_rv.tsv  — synonymous, non-metastatic arm (negative control)
  cnv_ddr.tsv                   — gene-level CNV counts (shown in B/C only)

Usage:
    python code/08_figures/fig1_ddr_oncoplot.py

Output:
    figures/fig1_ddr_oncoplot.{pdf,png}
"""

import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parents[1]
ONCOPLOT_DATA = PROJECT_DIR / "data" / "raw" / "oncoplot"
FIGURES_DIR = PROJECT_DIR / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

# ── Import shared constants and parsing functions ─────────────────────────────
sys.path.insert(0, str(SCRIPT_DIR))
from _oncoplot_helpers import (  # noqa: E402
    EXTRA_MULTI,
    EXTRA_PATHOGENIC,
    GENE_SORT,
    HRD_THRESHOLD,
    PALETTE,
    TILE_H,
    ZYGOSITY_BOOST,
    draw_panel_a,
    parse_ddr_tex,
    parse_tsv,
)

# ═════════════════════════════════════════════════════════════════════════════
#  Extended oncoplot panel — draw_panel_ext()
#  Extends the base oncoplot panel drawing with:
#    show_hrd   : bool  — draw the HRD colour track (hide for synonymous panels)
#    show_cnv   : bool  — draw the CNV side bar     (hide for synonymous panels)
#    top_label  : str   — y-axis label on the top count bar
#    right_title: str   — title above the patients-per-gene bar
# ═════════════════════════════════════════════════════════════════════════════


def draw_panel_ext(
    fig,
    gs_slot,
    data_path,
    cnv_path=None,
    arm_label="metastatic arm",
    bar_color=None,
    zygosity_dict=None,
    show_hrd=True,
    show_cnv=True,
    show_patient=True,
    top_label="DDR gnsRV\ncount",
    right_title="EPC variants",
):
    if bar_color is None:
        bar_color = PALETTE["bar_met"]
    if zygosity_dict is None:
        zygosity_dict = ZYGOSITY_BOOST

    # ── Parse data ────────────────────────────────────────────────────────────
    if str(data_path).endswith(".tex"):
        patients, hrd, genes_raw, matrix, cnv, totals = parse_ddr_tex(data_path)
    else:
        patients, hrd, genes_raw, matrix, cnv, totals = parse_tsv(data_path, cnv_path)
    n_pat = len(patients)

    for g, p in EXTRA_MULTI:
        if g in matrix and p in matrix[g]:
            matrix[g][p]["tags"].add("multi")
    for g, p in EXTRA_PATHOGENIC:
        if g in matrix and p in matrix[g]:
            matrix[g][p]["tags"].add("pathogenic")

    if GENE_SORT == "frequency":
        genes = sorted(genes_raw, key=lambda g: len(matrix[g]), reverse=True)
    elif GENE_SORT == "alphabetical":
        genes = sorted(genes_raw)
    else:
        genes = list(genes_raw)
    n_genes = len(genes)
    gene_freq = {g: len(matrix[g]) for g in genes}

    # ── Inner GridSpec ────────────────────────────────────────────────────────
    n_rows = 3 if show_hrd else 2
    n_cols = 3 if show_cnv else 2
    h_ratios = [0.8, 0.2, n_genes * 0.36] if show_hrd else [0.8, n_genes * 0.36]
    w_ratios = [n_pat * 0.36, 1.1, 1.1] if show_cnv else [n_pat * 0.36, 1.1]

    inner = GridSpecFromSubplotSpec(
        nrows=n_rows,
        ncols=n_cols,
        subplot_spec=gs_slot,
        height_ratios=h_ratios,
        width_ratios=w_ratios,
        hspace=0.06,
        wspace=0.10,
    )

    r_top = 0
    r_hrd = 1 if show_hrd else None
    r_main = 2 if show_hrd else 1

    ax_top = fig.add_subplot(inner[r_top, 0])
    ax_hrd = fig.add_subplot(inner[r_hrd, 0]) if show_hrd else None
    ax_main = fig.add_subplot(inner[r_main, 0])
    ax_right = fig.add_subplot(inner[r_main, 1])
    ax_cnv = fig.add_subplot(inner[r_main, 2]) if show_cnv else None

    # ── Main matrix ───────────────────────────────────────────────────────────
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
    if show_patient:
        ax_main.set_xticks(range(n_pat))
        ax_main.set_xticklabels(patients, rotation=90, fontsize=7.5, ha="center")
        ax_main.tick_params(axis="x", length=0)
        ax_main.tick_params(axis="y", length=0)
    else:
        ax_main.set_xticks([])
        ax_main.tick_params(axis="both", length=0)
    ax_main.set_xlabel(
        f"Patients ({arm_label}, n={n_pat}; sorted by decreasing HRD score)",
        fontsize=10.5,
        labelpad=6,
    )
    for sp in ax_main.spines.values():
        sp.set_visible(False)

    # ── Top bar (variant count per patient) ───────────────────────────────────
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
    ax_top.set_ylabel(top_label, fontsize=9.5)
    ax_top.set_xticks([])
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_top.spines["bottom"].set_visible(False)
    ax_top.tick_params(axis="y", labelsize=8.5)

    # ── HRD track (optional) ─────────────────────────────────────────────────
    if ax_hrd is not None:
        hrd_vals = [hrd.get(p, 0) for p in patients]
        hrd_max = max(hrd_vals) if hrd_vals else 50
        hrd_cmap = mcolors.LinearSegmentedColormap.from_list(
            "hrd", ["#95a5a6", "#f2f3f4", "#1e8449"]
        )  # gray → near-white → dark green
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

    # ── Right bar (patients per gene) ─────────────────────────────────────────
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
    ax_right.set_title(right_title, fontsize=9.5, pad=3)

    # ── CNV bar (optional) ────────────────────────────────────────────────────
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
        # ax_cnv.legend(fontsize=8.5, frameon=False, loc="upper right")
        ax_cnv.legend(
            fontsize=10,  # 8.5,
            frameon=False,
            loc="lower right",  # anchor point on the legend box itself
            bbox_to_anchor=(1.0, 1.05),  # (x, y) in axes coords — y > 1 = above axes
            bbox_transform=ax_cnv.transAxes,
        )
        ax_cnv.set_title("Somatic CNV\n (EPC)", fontsize=9.5, pad=3)


# ═════════════════════════════════════════════════════════════════════════════
#  Build the five-panel figure
# ═════════════════════════════════════════════════════════════════════════════


def build_figure(
    gnsrv_met_path=None,
    gnsrv_nonmet_path=None,
    gsrv_met_path=None,
    gsrv_nonmet_path=None,
    cnv_path=None,
    out_dir=FIGURES_DIR,
):
    # ── Resolve data paths ────────────────────────────────────────────────────
    def _resolve(path, default_name):
        if path is None:
            p = ONCOPLOT_DATA / default_name
            return p if p.exists() else None
        return Path(path)

    gnsrv_met_path = _resolve(gnsrv_met_path, "met_ddr_gnsrv.tsv")
    gnsrv_nonmet_path = _resolve(gnsrv_nonmet_path, "nonmet_ddr_gnsrv.tsv")
    # NOTE: filenames renamed from `met_ddr_gsrv.tsv` / `nonmet_ddr_gsrv.tsv`
    # to `met_ddr_synonymous_rv.tsv` / `nonmet_ddr_synonymous_rv.tsv` to
    # avoid one-letter visual confusion with the nonsynonymous gnsRV files
    # in `ls` listings. Internal `gsrv_*` variable names are kept for code
    # consistency (field-standard abbreviation paired with gnsRV).
    gsrv_met_path = _resolve(gsrv_met_path, "met_ddr_synonymous_rv.tsv")
    gsrv_nonmet_path = _resolve(gsrv_nonmet_path, "nonmet_ddr_synonymous_rv.tsv")

    if cnv_path is None:
        _cnv = ONCOPLOT_DATA / "cnv_ddr.tsv"
        cnv_path = _cnv if _cnv.exists() else None

    if gnsrv_met_path is None:
        # Fall back to legacy .tex if no TSV found
        gnsrv_met_path = ONCOPLOT_DATA / "ddr.tex"

    _has_nonmet = gnsrv_nonmet_path is not None and gnsrv_nonmet_path.exists()
    _has_gsrv_met = gsrv_met_path is not None and gsrv_met_path.exists()
    _has_gsrv_nonmet = gsrv_nonmet_path is not None and gsrv_nonmet_path.exists()

    # ── Figure and outer GridSpec ─────────────────────────────────────────────
    # 3 rows × 2 cols; row 0 spans both cols (Panel A)
    n_rows = 3 if _has_nonmet else 2
    h_ratios = [0.55, 2.6, 2.6][:n_rows]
    fig_w = 26 if _has_gsrv_met else 14
    fig_h = 20 if _has_nonmet else 13

    fig = plt.figure(figsize=(fig_w, fig_h))
    outer = GridSpec(
        nrows=n_rows,
        ncols=2,
        height_ratios=h_ratios,
        width_ratios=[1, 1],
        hspace=0.2,  # hspace=0.14,
        wspace=0.06,
        figure=fig,
    )

    # ── Panel A: bar chart + legend (spans both columns) ─────────────────────
    top_gs = GridSpecFromSubplotSpec(
        1,
        2,
        subplot_spec=outer[0, :],
        width_ratios=[1, 1.6],
        wspace=0.25,
    )
    ax_bar = fig.add_subplot(top_gs[0])
    ax_leg = fig.add_subplot(top_gs[1])

    draw_panel_a(ax_bar)
    ax_bar.text(
        -0.02,  # 0.0,  # -0.15,
        1.2,  # 08,
        "A",
        transform=ax_bar.transAxes,
        fontsize=20,
        fontweight="bold",
        va="top",
    )

    ax_leg.axis("off")
    legend_handles = [
        mpatches.Patch(
            facecolor=PALETTE["deleterious"], label="Deleterious (fathmm-MKL)"
        ),
        mpatches.Patch(facecolor=PALETTE["tolerated"], label="Tolerated (fathmm-MKL)"),
        mpatches.Patch(
            facecolor=PALETTE["bg_empty"],
            edgecolor=PALETTE["multi_border"],
            linewidth=2.0,
            label="Multiple variants in gene",
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
            facecolor=PALETTE["bar_nonmet"], alpha=0.9, label="Non-metastatic (n = 26)"
        ),
    ]
    ax_leg.legend(
        handles=legend_handles,
        loc="center",
        ncol=2,
        fontsize=20,  # 10.5,
        frameon=True,
        fancybox=True,
        edgecolor="#cccccc",
        title="Oncoplot legend",
        title_fontsize=11.5,
        handlelength=1.8,
        handleheight=1.2,
        columnspacing=1.5,
    )

    # ── Panel B: gnsRV metastatic (top-left) ─────────────────────────────────
    _label_panel(fig, outer[1, 0], "B")
    draw_panel_ext(
        fig,
        outer[1, 0],
        gnsrv_met_path,
        cnv_path=cnv_path,
        arm_label="metastatic arm",
        bar_color=PALETTE["bar_met"],
        zygosity_dict=ZYGOSITY_BOOST,
        show_hrd=True,
        show_cnv=True,
        top_label="DDR gnsRV\ncount",
        right_title="EPC gnsRV\n variants",
    )

    # ── Panel D: gsRV metastatic (top-right, negative control) ───────────────
    if _has_gsrv_met:
        _label_panel(fig, outer[1, 1], "D")
        draw_panel_ext(
            fig,
            outer[1, 1],
            gsrv_met_path,
            cnv_path=None,
            arm_label="metastatic arm",
            bar_color=PALETTE["bar_met"],
            zygosity_dict={},
            show_hrd=True,
            show_cnv=False,
            top_label="DDR gsRV\ncount",
            right_title="EPC gsRV\n variants",
        )

    # ── Panel C: gnsRV non-metastatic (bottom-left) ───────────────────────────
    if _has_nonmet:
        _label_panel(fig, outer[2, 0], "C")
        draw_panel_ext(
            fig,
            outer[2, 0],
            gnsrv_nonmet_path,
            cnv_path=cnv_path,
            arm_label="non-metastatic arm",
            bar_color=PALETTE["bar_nonmet"],
            zygosity_dict={},
            show_hrd=True,
            show_cnv=True,
            top_label="DDR gnsRV\ncount",
            right_title="EPC gnsRV\n variants",
        )

    # ── Panel E: gsRV non-metastatic (bottom-right, negative control) ─────────
    if _has_nonmet and _has_gsrv_nonmet:
        _label_panel(fig, outer[2, 1], "E")
        draw_panel_ext(
            fig,
            outer[2, 1],
            gsrv_nonmet_path,
            cnv_path=None,
            arm_label="non-metastatic arm",
            bar_color=PALETTE["bar_nonmet"],
            zygosity_dict={},
            show_hrd=True,
            show_cnv=False,
            top_label="DDR gsRV\ncount",
            right_title="EPC gsRV\n variants",
        )

    # ── Column headers ────────────────────────────────────────────────────────
    # Derive x-centres and y-top of the oncoplot row from the GridSpec geometry.
    # get_position() returns a Bbox in normalised figure coordinates (0–1).
    if _has_gsrv_met:
        left_bb = outer[1, 0].get_position(fig)
        right_bb = outer[1, 1].get_position(fig)
        x_left = (left_bb.x0 + left_bb.x1) / 2
        x_right = (right_bb.x0 + right_bb.x1) / 2
        y_hdr = left_bb.y1 + 0.012  # just above the top of the oncoplot rows

        fig.text(
            x_left - 0.05,
            y_hdr,
            "Nonsynonymous variants (gnsRV)",
            ha="center",
            va="bottom",
            fontsize=13.5,
            fontweight="bold",
            color="#2c3e50",
        )

        fig.text(
            x_right,
            y_hdr,
            "Synonymous variants (gsRV)  —  negative control",
            ha="center",
            va="bottom",
            fontsize=13.5,
            fontweight="bold",
            color="#7f8c8d",
        )  # muted: secondary/control panels

    # ── Save ──────────────────────────────────────────────────────────────────
    out = Path(out_dir)
    for ext in ("pdf", "png"):
        p = out / f"fig1_ddr_oncoplot.{ext}"
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(f"Saved: {p}")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
def _label_panel(fig, gs_slot, letter):
    """Add a bold panel letter (A, B, …) to the top-left of a GridSpec slot."""
    ax = fig.add_subplot(gs_slot)
    ax.axis("off")
    ax.text(
        -0.02,
        1.04,  # 1.02,
        letter,
        transform=ax.transAxes,
        fontsize=20,  # 16,
        fontweight="bold",
        va="top",
    )


# ═════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    build_figure()
