#!/usr/bin/env python3
"""
Shared helpers for the DDR / genome-wide oncoplot figures.

Parsing functions, colour palette, annotation constants, and panel-drawing
routines imported by `fig1_ddr_oncoplot.py` and `fig2_genomewide_oncoplot.py`.

This is a helper module, not a figure script — it is not meant to be run
directly (the `__main__` block at the end is a standalone demo only).

Annotation layers provided to the importing scripts:
  Pathogenic     — ClinVar pathogenic / likely pathogenic  → crimson tile
  Multi-variant  — ≥2 gnsRVs in same gene×patient          → black border
  Zygosity boost — dosage-dependent effect observed        → gold ▲ marker
"""

import re

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


# ═══════════════════════════════════════════════════════════════════════
#  SECTION 1 — Parse data from ddr.tex
# ═══════════════════════════════════════════════════════════════════════


def parse_ddr_tex(filepath):
    """Return (patients, hrd, genes_order, matrix, cnv, totals).

    matrix[gene][patient] = {'type': 'deleterious'|'tolerated',
                             'tags': set of 'pathogenic','multi'}
    cnv[gene] = (M_count, N_count)
    """
    text = Path(filepath).read_text()
    lines = text.strip().splitlines()

    patients = []
    hrd_scores = {}
    matrix = {}
    cnv = {}
    totals = {}
    genes_order = []

    for line in lines:
        line = line.strip()
        if (
            not line
            or "\\begin{" in line
            or line == "\\hline"
            or "\\multicolumn" in line
        ):
            continue

        row = re.sub(r"\\\\.*$", "", line)
        cells = [c.strip() for c in row.split("&")]
        if len(cells) < 3:
            continue

        label = cells[0].strip()

        if label == "Patient":
            for c in cells[1:]:
                c = c.strip()
                if c and c not in ("M", "N"):
                    patients.append(int(c))
            continue

        n_pat = len(patients)

        if label == "HRD":
            for i, cell in enumerate(cells[1 : n_pat + 1]):
                m = re.search(r"(\d+)", cell)
                if m:
                    hrd_scores[patients[i]] = int(m.group(1))
            continue

        if label == "Total":
            for i, cell in enumerate(cells[1 : n_pat + 1]):
                cell = cell.strip()
                if cell.isdigit():
                    totals[patients[i]] = int(cell)
            continue

        gene = label
        genes_order.append(gene)
        matrix[gene] = {}

        for i, cell in enumerate(cells[1 : n_pat + 1]):
            cell = cell.strip()
            if not cell:
                continue
            tags = set()
            if "\\cellcolor{yellow}" in cell:
                tags.add("multi")
            if "\\color{red}" in cell:
                tags.add("pathogenic")
            if "\\ding{51}" in cell:
                vtype = "deleterious"
            elif "\\ding{55}" in cell:
                vtype = "tolerated"
            else:
                continue
            matrix[gene][patients[i]] = {"type": vtype, "tags": tags}

        m_val = cells[n_pat + 1].strip() if len(cells) > n_pat + 1 else ""
        n_val = cells[n_pat + 2].strip() if len(cells) > n_pat + 2 else ""
        cnv[gene] = (
            int(m_val) if m_val.isdigit() else 0,
            int(n_val) if n_val.isdigit() else 0,
        )

    return patients, hrd_scores, genes_order, matrix, cnv, totals


# ═══════════════════════════════════════════════════════════════════════
#  SECTION 2 — User-configurable annotations
# ═══════════════════════════════════════════════════════════════════════
#
#  The parser auto-detects from LaTeX markup:
#    • 'multi'      ← \cellcolor{yellow}  (≥2 gnsRVs in gene for patient)
#    • 'pathogenic'  ← \color{red}         (ClinVar pathogenic/LP)
#
#  Edit the dicts below to ADD or OVERRIDE annotations.
#  Keys are (gene_name, patient_id) tuples.
# ───────────────────────────────────────────────────────────────────────

# Additional multi-variant entries.  Value = number of gnsRVs (for label).
EXTRA_MULTI = {
    # ("BRCA2", 66): 2,
}

# Additional pathogenic entries.
EXTRA_PATHOGENIC = {
    # ("CHEK2", 16): True,
}

# Zygosity / dosage-boost entries.  Gold ▲ marker on tile.
ZYGOSITY_BOOST = {
    # ("KDM6B", 9): True,    # not in DDR table; example only
}

# ── Gene sort order ────────────────────────────────────────────────────
# "frequency" = most-affected gene at top (standard oncoplot)
# "alphabetical" = A–Z (matches the original LaTeX)
# "original" = keep the row order from ddr.tex
GENE_SORT = "original"


# ═══════════════════════════════════════════════════════════════════════
#  SECTION 3 — Panel (A) bar chart data
# ═══════════════════════════════════════════════════════════════════════
#  From manuscript results: arm-exclusive variant counts in 25 DDR genes.
#  Update these if the underlying numbers change.

BAR_DATA = {
    "gnsRV_met": 81,  # nonsynonymous, metastatic arm
    "gnsRV_nonmet": 21,  # nonsynonymous, non-metastatic arm
    "gsRV_met": 15,  # synonymous, metastatic arm
    "gsRV_nonmet": 17,  # synonymous, non-metastatic arm
    "p_gnsRV": "4.57×10⁻⁶",
    "p_gsRV": "n.s.",
}


# ═══════════════════════════════════════════════════════════════════════
#  SECTION 4 — Color palette & layout constants
# ═══════════════════════════════════════════════════════════════════════

PALETTE = {
    "deleterious": "#2d6a4f",  # forest green
    "tolerated": "#8aadbd",  # steel blue
    "pathogenic": "#c1272d",  # crimson
    "multi_border": "#1a1a1a",  # black border
    "zygosity_marker": "#f1c40f",  # gold triangle
    "bg_empty": "#eaeaea",  # cell background (no variant)
    "hrd_high": "#27ae60",  # HRD ≥ threshold
    "hrd_low": "#f0f0f0",  # HRD < threshold
    "bar_met": "#b03a2e",  # metastatic
    "bar_nonmet": "#2e86c1",  # non-metastatic
    "bar_top": "#5b2c6f",  # top bar (total gnsRV count)
}

HRD_THRESHOLD = 21

# Tile height as fraction of cell (1.0 = full cell)
TILE_H = {
    "deleterious": 0.72,
    "tolerated": 0.42,
    "pathogenic": 0.88,
}


# ═══════════════════════════════════════════════════════════════════════
#  SECTION 5 — Draw panel (A): bar chart
# ═══════════════════════════════════════════════════════════════════════


def draw_panel_a(ax):
    """gnsRV vs gsRV enrichment bar chart."""
    met_counts = [BAR_DATA["gnsRV_met"], BAR_DATA["gsRV_met"]]
    nonmet_counts = [BAR_DATA["gnsRV_nonmet"], BAR_DATA["gsRV_nonmet"]]
    labels = ["Nonsynonymous\n(gnsRVs)", "Synonymous\n(gsRVs)"]

    x = np.array([0, 1])
    width = 0.30

    bars_m = ax.bar(
        x - width / 2,
        met_counts,
        width,
        color=PALETTE["bar_met"],
        alpha=0.90,
        edgecolor="black",
        linewidth=0.7,
        label="Metastatic (n=26)",
    )
    bars_n = ax.bar(
        x + width / 2,
        nonmet_counts,
        width,
        color=PALETTE["bar_nonmet"],
        alpha=0.90,
        edgecolor="black",
        linewidth=0.7,
        label="Non-metastatic (n=26)",
    )

    # Value labels
    for bar in list(bars_m) + list(bars_n):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 1.0,
            str(int(bar.get_height())),
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
        )

    # Significance bracket — nonsynonymous
    sig_y = 77
    ax.plot(
        [x[0] - width / 2, x[0] + width / 2], [sig_y + 12, sig_y + 12], "k-", lw=1.2
    )
    ax.text(
        x[0],
        sig_y + 15,  # 1.5,
        f"p = {BAR_DATA['p_gnsRV']}",
        ha="center",
        va="bottom",
        fontsize=9,
        style="italic",
    )

    # n.s. — synonymous
    ax.text(
        x[1],
        max(met_counts[1], nonmet_counts[1]) + 2.5,
        BAR_DATA["p_gsRV"],
        ha="center",
        va="bottom",
        fontsize=9,
        color="#555555",
    )

    ax.set_ylabel("Arm-exclusive variants\nin 25 DDR genes", fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylim(0, 92)
    ax.set_xlim(-0.55, 1.55)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", length=3)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")


# ═══════════════════════════════════════════════════════════════════════
#  SECTION 5b — TSV parser (alternative to ddr.tex)
# ═══════════════════════════════════════════════════════════════════════
#
#  Expected TSV layout (one file per arm):
#    Row index (col 0) : gene name  |  special rows: HRD, TOTAL
#    Column headers    : patient IDs (integers)
#
#  Cell codes (pipe-separated flags):
#    D = deleterious (fathmm-MKL)   T = tolerated    P = pathogenic (ClinVar)
#    M = multiple gnsRVs            Z = zygosity / dosage boost
#    . or empty = no variant
#  Combinations: D|M   D|Z   D|M|Z   T|M   P|M   etc.
#
#  Optional separate CNV file (gene-level, shared across arms):
#    Columns: GENE  MET  NONMET

import pandas as _pd  # local alias to avoid polluting namespace


def parse_tsv(filepath, cnv_path=None):
    """Return (patients, hrd, genes_order, matrix, cnv, totals) — same
    tuple as parse_ddr_tex() so draw_panel_b() accepts both formats."""
    df = _pd.read_csv(filepath, sep="\t", index_col=0, dtype=str).fillna(".")
    df.columns = [c.strip() for c in df.columns]

    patients = []
    for c in df.columns:
        try:
            patients.append(int(c))
        except ValueError:
            pass

    hrd_scores = {}
    totals = {}
    matrix = {}
    genes_order = []

    for gene, row in df.iterrows():
        gene = str(gene).strip()

        if gene.upper() == "HRD":
            for pat in patients:
                val = row.get(str(pat), ".").strip()
                if val not in (".", ""):
                    try:
                        hrd_scores[pat] = int(val)
                    except ValueError:
                        pass
            continue

        if gene.upper() == "TOTAL":
            for pat in patients:
                val = row.get(str(pat), ".").strip()
                if val not in (".", ""):
                    try:
                        totals[pat] = int(val)
                    except ValueError:
                        pass
            continue

        genes_order.append(gene)
        matrix[gene] = {}

        for pat in patients:
            cell = row.get(str(pat), ".").strip().upper()
            if not cell or cell == ".":
                continue

            flags = set(cell.split("|"))

            if "T" in flags:
                vtype = "tolerated"
            elif "D" in flags or "P" in flags:  # P alone = pathogenic+deleterious
                vtype = "deleterious"
            else:
                continue  # unrecognised — skip

            tags = set()
            if "M" in flags:
                tags.add("multi")
            if "Z" in flags:
                tags.add("zygosity")
            if "P" in flags:
                tags.add("pathogenic")

            matrix[gene][pat] = {"type": vtype, "tags": tags}

    # CNV (optional)
    cnv = {}
    if cnv_path and Path(cnv_path).exists():
        cnv_df = _pd.read_csv(cnv_path, sep="\t", index_col=0, dtype=str).fillna("0")
        for g in genes_order:
            if g in cnv_df.index:
                m = int(cnv_df.loc[g, "MET"]) if "MET" in cnv_df.columns else 0
                n = int(cnv_df.loc[g, "NONMET"]) if "NONMET" in cnv_df.columns else 0
                cnv[g] = (m, n)
            else:
                cnv[g] = (0, 0)

    return patients, hrd_scores, genes_order, matrix, cnv, totals


# ═══════════════════════════════════════════════════════════════════════
#  SECTION 6 — Draw panel (B / C): oncoplot
# ═══════════════════════════════════════════════════════════════════════


def draw_panel_b(
    fig,
    gs_slot,
    data_path,
    cnv_path=None,
    arm_label="metastatic arm",
    bar_color=None,
    zygosity_dict=None,
):
    """Build the oncoplot inside the given GridSpec slot.

    data_path    : path to .tex (legacy) or .tsv (new format)
    cnv_path     : optional path to cnv_ddr.tsv (gene-level CNV counts)
    arm_label    : x-axis label suffix, e.g. 'metastatic arm'
    bar_color    : colour for the right-bar (patients per gene)
    zygosity_dict: {(gene, patient_id): True} for .tex files; ignored for .tsv
                   (TSV encodes zygosity via the Z flag in the cell)
    """
    if bar_color is None:
        bar_color = PALETTE["bar_met"]
    if zygosity_dict is None:
        zygosity_dict = ZYGOSITY_BOOST

    # Auto-detect format
    if str(data_path).endswith(".tex"):
        patients, hrd, genes_raw, matrix, cnv, totals = parse_ddr_tex(data_path)
    else:
        patients, hrd, genes_raw, matrix, cnv, totals = parse_tsv(data_path, cnv_path)
    n_pat = len(patients)

    # Apply user annotations
    for g, p in EXTRA_MULTI:
        if g in matrix and p in matrix[g]:
            matrix[g][p]["tags"].add("multi")
    for g, p in EXTRA_PATHOGENIC:
        if g in matrix and p in matrix[g]:
            matrix[g][p]["tags"].add("pathogenic")

    # Sort genes
    if GENE_SORT == "frequency":
        genes = sorted(genes_raw, key=lambda g: len(matrix[g]), reverse=True)
    elif GENE_SORT == "alphabetical":
        genes = sorted(genes_raw)
    else:
        genes = list(genes_raw)
    n_genes = len(genes)

    gene_freq = {g: len(matrix[g]) for g in genes}

    # ── Sub-gridspec inside the panel B slot ──────────────────────
    inner = GridSpecFromSubplotSpec(
        nrows=3,
        ncols=3,
        subplot_spec=gs_slot,
        height_ratios=[0.8, 0.2, n_genes * 0.36],
        width_ratios=[n_pat * 0.36, 1.8, 1.8],
        hspace=0.06,
        wspace=0.10,
    )

    ax_top = fig.add_subplot(inner[0, 0])
    ax_hrd = fig.add_subplot(inner[1, 0])
    ax_main = fig.add_subplot(inner[2, 0])
    ax_right = fig.add_subplot(inner[2, 1])
    ax_cnv = fig.add_subplot(inner[2, 2])

    # ── Main matrix ───────────────────────────────────────────────
    cell_w = 0.90

    for j in range(n_pat):
        pat = patients[j]
        for i, gene in enumerate(genes):
            row_y = n_genes - 1 - i

            bg = mpatches.FancyBboxPatch(
                (j - cell_w / 2, row_y - 0.45),
                cell_w,
                0.90,
                boxstyle="round,pad=0.02",
                facecolor=PALETTE["bg_empty"],
                edgecolor="white",
                linewidth=0.6,
            )
            ax_main.add_patch(bg)

            if pat not in matrix[gene]:
                continue

            entry = matrix[gene][pat]
            tags = entry["tags"]

            color = PALETTE[entry["type"]]
            h = TILE_H[entry["type"]]

            tile = mpatches.FancyBboxPatch(
                (j - cell_w / 2 + 0.03, row_y - h / 2),
                cell_w - 0.06,
                h,
                boxstyle="round,pad=0.01",
                facecolor=color,
                edgecolor="none",
            )
            ax_main.add_patch(tile)

            if "multi" in tags:
                border = mpatches.FancyBboxPatch(
                    (j - cell_w / 2 + 0.03, row_y - h / 2),
                    cell_w - 0.06,
                    h,
                    boxstyle="round,pad=0.01",
                    facecolor="none",
                    edgecolor=PALETTE["multi_border"],
                    linewidth=2.0,
                )
                ax_main.add_patch(border)

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
        fontsize=7.5,
        fontstyle="italic",
    )
    ax_main.set_xticks([])
    ax_main.tick_params(axis="both", length=0)
    ax_main.set_xlabel(
        f"Patients ({arm_label}, n={n_pat}; sorted by decreasing HRD score)",
        fontsize=8.5,
        labelpad=6,
    )
    for sp in ax_main.spines.values():
        sp.set_visible(False)

    # ── Top bar ───────────────────────────────────────────────────
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
    ax_top.set_ylabel("DDR gnsRV\ncount", fontsize=7.5)
    ax_top.set_xticks([])
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_top.spines["bottom"].set_visible(False)
    ax_top.tick_params(axis="y", labelsize=6.5)

    # ── HRD track ─────────────────────────────────────────────────
    hrd_vals = [hrd.get(p, 0) for p in patients]
    hrd_max = max(hrd_vals) if hrd_vals else 50

    # Diverging colormap centered at HRD_THRESHOLD:
    #   blue gradient (below threshold) → pale midpoint → dark green (above)
    hrd_cmap = mcolors.LinearSegmentedColormap.from_list(
        "hrd",
        ["#95a5a6", "#f2f3f4", "#1e8449"],  # gray → near-white → dark green
    )
    hrd_norm = mcolors.TwoSlopeNorm(
        vmin=0,
        vcenter=HRD_THRESHOLD,
        vmax=max(hrd_max, HRD_THRESHOLD + 1),
    )

    for j, (pat, val) in enumerate(zip(patients, hrd_vals)):
        fc = hrd_cmap(hrd_norm(val)) if hrd_max > 0 else "#f2f3f4"
        rect = mpatches.Rectangle(
            (j - 0.45, 0),
            0.9,
            1,
            facecolor=fc,
            edgecolor="white",
            linewidth=0.5,
        )
        ax_hrd.add_patch(rect)
        txt_col = "white" if val >= HRD_THRESHOLD else "#333333"
        ax_hrd.text(
            j,
            0.5,
            str(val),
            ha="center",
            va="center",
            fontsize=5,
            fontweight="bold",
            color=txt_col,
        )

    ax_hrd.set_xlim(-0.5, n_pat - 0.5)
    ax_hrd.set_ylim(0, 1)
    ax_hrd.set_yticks([0.5])
    ax_hrd.set_yticklabels(["HRD"], fontsize=7.5, fontweight="bold")
    ax_hrd.set_xticks([])
    ax_hrd.tick_params(axis="both", length=0)
    for sp in ax_hrd.spines.values():
        sp.set_visible(False)

    # ── Right bar — patients per gene ─────────────────────────────
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
                fontsize=6,
            )
    ax_right.set_ylim(-0.5, n_genes - 0.5)
    ax_right.set_xlabel("# patients", fontsize=7.5)
    ax_right.set_yticks([])
    ax_right.spines["top"].set_visible(False)
    ax_right.spines["right"].set_visible(False)
    ax_right.tick_params(axis="x", labelsize=6.5)
    ax_right.set_title("EPC variants", fontsize=7.5, pad=3)

    # ── CNV bar ───────────────────────────────────────────────────
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
    ax_cnv.set_xlabel("Somatic CNV count", fontsize=7.5)
    ax_cnv.set_yticks([])
    ax_cnv.spines["top"].set_visible(False)
    ax_cnv.spines["right"].set_visible(False)
    ax_cnv.tick_params(axis="x", labelsize=6.5)
    ax_cnv.legend(fontsize=6.5, frameon=False, loc="lower right")
    ax_cnv.set_title("Somatic CNV (EPC)", fontsize=7.5, pad=3)


# ═══════════════════════════════════════════════════════════════════════
#  SECTION 7 — Compose the full figure
# ═══════════════════════════════════════════════════════════════════════


def build_figure(
    tex_path=ONCOPLOT_DATA / "ddr.tex",
    met_path=None,
    nonmet_path=None,
    cnv_path=None,
    out_dir=FIGURES_DIR,
):
    # Resolve data paths — prefer new TSV format, fall back to legacy .tex
    if met_path is None:
        _tsv = ONCOPLOT_DATA / "met_ddr_gnsrv.tsv"
        met_path = _tsv if _tsv.exists() else tex_path

    if nonmet_path is None:
        nonmet_path = ONCOPLOT_DATA / "nonmet_ddr_gnsrv.tsv"

    if cnv_path is None:
        _cnv = ONCOPLOT_DATA / "cnv_ddr.tsv"
        cnv_path = _cnv if _cnv.exists() else None

    _has_nonmet = Path(nonmet_path).exists()

    # Figure height grows when Panel C is included
    _fig_h = 20 if _has_nonmet else 13
    fig = plt.figure(figsize=(15, _fig_h))

    # Outer grid: A (top) | B (met oncoplot) | C (non-met oncoplot, optional)
    _hr = [1, 2.6, 2.6] if _has_nonmet else [1, 2.6]
    outer = GridSpec(
        nrows=3 if _has_nonmet else 2,
        ncols=1,
        height_ratios=_hr,
        hspace=0.14,
        figure=fig,
    )

    # ── Panel (A): bar chart + legend ─────────────────────────────
    top_gs = GridSpecFromSubplotSpec(
        1,
        2,
        subplot_spec=outer[0],
        width_ratios=[1, 1.6],
        wspace=0.25,
    )
    ax_bar = fig.add_subplot(top_gs[0])
    ax_leg = fig.add_subplot(top_gs[1])

    draw_panel_a(ax_bar)

    # Panel label
    ax_bar.text(
        -0.15,
        1.08,
        "A",
        transform=ax_bar.transAxes,
        fontsize=16,
        fontweight="bold",
        va="top",
    )

    # ── Legend (shared for both panels) ───────────────────────────
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
            facecolor=PALETTE["bar_nonmet"], alpha=0.9, label="Non-metastatic (n = 26)"
        ),
    ]
    ax_leg.legend(
        handles=legend_handles,
        loc="center",
        ncol=3,
        fontsize=9,
        frameon=True,
        fancybox=True,
        edgecolor="#cccccc",
        title="Oncoplot legend",
        title_fontsize=10,
        handlelength=1.8,
        handleheight=1.2,
        columnspacing=1.5,
    )

    # ── Panel (B): metastatic oncoplot ───────────────────────────
    ax_b_label = fig.add_subplot(outer[1])
    ax_b_label.axis("off")
    ax_b_label.text(
        -0.02,
        1.02,
        "B",
        transform=ax_b_label.transAxes,
        fontsize=16,
        fontweight="bold",
        va="top",
    )

    draw_panel_b(
        fig,
        outer[1],
        met_path,
        cnv_path=cnv_path,
        arm_label="metastatic arm",
        bar_color=PALETTE["bar_met"],
        zygosity_dict=ZYGOSITY_BOOST,
    )

    # ── Panel (C): non-metastatic oncoplot (only if data file exists) ─
    if _has_nonmet:
        ax_c_label = fig.add_subplot(outer[2])
        ax_c_label.axis("off")
        ax_c_label.text(
            -0.02,
            1.02,
            "C",
            transform=ax_c_label.transAxes,
            fontsize=16,
            fontweight="bold",
            va="top",
        )
        draw_panel_b(
            fig,
            outer[2],
            nonmet_path,
            cnv_path=cnv_path,
            arm_label="non-metastatic arm",
            bar_color=PALETTE["bar_nonmet"],
            zygosity_dict={},
        )

    # ── Save ──────────────────────────────────────────────────────
    out = Path(out_dir)
    fig.savefig(out / "fig_ddr_oncoplot.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(out / "fig_ddr_oncoplot.png", dpi=300, bbox_inches="tight")
    print(f"Saved: {out / 'fig_ddr_oncoplot.pdf'}")
    print(f"Saved: {out / 'fig_ddr_oncoplot.png'}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    build_figure()
