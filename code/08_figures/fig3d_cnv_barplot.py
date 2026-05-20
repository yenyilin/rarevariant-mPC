"""
fig3d_cnv_barplot.py — Fig. 3D cBioPortal CNV alteration frequency (stacked bar plot)
-------------------------------------------------------------------------------------
Displays homdel and amp percentage for each cohort.
Source: data/raw/cbioportal/fig3d_cnv_frequency.tsv
Output: figures/fig3d_cnv_barplot.{pdf,png}
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]

# ── Data ──────────────────────────────────────────────────────────────────────
df = pd.read_csv(
    REPO / "data" / "raw" / "cbioportal" / "fig3d_cnv_frequency.tsv", sep="\t"
)
df.columns = df.columns.str.strip()

# Keep only homdel and amp
df = df[df["Alteration Type"].isin(["homdel", "amp"])].copy()

# Derive sample size from count / frequency * 100
df["n"] = (
    (df["Alteration Count"] / df["Alteration Frequency"] * 100).round().astype(int)
)

# Short cohort labels (derived from study name)
label_map = {
    "Metastatic Prostate Adenocarcinoma (SU2C/PCF Dream Team, PNAS 2019)": "SU2C/PCF\nPNAS 2019\n(mPCa, n=429)",
    "Metastatic Prostate Cancer (SU2C/PCF Dream Team, Cell 2015)": "SU2C/PCF\nCell 2015\n(mPCa, n=150)",
    "Prostate Adenocarcinoma (TCGA, PanCancer Atlas)": "TCGA\nPanCancer Atlas\n(PCa, n=494)",
    "Prostate Adenocarcinoma (MSK, Cancer Cell 2010)": "MSK\nCancer Cell 2010\n(PCa, n=205)",
}
df["Label"] = df["Cancer Study"].map(label_map)

# Pivot to wide format: rows = cohorts, cols = alteration type
pivot = df.pivot_table(
    index="Label",
    columns="Alteration Type",
    values="Alteration Frequency",
    aggfunc="first",
).fillna(0)

# Preserve clinical order (metastatic first, then localised)
order = [
    "SU2C/PCF\nPNAS 2019\n(mPCa, n=429)",
    "SU2C/PCF\nCell 2015\n(mPCa, n=150)",
    "TCGA\nPanCancer Atlas\n(PCa, n=494)",
    "MSK\nCancer Cell 2010\n(PCa, n=205)",
]
pivot = pivot.reindex(order)

# ── Plot ──────────────────────────────────────────────────────────────────────
# HOMDEL_COLOR = "#C0392B"  # red — deletion
# AMP_COLOR = "#2980B9"  # blue — amplification
HOMDEL_COLOR = "#2980B9"  # blue — deletion
AMP_COLOR = "#C0392B"  # red  — amplification
EDGE_COLOR = "#444444"

fig, ax = plt.subplots(figsize=(5.5, 5.2))

x = np.arange(len(pivot))
bar_width = 0.6

# ── Background shading by cohort group ───────────────────────────────────────
# SU2C (metastatic): warm amber tint; TCGA+MSK (adenocarcinoma): cool gray tint
ax.axvspan(-0.5, 1.5, color="#FEF5E7", zorder=0)  # SU2C group (cols 0-1)
ax.axvspan(1.5, 3.5, color="#EAF2F8", zorder=0)  # TCGA+MSK group (cols 2-3)

bars_homdel = ax.bar(
    x,
    pivot["homdel"],
    width=bar_width,
    color=HOMDEL_COLOR,
    edgecolor=EDGE_COLOR,
    linewidth=0.6,
    label="Homozygous deletion (homdel)",
    zorder=3,
)
bars_amp = ax.bar(
    x,
    pivot["amp"],
    width=bar_width,
    bottom=pivot["homdel"],
    color=AMP_COLOR,
    edgecolor=EDGE_COLOR,
    linewidth=0.6,
    label="Amplification (amp)",
    zorder=3,
)

# ── Value labels ──────────────────────────────────────────────────────────────
for i, (hd, am) in enumerate(zip(pivot["homdel"], pivot["amp"])):
    if hd > 0:
        ax.text(
            i,
            hd / 2,
            f"{hd:.1f}%",
            ha="center",
            va="center",
            fontsize=9,
            color="white",
            fontweight="bold",
        )
    if am > 0:
        ax.text(
            i,
            hd + am / 2,
            f"{am:.1f}%",
            ha="center",
            va="center",
            fontsize=9,
            color="white",
            fontweight="bold",
        )

# ── Group bracket annotations ────────────────────────────────────────────────
GROUP_MET = "#B7770D"  # amber-brown for metastatic group label
GROUP_LOC = "#1A5276"  # navy for localised group label
bracket_y = 61

ax.annotate(
    "",
    xy=(1.45, bracket_y),
    xytext=(-0.45, bracket_y),
    arrowprops=dict(arrowstyle="-", color=GROUP_MET, lw=1.5),
)
ax.text(
    0.5,
    bracket_y + 1.2,
    "Metastatic",
    ha="center",
    va="bottom",
    fontsize=8.5,
    color=GROUP_MET,
    fontweight="bold",
)

ax.annotate(
    "",
    xy=(3.45, bracket_y),
    xytext=(1.55, bracket_y),
    arrowprops=dict(arrowstyle="-", color=GROUP_LOC, lw=1.5),
)
ax.text(
    2.5,
    bracket_y + 1.2,
    "Localised adenocarcinoma",
    ha="center",
    va="bottom",
    fontsize=8.5,
    color=GROUP_LOC,
    fontweight="bold",
)

# ── Axes formatting ───────────────────────────────────────────────────────────
ax.set_xticks(x)
ax.set_xticklabels(pivot.index, fontsize=8.5)

# Color x-tick labels by group
for tick, color in zip(
    ax.get_xticklabels(), [GROUP_MET, GROUP_MET, GROUP_LOC, GROUP_LOC]
):
    tick.set_color(color)

ax.set_ylabel("Alteration frequency (%)", fontsize=10)
ax.set_xlim(-0.5, 3.5)
ax.set_ylim(0, 68)
ax.yaxis.grid(True, linestyle="--", linewidth=0.5, alpha=0.6, zorder=0)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.set_title(
    "Somatic CNV frequency across PCa cohorts\n(cBioPortal)",
    fontsize=10,
    fontweight="bold",
    pad=10,
)

# ── Legend ────────────────────────────────────────────────────────────────────
patches = [
    mpatches.Patch(
        facecolor=HOMDEL_COLOR,
        edgecolor=EDGE_COLOR,
        linewidth=0.6,
        label="Homozygous deletion (homdel)",
    ),
    mpatches.Patch(
        facecolor=AMP_COLOR,
        edgecolor=EDGE_COLOR,
        linewidth=0.6,
        label="Amplification (amp)",
    ),
]
# ax.legend(handles=patches, fontsize=8.5, frameon=False, loc="center right")
ax.legend(
    handles=patches,
    fontsize=12,
    frameon=False,
    loc="center right",
    bbox_to_anchor=(1.0, 0.7),
)
# ── Footnote ─────────────────────────────────────────────────────────────────
fig.text(
    0.05,
    0.05,
    "Source: cBioPortal. mPCa = metastatic prostate cancer; PCa = localised prostate adenocarcinoma.",
    fontsize=7,
    color="#555",
    style="italic",
)

plt.tight_layout(rect=[0, 0.04, 1, 1])
plt.savefig(REPO / "figures" / "fig3d_cnv_barplot.pdf", dpi=300, bbox_inches="tight")
plt.savefig(REPO / "figures" / "fig3d_cnv_barplot.png", dpi=300, bbox_inches="tight")
print("Saved: figures/fig3d_cnv_barplot.pdf / .png")
