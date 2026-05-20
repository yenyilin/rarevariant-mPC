"""
fig3c_cnv_oncoprint.py — Fig. 3C CNV oncoprint for 12 gnsRV-harboring genes.
Source: data/raw/cbioportal/fig3c_cnv_oncoprint.tsv (cBioPortal export)
Output: figures/fig3c_cnv_oncoprint.{pdf,png}

Cohorts ordered by disease aggressiveness:
  TCGA (localized) → MSK (localized) → SU2C 2015 (mCRPC) → SU2C 2019 (metastatic)
"""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[2]

# ── Constants ─────────────────────────────────────────────────────────────────
GENE_ROWS = list(range(5, 17))  # 0-indexed rows for 12 CNA genes
COLORS = {
    "Amplification": "#D73027",
    "Deep Deletion": "#4575B4",
    "none": "#F2F2F2",
}

COHORT_COLORS = {
    "TCGA": "#7FB3D3",  # localized — light blue
    "MSK": "#2E86C1",  # localized — mid blue
    "SU2C2019": "#E59866",  # metastatic — orange
    "SU2C2015": "#C0392B",  # mCRPC — red
}

COHORT_LABELS = {
    "TCGA": "PCa\n(TCGA, n=494)",
    "MSK": "PCa\n(MSK, n=238)",
    "SU2C2019": "Met PCa\n(SU2C 2019, n=454)",
    "SU2C2015": "mCRPC\n(SU2C 2015, n=125)",
}

COHORT_ORDER = ["TCGA", "MSK", "SU2C2019", "SU2C2015"]


def cohort_key(study):
    if "TCGA" in study:
        return "TCGA"
    if "MSK" in study:
        return "MSK"
    if "PNAS 2019" in study:
        return "SU2C2019"
    if "Cell 2015" in study:
        return "SU2C2015"
    return "TCGA"


# ── Load data ─────────────────────────────────────────────────────────────────
with open(REPO / "data" / "raw" / "cbioportal" / "fig3c_cnv_oncoprint.tsv") as f:
    rows = list(csv.reader(f, delimiter="\t"))

header = rows[0]  # patient IDs from col 2 onward
study_row = rows[1]  # Study of origin
patient_ids = header[2:]
n_patients = len(patient_ids)

# Gene names and CNA matrix
genes = [rows[i][0] for i in GENE_ROWS]
n_genes = len(genes)

# matrix[gene_idx][patient_idx] = alteration string
matrix = []
for i in GENE_ROWS:
    matrix.append(
        [rows[i][j + 2] if j + 2 < len(rows[i]) else "" for j in range(n_patients)]
    )

# ── Assign cohort and sort patients ───────────────────────────────────────────
cohorts = [cohort_key(study_row[j + 2]) for j in range(n_patients)]


# sort: cohort order first, then by alteration count descending
def alteration_count(pidx):
    return sum(1 for g in range(n_genes) if matrix[g][pidx].strip() != "")


order = sorted(
    range(n_patients),
    key=lambda p: (COHORT_ORDER.index(cohorts[p]), -alteration_count(p)),
)

# reorder everything
matrix_sorted = [[matrix[g][p] for p in order] for g in range(n_genes)]
cohorts_sorted = [cohorts[p] for p in order]

# cohort boundary indices
boundaries = []
prev = None
for i, c in enumerate(cohorts_sorted):
    if c != prev:
        boundaries.append(i)
        prev = c
boundaries.append(n_patients)


# ── Build numeric array for rendering ─────────────────────────────────────────
# 0 = none, 1 = Amplification, 2 = Deep Deletion
def encode(val):
    v = val.strip()
    if v == "Amplification":
        return 1
    if v == "Deep Deletion":
        return 2
    return 0


num = np.array(
    [[encode(matrix_sorted[g][p]) for p in range(n_patients)] for g in range(n_genes)],
    dtype=np.int8,
)

# ── Figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 6))
gs = gridspec.GridSpec(
    3,
    2,
    height_ratios=[0.10, 1, 0.04],
    width_ratios=[1, 0.13],
    hspace=0.04,
    wspace=0.02,
)

ax_cohort = fig.add_subplot(gs[0, 0])  # cohort annotation bar
ax_main = fig.add_subplot(gs[1, 0])  # oncoprint
ax_freq = fig.add_subplot(gs[1, 1])  # alteration frequency
ax_leg = fig.add_subplot(gs[2, 0])  # legend space (invisible)
ax_leg.axis("off")

# ── Cohort annotation bar ─────────────────────────────────────────────────────
cohort_colors_arr = np.array(
    [list(matplotlib.colors.to_rgb(COHORT_COLORS[c])) for c in cohorts_sorted]
).reshape(1, n_patients, 3)

ax_cohort.imshow(cohort_colors_arr, aspect="auto", interpolation="none")
ax_cohort.set_xlim(-0.5, n_patients - 0.5)
ax_cohort.axis("off")

# cohort midpoint labels
for idx in range(len(COHORT_ORDER)):
    coh = COHORT_ORDER[idx]
    start = boundaries[idx]
    end = boundaries[idx + 1] if idx + 1 < len(boundaries) else n_patients
    mid = (start + end) / 2
    ax_cohort.text(
        mid,
        0.5,
        COHORT_LABELS[coh],
        ha="center",
        va="center",
        fontsize=7.5,
        color="white",
        fontweight="bold",
        transform=ax_cohort.get_xaxis_transform(),
    )

# ── Oncoprint main panel ──────────────────────────────────────────────────────
bg = np.ones((n_genes, n_patients, 3))
for g in range(n_genes):
    c = matplotlib.colors.to_rgb(COLORS["none"])
    bg[g, :] = c

ax_main.imshow(
    bg,
    aspect="auto",
    interpolation="none",
    extent=[-0.5, n_patients - 0.5, -0.5, n_genes - 0.5],
)

# draw colored rectangles for alterations
bar_h = 0.72
for g in range(n_genes):
    for p in range(n_patients):
        v = num[g, p]
        if v == 0:
            continue
        color = COLORS["Amplification"] if v == 1 else COLORS["Deep Deletion"]
        ax_main.add_patch(
            mpatches.Rectangle(
                (p - 0.5, (n_genes - 1 - g) - bar_h / 2),
                1,
                bar_h,
                color=color,
                linewidth=0,
            )
        )

# cohort dividers
for b in boundaries[1:-1]:
    ax_main.axvline(b - 0.5, color="white", linewidth=1.5, zorder=5)

ax_main.set_xlim(-0.5, n_patients - 0.5)
ax_main.set_ylim(-0.5, n_genes - 0.5)
ax_main.set_yticks(range(n_genes))
ax_main.set_yticklabels(genes[::-1], fontsize=9, fontstyle="italic")
ax_main.set_xticks([])
ax_main.tick_params(left=False, pad=3)
for spine in ax_main.spines.values():
    spine.set_visible(False)

# ── Alteration frequency bar (diverging: AMP right, DEL left) ────────────────
amp_freq = [(num[g] == 1).sum() / n_patients * 100 for g in range(n_genes)]
del_freq = [(num[g] == 2).sum() / n_patients * 100 for g in range(n_genes)]
y_pos = np.arange(n_genes)[::-1]

ax_freq.barh(
    y_pos,
    amp_freq,
    color=COLORS["Amplification"],
    height=0.65,
    edgecolor="white",
    linewidth=0.3,
)
ax_freq.barh(
    y_pos,
    [-d for d in del_freq],
    color=COLORS["Deep Deletion"],
    height=0.65,
    edgecolor="white",
    linewidth=0.3,
)
ax_freq.axvline(0, color="#888888", linewidth=0.6, zorder=5)

max_side = max(max(amp_freq), max(del_freq))
pad = max_side * 0.45
ax_freq.set_xlim(-max_side - pad, max_side + pad)
ax_freq.set_ylim(-0.5, n_genes - 0.5)
ax_freq.set_yticks([])
ax_freq.set_xlabel("% altered", fontsize=8)
ax_freq.tick_params(labelsize=7.5)
ax_freq.spines["top"].set_visible(False)
ax_freq.spines["right"].set_visible(False)
ax_freq.spines["left"].set_visible(False)
ax_freq.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{abs(v):.0f}"))
ax_freq.xaxis.set_major_locator(plt.MaxNLocator(5, integer=True, symmetric=True))

# percentage labels: AMP to the right, DEL to the left
for i in range(n_genes):
    if amp_freq[i] > 0:
        ax_freq.text(
            amp_freq[i] + max_side * 0.06,
            y_pos[i],
            f"{amp_freq[i]:.1f}%",
            va="center",
            fontsize=6.5,
            color=COLORS["Amplification"],
        )
    if del_freq[i] > 0:
        ax_freq.text(
            -del_freq[i] - max_side * 0.06,
            y_pos[i],
            f"{del_freq[i]:.1f}%",
            va="center",
            ha="right",
            fontsize=6.5,
            color=COLORS["Deep Deletion"],
        )

# ── Legend ────────────────────────────────────────────────────────────────────
legend_elements = [
    mpatches.Patch(facecolor=COLORS["Amplification"], label="Amplification"),
    mpatches.Patch(facecolor=COLORS["Deep Deletion"], label="Deep Deletion"),
    mpatches.Patch(
        facecolor=COLORS["none"],
        edgecolor="#AAAAAA",
        linewidth=0.5,
        label="No alteration",
    ),
]
fig.legend(
    handles=legend_elements,
    loc="lower center",
    ncol=3,
    fontsize=9,
    frameon=False,
    bbox_to_anchor=(0.44, 0.1),  # -0.01),
)

plt.suptitle(
    "Somatic Copy Number Alterations in gnsRV-Harboring Genes Across Public Prostate Cancer Cohorts",
    fontsize=11,
    fontweight="bold",
    y=0.92,
)

plt.savefig(REPO / "figures" / "fig3c_cnv_oncoprint.pdf", dpi=300, bbox_inches="tight")
plt.savefig(REPO / "figures" / "fig3c_cnv_oncoprint.png", dpi=300, bbox_inches="tight")
print("Saved: figures/fig3c_cnv_oncoprint.pdf / .png")

# ── stdout summary ────────────────────────────────────────────────────────────
print(f"\nTotal patients: {n_patients}")
for coh in COHORT_ORDER:
    n = cohorts.count(coh)
    print(f"  {COHORT_LABELS[coh]}")
print(f"\nAlteration frequencies (AMP / DEL):")
for g, a, d in zip(genes, amp_freq, del_freq):
    dom = "DEL" if d > a else "AMP"
    print(f"  {g:<12}  AMP {a:5.1f}%  DEL {d:5.1f}%  [{dom}]")
