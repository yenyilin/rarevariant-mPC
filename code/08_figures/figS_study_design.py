"""
Supplementary Figure S1 — study design schematic.

Shows the two-stage discovery-replication design, with visual separation
of prespecified (DDR) vs. exploratory (genome-wide) analyses.

Key messages this figure communicates:
  1. EPC design rationale: extreme phenotype (met vs. non-met), matched, long follow-up
  2. Two-stage analytic structure: EPC discovery → independent replication
  3. Prespecified (DDR) vs. exploratory (genome-wide) — visually distinct
  4. Functional characterization is illustrative, not the primary evidence stream

Run with: python code/08_figures/figS_study_design.py
Output:   figures/figS_study_design.{pdf,png}
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

# ── Color palette ─────────────────────────────────────────────────────────────
C_MET = "#B03A2E"  # red — metastatic
C_NONMET = "#2E86C1"  # blue — non-metastatic
C_PRESPEC = "#1A5276"  # dark blue — prespecified
C_EXPLOR = "#6C3483"  # purple — exploratory
C_REPLIC = "#1E8449"  # green — replication
C_FUNC = "#B7770D"  # amber — functional
C_BOX_LIGHT = "#F4F6F7"  # light grey fill
C_BOX_MID = "#EBF5FB"  # light blue fill
C_BORDER = "#717D7E"  # grey border

fig, ax = plt.subplots(figsize=(11, 8.5))
ax.set_xlim(0, 11)
ax.set_ylim(0, 8.5)
ax.axis("off")


# ── Helper functions ──────────────────────────────────────────────────────────
def box(
    ax,
    x,
    y,
    w,
    h,
    text,
    fontsize=12,  # 9,
    facecolor=C_BOX_LIGHT,
    edgecolor=C_BORDER,
    lw=1.2,
    bold=False,
    textcolor="black",
    style="round,pad=0.1",
    valign="center",
):
    rect = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=style,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=lw,
        zorder=3,
    )
    ax.add_patch(rect)
    weight = "bold" if bold else "normal"
    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        ha="center",
        va=valign,
        fontsize=fontsize,
        color=textcolor,
        weight=weight,
        multialignment="center",
        zorder=4,
    )


def arrow(ax, x1, y1, x2, y2, color="#444444", lw=1.5, style="->", head=12):
    ax.annotate(
        "",
        xy=(x2, y2),
        xytext=(x1, y1),
        arrowprops=dict(arrowstyle=style, color=color, lw=lw, mutation_scale=head),
        zorder=5,
    )


def label_badge(ax, x, y, text, color):
    ax.text(
        x,
        y,
        text,
        fontsize=7.5,
        color="white",
        weight="bold",
        ha="center",
        va="center",
        zorder=6,
        bbox=dict(
            boxstyle="round,pad=0.3", facecolor=color, edgecolor="none", alpha=0.9
        ),
    )


# ══════════════════════════════════════════════════════════════════════════════
# ROW 1 — EPC cohort (discovery)
# ══════════════════════════════════════════════════════════════════════════════
box(
    ax,
    0.3,
    6.7,
    10.4,
    1.3,
    "\n\nExtreme Phenotype Cohort (EPC): n = 52\n"
    "High-grade, localized, treatment-naive PCa with ≥ 7 years follow-up\n"
    "Whole-exome sequencing (tumor + matched benign tissue)  gnsRVs: MAF ≤2%\n\n\n\n\n\n",
    fontsize=12,
    facecolor="#EBF5FB",
    edgecolor="#2E86C1",
    lw=2.5,
    bold=False,
)  # lw=1.8, bold=False)

# Split EPC into two arms
box(
    ax,
    2,
    7,
    2.3,
    0.25,  # 0.45, 7.15, 2.3, 0.75,
    "Metastatic\n(n = 26)",
    fontsize=12,  # 9,
    facecolor="#FADBD8",
    edgecolor=C_MET,
    lw=1.4,
    textcolor=C_MET,
)
box(
    ax,
    7,
    7,
    2.3,
    0.25,  # 3.0, 7.15, 2.6, 0.75,
    "Non-metastatic\n(n = 26)",
    fontsize=12,  # 9,
    facecolor="#D6EAF8",
    edgecolor=C_NONMET,
    lw=1.4,
    textcolor=C_NONMET,
)

ax.text(
    4.6,  # 5.85,
    7.2,  # 7.52,
    "Matched: ancestry, GS, PSA,\nstage, margins, follow-up",
    fontsize=10,  # 7.8,
    ha="left",
    va="center",
    color="#5D6D7E",
    style="italic",
)

# ── Bracket connecting two arms ───────────────────────────────────────────────
ax.annotate(
    "",
    xy=(3.15, 6.95),  # (1.6, 7.15),
    xytext=(8.15, 6.95),  # (4.3, 7.15),
    arrowprops=dict(
        arrowstyle="-",
        color=C_BORDER,
        lw=1.5,
        connectionstyle="bar,fraction=-0.12",  # lw=1.0, connectionstyle="bar,fraction=-0.3"
    ),
    zorder=10,  # 2,
)

# ══════════════════════════════════════════════════════════════════════════════
# ROW 2 — Two analysis branches
# ══════════════════════════════════════════════════════════════════════════════
arrow(ax, 5.5, 6.9, 5.5, 6.35, color=C_BORDER)

ax.text(
    5.5,
    6.3,  # 6.22,
    "Variant discovery",
    ha="center",
    va="top",
    fontsize=12,  # 8.5,
    color="#444",
    weight="bold",
)

# Prespecified branch (left)
box(
    ax,
    0.3,
    5.0,  # 4.8,
    4.6,
    0.8,  # 1.1,
    "PRIMARY HYPOTHESIS TEST  (prespecified)\n"
    "25 canonical DDR genes\n"
    "Bootstrapping + Wilcoxon rank-sum",
    fontsize=11,  # 8.8,
    facecolor="#D6EAF8",
    edgecolor=C_PRESPEC,
    lw=1.8,
    textcolor="#1A3A4A",
)

label_badge(
    ax, 1.35, 5.75, "PRESPECIFIED", C_PRESPEC
)  # 1.15, 5.75, "PRESPECIFIED", C_PRESPEC)

# Result box under prespecified
box(
    ax,
    0.75,  # 0.55,
    4.05,
    3.6,  # 4.1,
    0.5,  # 0.65,
    "81 met-exclusive gnsRVs in DDR genes  (p = 4.57×10⁻⁶)\n"
    "Synonymous variants: no difference (internal control)",
    fontsize=9,  # 8,
    facecolor="#EBF5FB",
    edgecolor=C_PRESPEC,
    lw=1.0,
)
arrow(
    ax, 2.6, 4.9, 2.6, 4.6, color=C_PRESPEC, lw=1.2
)  # 2.6, 4.8, 2.6, 4.7, color=C_PRESPEC, lw=1.2)

# Sensitivity analysis callout — compact box in center gap, connected to DDR result
box(
    ax,
    4.75,
    4.02,
    1.48,  # 1.35,
    0.65,
    "Sensitivity analyses\n"
    "Firth LR(adj. age+ancestry):\n OR=25.0, p<0.001 \n"
    "Follow-up MWU: p<0.001\n Burden ρ=−0.12, p=0.55",
    fontsize=8,  # 6.3,
    facecolor="#F2F3F4",
    edgecolor="#85929E",
    lw=0.8,
)
# ax.annotate(
#     "",
#     xy=(4.75, 4.35),
#     xytext=(4.65, 4.35),
#     arrowprops=dict(arrowstyle="-", color="#85929E", lw=0.8),
#     zorder=4,
# )

# Exploratory branch (right)
box(
    ax,
    5.9,
    5.0,  # 4.8,
    4.8,
    0.8,  # 1.1,
    "SECONDARY EXPLORATORY SCREEN (genome-wide)\n"
    "gnsRVs exclusive to ≥3 metastatic patients\n"
    "CNV + survival analysis in public cohorts",
    fontsize=11,  # 8.8,
    facecolor="#F5EEF8",
    edgecolor=C_EXPLOR,
    lw=1.8,
    textcolor="#3B1A5A",
)

label_badge(ax, 7.1, 5.75, "EXPLORATORY", C_EXPLOR)

# Result box under exploratory
box(
    ax,
    6.75,  # 6.15,
    4.05,
    3.6,  # 4.3,
    0.5,  # 0.65,
    "56 gnsRVs in 53 genes  (p = 2.98×10⁻⁸)\n"
    "Somatic selection & survival signal in public PCa cohorts",
    fontsize=9,  # 8,
    facecolor="#F9F3FD",
    edgecolor=C_EXPLOR,
    lw=1.0,
)
(
    arrow(ax, 8.3, 4.9, 8.3, 4.6, color=C_EXPLOR, lw=1.2),
)  # 8.3, 4.8, 8.3, 4.7, color=C_EXPLOR, lw=1.2)

# Arrow from EPC down to two branches
arrow(ax, 5.5, 6.15, 2.6, 5.9, color=C_BORDER, lw=1.3)
arrow(ax, 5.5, 6.15, 8.3, 5.9, color=C_BORDER, lw=1.3)

# ══════════════════════════════════════════════════════════════════════════════
# ROW 3 — Replication (independent cohorts)
# ══════════════════════════════════════════════════════════════════════════════
ax.text(
    5.5,
    3.8,
    "Independent replication",
    ha="center",
    va="top",
    fontsize=12,  # 8.5,
    color="#444",
    weight="bold",
)

arrow(ax, 2.6, 4.05, 3.5, 3.35, color=C_REPLIC, lw=1.5)
arrow(ax, 8.3, 4.05, 7.5, 3.35, color=C_REPLIC, lw=1.5)

# Australian EPC box
box(
    ax,
    0.3,
    2.35,
    4.7,  # 4.6,
    1.1,
    "Australian EPC  (n = 53)\n"
    "Carrier frequency > expected (binomial)\n"
    "p = 0.03–0.05 (metastatic)  vs.  p = 0.45 (non-metastatic)",
    fontsize=11,  # 8.8,
    facecolor="#EAFAF1",
    edgecolor=C_REPLIC,
    lw=1.8,
)

label_badge(ax, 1.35, 3.3, "REPLICATION 1", C_REPLIC)

# PPCG box
box(
    ax,
    5.9,
    2.35,
    4.8,
    1.1,
    "PPCG cohort  (n = 976,  including 200 metastatic)\n"
    "Aggregate gnsRV score: high-risk disease  (SKAT p < 0.03)\n"
    "Gene-level: CHEK2, FOCAD, SPATA9, ZSWIM4, DNAJC10",
    fontsize=11,  # 8.8,
    facecolor="#EAFAF1",
    edgecolor=C_REPLIC,
    lw=1.8,
)

label_badge(ax, 7.15, 3.3, "REPLICATION 2", C_REPLIC)

# ══════════════════════════════════════════════════════════════════════════════
# ROW 4 — Functional characterization (side branch, illustrative)
# ══════════════════════════════════════════════════════════════════════════════
# NOTE: functional arrow NO LONGER originates from Australian EPC box.
# Routing arrows from discovery result boxes are added below (after FUNCTIONAL badge).

box(
    ax,
    0.3,
    0.65,
    10.4,
    1.0,
    "Functional characterization  (2 discovery-stage candidates selected prior to and independent of replication)\n"
    "KDM6B K973Q  [genome-wide screen, ≥3 metastatic carriers] → migration & invasion in LNCaP\n"
    # "BRCA2 I1962T  [prespecified DDR candidate] → Olaparib sensitivity in HEK293T",
    "BRCA2 I1962T  [prespecified DDR candidate; ultrarare] → Olaparib sensitivity in HEK293T",
    fontsize=12,  # 8.8,
    facecolor="#FEF9E7",
    edgecolor=C_FUNC,
    lw=1.8,
)

# label_badge(ax, 1.6, 1.5, "FUNCTIONAL", C_FUNC)
label_badge(ax, 1.35, 1.6, "FUNCTIONAL", C_FUNC)

# ── Dashed routing arrows: discovery → functional (bypasses replication row) ──
# Encodes the correct selection logic: candidates chosen at discovery stage,
# not after Australian EPC replication.

# LEFT route: DDR result box left edge → left margin → functional top-left
#   Represents BRCA2 I1962T (prespecified DDR candidate)
lx = [0.55, 0.08, 0.08, 0.30]
ly = [4.37, 4.37, 1.65, 1.65]
ax.plot(lx, ly, color=C_FUNC, lw=1.3, linestyle="--", zorder=4)
ax.annotate(
    "",
    xy=(0.30, 1.65),
    xytext=(0.14, 1.65),
    arrowprops=dict(arrowstyle="->", color=C_FUNC, lw=1.3, mutation_scale=10),
    zorder=5,
)
ax.text(
    0.22,
    2.85,
    "DDR\ncandidate",
    ha="center",
    va="center",
    fontsize=7.0,  # 6.0,
    color=C_FUNC,
    style="italic",
    rotation=90,
)

# RIGHT route: exploratory result box right edge → right margin → functional top-right
#   Represents KDM6B K973Q (genome-wide screen, ≥3 metastatic carriers)
rx = [10.45, 10.92, 10.92, 10.70]
ry = [4.37, 4.37, 1.65, 1.65]
ax.plot(rx, ry, color=C_FUNC, lw=1.3, linestyle="--", zorder=4)
ax.annotate(
    "",
    xy=(10.70, 1.65),
    xytext=(10.86, 1.65),
    arrowprops=dict(arrowstyle="->", color=C_FUNC, lw=1.3, mutation_scale=10),
    zorder=5,
)
ax.text(
    10.78,
    2.85,
    "DDR and Genome-wide\ncandidate",
    ha="center",
    va="center",
    fontsize=7.0,  # 6.0,
    color=C_FUNC,
    style="italic",
    rotation=90,
)

# Floating 3-criterion selection label between replication and functional rows
ax.text(
    5.5,
    1.92,
    "3-criterion selection: exclusive to metastatic arm  ·  PCa biological relevance  ·  cell model availability",
    ha="center",
    va="center",
    fontsize=11.0,  # 7.0,
    color=C_FUNC,
    style="italic",
    zorder=6,
)

# ══════════════════════════════════════════════════════════════════════════════
# Stage labels (left margin)
# ══════════════════════════════════════════════════════════════════════════════
for y, txt, color in [
    (7.5, "COHORT", "#2E86C1"),
    (5.35, "DISCOVERY", "#444444"),
    (2.9, "REPLICATION", C_REPLIC),
    (1.2, "FUNCTIONAL", C_FUNC),
]:
    ax.text(
        -0.05,
        y,
        txt,
        ha="right",
        va="center",
        fontsize=11,  # 7.5,
        color=color,
        weight="bold",
        rotation=90,
    )

# ══════════════════════════════════════════════════════════════════════════════
# Title
# ══════════════════════════════════════════════════════════════════════════════
ax.text(
    5.5,
    8.4,
    "Study Design Overview",
    ha="center",
    va="top",
    fontsize=12,
    weight="bold",
    color="#1C2833",
)

plt.tight_layout(pad=0.2)
from pathlib import Path as _Path

_FIG = _Path(__file__).resolve().parents[2] / "figures"
_FIG.mkdir(exist_ok=True)
plt.savefig(_FIG / "figS_study_design.pdf", dpi=300, bbox_inches="tight")
plt.savefig(_FIG / "figS_study_design.png", dpi=300, bbox_inches="tight")
print("Saved: figS_study_design.pdf / .png")
# plt.show()  # disabled for non-interactive rendering
