"""
Figure 3B — PPCG cohort (n=976) gene-level associations.

Single-variant Wald and gene-level Burden p-values for candidate gnsRVs,
plus the aggregate SKAT result. The Burden test is the per-gene statistic
(column marked †); SKAT is the single aggregate test in the final row.

Output: figures/fig3b_ppcg_associations.{pdf,png}
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

REPO = Path(__file__).resolve().parents[2]

# ── Author-confirmed values (2026-04-16) ───────────────────────────────────────
data = {
    "Gene": [
        "FOCAD",
        "CHEK2",
        "SIVA1",
        "ACSF3",
        "HPGDS",
        "SPATA9",
        "ZSWIM4",
        "DNAJC10",
    ],
    "gnsRVs\ntested": [
        "rs117863779",
        "rs17879961",
        "rs8010264,\nrs141296448",
        "rs150487794,\nrs141090143",
        "rs76328980,\nrs34124298,\nrs61752528",
        "rs55796768",
        "rs76549352,\nrs56177954",
        "rs13414223",
    ],
    "PPCG\ncarriers": ["1", "9", "2,6", "2,6", "13,9,8", "32", "56,56", "23"],
    "Phenotype": [
        "Age at RP",
        "Age at RP",
        "PSA at RP",
        "PSA at RP",
        "PSA at RP",
        "Metastasis",
        "High-risk",
        "High-risk",
    ],
    "Wald\np-value": [
        "<0.001",
        "<0.001",
        "<0.001",
        "<0.001",
        "<0.001",
        "0.005",
        "0.001",
        "0.001",
    ],
    "p-value\n†": [
        "<0.01",
        "<0.001",
        "<0.001",
        "<0.001",
        "<0.01",
        "<0.001",
        "<0.001",
        "<0.001",
    ],
    "Single-allele\ndriven?": ["Yes", "Yes", "No", "No", "No", "Yes", "No", "Yes"],
    "Discovery\ntier": [
        "Recurrent\nscreen",
        "DDR panel",
        "Extended",
        "Extended",
        "Extended",
        "Extended",
        "Recurrent\nscreen",
        "Recurrent\nscreen",
    ],
}

skat_row = {
    "Gene": "ALL 289 genes (aggregate)",
    "gnsRVs\ntested": "352 detected",
    "PPCG\ncarriers": "-",
    "Phenotype": "High-risk disease",
    "Wald\np-value": "—",
    "p-value\n†": "0.03",
    "Single-allele\ndriven?": "N/A",
    "Discovery\ntier": "Mixed",
}

df_full = pd.concat([pd.DataFrame(data), pd.DataFrame([skat_row])], ignore_index=True)
df_r = df_full.copy().astype(str)

# ── Render ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(14, 4.2))
ax.axis("off")

col_widths = [0.16, 0.10, 0.09, 0.12, 0.09, 0.09, 0.08, 0.08]

tbl = ax.table(
    cellText=df_r.values,
    colLabels=df_r.columns,
    cellLoc="center",
    loc="center",
    colWidths=col_widths,
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(9)
tbl.scale(1, 2.1)

# Taller rows for genes listing multiple rsIDs (row index → number of text lines)
line_counts = {1: 1, 2: 1, 3: 2, 4: 2, 5: 3, 6: 1, 7: 2, 8: 1, 9: 1}
base_height = tbl[1, 0].get_height()
for row_idx, n_lines in line_counts.items():
    if n_lines > 1:
        for j in range(len(df_r.columns)):
            tbl[row_idx, j].set_height(base_height * (n_lines / 2))

# ── Styling ────────────────────────────────────────────────────────────────────
BLUE = "#1A5276"
AMBER = "#FEF9E7"

for j in range(len(df_r.columns)):
    tbl[0, j].set_facecolor(BLUE)
    tbl[0, j].set_text_props(color="white", weight="bold", fontsize=8.5)

wald_col = list(df_r.columns).index("Wald\np-value")
burden_col = list(df_r.columns).index("p-value\n†")

for i in range(1, len(df_r) + 1):
    row = df_full.iloc[i - 1]
    is_skat = row["Gene"] == "ALL 289 genes (aggregate)"
    is_ddr = row["Discovery\ntier"] == "DDR panel"

    bg = AMBER if is_skat else ("#D6EAF8" if is_ddr else "#F2F3F4")
    for j in range(len(df_r.columns)):
        tbl[i, j].set_facecolor(bg)
        if is_skat:
            tbl[i, j].set_text_props(weight="bold")

    for col_idx, col_key in [(wald_col, "Wald\np-value"), (burden_col, "p-value\n†")]:
        if str(row[col_key]) not in ["[?]", "—", "N/A"]:
            tbl[i, col_idx].set_text_props(weight="bold", color=BLUE)

# ── Subtitle and footnote (no panel label — add in Inkscape) ──────────────────
# ax.set_title(
#    "PPCG cohort (n=976) — single-variant Wald and gene-level SKAT associations",
#    fontsize=9.5,
#    pad=10,
#    weight="bold",
#    loc="left",
# )

# ── Title ──────────────────────────────────────────────────────────────────────
ax.text(
    x=0.432,  # 0.5
    y=1.08,  # 1.02,
    s="PPCG cohort (n=976): single-variant Wald, gene-level Burden, and aggregated SKAT associations",
    transform=ax.transAxes,
    ha="center",
    va="bottom",
    fontsize=12,
    fontweight="bold",
)

# # ── Footnote ───────────────────────────────────────────────────────────────────
# fig.text(
#     0.11,
#     0.06,
#     "† Gene-level Burden test p-value for individual gene associations; aggregate SKAT "
#     "p-value for the ALL 289 genes row.\n"
#     "PPCG: 531 high-risk (including 200 metastatic), 305 intermediate, 140 low-risk. "
#     "Discovery tiers: DDR panel = 25 prespecified DDR genes from Figure 1; "
#     "Recurrent screen = 56 metastatic-exclusive gnsRVs (≥3 \npatients) from Figure 2; "
#     "Extended = 343 gnsRVs from 206 metastatic-exclusive genes. "
#     '"Single-allele driven?" flags associations attributable to one recurrent allele.',
#     fontsize=7.5,
#     color="#555",
#     style="italic",
# )


# ── Legend ─────────────────────────────────────────────────────────────────────
patches = [
    mpatches.Patch(
        facecolor="#D6EAF8",
        edgecolor="#888888",
        linewidth=0.8,
        label="Prespecified\nDDR gene",
    ),
    mpatches.Patch(
        facecolor="#F2F3F4",
        edgecolor="#888888",
        linewidth=0.8,
        label="Exploratory\ncandidate gene",
    ),
    mpatches.Patch(
        facecolor=AMBER, edgecolor="#888888", linewidth=0.8, label="SKAT aggregate"
    ),
]
ax.legend(
    handles=patches,
    loc="lower right",
    fontsize=8,
    frameon=False,
    bbox_to_anchor=(1.02, 0.0),
    handleheight=1.2,
    handlelength=1.5,
)

plt.tight_layout(rect=[0, 0.1, 1, 1])
out = REPO / "figures" / "fig3b_ppcg_associations"
plt.savefig(f"{out}.pdf", dpi=300, bbox_inches="tight")
plt.savefig(f"{out}.png", dpi=300, bbox_inches="tight")
print(f"Saved: {out}.pdf / .png")
