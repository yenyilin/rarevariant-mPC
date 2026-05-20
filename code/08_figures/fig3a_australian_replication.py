"""
Figure 3A — Australian EPC replication.

Carrier frequency of the 56 genome-wide-screen gnsRVs across phenotype
subgroups of the Australian EPC, tested against the gnomAD-derived
population expectation with a one-sided binomial test.

p_carrier = 0.3115, 1 - prod(1 - MAF_i) over the 56 gnomAD MAFs;
this is the value computed by code/06_replication/australian_binomial.py.

Output: figures/fig3a_australian_replication.{pdf,png}
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import binom

REPO = Path(__file__).resolve().parents[2]

# ── Author-confirmed values ───────────────────────────────────────
P_CARRIER = 0.3115

groups_raw = [
    ("Metastatic (bone/visceral)", 17, 9),
    ("Metastatic + BCR", 29, 14),
    ("Metastatic + BCR + nodal", 33, 16),
    ("Non-metastatic (control)", 14, 5),
]

subgroups, ns, carriers, obs_freqs, exp_carriers, pvals, directions = (
    [],
    [],
    [],
    [],
    [],
    [],
    [],
)
for name, n, k in groups_raw:
    p_val = 1 - binom.cdf(k - 1, n, P_CARRIER)
    subgroups.append(name)
    ns.append(str(n))
    carriers.append(str(k))
    obs_freqs.append(f"{k / n * 100:.1f}%")
    exp_carriers.append(f"{n * P_CARRIER:.1f}")
    pvals.append(f"{p_val:.3f}")
    directions.append("↑ Enriched" if "Non-met" not in name else "— No enrichment")

data = {
    "Phenotype subgroup": subgroups,
    "n": ns,
    "Observed\ncarriers": carriers,
    "Expected\ncarriers": exp_carriers,
    "Carrier freq.\n(observed)": obs_freqs,
    "p-value\n(one-sided binomial)": pvals,
    "Direction": directions,
    "Variants\ntested": ["56 gnsRVs\n(genome-wide screen)"] * 4,
}
df = pd.DataFrame(data)

# ── Render ─────────────────────────────────────────────────────────────────────
df_render = df.copy()
for col in df_render.columns:
    df_render[col] = df_render[col].str.replace("\n", " ")

fig, ax = plt.subplots(figsize=(14, 3.2))
ax.axis("off")

col_widths = [0.18, 0.04, 0.08, 0.08, 0.08, 0.12, 0.1, 0.18]

tbl = ax.table(
    cellText=df_render.values,
    colLabels=df_render.columns,
    cellLoc="center",
    loc="center",
    colWidths=col_widths,
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(9)
tbl.scale(1, 2.4)

GREEN = "#1E8449"
ROW_COLS = ["#EAFAF1", "#D5F5E3", "#A9DFBF", "#F2F3F4"]

for j in range(len(df_render.columns)):
    tbl[0, j].set_facecolor(GREEN)
    tbl[0, j].set_text_props(color="white", weight="bold", fontsize=8.5)

p_col = [i for i, c in enumerate(df_render.columns) if "p-value" in c][0]

for i in range(1, len(df_render) + 1):
    is_ctrl = i == 4  # non-metastatic control row
    for j in range(len(df_render.columns)):
        tbl[i, j].set_facecolor(ROW_COLS[i - 1])
        if is_ctrl:
            tbl[i, j].set_text_props(color="#888888", fontweight="bold")
        else:
            tbl[i, j].set_text_props(fontweight="bold")
    if not is_ctrl:
        tbl[i, p_col].set_text_props(weight="bold", color=GREEN)

# ── Subtitle and footnote (no panel label — add in Inkscape) ──────────────────
# ax.set_title(
#    "Australian EPC (n=53) — carrier frequency of 56 genome-wide gnsRVs by phenotype subgroup",
#    fontsize=9.5,
#    pad=0,  # 8,
#    weight="bold",
#    loc="left",
# )

ax.text(
    x=0.4,  # 0.5
    y=0.95,  # 1.02,  # in axes coordinates (0–1); >1 is above the axes
    s="Australian EPC (n=53) — carrier frequency of 56 genome-wide gnsRVs by phenotype subgroup",
    transform=ax.transAxes,
    ha="center",
    va="bottom",
    fontsize=12,
    fontweight="bold",
)

# fig.text(
#     0.09,  # 0.01,
#     0.14,  # 0.01,
#     "BCR, biochemical recurrence. Carrier defined as ≥1 of the 56 genome-wide screen gnsRVs. "
#     "Expected frequency from gnomAD v2.1.1 via one-sided binomial test "
#     "(p_carrier = 0.310).",
#     fontsize=7.5,
#     color="#555",
#     style="italic",
# )
#
plt.tight_layout(rect=[0, 0.08, 1, 1])
out = REPO / "figures" / "fig3a_australian_replication"
plt.savefig(f"{out}.pdf", dpi=300, bbox_inches="tight")
plt.savefig(f"{out}.png", dpi=300, bbox_inches="tight")
print(f"Saved: {out}.pdf / .png")
