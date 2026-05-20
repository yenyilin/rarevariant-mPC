"""
fig4a_proliferation.py — CCK-8 proliferation: LNCaP parental vs KDM6B K973Q
knock-in clones C3 (homozygous) and C6 (heterozygous).

Source: data/raw/functional_assays/fig4a_cck8_proliferation.tsv
Output: figures/fig4a_proliferation.{pdf,png}

Data structure: 3 cell lines × 3 biological replicates × 3 technical replicates.

Toggle:
  USE_BR_MEANS = False  → 9 TR values per clone (default)
  USE_BR_MEANS = True   → 3 BR means per clone
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import csv
from itertools import combinations

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
DATA_FILE = REPO / "data" / "raw" / "functional_assays" / "fig4a_cck8_proliferation.tsv"
FIGURES_DIR = REPO / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

# ── Toggle ───────────────────────────────────────────────────────────────────
USE_BR_MEANS = False  # False = 9 TR values (default) | True = 3 BR means

# ── Style 2 editorial palette ─────────────────────────────────────────────────
COLORS = {"LNCaP": "#C0513A", "C3": "#2E7D6A", "C6": "#4B5EA8"}
COLORS_LITE = {"LNCaP": "#E8A898", "C3": "#82BFB0", "C6": "#A8B4DC"}
FONT = "Arial"

# ── Load data ─────────────────────────────────────────────────────────────────
rows = []
with open(DATA_FILE) as f:
    reader = csv.reader(f, delimiter="\t")
    next(reader)
    for row in reader:
        if row and row[0].strip():
            rows.append(row)


def parse_row(row):
    name = row[0].strip()
    vals = [float(x) for x in row[1:] if x.strip()]
    return name, vals[0:3], vals[3:6], vals[6:9]


cells = {}
for row in rows:
    name, tr1, tr2, tr3 = parse_row(row)
    br_means = [np.mean(tr1), np.mean(tr2), np.mean(tr3)]
    cells[name] = {"br_means": br_means, "all": tr1 + tr2 + tr3}

labels = ["LNCaP", "C3", "C6"]
sublabel = ["Parental\n(LNCaP)", "Homozygous\n(C3)", "Heterozygous\n(C6)"]

# ── Select values for statistics and scatter ──────────────────────────────────
key = "br_means" if USE_BR_MEANS else "all"
vals = {l: np.array(cells[l][key]) for l in labels}

f_stat, p_anova = stats.f_oneway(*[vals[l] for l in labels])

pairs = list(combinations(labels, 2))
raw_p = {}
bonf_p = {}
for a, b in pairs:
    _, p = stats.ttest_ind(vals[a], vals[b], equal_var=False)
    raw_p[(a, b)] = float(p)
    bonf_p[(a, b)] = min(float(p) * len(pairs), 1.0)


def fmt_stars(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def fmt_p_num(p):
    if p >= 0.001:
        return f"p = {p:.3f}"
    exp = int(np.floor(np.log10(p)))
    coef = p / 10**exp
    return f"p = {coef:.2f}×10$^{{{exp}}}$"


# ── stdout summary ────────────────────────────────────────────────────────────
SEP = "─" * 70
print(SEP)
print("Proliferation rate — Statistical Summary")
print(SEP)
mode_str = "3 BR means" if USE_BR_MEANS else "9 TR values"
print(f"Mode    : USE_BR_MEANS={USE_BR_MEANS}  →  {mode_str} per clone")
print("Test    : One-way ANOVA (omnibus) + pairwise Welch t-test, Bonferroni ×3")
print()
print("Values used (confluency %):")
for l in labels:
    v = vals[l]
    print(
        f"  {l:<6}: {np.array2string(v, precision=2, separator='  ')}"
        f"  mean={v.mean():.2f}  SEM={stats.sem(v):.2f}"
    )
print()
print(f"One-way ANOVA: F={f_stat:.3f}, {fmt_p_num(p_anova)}")
print()
print(f"{'Comparison':<18} {'mean diff':>10}  {'raw p':>8}  {'Bonf p':>8}  sig")
print("─" * 60)
for a, b in pairs:
    diff = vals[a].mean() - vals[b].mean()
    rp = raw_p[(a, b)]
    bp = bonf_p[(a, b)]
    print(
        f"  {a} vs {b:<10}  {diff:>+8.2f}    {rp:>8.4f}    {bp:>8.4f}  {fmt_stars(bp)}"
    )
print(SEP)

# ── Plot ──────────────────────────────────────────────────────────────────────
plt.rcParams.update({"font.family": FONT, "font.size": 11})

fig, ax = plt.subplots(figsize=(5.0, 5.5))

x = np.arange(len(labels))
bar_width = 0.52
means = [vals[l].mean() for l in labels]
sems = [stats.sem(vals[l]) for l in labels]

# bars
ax.bar(
    x,
    means,
    width=bar_width,
    color=[COLORS[l] for l in labels],
    edgecolor="white",
    linewidth=0.8,
    alpha=0.88,
    zorder=3,
)

# error bars
ax.errorbar(
    x,
    means,
    yerr=sems,
    fmt="none",
    color="#222222",
    capsize=4,
    capthick=1.1,
    linewidth=1.1,
    zorder=6,
)

# scatter — BR means with lighter tint
rng = np.random.default_rng(42)
for i, lab in enumerate(labels):
    v = vals[lab]
    jitter = rng.uniform(-0.10, 0.10, len(v))
    ax.scatter(
        x[i] + jitter,
        v,
        color=COLORS_LITE[lab],
        s=40,
        zorder=5,
        edgecolors=COLORS[lab],
        linewidths=0.7,
    )

# ── Significance brackets ─────────────────────────────────────────────────────
y_top = max(m + s for m, s in zip(means, sems))
h = y_top * 0.022
bracket_pairs = [
    ("LNCaP", "C3", y_top * 1.10),
    ("LNCaP", "C6", y_top * 1.22),
    ("C3", "C6", y_top * 1.34),
]
label_idx = {l: i for i, l in enumerate(labels)}

for a, b, yb in bracket_pairs:
    xa, xb = x[label_idx[a]], x[label_idx[b]]
    bp = bonf_p[(a, b)]
    ax.plot([xa, xa, xb, xb], [yb, yb + h, yb + h, yb], color="#444444", linewidth=0.9)
    ax.text(
        (xa + xb) / 2,
        yb + h * 1.4,
        f"{fmt_stars(bp)}  ({fmt_p_num(bp)})",
        ha="center",
        va="bottom",
        fontsize=14,  # fontsize=9,
        fontweight="bold",
        color="#222222",
    )

# ── Axes ──────────────────────────────────────────────────────────────────────
ax.set_xticks(x)
ax.set_xticklabels(sublabel, fontsize=10)
ax.set_ylabel("Cell Confluency (%)", fontsize=12)
ax.set_ylim(0, y_top * 1.58)
ax.yaxis.set_major_locator(mticker.MultipleLocator(20))
ax.yaxis.grid(True, linestyle="--", linewidth=0.45, alpha=0.5, zorder=0)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="y", labelsize=10)

ax.set_title(
    "KDM6B K973Q: Proliferation rate\n(LNCaP parental vs knock-in clones)",
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
    fontsize=14,  # 9,
    color="#666666",
    style="italic",
)

plt.tight_layout()
plt.savefig(FIGURES_DIR / "fig4a_proliferation.pdf", dpi=300, bbox_inches="tight")
plt.savefig(FIGURES_DIR / "fig4a_proliferation.png", dpi=300, bbox_inches="tight")
print(f"\nSaved: {FIGURES_DIR / 'fig4a_proliferation.pdf'}")
print(f"Saved: {FIGURES_DIR / 'fig4a_proliferation.png'}")
